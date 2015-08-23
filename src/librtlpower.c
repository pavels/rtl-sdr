/*
 * rtl-sdr, turns your Realtek RTL2832 based DVB dongle into a SDR receiver
 * Copyright (C) 2012-2014 by Steve Markgraf <steve@steve-m.de>
 * Copyright (C) 2012 by Dimitri Stolnikov <horiz0n@gmx.net>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <stdlib.h>
#include <stdio.h>
#ifndef _WIN32
#include <unistd.h>
#define min(a, b) (((a) < (b)) ? (a) : (b))
#endif

#include <math.h>
#include "rtl-sdr.h"
#include "rtl-power.h"

/* {length, coef, coef, coef}  and scaled by 2^15
   for now, only length 9, optimal way to get +85% bandwidth */
#define CIC_TABLE_MAX 10
int cic_9_tables[][10] = {
	{0,},
	{9, -156,  -97, 2798, -15489, 61019, -15489, 2798,  -97, -156},
	{9, -128, -568, 5593, -24125, 74126, -24125, 5593, -568, -128},
	{9, -129, -639, 6187, -26281, 77511, -26281, 6187, -639, -129},
	{9, -122, -612, 6082, -26353, 77818, -26353, 6082, -612, -122},
	{9, -120, -602, 6015, -26269, 77757, -26269, 6015, -602, -120},
	{9, -120, -582, 5951, -26128, 77542, -26128, 5951, -582, -120},
	{9, -119, -580, 5931, -26094, 77505, -26094, 5931, -580, -119},
	{9, -119, -578, 5921, -26077, 77484, -26077, 5921, -578, -119},
	{9, -119, -577, 5917, -26067, 77473, -26067, 5917, -577, -119},
	{9, -199, -362, 5303, -25505, 77489, -25505, 5303, -362, -199},
};

#if defined(_MSC_VER) && _MSC_VER < 1800
double log2(double n)
{
	return log(n) / log(2.0);
}
#endif

/* FFT based on fix_fft.c by Roberts, Slaney and Bouras
   http://www.jjj.de/fft/fftpage.html
   16 bit ints for everything
   -32768..+32768 maps to -1.0..+1.0
*/

void sine_table(struct sine_table *s_tables, int size)
{
	int i;
	double d;
	struct sine_table *sine;
	if (size > (MAXIMUM_FFT-1)) {
		exit(1);
	}
	sine = &s_tables[size];
	if (sine->LOG2_N_WAVE == size) {
		return;}
	sine->LOG2_N_WAVE = size;
	sine->N_WAVE = 1 << sine->LOG2_N_WAVE;
	sine->Sinewave = malloc(sizeof(int16_t) * sine->N_WAVE*3/4);
	for (i=0; i<sine->N_WAVE*3/4; i++)
	{
		d = (double)i * 2.0 * M_PI / sine->N_WAVE;
		sine->Sinewave[i] = (int)round(32767*sin(d));
		//printf("%i\n", sine->Sinewave[i]);
	}
}

void generate_sine_tables(struct sine_table *s_tables,struct tuning_state *tunes, int tune_count)
{
	struct tuning_state *ts;
	int i;
	for (i=0; i < tune_count; i++) {
		ts = &tunes[i];
		sine_table(s_tables,ts->bin_e);
		ts->sine = &s_tables[ts->bin_e];
		ts->fft_buf = malloc(ts->buf_len * sizeof(int16_t));
	}
}

inline int16_t FIX_MPY(int16_t a, int16_t b)
/* fixed point multiply and scale */
{
	int c = ((int)a * (int)b) >> 14;
	b = c & 0x01;
	return (c >> 1) + b;
}

int fix_fft(int16_t iq[], int m, struct sine_table *sine)
/* interleaved iq[], 0 <= n < 2**m, changes in place */
{
	int mr, nn, i, j, l, k, istep, n, shift;
	int16_t qr, qi, tr, ti, wr, wi;
	n = 1 << m;
	if (n > sine->N_WAVE)
		{return -1;}
	mr = 0;
	nn = n - 1;
	/* decimation in time - re-order data */
	for (m=1; m<=nn; ++m) {
		l = n;
		do
			{l >>= 1;}
		while (mr+l > nn);
		mr = (mr & (l-1)) + l;
		if (mr <= m)
			{continue;}
		// real = 2*m, imag = 2*m+1
		tr = iq[2*m];
		iq[2*m] = iq[2*mr];
		iq[2*mr] = tr;
		ti = iq[2*m+1];
		iq[2*m+1] = iq[2*mr+1];
		iq[2*mr+1] = ti;
	}
	l = 1;
	k = sine->LOG2_N_WAVE-1;
	while (l < n) {
		shift = 1;
		istep = l << 1;
		for (m=0; m<l; ++m) {
			j = m << k;
			wr =  sine->Sinewave[j+sine->N_WAVE/4];
			wi = -sine->Sinewave[j];
			if (shift) {
				wr >>= 1; wi >>= 1;}
			for (i=m; i<n; i+=istep) {
				j = i + l;
				tr = FIX_MPY(wr,iq[2*j]) - FIX_MPY(wi,iq[2*j+1]);
				ti = FIX_MPY(wr,iq[2*j+1]) + FIX_MPY(wi,iq[2*j]);
				qr = iq[2*i];
				qi = iq[2*i+1];
				if (shift) {
					qr >>= 1; qi >>= 1;}
				iq[2*j] = qr - tr;
				iq[2*j+1] = qi - ti;
				iq[2*i] = qr + tr;
				iq[2*i+1] = qi + ti;
			}
		}
		--k;
		l = istep;
	}
	return 0;
}

double rectangle(int i, int length)
{
	return 1.0;
}

double hamming(int i, int length)
{
	double a, b, w, N1;
	a = 25.0/46.0;
	b = 21.0/46.0;
	N1 = (double)(length-1);
	w = a - b*cos(2*i*M_PI/N1);
	return w;
}

double blackman(int i, int length)
{
	double a0, a1, a2, w, N1;
	a0 = 7938.0/18608.0;
	a1 = 9240.0/18608.0;
	a2 = 1430.0/18608.0;
	N1 = (double)(length-1);
	w = a0 - a1*cos(2*i*M_PI/N1) + a2*cos(4*i*M_PI/N1);
	return w;
}

double blackman_harris(int i, int length)
{
	double a0, a1, a2, a3, w, N1;
	a0 = 0.35875;
	a1 = 0.48829;
	a2 = 0.14128;
	a3 = 0.01168;
	N1 = (double)(length-1);
	w = a0 - a1*cos(2*i*M_PI/N1) + a2*cos(4*i*M_PI/N1) - a3*cos(6*i*M_PI/N1);
	return w;
}

double hann_poisson(int i, int length)
{
	double a, N1, w;
	a = 2.0;
	N1 = (double)(length-1);
	w = 0.5 * (1 - cos(2*M_PI*i/N1)) * \
	    pow(M_E, (-a*(double)abs((int)(N1-1-2*i)))/N1);
	return w;
}

double youssef(int i, int length)
/* really a blackman-harris-poisson window, but that is a mouthful */
{
	double a, a0, a1, a2, a3, w, N1;
	a0 = 0.35875;
	a1 = 0.48829;
	a2 = 0.14128;
	a3 = 0.01168;
	N1 = (double)(length-1);
	w = a0 - a1*cos(2*i*M_PI/N1) + a2*cos(4*i*M_PI/N1) - a3*cos(6*i*M_PI/N1);
	a = 0.0025;
	w *= pow(M_E, (-a*(double)abs((int)(N1-1-2*i)))/N1);
	return w;
}

double kaiser(int i, int length)
// todo, become more smart
{
	return 1.0;
}

double bartlett(int i, int length)
{
	double N1, L, w;
	L = (double)length;
	N1 = L - 1;
	w = (i - N1/2) / (L/2);
	if (w < 0) {
		w = -w;}
	w = 1 - w;
	return w;
}

void rms_power(struct tuning_state *ts)
/* for bins between 1MHz and 2MHz */
{
	int i, s;
	uint8_t *buf = ts->buf8;
	int buf_len = ts->buf_len;
	int64_t p, t;
	double dc, err;

	p = t = 0L;
	for (i=0; i<buf_len; i++) {
		s = (int)buf[i] - 127;
		t += (int64_t)s;
		p += (int64_t)(s * s);
	}
	/* correct for dc offset in squares */
	dc = (double)t / (double)buf_len;
	err = t * 2 * dc - dc * dc * buf_len;
	p -= (int64_t)round(err);

	if (!ts->peak_hold) {
		ts->avg[0] += p;
	} else {
		ts->avg[0] = MAX(ts->avg[0], p);
	}
	ts->samples += 1;
}

int solve_giant_bins(struct channel_solve *c)
{
	c->bw_wanted = c->bin_spec;
	c->bw_needed = c->bin_spec;
	c->hops = (c->upper - c->lower) / c->bin_spec;
	c->bin_e = 0;
	c->crop_tmp = 0;
	return 0;
}

int solve_downsample(struct channel_solve *c, int target_rate, int boxcar)
{
	int scan_size, bins_wanted, bins_needed, ds_next, bw;

	scan_size = c->upper - c->lower;
	c->hops = 1;
	c->bw_wanted = scan_size;

	bins_wanted = (int)ceil((double)scan_size / (double)c->bin_spec);
	c->bin_e = (int)ceil(log2(bins_wanted));
	while (1) {
		bins_needed = 1 << c->bin_e;
		c->crop_tmp = (double)(bins_needed - bins_wanted) / (double)bins_needed;
		if (c->crop_tmp >= c->crop) {
			break;}
		c->bin_e++;
	}

	c->downsample = 1;
	c->downsample_passes = 0;
	while (1) {
		bw = (int)((double)scan_size / (1.0 - c->crop_tmp));
		c->bw_needed = bw * c->downsample;

		if (boxcar) {
			ds_next = c->downsample + 1;
		} else {
			ds_next = c->downsample * 2;
		}
		if ((bw * ds_next) > target_rate) {
			break;}

		c->downsample = ds_next;
		if (!boxcar) {
			c->downsample_passes++;}
	}

	return 0;
}

int solve_single(struct channel_solve *c, int target_rate)
{
	int i, scan_size, bins_all, bins_crop, bin_e, bins_2, bw_needed;
	scan_size = c->upper - c->lower;
	bins_all = scan_size / c->bin_spec;
	bins_crop = (int)ceil((double)bins_all * (1.0 + c->crop));
	bin_e = (int)ceil(log2(bins_crop));
	bins_2 = 1 << bin_e;
	bw_needed = bins_2 * c->bin_spec;

	if (bw_needed > target_rate) {
		/* actually multi-hop */
		return 1;}

	c->bw_wanted = scan_size;
	c->bw_needed = bw_needed;
	c->hops = 1;
	c->bin_e = bin_e;
	/* crop will always be bigger than specified crop */
	c->crop_tmp = (double)(bins_2 - bins_all) / (double)bins_2;
	return 0;
}

int solve_hopping(struct channel_solve *c, int target_rate)
{
	int i, scan_size, bins_all, bins_sub, bins_crop, bins_2, min_hops;
	scan_size = c->upper - c->lower;
	min_hops = scan_size / MAXIMUM_RATE - 1;
	if (min_hops < 1) {
		min_hops = 1;}
	/* evenly sized ranges, as close to target_rate as possible */
	for (i=min_hops; i<MAX_TUNES; i++) {
		c->bw_wanted = scan_size / i;
		bins_all = scan_size / c->bin_spec;
		bins_sub = (int)ceil((double)bins_all / (double)i);
		bins_crop = (int)ceil((double)bins_sub * (1.0 + c->crop));
		c->bin_e = (int)ceil(log2(bins_crop));
		bins_2 = 1 << c->bin_e;
		c->bw_needed = bins_2 * c->bin_spec;
		c->crop_tmp = (double)(bins_2 - bins_sub) / (double)bins_2;
		if (c->bw_needed > target_rate) {
			continue;}
		if (c->crop_tmp < c->crop) {
			continue;}
		c->hops = i;
		break;
	}
	return 0;
}

int frequency_range(struct misc_settings *ms, struct tuning_state *tunes, struct channel_solve *c,  int tune_count)
/* flesh out the tunes[] for scanning */
{
	struct tuning_state *ts;
	int r, i, j, buf_len, length, hop_bins, logged_bins, planned_bins;
	int lower_edge, actual_bw, upper_perfect, remainder;

	c->downsample = 1;
	c->downsample_passes = 0;
	c->crop = ms->crop;

	if (ms->target_rate < 2 * MINIMUM_RATE) {
		ms->target_rate = 2 * MINIMUM_RATE;
	}
	if (ms->target_rate > MAXIMUM_RATE) {
		ms->target_rate = MAXIMUM_RATE;
	}
	if ((ms->crop < 0.0) || (ms->crop > 1.0)) {
		exit(1);
	}

	r = -1;
	if (c->bin_spec >= MINIMUM_RATE) {
		fprintf(stderr, "Mode: rms power\n");
		solve_giant_bins(c);
	} else if ((c->upper - c->lower) < MINIMUM_RATE) {
		fprintf(stderr, "Mode: downsampling\n");
		solve_downsample(c, ms->target_rate, ms->boxcar);
	} else if ((c->upper - c->lower) < MAXIMUM_RATE) {
		r = solve_single(c, ms->target_rate);
	} else {
		fprintf(stderr, "Mode: hopping\n");
		solve_hopping(c, ms->target_rate);
	}

	if (r == 0) {
		fprintf(stderr, "Mode: single\n");
	} else if (r == 1) {
		fprintf(stderr, "Mode: hopping\n");
		solve_hopping(c, ms->target_rate);
	}
	c->crop = c->crop_tmp;

	if ((tune_count+c->hops) > MAX_TUNES) {
		fprintf(stderr, "Error: bandwidth too wide.\n");
		exit(1);
	}
	buf_len = 2 * (1<<c->bin_e) * c->downsample;
	if (buf_len < DEFAULT_BUF_LENGTH) {
		buf_len = DEFAULT_BUF_LENGTH;
	}
	/* build the array */
	logged_bins = 0;
	lower_edge = c->lower;
	planned_bins = (c->upper - c->lower) / c->bin_spec;
	for (i=0; i < c->hops; i++) {
		ts = &tunes[tune_count + i];
		/* copy common values */
		ts->rate = c->bw_needed;
		ts->gain = ms->gain;
		ts->bin_e = c->bin_e;
		ts->samples = 0;
		ts->bin_spec = c->bin_spec;
		ts->crop = c->crop;
		ts->downsample = c->downsample;
		ts->downsample_passes = c->downsample_passes;
		ts->comp_fir_size = ms->comp_fir_size;
		ts->peak_hold = ms->peak_hold;
		ts->linear = ms->linear;
		ts->avg = (int64_t*)malloc((1<<c->bin_e) * sizeof(int64_t));
		if (!ts->avg) {
			fprintf(stderr, "Error: malloc->\n");
			exit(1);
		}
		for (j=0; j<(1<<c->bin_e); j++) {
			ts->avg[j] = 0L;
		}
		ts->buf8 = (uint8_t*)malloc(buf_len * sizeof(uint8_t));
		if (!ts->buf8) {
			fprintf(stderr, "Error: malloc->\n");
			exit(1);
		}
		ts->buf_len = buf_len;
		length = 1 << c->bin_e;
		ts->window_coefs = malloc(length * sizeof(int));
		for (j=0; j<length; j++) {
			ts->window_coefs[j] = (int)(256*ms->window_fn(j, length));
		}
		/* calculate unique values */
		ts->freq_low = lower_edge;
		hop_bins = c->bw_wanted / c->bin_spec;
		actual_bw = hop_bins * c->bin_spec;
		ts->freq_high = lower_edge + actual_bw;
		upper_perfect = c->lower + (i+1) * c->bw_wanted;
		if (ts->freq_high + c->bin_spec <= upper_perfect) {
			hop_bins += 1;
			actual_bw = hop_bins * c->bin_spec;
			ts->freq_high = lower_edge + actual_bw;
		}
		remainder = planned_bins - logged_bins - hop_bins;
		if (i == c->hops-1 && remainder > 0) {
			hop_bins += remainder;
			actual_bw = hop_bins * c->bin_spec;
			ts->freq_high = lower_edge + actual_bw;
		}
		logged_bins += hop_bins;
		ts->crop_i1 = (length - hop_bins) / 2;
		ts->crop_i2 = ts->crop_i1 + hop_bins - 1;
		ts->dbm = malloc((ts->crop_i2 - ts-> crop_i1 + 1) * sizeof(double));
		ts->freq = (lower_edge - ts->crop_i1 * c->bin_spec) + c->bw_needed/(2*c->downsample);
		/* prep for next hop */
		lower_edge = ts->freq_high;
	}
	return tune_count + c->hops;
}

void retune(rtlsdr_dev_t *dev, int freq)
{
	uint8_t dump[BUFFER_DUMP];
	int f, n_read;
	f = (int)rtlsdr_get_center_freq(dev);
	if (f == freq) {
		return;}
	rtlsdr_set_center_freq(dev, (uint32_t)freq);
	/* wait for settling and flush buffer */
	usleep(5000);
	rtlsdr_read_sync(dev, &dump, BUFFER_DUMP, &n_read);
	if (n_read != BUFFER_DUMP) {
	 	fprintf(stderr, "Error: bad retune.\n");}
}

void rerate(rtlsdr_dev_t *dev, int rate)
{
	uint32_t r;
	r = rtlsdr_get_sample_rate(dev);
	if ((int)r == rate) {
		return;}
	rtlsdr_set_sample_rate(dev, (uint32_t)rate);
}

void regain(rtlsdr_dev_t *dev, int gain)
{
	int g;
	g = rtlsdr_get_tuner_gain(dev);
	if (g == gain) {
		return;}
	if (gain == AUTO_GAIN) {
		/* switch to auto */
		rtlsdr_set_tuner_gain_mode(dev, 0);
		return;
	}
	if (g == AUTO_GAIN) {
		/* switch to manual */
		rtlsdr_set_tuner_gain_mode(dev, 1);
	}
	rtlsdr_set_tuner_gain(dev, gain);
}

void fifth_order(int16_t *data, int length)
/* for half of interleaved data */
{
	int i;
	int a, b, c, d, e, f;
	a = data[0];
	b = data[2];
	c = data[4];
	d = data[6];
	e = data[8];
	f = data[10];
	/* a downsample should improve resolution, so don't fully shift */
	/* ease in instead of being stateful */
	data[0] = ((a+b)*10 + (c+d)*5 + d + f) >> 4;
	data[2] = ((b+c)*10 + (a+d)*5 + e + f) >> 4;
	data[4] = (a + (b+e)*5 + (c+d)*10 + f) >> 4;
	for (i=12; i<length; i+=4) {
		a = c;
		b = d;
		c = e;
		d = f;
		e = data[i-2];
		f = data[i];
		data[i/2] = (a + (b+e)*5 + (c+d)*10 + f) >> 4;
	}
}

void remove_dc(int16_t *data, int length)
/* works on interleaved data */
{
	int i;
	int16_t ave;
	int64_t sum = 0L;
	for (i=0; i < length; i+=2) {
		sum += data[i];
	}
	ave = (int16_t)(sum / (int64_t)(length));
	if (ave == 0) {
		return;}
	for (i=0; i < length; i+=2) {
		data[i] -= ave;
	}
}

void generic_fir(int16_t *data, int length, int *fir)
/* Okay, not at all generic->  Assumes length 9, fix that eventually. */
{
	int d, temp, sum;
	int hist[9] = {0,};
	/* cheat on the beginning, let it go unfiltered */
	for (d=0; d<18; d+=2) {
		hist[d/2] = data[d];
	}
	for (d=18; d<length; d+=2) {
		temp = data[d];
		sum = 0;
		sum += (hist[0] + hist[8]) * fir[1];
		sum += (hist[1] + hist[7]) * fir[2];
		sum += (hist[2] + hist[6]) * fir[3];
		sum += (hist[3] + hist[5]) * fir[4];
		sum +=            hist[4]  * fir[5];
		data[d] = (int16_t)(sum >> 15) ;
		hist[0] = hist[1];
		hist[1] = hist[2];
		hist[2] = hist[3];
		hist[3] = hist[4];
		hist[4] = hist[5];
		hist[5] = hist[6];
		hist[6] = hist[7];
		hist[7] = hist[8];
		hist[8] = temp;
	}
}

void downsample_iq(int16_t *data, int length)
{
	fifth_order(data, length);
	//remove_dc(data, length);
	fifth_order(data+1, length-1);
	//remove_dc(data+1, length-1);
}

int64_t real_conj(int16_t real, int16_t imag)
/* real(n * conj(n)) */
{
	return ((int64_t)real*(int64_t)real + (int64_t)imag*(int64_t)imag);
}

void scan_tune(rtlsdr_dev_t *dev, struct tuning_state *ts)
{
	int i, j, j2, n_read, offset, bin_e, bin_len, buf_len, ds, ds_p;
	int64_t tmp;
	int32_t w;
	int16_t *fft_buf;
	double dbm;
	bin_e = ts->bin_e;
	bin_len = 1 << bin_e;

	for (i=0; i<bin_len; i++) {
		ts->avg[i] = 0L;
	}
	ts->samples = 0;

	buf_len = ts->buf_len;
	fft_buf = ts->fft_buf;
	regain(dev, ts->gain);
	rerate(dev, ts->rate);
	retune(dev, ts->freq);
	rtlsdr_read_sync(dev, ts->buf8, buf_len, &n_read);
	if (n_read != buf_len) {
		fprintf(stderr, "Error: dropped samples.\n");}
	/* rms */
	if (bin_len == 1) {
		rms_power(ts);
		return;
	}
	/* prep for fft */
	for (j=0; j<buf_len; j++) {
		fft_buf[j] = (int16_t)ts->buf8[j] - 127;
	}
	ds = ts->downsample;
	ds_p = ts->downsample_passes;
	if (ds_p) {  /* recursive */
		for (j=0; j < ds_p; j++) {
			downsample_iq(fft_buf, buf_len >> j);
		}
		/* droop compensation */
		if (ts->comp_fir_size == 9 && ds_p <= CIC_TABLE_MAX) {
			generic_fir(fft_buf, buf_len >> j, cic_9_tables[ds_p]);
			generic_fir(fft_buf+1, (buf_len >> j)-1, cic_9_tables[ds_p]);
		}
	} else if (ds > 1) {  /* boxcar */
		j=2, j2=0;
		while (j < buf_len) {
			fft_buf[j2]   += fft_buf[j];
			fft_buf[j2+1] += fft_buf[j+1];
			fft_buf[j] = 0;
			fft_buf[j+1] = 0;
			j += 2;
			if (j % (ds*2) == 0) {
				j2 += 2;}
		}
	}
	remove_dc(fft_buf, buf_len / ds);
	remove_dc(fft_buf+1, (buf_len / ds) - 1);
	/* window function and fft */
	for (offset=0; offset<(buf_len/ds); offset+=(2*bin_len)) {
		// todo, let rect skip this
		for (j=0; j<bin_len; j++) {
			w =  (int32_t)fft_buf[offset+j*2];
			w *= (int32_t)(ts->window_coefs[j]);
			//w /= (int32_t)(ds);
			fft_buf[offset+j*2]   = (int16_t)w;
			w =  (int32_t)fft_buf[offset+j*2+1];
			w *= (int32_t)(ts->window_coefs[j]);
			//w /= (int32_t)(ds);
			fft_buf[offset+j*2+1] = (int16_t)w;
		}
		fix_fft(fft_buf+offset, bin_e, ts->sine);
		if (!ts->peak_hold) {
			for (j=0; j<bin_len; j++) {
				ts->avg[j] += real_conj(fft_buf[offset+j*2], fft_buf[offset+j*2+1]);
			}
		} else {
			for (j=0; j<bin_len; j++) {
				ts->avg[j] = MAX(real_conj(fft_buf[offset+j*2], fft_buf[offset+j*2+1]), ts->avg[j]);
			}
		}
		ts->samples += ds;
	}

	/* fix FFT stuff quirks */
	if (ts->bin_e > 0) {
		/* nuke DC component (not effective for all windows) */
		ts->avg[0] = ts->avg[1];
		/* FFT is translated by 180 degrees */
		for (i=0; i<bin_len/2; i++) {
			tmp = ts->avg[i];
			ts->avg[i] = ts->avg[i+bin_len/2];
			ts->avg[i+bin_len/2] = tmp;
		}
	}

	// something seems off with the dbm math
	for (i=ts->crop_i1; i<=ts->crop_i2; i++) {
		dbm  = (double)ts->avg[i];
		dbm /= (double)ts->rate;
		if (!ts->peak_hold) {
			dbm /= (double)ts->samples;
		}
		if (!ts->linear) {
			dbm = 10 * log10(dbm);
		}
		ts->dbm[i - ts->crop_i1] = dbm;
	}

}
