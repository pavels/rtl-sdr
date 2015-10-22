/*
 * rtl-sdr, turns your Realtek RTL2832 based DVB dongle into a SDR receiver
 * Copyright (C) 2012-2013 by Steve Markgraf <steve@steve-m.de>
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

#ifndef __RTL_POWER_H
#define __RTL_POWER_H

#include <rtl-power_export.h>

#ifdef __cplusplus
extern "C" {
#endif

#define MAXIMUM_RATE			3200000
#define MINIMUM_RATE			1000000
#define DEFAULT_TARGET			2400000

#define MAXIMUM_FFT			32

#define DEFAULT_BUF_LENGTH		(1 * 16384)
#define AUTO_GAIN			-100
#define BUFFER_DUMP			(1<<12)

#define MAX_TUNES	4000

#define MAX(x, y) (((x) > (y)) ? (x) : (y))

struct sine_table
{
	int16_t* Sinewave;
	int N_WAVE;
	int LOG2_N_WAVE;
};

struct tuning_state
/* one per tuning range */
{
	int freq;
	int rate;
	int gain;
	int bin_e;
	int16_t *fft_buf;
	int64_t *avg;  /* length == 2^bin_e */
	double *dbm; /* lebgth == crop_i2 - crop_i1 + 1 */
	int samples;
	int downsample;
	int downsample_passes;  /* for the recursive filter */
	int comp_fir_size;
	int peak_hold;
	int linear;
	int bin_spec;
	double crop;
	int crop_i1, crop_i2;
	int freq_low, freq_high;
	//pthread_rwlock_t avg_lock;
	//pthread_mutex_t avg_mutex;
	/* having the iq buffer here is wasteful, but will avoid contention */
	uint8_t *buf8;
	int buf_len;
	//int *comp_fir;
	//pthread_rwlock_t buf_lock;
	//pthread_mutex_t buf_mutex;
	int *window_coefs;
	struct sine_table *sine;  /* points to an element of s_tables */	
};

struct channel_solve
/* details required to find optimal tuning */
{
	int upper, lower, bin_spec;
	int hops, bw_wanted, bw_needed;
	int bin_e, downsample, downsample_passes;
	double crop, crop_tmp;
};

struct misc_settings
{
	int boxcar;
	int comp_fir_size;
	int peak_hold;
	int linear;
	int target_rate;
	double crop;
	int gain;
	double (*window_fn)(int, int);
	int smoothing;	
};

RTLPOWER_API int frequency_range(struct misc_settings *ms, struct tuning_state *tunes, struct channel_solve *c,  int tune_count);
RTLPOWER_API void free_frequency_range(struct tuning_state *tunes, int tune_count);
RTLPOWER_API void scan_tune(rtlsdr_dev_t *dev,struct tuning_state *ts);
RTLPOWER_API void generate_sine_tables(struct sine_table *s_tables,struct tuning_state *tunes, int tune_count);

RTLPOWER_API double rectangle(int i, int length);
RTLPOWER_API double hamming(int i, int length);
RTLPOWER_API double blackman(int i, int length);
RTLPOWER_API double blackman_harris(int i, int length);
RTLPOWER_API double hann_poisson(int i, int length);
RTLPOWER_API double youssef(int i, int length);
RTLPOWER_API double kaiser(int i, int length);
RTLPOWER_API double bartlett(int i, int length);

#endif
