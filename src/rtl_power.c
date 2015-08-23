/*
 * rtl-sdr, turns your Realtek RTL2832 based DVB dongle into a SDR receiver
 * Copyright (C) 2012 by Steve Markgraf <steve@steve-m.de>
 * Copyright (C) 2012 by Hoernchen <la@tfc-server.de>
 * Copyright (C) 2012 by Kyle Keen <keenerd@gmail.com>
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


/*
 * rtl_power: general purpose FFT integrator
 * -f low_freq:high_freq:max_bin_size
 * -i seconds
 * outputs CSV
 * time, low, high, step, db, db, db ...
 * db optional?  raw output might be better for noise correction
 * todo:
 *	threading
 *	randomized hopping
 *	noise correction
 *	continuous IIR
 *	general astronomy usefulness
 *	multiple dongles
 *	multiple FFT workers
 *	check edge cropping for off-by-one and rounding errors
 *	1.8MS/s for hiding xtal harmonics
 */

#include <errno.h>
#include <signal.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#ifndef _WIN32
#include <unistd.h>
#else
#include <windows.h>
#include <fcntl.h>
#include <io.h>
#include "getopt/getopt.h"
#define usleep(x) Sleep(x/1000)
#if defined(_MSC_VER) && _MSC_VER < 1800
#define round(x) (x > 0.0 ? floor(x + 0.5): ceil(x - 0.5))
#endif
#define _USE_MATH_DEFINES
#endif

#include <math.h>
#include <pthread.h>
#include <libusb.h>

#include "rtl-sdr.h"
#include "rtl-power.h"
#include "convenience/convenience.h"

#define MAX(x, y) (((x) > (y)) ? (x) : (y))

#define DEFAULT_BUF_LENGTH		(1 * 16384)
#define AUTO_GAIN			-100
#define BUFFER_DUMP			(1<<12)

#define MAXIMUM_RATE			3200000
#define MINIMUM_RATE			1000000
#define DEFAULT_TARGET			2400000

#define MAXIMUM_FFT			32

enum time_modes { VERBOSE_TIME, EPOCH_TIME };
enum time_modes time_mode;

static volatile int do_exit = 0;
static rtlsdr_dev_t *dev = NULL;
FILE *file;

struct sine_table s_tables[MAXIMUM_FFT];

struct tuning_state tunes[MAX_TUNES];
int tune_count = 0;

void usage(void)
{
	fprintf(stderr,
		"rtl_power, a simple FFT logger for RTL2832 based DVB-T receivers\n\n"
		"Use:\trtl_power -f freq_range [-options] [-f freq2 -opts2] [filename]\n"
		"\t-f lower:upper:bin_size [Hz]\n"
		"\t  valid range for bin_size is 1Hz - 2.8MHz\n"
		"\t  multiple frequency ranges are supported\n"
		"\t[-i integration_interval (default: 10 seconds)]\n"
		"\t  buggy if a full sweep takes longer than the interval\n"
		"\t[-1 enables single-shot mode (default: off)]\n"
		"\t[-e exit_timer (default: off/0)]\n"
		//"\t[-s avg/iir smoothing (default: avg)]\n"
		//"\t[-t threads (default: 1)]\n"
		"\t[-d device_index (default: 0)]\n"
		"\t[-g tuner_gain (default: automatic)]\n"
		"\t[-p ppm_error (default: 0)]\n"
		"\tfilename (a '-' dumps samples to stdout)\n"
		"\t  omitting the filename also uses stdout\n"
		"\n"
		"Experimental options:\n"
		"\t[-w window (default: rectangle)]\n"
		"\t  hamming, blackman, blackman-harris, hann-poisson, bartlett, youssef\n"
		// kaiser
		"\t[-c crop_percent (default: 0%% suggested: 20%%)]\n"
		"\t  discards data at the edges, 100%% discards everything\n"
		"\t  has no effect for bins larger than 1MHz\n"
		"\t  this value is a minimum crop size, more may be discarded\n"
		"\t[-F fir_size (default: disabled)]\n"
		"\t  enables low-leakage downsample filter,\n"
		"\t  fir_size can be 0 or 9.  0 has bad roll off,\n"
		"\t  try -F 0 with '-c 50%%' to hide the roll off\n"
		"\t[-r max_sample_rate (default: 2.4M)]\n"
		"\t  possible values are 2M to 3.2M\n"
		"\t[-E enables epoch timestamps (default: off/verbose)]\n"
		"\t[-P enables peak hold (default: off/averaging)]\n"
		"\t[-L enable linear output (default: off/dB)]\n"
		"\t[-D direct_sampling_mode, 0 (default/off), 1 (I), 2 (Q), 3 (no-mod)]\n"
		"\t[-O enable offset tuning (default: off)]\n"
		"\n"
		"CSV FFT output columns:\n"
		"\tdate, time, Hz low, Hz high, Hz step, samples, dbm, dbm, ...\n\n"
		"Examples:\n"
		"\trtl_power -f 88M:108M:125k fm_stations.csv\n"
		"\t  creates 160 bins across the FM band,\n"
		"\t  individual stations should be visible\n"
		"\trtl_power -f 100M:1G:1M -i 5m -1 survey.csv\n"
		"\t  a five minute low res scan of nearly everything\n"
		"\trtl_power -f ... -i 15m -1 log.csv\n"
		"\t  integrate for 15 minutes and exit afterwards\n"
		"\trtl_power -f ... -e 1h | gzip > log.csv.gz\n"
		"\t  collect data for one hour and compress it on the fly\n\n"
		"\tIf you have issues writing +2GB logs on a 32bit platform\n"
		"\tuse redirection (rtl_power ... > filename.csv) instead\n\n"
		"Convert CSV to a waterfall graphic with:\n"
		"  https://github.com/keenerd/rtl-sdr-misc/blob/master/heatmap/heatmap.py \n"
		"More examples at http://kmkeen.com/rtl-power/\n");
	exit(1);
}

void multi_bail(void)
{
	if (do_exit == 1)
	{
		fprintf(stderr, "Signal caught, finishing scan pass.\n");
	}
	if (do_exit >= 2)
	{
		fprintf(stderr, "Signal caught, aborting immediately.\n");
	}
}

#ifdef _WIN32
BOOL WINAPI
sighandler(int signum)
{
	if (CTRL_C_EVENT == signum) {
		do_exit++;
		multi_bail();
		return TRUE;
	}
	return FALSE;
}
#else
static void sighandler(int signum)
{
	do_exit++;
	multi_bail();
}
#endif

/* more cond dumbness */
#define safe_cond_signal(n, m) pthread_mutex_lock(m); pthread_cond_signal(n); pthread_mutex_unlock(m)
#define safe_cond_wait(n, m) pthread_mutex_lock(m); pthread_cond_wait(n, m); pthread_mutex_unlock(m)

/* todo, add errors to parse_freq, solve_foo */

int parse_frequency(char *arg, struct channel_solve *c)
{
	char *start, *stop, *step;
	/* hacky string parsing */
	start = arg;
	stop = strchr(start, ':') + 1;
	stop[-1] = '\0';
	step = strchr(stop, ':') + 1;
	step[-1] = '\0';
	c->lower = (int)atofs(start);
	c->upper = (int)atofs(stop);
	c->bin_spec = (int)atofs(step);
	stop[-1] = ':';
	step[-1] = ':';
	return 0;
}

void scanner(void)
{
	int i;
	for (i=0; i<tune_count; i++) {
		if (do_exit >= 2)
			{return;}
		scan_tune(dev,&tunes[i]);
	}
}

void csv_dbm(struct tuning_state *ts)
{
	int i;
	char *sep = ", ";

	/* Hz low, Hz high, Hz step, samples, dbm, dbm, ... */
	fprintf(file, "%i, %i, %i, %i, ", ts->freq_low, ts->freq_high,
		ts->bin_spec, ts->samples);

	for (i=0; i<=(ts->crop_i2 - ts->crop_i1); i++) {
		if (i == ts->crop_i2) {
			sep = "\n";
		}
		if (ts->linear) {
			fprintf(file, "%.5g%s", ts->dbm[i], sep);
		} else {
			fprintf(file, "%.2f%s", ts->dbm[i], sep);
		}
	}
}

void init_misc(struct misc_settings *ms)
{
	ms->target_rate = DEFAULT_TARGET;
	ms->boxcar = 1;
	ms->comp_fir_size = 0;
	ms->crop = 0.0;
	ms->gain = AUTO_GAIN;
	ms->window_fn = rectangle;
	ms->smoothing = 0;
	ms->peak_hold = 0;
	ms->linear = 0;
	time_mode = VERBOSE_TIME;
}

int main(int argc, char **argv)
{
#ifndef _WIN32
	struct sigaction sigact;
#endif
	char *filename = NULL;
	int i, r, opt;
	int f_set = 0;
	int dev_index = 0;
	char dev_label[255];
	int dev_given = 0;
	int ppm_error = 0;
	int custom_ppm = 0;
	int interval = 10;
	int fft_threads = 1;
	int single = 0;
	int direct_sampling = 0;
	int offset_tuning = 0;
	char *freq_optarg;
	time_t next_tick;
	time_t time_now;
	time_t exit_time = 0;
	char t_str[50];
	struct tm *cal_time;
	struct misc_settings ms;
	struct channel_solve c;
	freq_optarg = "";
	init_misc(&ms);
	strcpy(dev_label, "DEFAULT");

	while ((opt = getopt(argc, argv, "f:i:s:r:t:d:g:p:e:w:c:F:1EPLD:Oh")) != -1) {
		switch (opt) {
		case 'f': // lower:upper:bin_size
			if (f_set) {
				parse_frequency(freq_optarg, &c);
				tune_count = frequency_range(&ms, tunes, &c, tune_count);
			}
			freq_optarg = strdup(optarg);
			f_set = 1;
			break;
		case 'd':
			dev_index = verbose_device_search(optarg);
			strncpy(dev_label, optarg, 255);
			dev_given = 1;
			break;
		case 'g':
			ms.gain = (int)(atof(optarg) * 10);
			break;
		case 'c':
			ms.crop = atofp(optarg);
			break;
		case 'i':
			interval = (int)round(atoft(optarg));
			break;
		case 'e':
			exit_time = (time_t)((int)round(atoft(optarg)));
			break;
		case 's':
			if (strcmp("avg",  optarg) == 0) {
				ms.smoothing = 0;}
			if (strcmp("iir",  optarg) == 0) {
				ms.smoothing = 1;}
			break;
		case 'w':
			if (strcmp("rectangle",  optarg) == 0) {
				ms.window_fn = rectangle;}
			if (strcmp("hamming",  optarg) == 0) {
				ms.window_fn = hamming;}
			if (strcmp("blackman",  optarg) == 0) {
				ms.window_fn = blackman;}
			if (strcmp("blackman-harris",  optarg) == 0) {
				ms.window_fn = blackman_harris;}
			if (strcmp("hann-poisson",  optarg) == 0) {
				ms.window_fn = hann_poisson;}
			if (strcmp("youssef",  optarg) == 0) {
				ms.window_fn = youssef;}
			if (strcmp("kaiser",  optarg) == 0) {
				ms.window_fn = kaiser;}
			if (strcmp("bartlett",  optarg) == 0) {
				ms.window_fn = bartlett;}
			break;
		case 't':
			fft_threads = atoi(optarg);
			break;
		case 'p':
			ppm_error = atoi(optarg);
			custom_ppm = 1;
			break;
		case 'r':
			ms.target_rate = (int)atofs(optarg);
			break;
		case '1':
			single = 1;
			break;
		case 'E':
			time_mode = EPOCH_TIME;
			break;
		case 'P':
			ms.peak_hold = 1;
			break;
		case 'L':
			ms.linear = 1;
			break;
		case 'D':
			direct_sampling = atoi(optarg);
			break;
		case 'O':
			offset_tuning = 1;
			break;
		case 'F':
			ms.boxcar = 0;
			ms.comp_fir_size = atoi(optarg);
			break;
		case 'h':
		default:
			usage();
			break;
		}
	}

	if (!f_set) {
		fprintf(stderr, "No frequency range provided.\n");
		exit(1);
	}

	parse_frequency(freq_optarg, &c);
	tune_count = frequency_range(&ms, tunes, &c, tune_count);

	if (tune_count == 0) {
		usage();}

	if (argc <= optind) {
		filename = "-";
	} else {
		filename = argv[optind];
	}

	if (interval < 1) {
		interval = 1;}

	fprintf(stderr, "Reporting every %i seconds\n", interval);

	if (!dev_given) {
		dev_index = verbose_device_search("0");
	}

	if (dev_index < 0) {
		exit(1);
	}

	r = rtlsdr_open(&dev, (uint32_t)dev_index);
	if (r < 0) {
		fprintf(stderr, "Failed to open rtlsdr device #%d.\n", dev_index);
		exit(1);
	}
#ifndef _WIN32
	sigact.sa_handler = sighandler;
	sigemptyset(&sigact.sa_mask);
	sigact.sa_flags = 0;
	sigaction(SIGINT, &sigact, NULL);
	sigaction(SIGTERM, &sigact, NULL);
	sigaction(SIGQUIT, &sigact, NULL);
	sigaction(SIGPIPE, &sigact, NULL);
	signal(SIGPIPE, SIG_IGN);
#else
	SetConsoleCtrlHandler( (PHANDLER_ROUTINE) sighandler, TRUE );
#endif

	if (direct_sampling) {
		verbose_direct_sampling(dev, direct_sampling);
	}

	if (offset_tuning) {
		verbose_offset_tuning(dev);
	}

	/* Set the tuner gain */
	for (i=0; i<tune_count; i++) {
		if (tunes[i].gain == AUTO_GAIN) {
			continue;}
		tunes[i].gain = nearest_gain(dev, tunes[i].gain);
	}
	if (ms.gain == AUTO_GAIN) {
		verbose_auto_gain(dev);
	} else {
		ms.gain = nearest_gain(dev, ms.gain);
		verbose_gain_set(dev, ms.gain);
	}

	if (!custom_ppm) {
		verbose_ppm_eeprom(dev, &ppm_error);
	}
	verbose_ppm_set(dev, ppm_error);

	if (strcmp(filename, "-") == 0) { /* Write log to stdout */
		file = stdout;
#ifdef _WIN32
		// Is this necessary?  Output is ascii.
		_setmode(_fileno(file), _O_BINARY);
#endif
	} else {
		file = fopen(filename, "wb");
		if (!file) {
			fprintf(stderr, "Failed to open %s\n", filename);
			exit(1);
		}
	}
	generate_sine_tables(s_tables, tunes, tune_count);

	/* Reset endpoint before we start reading from it (mandatory) */
	verbose_reset_buffer(dev);

	/* actually do stuff */
	rtlsdr_set_sample_rate(dev, (uint32_t)tunes[0].rate);
	next_tick = time(NULL) + interval;
	if (exit_time) {
		exit_time = time(NULL) + exit_time;}
	while (!do_exit) {
		scanner();
		time_now = time(NULL);
		if (time_now < next_tick) {
			continue;}
		// time, Hz low, Hz high, Hz step, samples, dbm, dbm, ...
		cal_time = localtime(&time_now);
		if (time_mode == VERBOSE_TIME) {
			strftime(t_str, 50, "%Y-%m-%d, %H:%M:%S", cal_time);
		}
		if (time_mode == EPOCH_TIME) {
			snprintf(t_str, 50, "%u, %s", (unsigned)time_now, dev_label);
		}
		for (i=0; i<tune_count; i++) {
			fprintf(file, "%s, ", t_str);
			csv_dbm(&tunes[i]);
		}
		fflush(file);
		while (time(NULL) >= next_tick) {
			next_tick += interval;}
		if (single) {
			do_exit = 1;}
		if (exit_time && time(NULL) >= exit_time) {
			do_exit = 1;}
	}

	/* clean up */

	if (do_exit) {
		fprintf(stderr, "\nUser cancel, exiting...\n");}
	else {
		fprintf(stderr, "\nLibrary error %d, exiting...\n", r);}

	if (file != stdout) {
		fclose(file);}

	rtlsdr_close(dev);
	//free(fft_buf);
	//for (i=0; i<tune_count; i++) {
	//	free(tunes[i].avg);
	//	free(tunes[i].buf8);
	//}
	return r >= 0 ? r : -r;
}

// vim: tabstop=8:softtabstop=8:shiftwidth=8:noexpandtab
