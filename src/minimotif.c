/*
 *   minimotif: A small super-fast DNA/RNA motif scanner
 *   Copyright (C) 2022  Benjamin Jean-Marie Tremblay
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <locale.h>
#include <getopt.h>
#include <math.h>
#include <limits.h>
#include <time.h>
#include <pthread.h>
#include <zlib.h>
#include "kseq.h"

KSEQ_INIT(gzFile, gzread)

#define MINIMOTIF_VERSION                  "1.0"
#define MINIMOTIF_YEAR                      2022

/* These defaults can be safely changed. The only effects of doing so will be
 * on performance. Depending on whether your motifs are extremely large, or
 * you are working with extreme imbalances in your background, changing these
 * might be a good idea to get around internal size limits. Do be careful that
 * any changes don't lead to int overflows or going out of bounds.
 */

/* Max stored size of motif names.
 */
#define MAX_NAME_SIZE             ((size_t) 256)

/* The motif cannot be larger than 50 positions. This ensures no integer
 * overflow occurs if there are too many non-standard letters (as the current
 * solution to dealing with non-standard letters is to assign them a score of
 * -10,000,000; for 50 positions, this makes for a possible min score of
 * -500,000,000, or approx 1/4 of the way until INT_MIN [-2,147,483,648]).
 * This value is also used when performing memory allocation for new motifs,
 * regardless of the actual size of the motif (thus sticking to lower values
 * may be better for performance).
 * Note: Motif size cannot exceed INT_MAX, since it has to be casted to an int
 * in order to print the match (see score_seq). But realistically having such
 * a big motif will cause the max score to overflow long before then.
 */
#define MAX_MOTIF_SIZE            ((size_t) 250)    /* 5 rows per position */
#define AMBIGUITY_SCORE                -10000000

/* No bkg prob can be smaller than 0.001, to allow for a relatively small
 * max CDF size. (PWM scores are multiplied by 1000 and used as ints.)
 *     max score: (int) 1000*log2(1/0.001)      =>   9,965
 *     min score: (int) 1000*log2(0.001/0.997)  =>  -9,961
 *     cdf size:        (9965+9961)*50          => 996,300
 */
#define MIN_BKG_VALUE                      0.001
#define MAX_CDF_SIZE          ((size_t) 2097152)
#define PWM_INT_MULTIPLIER                1000.0    /* Needs to be a double */

/* Max size of the parsed -b char array.
 */
#define USER_BKG_MAX_SIZE         ((size_t) 256)

/* Max size of the parsed MEME background probabilities.
 */
#define MEME_BKG_MAX_SIZE         ((size_t) 256)

/* Max size of PCM/PPM values in parsed motifs.
 */
#define MOTIF_VALUE_MAX_CHAR      ((size_t) 256)

/* Max size of sequence names.
 */
#define SEQ_NAME_MAX_CHAR         ((size_t) 512)

/* Front-facing defaults.
 */
#define DEFAULT_NSITES                      1000
#define DEFAULT_PVALUE                   0.00001
#define DEFAULT_PSEUDOCOUNT                    1

/* Minimum amount of additional memory to request when reading sequences and
 * more memory is needed. Current default: request memory in 1/2 MB chunks.
 */
#define SEQ_REALLOC_SIZE                  524288

#define VEC_ADD(VEC, X, VEC_LEN)                                \
  do {                                                          \
    for (size_t Xi = 0; Xi < VEC_LEN; Xi++) VEC[Xi] += X;       \
  } while (0)

#define VEC_DIV(VEC, X, VEC_LEN)                                \
  do {                                                          \
    for (size_t Xi = 0; Xi < VEC_LEN; Xi++) VEC[Xi] /= X;       \
  } while (0)

#define VEC_SUM(VEC, SUM_RES, VEC_LEN)                          \
  do {                                                          \
    SUM_RES = 0;                                                \
    for (size_t Xi = 0; Xi < VEC_LEN; Xi++) SUM_RES += VEC[Xi]; \
  } while (0)

#define VEC_MIN(VEC, MIN_RES, VEC_LEN)                          \
  do {                                                          \
    MIN_RES = VEC[0];                                           \
    for (size_t Xi = 1; Xi < VEC_LEN; Xi++) {                   \
      if (VEC[Xi] < MIN_RES) MIN_RES = VEC[Xi];                 \
    }                                                           \
  } while (0)

#define ERASE_ARRAY(ARR, LEN) memset(ARR, 0, sizeof(ARR[0]) * (LEN))

#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define MAX(x, y) (((x) > (y)) ? (x) : (y))

/* Size of progress bar.
 */
#define PROGRESS_BAR_WIDTH                    60
#define PROGRESS_BAR_STRING                                     \
  "============================================================"

void usage(void) {
  printf(
    "minimotif v%s  Copyright (C) %d  Benjamin Jean-Marie Tremblay              \n"
    "                                                                              \n"
    "Usage:  minimotif [options] [ -m motifs.txt | -1 CONSENSUS ] -s sequences.fa  \n"
    "                                                                              \n"
    " -m <str>   Filename of text file containing motifs. Acceptable formats: MEME,\n"
    "            JASPAR, HOMER, HOCOMOCO (PCM). Must be 1-%zu bases wide.           \n"
    " -1 <str>   Instead of -m, scan a single consensus sequence. Ambiguity letters\n"
    "            are allowed. Must be 1-%zu bases wide. The -b, -t, -0, -p, and -n \n"
    "            flags are unused.                                                 \n"
    " -s <str>   Filename of fast(a|q)-formatted file containing DNA/RNA sequences \n"
    "            to scan. Can be gzipped. Use '-' for stdin. Omitting -s will cause\n"
    "            minimotif to print the parsed motifs instead of scanning.         \n"
    "            Alternatively, solely providing -s and not -m/-1 will cause       \n"
    "            minimotif to return sequence stats. Non-standard characters (i.e. \n"
    "            other than ACGTU) will be read but are treated as gaps during     \n"
    "            scanning.                                                         \n"
    " -o <str>   Filename to output results. By default output goes to stdout.     \n"
    " -b <dbl>   Comma-separated background probabilities for A,C,G,T|U. By default\n"
    "            the background probability values from the motif file (MEME only) \n"
    "            are used, or a uniform background is assumed. Used in PWM         \n"
    "            generation.                                                       \n"
    " -f         Only scan the forward strand.                                     \n"
    " -t <dbl>   Threshold P-value. Default: %g.                         \n"
    " -0         Instead of using a threshold, simply report all hits with a score \n"
    "            of zero or greater. Useful for manual filtering.                  \n"
    " -p <int>   Pseudocount for PWM generation. Default: %d. Must be a positive    \n"
    "            integer.                                                          \n"
    " -n <int>   Number of motif sites used in PWM generation. Default: %d.         \n"
    " -d         Deduplicate motif/sequence names. Default: abort. Duplicates will \n"
    "            have the motif/sequence numbers appended.                         \n"
    " -r         Don't trim motif (HOCOMOCO/JASPAR only) and sequence names to the \n"
    "            first word.                                                       \n"
    " -l         Deactivate low memory mode. Normally only a single sequence is    \n"
    "            stored in memory at a time. Setting this flag allows the program  \n"
    "            to store the entire input into memory at a time, which can help   \n"
    "            with performance when scanning multiple motifs. Note that this    \n"
    "            flag is automatically set when reading sequences from stdin, and  \n"
    "            when multithreading is enabled.                                   \n"
    " -j <int>   Number of threads minimotif can use to scan. Default: 1. Note that\n"
    "            increasing this number will also increase memory usage slightly.  \n"
    "            The number of threads is limited by the number of motifs being    \n"
    "            scanned.                                                          \n"
    " -g         Print a progress bar during scanning. This turns off some of the  \n"
    "            messages printed by -w. Note that it's only useful if there is    \n"
    "            more than one input motif.                                        \n"
    " -v         Verbose mode.                                                     \n"
    " -w         Very verbose mode.                                                \n"
    " -h         Print this help message.                                          \n"
    , MINIMOTIF_VERSION, MINIMOTIF_YEAR, MAX_MOTIF_SIZE / 5, MAX_MOTIF_SIZE / 5,
      DEFAULT_PVALUE, DEFAULT_PSEUDOCOUNT, DEFAULT_NSITES
  );
}

const unsigned char char2index[] = {
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
  4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
  4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4 
};

size_t char_counts[256];

const double consensus2probs[] = {
  1.0,   0.0,   0.0,   0.0,        /*  0. A */
  0.0,   1.0,   0.0,   0.0,        /*  1. C */
  0.0,   0.0,   1.0,   0.0,        /*  2. G */
  0.0,   0.0,   0.0,   1.0,        /*  3. T */
  0.0,   0.5,   0.0,   0.5,        /*  4. Y */
  0.5,   0.0,   0.5,   0.0,        /*  5. R */
  0.5,   0.0,   0.0,   0.5,        /*  6. W */
  0.0,   0.5,   0.5,   0.0,        /*  7. S */
  0.0,   0.0,   0.5,   0.5,        /*  8. K */
  0.5,   0.5,   0.0,   0.0,        /*  9. M */
  0.333, 0.0,   0.333, 0.333,      /* 10. D */
  0.333, 0.333, 0.333, 0.0,        /* 11. V */
  0.333, 0.333, 0.0,   0.333,      /* 12. H */
  0.0,   0.333, 0.333, 0.333,      /* 13. B */
  0.25,  0.25,  0.25,  0.25        /* 14. N */
};

const int consensus2index[] = {
  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  -1,  0, 13,  1, 10, -1, -1,  2, 12, -1, -1,  8, -1,  9, 14, -1,
  -1, -1,  5,  7,  3,  3, 11,  6, -1,  4, -1, -1, -1, -1, -1, -1,
  -1,  0, 13,  1, 10, -1, -1,  2, 12, -1, -1,  8, -1,  9, 14, -1,
  -1, -1,  5,  7,  3,  3, 11,  6, -1,  4, -1, -1, -1, -1, -1, -1,
  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1
};

enum MOTIF_FMT {
  FMT_MEME     = 1,
  FMT_HOMER    = 2,
  FMT_JASPAR   = 3,
  FMT_HOCOMOCO = 4,
  FMT_UNKNOWN  = 5
};

typedef struct args_t {
  double   bkg[4];
  double   pvalue;
  int      nsites;
  int      pseudocount; 
  int      scan_rc;
  int      dedup;
  int      trim_names;
  int      use_user_bkg;
  int      low_mem;
  int      nthreads;
  int      thresh0;
  int      progress;
  int      v;
  int      w;
} args_t;

args_t args = {
  .bkg             = {0.25, 0.25, 0.25, 0.25},
  .pvalue          = DEFAULT_PVALUE,
  .nsites          = DEFAULT_NSITES,
  .pseudocount     = DEFAULT_PSEUDOCOUNT,
  .scan_rc         = 1,
  .dedup           = 0,
  .trim_names      = 1,
  .use_user_bkg    = 0,
  .low_mem         = 1,
  .nthreads        = 1,
  .thresh0         = 0,
  .progress        = 0,
  .v               = 0,
  .w               = 0
};

typedef struct motif_t {
  size_t    size;
  size_t    cdf_size;
  size_t    thread;
  size_t    file_line_num;
  int       threshold;
  int       min;                         /* Smallest single PWM score */
  int       max;                         /* Largest single PWM score  */
  int       max_score;                   /* Largest total PWM score   */
  int       min_score;                   /* Smallest total PWM score  */
  int       cdf_max;
  int       cdf_offset;
  char      name[MAX_NAME_SIZE];
  int       pwm[MAX_MOTIF_SIZE];
  int       pwm_rc[MAX_MOTIF_SIZE];
  double   *cdf;
  double   *tmp_pdf;
} motif_t;

motif_t **motifs;

typedef struct motif_info_t {
  int     is_consensus;
  int     fmt;
  size_t  n;
} motif_info_t;

motif_info_t motif_info = {
  .is_consensus = 0,
  .fmt = 0,
  .n = 0
};

size_t   *cdf_real_size;
double  **cdf;
double  **tmp_pdf;

int alloc_cdf(void) {
  cdf_real_size = malloc(sizeof(size_t) * args.nthreads);
  if (cdf_real_size == NULL) {
    fprintf(stderr, "Error: Failed to allocate memory for real CDF sizes.");
    return 1;
  }
  for (size_t i = 0; i < args.nthreads; i++) {
    cdf_real_size[i] = 1;
  }
  cdf = malloc(sizeof(double *) * args.nthreads);
  if (cdf == NULL) {
    fprintf(stderr, "Error: Failed to allocate memory for CDFs.");
    return 1;
  }
  tmp_pdf = malloc(sizeof(double *) * args.nthreads);
  if (tmp_pdf == NULL) {
    fprintf(stderr, "Error: Failed to allocate memory for temporary PDFs.");
    return 1;
  }
  for (size_t i = 0; i < args.nthreads; i++) {
    cdf[i] = malloc(sizeof(double));
    if (cdf[i] == NULL) {
      fprintf(stderr, "Error: Failed to allocate memory for CDF (#%zu).", i);
      return 1;
    }
    tmp_pdf[i] = malloc(sizeof(double));
    if (tmp_pdf[i] == NULL) {
      fprintf(stderr, "Error: Failed to allocate memory for temporary PDF (#%zu).", i);
      return 1;
    }
  }
  return 0;
}

typedef struct seq_info_t {
  size_t     n;
  size_t     total_bases;
  size_t     unknowns;
  double     gc_pct;
} seq_info_t;

seq_info_t seq_info = {
  .n = 0,
  .total_bases = 0,
  .unknowns = 0,
  .gc_pct = 0.0
};

char            **seq_names;
unsigned char   **seqs;
size_t           *seq_sizes;

void free_seqs(void) {
  for (size_t i = 0; i < seq_info.n; i++) {
    free(seq_names[i]);
    if (!args.low_mem) free(seqs[i]);
  }
  free(seq_names);
  free(seq_sizes);
  free(seqs);
}

void free_motifs(void) {
  for (size_t i = 0; i < motif_info.n; i++) {
    free(motifs[i]);
  }
  free(motifs);
}

void free_cdf(void) {
  for (size_t i = 0; i < args.nthreads; i++) {
    free(cdf[i]);
    free(tmp_pdf[i]);
  }
  free(cdf);
  free(tmp_pdf);
  free(cdf_real_size);
}

pthread_t         *threads;
pthread_mutex_t    pb_lock = PTHREAD_MUTEX_INITIALIZER;
size_t             pb_counter = 0;

typedef struct files_t {
  FILE     *m;
  gzFile    s;
  FILE     *o;
  int       m_open;
  int       s_open;
  int       o_open;
} files_t;

files_t files = {
  .m_open = 0,
  .s_open = 0,
  .o_open = 0
};

void close_files(void) {
  if (files.m_open) fclose(files.m);
  if (files.s_open) gzclose(files.s);
  if (files.o_open) fclose(files.o);
}

void init_motif(motif_t *motif) {
  ERASE_ARRAY(motif->name, MAX_NAME_SIZE);
  motif->name[0] = 'm';
  motif->name[1] = 'o';
  motif->name[2] = 't';
  motif->name[3] = 'i';
  motif->name[4] = 'f';
  motif->name[5] = '\0';
  motif->size = 0;
  motif->threshold = 0;
  motif->max_score = 0;
  motif->min = 0;
  motif->max = 0;
  motif->cdf_size = 0;
  motif->file_line_num = 0;
  motif->min_score = 0;
  motif->cdf_max = 0;
  motif->thread = 0;
  for (size_t i = 0; i < MAX_MOTIF_SIZE; i++) {
    motif->pwm[i] = 0;
    motif->pwm_rc[i] = 0;
  }
  for (size_t i = 4; i < MAX_MOTIF_SIZE; i += 5) {
    motif->pwm[i] = AMBIGUITY_SCORE;
    motif->pwm_rc[i] = AMBIGUITY_SCORE;
  }
}

static inline void set_score(motif_t *motif, const unsigned char let, const size_t pos, const int score) {
  motif->pwm[char2index[let] + pos * 5] = score;
}

static inline int get_score(const motif_t *motif, const unsigned char let, const size_t pos) {
  return motif->pwm[char2index[let] + pos * 5];
}

static inline int get_score_i(const motif_t *motif, const int i, const size_t pos) {
  return motif->pwm[i + pos * 5];
}

static inline void set_score_rc(motif_t *motif, const unsigned char let, const size_t pos, const int score) {
  motif->pwm_rc[char2index[let] + pos * 5] = score;
}

static inline int get_score_rc(const motif_t *motif, const unsigned char let, const size_t pos) {
  return motif->pwm_rc[char2index[let] + pos * 5];
}

void badexit(const char *msg) {
  fprintf(stderr, "%s\nRun minimotif -h to see usage.\n", msg);
  free(threads);
  free_motifs();
  free_seqs();
  close_files();
  exit(EXIT_FAILURE);
}

/* For the motif half of minimotif, this function is (by far) where it spends
 * most of its time.
 */
void fill_cdf(motif_t *motif) {
  size_t max_step, s;
  double pdf_sum = 0.0;
  if (args.w && args.nthreads == 1 && !args.progress) {
    fprintf(stderr, "        Generating CDF for [%s] (n=%'zu) ... ",
      motif->name, motif->cdf_size);
  } else if (args.w && !args.progress) {
    fprintf(stderr, "        Generating CDF for [%s] (n=%'zu)\n",
      motif->name, motif->cdf_size);
  }
  if (motif->cdf_size > MAX_CDF_SIZE) {
    if (args.w&& args.nthreads == 1 && !args.progress) fprintf(stderr, "\n");
    fprintf(stderr,
        "Internal error: Requested CDF size for [%s] is too large (%'zu>%'zu).\n",
        motif->name, motif->cdf_size, MAX_CDF_SIZE);
    fprintf(stderr, "    Make sure no background values are below %f.",
        MIN_BKG_VALUE);
    badexit("");
  }
  /* Instead of allocating and freeing a CDF for every motif, share a
   * single one for all motifs -- just reset it every time and realloc to a
   * larger size if needed.
   */
  if (cdf_real_size[motif->thread] < motif->cdf_size) {
    double *cdf_rl = realloc(cdf[motif->thread], motif->cdf_size * sizeof(double));
    if (cdf_rl == NULL) {
      badexit("Error: Memory re-allocation for motif CDF failed.");
    }
    cdf[motif->thread] = cdf_rl;
    double *tmp_pdf_rl = realloc(tmp_pdf[motif->thread], motif->cdf_size * sizeof(double));
    if (cdf_rl == NULL) {
      badexit("Error: Memory re-allocation for temporary motif PDF failed.");
    }
    tmp_pdf[motif->thread] = tmp_pdf_rl;
    cdf_real_size[motif->thread] = motif->cdf_size;
  }
  motif->cdf = cdf[motif->thread];
  motif->tmp_pdf = tmp_pdf[motif->thread];
  for (size_t i = 0; i < motif->cdf_size; i++) motif->cdf[i] = 1.0;
  for (size_t i = 0; i < motif->size; i++) {
    max_step = i * motif->cdf_max;
    for (size_t j = 0; j < motif->cdf_size; j++) {
      motif->tmp_pdf[j] = motif->cdf[j];
    }
    ERASE_ARRAY(motif->cdf, max_step + motif->cdf_max + 1);
    for (int j = 0; j < 4; j++) {
      s = get_score_i(motif, j, i) - motif->min;
      /* This loop is where the majority of time is spent for motif-related code. */
      for (size_t k = 0; k <= max_step; k++) {
        motif->cdf[k+s] += motif->tmp_pdf[k] * args.bkg[j];
      }
    }
  }
  for (size_t i = 0; i < motif->cdf_size; i++) pdf_sum += motif->cdf[i];
  if (fabs(pdf_sum - 1.0) > 0.0001) {
    if (args.w && args.nthreads == 1 && !args.progress) {
      fprintf(stderr, "Internal warning: sum(PDF)!= 1.0 for [%s] (sum=%.2g)\n",
          motif->name, pdf_sum);
    }
    for (size_t i = 0; i < motif->cdf_size; i++) {
      motif->cdf[i] /= pdf_sum;
    }
  }
  for (size_t i = motif->cdf_size - 2; i >= 0 && i < -1; i--) {
    motif->cdf[i] += motif->cdf[i + 1];
  }
  if (args.w && args.nthreads == 1 && !args.progress) fprintf(stderr, "done.\n");
}

static inline double score2pval(const motif_t *motif, const int score) {
  return motif->cdf[score - motif->cdf_offset];
}

void set_threshold(motif_t *motif) {
  size_t threshold_i = motif->cdf_size;
  for (size_t i = 0; i < motif->cdf_size; i++) {
    if (motif->cdf[i] < args.pvalue) {
      threshold_i = i;
      break;
    }
  }
  motif->threshold -= motif->min;
  motif->threshold *= motif->size;
  motif->threshold = threshold_i - motif->threshold;
  for (size_t i = 0; i < motif->size; i++) {
    int max_pos = get_score_i(motif, 0, i);
    int min_pos = max_pos;
    for (int j = 1; j < 4; j++) {
      int tmp_pos = get_score_i(motif, j, i);
      if (tmp_pos > max_pos) max_pos = tmp_pos;
      if (tmp_pos < min_pos) min_pos = tmp_pos;
    }
    motif->max_score += max_pos;
    motif->min_score += min_pos;
  }
  double min_pvalue = score2pval(motif, motif->max_score);
  if (min_pvalue / args.pvalue > 1.0001) {
    if (args.w && !args.progress) {
      fprintf(stderr,
        "Warning: Min possible pvalue for [%s] is greater than the threshold,\n",
        motif->name);
      fprintf(stderr, "  motif will not be scored (%g>%g).\n",
        min_pvalue, args.pvalue);
    }
    motif->threshold = INT_MAX;
  }
  if (args.thresh0) {
    motif->threshold = 0;
  } else if (motif_info.is_consensus) {
    motif->threshold = motif->max_score;
  }
}

int check_and_load_bkg(double *bkg) {
  if (bkg[0] == -1.0 || bkg[1] == -1.0 || bkg[2] == -1.0 || bkg[3] == -1.0) {
    fprintf(stderr, "Error: Too few background values found (need 4)."); return 1;
  }
  double min = 0; VEC_MIN(bkg, min, 4);
  if (min < MIN_BKG_VALUE) {
    if (args.v) {
      fprintf(stderr,
        "Warning: Detected background values smaller than allowed min,\n");
      fprintf(stderr, "    adjusting (%.2g<%.2g).\n", min, MIN_BKG_VALUE);
    }
    VEC_ADD(bkg, MIN_BKG_VALUE, 4);
  }
  double sum = bkg[0] + bkg[1] + bkg[2] + bkg[3];
  if (fabs(sum - 1.0) > 0.001 && args.v) {
    fprintf(stderr,
      "Warning: Background values don't add up to 1.0, adjusting (sum=%.3g).\n",
      sum);
  }
  VEC_DIV(bkg, sum, 4);
  args.bkg[0] = bkg[0]; args.bkg[1] = bkg[1]; args.bkg[2] = bkg[2]; args.bkg[3] = bkg[3];
  return 0;
}

void parse_user_bkg(const char *bkg_usr) {
  size_t i = 0, j = 0, bi = 0;
  char bc[USER_BKG_MAX_SIZE];
  double b[] = {-1.0, -1.0, -1.0, -1.0};
  ERASE_ARRAY(bc, USER_BKG_MAX_SIZE);
  while (bkg_usr[i] != '\0') {
    if (bkg_usr[i] != ',' && bkg_usr[i] != ' ') {
      bc[j] = bkg_usr[i];
      j++;
    } else if (bkg_usr[i] == ',') {
      if (bi > 2) {
        badexit("Error: Too many background values provided (need 4).");
      }
      b[bi] = atof(bc);
      ERASE_ARRAY(bc, USER_BKG_MAX_SIZE);
      bi++;
      j = 0;
    }
    i++;
  }
  b[3] = atof(bc);
  if (check_and_load_bkg(b)) badexit("");
  if (args.w) {
    fprintf(stderr, "Using new background values:\n");
    fprintf(stderr, "    A=%.3g", args.bkg[0]);
    fprintf(stderr, "    C=%.3g\n", args.bkg[1]);
    fprintf(stderr, "    G=%.3g", args.bkg[2]);
    fprintf(stderr, "    T=%.3g\n", args.bkg[3]);
  }
}

static inline double b2mb(const size_t x) {
  return (x / 1024.0) / 1024.0;
}

int check_line_contains(const char *line, const char *substring) {
  size_t ss_len = strlen(substring);
  if (strlen(line) < ss_len) return 0;
  for (size_t i = 0; i < ss_len; i++) {
    if (line[i] != substring[i]) return 0;
  }
  return 1;
}

size_t count_nonempty_chars(const char *line) {
  /* Other whitespace characters: \v, \f. Assuming these would never show up. */
  size_t total_chars = 0, i = 0;
  for (;;) {
    switch (line[i]) {
      case ' ':
      case '\t':
      case '\r':
      case '\n': break;
      case '\0': return total_chars;
      default: total_chars++;
    }
    i++;
  }
  return total_chars;
}

int check_char_is_one_of(const char c, const char *list) {
  for (size_t i = 0; i < strlen(list); i++) {
    if (list[i] == c) return 1;
  }
  return 0;
}

int detect_motif_fmt(void) {
  int jaspar_or_hocomoco = 0, file_fmt = 0, has_tabs = 0;
  char *line = NULL;
  size_t len = 0;
  ssize_t read;
  while ((read = getline(&line, &len, files.m)) != -1) {
    if (!count_nonempty_chars(line)) continue;
    if (check_line_contains(line, "MEME version \0")) {
      if (args.w) {
        fprintf(stderr, "Detected MEME format (version %c).\n", line[13]);
      }
      file_fmt = FMT_MEME;
      break;
    }
    if (jaspar_or_hocomoco) {
      if (line[0] == 'A' &&
          check_char_is_one_of('[', line) &&
          check_char_is_one_of(']', line)) {
        file_fmt = FMT_JASPAR;
        if (args.w) fprintf(stderr, "Detected JASPAR format.\n");
        break;
      } else {
        if (line[0] == 'A'|| 
            check_char_is_one_of('[', line) ||
            check_char_is_one_of(']', line)) {
          badexit("Error: Detected malformed JASPAR format.");
        }
        if (has_tabs) {
          file_fmt = FMT_HOMER;
          if (args.w) fprintf(stderr, "Detected HOMER format.\n");
          break;
        } else {
          if (check_char_is_one_of('-', line)) {
            badexit("Error: minimotif cannot read HOCOMOCO PWMs.");
          }
          file_fmt = FMT_HOCOMOCO;
          if (args.w) fprintf(stderr, "Detected HOCOMOCO format.\n");
          break;
        }
      }
    } else if (line[0] == '>') {
      if (check_char_is_one_of('\t', line)) {
        has_tabs = 1;
      }
      jaspar_or_hocomoco = 1;
    }
  }
  rewind(files.m);
  free(line);
  if (!file_fmt) file_fmt = FMT_UNKNOWN;
  return file_fmt;
}

int add_motif(void) {
  motif_t **tmp_ptr = realloc(motifs, sizeof(*motifs) * (motif_info.n + 1));
  if (tmp_ptr == NULL) {
    fprintf(stderr, "Error: Failed to allocate memory for motifs.");
    return 1;
  } else {
    motifs = tmp_ptr;
  }
  motifs[motif_info.n] = malloc(sizeof(motif_t));
  if (motifs[motif_info.n] == NULL) {
    fprintf(stderr, "Error: Failed to allocate memory for motifs.");
    return 1;
  }
  motif_info.n++;
  init_motif(motifs[motif_info.n - 1]);
  return 0;
}

int calc_score(const double prob_i, const double bkg_i) {
  double x;
  x = prob_i * args.nsites;
  x += ((double) args.pseudocount) / 4.0;
  x /= (double) (args.nsites + args.pseudocount);
  return (int) (log2(x / bkg_i) * PWM_INT_MULTIPLIER);
}

int normalize_probs(double *probs, const char *name) {
  double sum = probs[0] + probs[1] + probs[2] + probs[3];
  if (fabs(sum - 1.0) > 0.1) {
    if (args.w) fprintf(stderr, "\n");
    fprintf(stderr,
      "Error: Position for [%s] does not add up to 1 (sum=%.3g)",
      name, sum);
    return 1;
  }
  if (fabs(sum - 1.0) > 0.02) {
    if (args.w) {
      fprintf(stderr,
        "\nWarning: Position for [%s] does not add up to 1, adjusting (sum=%.3g) ",
        name, sum);
    }
    VEC_DIV(probs, sum, 4);
  }
  return 0;
}

int get_line_probs(const motif_t *motif, const char *line, double *probs, const size_t n) {
  size_t i = 0, j = 0, which_i = -1;
  int prev_line_was_space = 1;
  char pos_i[MOTIF_VALUE_MAX_CHAR];
  ERASE_ARRAY(pos_i, MOTIF_VALUE_MAX_CHAR);
  for (;;) {
    if (line[i] != ' ' && line[i] != '\t') break;
    i++;
  }
  while (line[i] != '\0' && line[i] != '\r' && line[i] != '\n') {
    if (line[i] != ' ' && line[i] != '\t') {
      pos_i[j] = line[i];
      j++;
      prev_line_was_space = 0;
    } else {
      if (!prev_line_was_space) {
        which_i++; 
        if (which_i > n - 1) {
          if (args.w) fprintf(stderr, "\n");
          fprintf(stderr,
            "Error: Motif [%s] has too many columns (need %zu).",
            motif->name, n); return 1;
        }
        probs[which_i] = atof(pos_i);
        ERASE_ARRAY(pos_i, MOTIF_VALUE_MAX_CHAR);
        j = 0;
      }
      prev_line_was_space = 1;
    }
    i++;
  }

  if (!prev_line_was_space) {
    which_i++; 
    if (which_i > n - 1) {
      if (args.w) fprintf(stderr, "\n");
      fprintf(stderr,
        "Error: Motif [%s] has too many columns (need %zu).",
        motif->name, n); return 1;
    }
    probs[which_i] = atof(pos_i);
  }

  if (which_i == -1) {
    if (args.w) fprintf(stderr, "\n");
    fprintf(stderr, "Error: Motif [%s] has an empty row.",
      motif->name); return 1;
  }

  if (which_i < n - 1) {
    if (args.w) fprintf(stderr, "\n");
    fprintf(stderr, "Error: Motif [%s] has too few columns (need %zu).",
      motif->name, n); return 1;
  }

  return 0;
}

int add_motif_ppm_column(motif_t *motif, const char *line, const size_t pos) {
  double probs[] = {-1.0, -1.0, -1.0, -1.0};
  if (get_line_probs(motif, line, probs, 4)) return 1;
  if (normalize_probs(probs, motif->name)) return 1;
  set_score(motif, 'A', pos, calc_score(probs[0], args.bkg[0]));
  set_score(motif, 'C', pos, calc_score(probs[1], args.bkg[1]));
  set_score(motif, 'G', pos, calc_score(probs[2], args.bkg[2]));
  set_score(motif, 'T', pos, calc_score(probs[3], args.bkg[3]));
  return 0;
}

int check_meme_alph(const char *line, const size_t line_num) {
  if (check_line_contains(line, "ALPHABET= ACDEFGHIKLMNPQRSTVWY\0")) {
    fprintf(stderr, "Error: Detected protein alphabet (L%zu).", line_num);
    return 1;
  }
  return 0;
}

int check_meme_strand(const char *line, const size_t line_num) {
  size_t scan_fwd = 0, scan_rev = 0, i = 0;
  for (;;) {
    if (line[i] == '\0') break;
    if (line[i] == '+') scan_fwd++;
    if (line[i] == '-') scan_rev++;
    i++;
  }
  if (((scan_fwd > 1 || scan_rev > 1) || (!scan_fwd && !scan_rev)) && args.v) {
    fprintf(stderr, "Warning: Possible malformed strand field (L%zu).\n", line_num);
  }
  if (args.scan_rc && scan_fwd && !scan_rev && args.v) {
    fprintf(stderr, "Warning: MEME motifs are only for the forward strand (L%zu).\n",
      line_num);
  }
  if (!scan_fwd && scan_rev && args.v) {
    fprintf(stderr, "Warning: MEME motifs are only for the reverse strand (L%zu).\n",
      line_num);
  }
  if (!args.scan_rc && scan_fwd && scan_rev && args.v) {
    fprintf(stderr, "Warning: MEME motifs are for both strands (L%zu).\n",
      line_num);
  }
  return 0;
}

int get_meme_bkg(const char *line, const size_t line_num) {
  if (args.use_user_bkg) return 0;
  double bkg_probs[] = {-1.0, -1.0, -1.0, -1.0};
  size_t i = 1, let_i = 0, j = 0, empty = 0;
  char bkg_char[MEME_BKG_MAX_SIZE];
  ERASE_ARRAY(bkg_char, MEME_BKG_MAX_SIZE);
  if (line[0] != 'A') {
    fprintf(stderr, "Error: Expected first character of background line to be 'A' (L%zu).",
      line_num); return 1;
  }
  while (line[i] != '\0' && line[i] != '\n' && line[i] != '\r') {
    if (let_i > 3) {
      fprintf(stderr, "Error: Parsed too many background values in MEME file (L%zu).",
        line_num); return 1;
    }
    if (line[i] != ' ' && line[i] != '\t') {
      if (line[i] == 'C') {
        if (!empty) {
          fprintf(stderr, "Error: Expected whitespace before 'C' character (L%zu).",
            line_num); return 1;
        }
        if (let_i != 0) {
          fprintf(stderr,
            "Error: Expected 'C' to be second letter in MEME background (L%zu).",
            line_num); return 1;
        }
        bkg_probs[let_i] = atof(bkg_char);
        ERASE_ARRAY(bkg_char, MEME_BKG_MAX_SIZE);
        let_i = 1; j = 0;
      } else if (line[i] == 'G') {
        if (!empty) {
          fprintf(stderr, "Error: Expected whitespace before 'C' character (L%zu).",
            line_num); return 1;
        }
        if (let_i != 1) {
          fprintf(stderr,
            "Error: Expected 'G' to be third letter in MEME background (L%zu).",
            line_num); return 1;
        }
        bkg_probs[let_i] = atof(bkg_char);
        ERASE_ARRAY(bkg_char, MEME_BKG_MAX_SIZE);
        let_i = 2; j = 0;
      } else if (line[i] == 'T' || line[i] == 'U') {
        if (!empty) {
          fprintf(stderr, "Error: Expected whitespace before 'C' character (L%zu).",
            line_num); return 1;
        }
        if (let_i != 2) {
          fprintf(stderr,
            "Error: Expected 'T/U' to be fourth letter in MEME background (L%zu).",
            line_num); return 1;
        }
        bkg_probs[let_i] = atof(bkg_char);
        ERASE_ARRAY(bkg_char, MEME_BKG_MAX_SIZE);
        let_i = 3; j = 0;
      } else if (check_char_is_one_of(line[i], "0123456789.\0")) {
        bkg_char[j] = line[i]; 
        j++;
      } else {
        fprintf(stderr,
          "Error: Encountered unexpected character (%c) in MEME background (L%zu).",
          line[i], line_num); return 1;
      }
      empty = 0;
    } else {
      empty = 1;
    }
    i++;
  }
  if (bkg_char[0] != '\0') bkg_probs[let_i] = atof(bkg_char);
  if (check_and_load_bkg(bkg_probs)) return 1;
  if (args.w) {
    fprintf(stderr, "Found MEME background values:\n");
    fprintf(stderr, "    A=%.3g", args.bkg[0]);
    fprintf(stderr, "    C=%.3g\n", args.bkg[1]);
    fprintf(stderr, "    G=%.3g", args.bkg[2]);
    fprintf(stderr, "    T=%.3g\n", args.bkg[3]);
  }
  return 0;
}

void parse_meme_name(const char *line, const size_t motif_i) {
  size_t i = 5, j = 0, name_read = 0;
  while (line[i] != '\0' && line[i] != '\r' && line[i] != '\n') {
    if (line[i] == ' ' && name_read) break;
    else if (line[i] == ' ') {
      i++;
      continue;
    }
    name_read = 1;
    motifs[motif_i]->name[j] = line[i];
    j++; i++;
  }
  motifs[motif_i]->name[j] = '\0';
  if (args.w) fprintf(stderr, "    Found motif: %s (size=", motifs[motif_i]->name);
}

void read_meme(void) {
  motif_info.fmt = FMT_MEME;
  char *line = NULL;
  size_t len = 0;
  ssize_t read;
  size_t line_num = 0, l_p_m_L = 0, bkg_let_freqs_L = 0, motif_i = -1, pos_i = -1;
  int alph_detected = 0, strand_detected = 0, live_motif = 0;
  while ((read = getline(&line, &len, files.m)) != -1) {
    line_num++;
    if (check_line_contains(line, "Background letter frequencies\0")) {
      if (bkg_let_freqs_L) {
        free(line);
        fprintf(stderr,
          "Error: Detected multiple background definition lines in MEME file (L%zu).",
          line_num);
      } else {
        if (motif_i < -1) {
          free(line);
          fprintf(stderr, "Error: Found background definition line after motifs (L%zu).",
            line_num);
          badexit("");
        }
        bkg_let_freqs_L = line_num;
      }
    } else if (bkg_let_freqs_L && bkg_let_freqs_L == line_num - 1) {
      if (get_meme_bkg(line, line_num)) {
        free(line);
        badexit("");
      }
    } else if (check_line_contains(line, "ALPHABET\0")) {
      if (alph_detected) {
        free(line);
        fprintf(stderr,
          "Error: Detected multiple alphabet definition lines in MEME file (L%zu).",
          line_num);
        badexit("");
      }
      if (motif_i < -1) {
        free(line);
        fprintf(stderr, "Error: Found alphabet definition line after motifs (L%zu).",
          line_num);
        badexit("");
      }
      if (check_meme_alph(line, line_num)) {
        free(line);
        badexit("");
      }
      alph_detected = 1;
    } else if (check_line_contains(line, "strands:\0")) {
      if (strand_detected) {
        free(line);
        fprintf(stderr,
          "Error: Detected multiple strand information lines in MEME file (L%zu).",
          line_num);
        badexit("");
      }
      if (motif_i < -1) {
        free(line);
        fprintf(stderr, "Error: Found strand information line after motifs (L%zu).",
          line_num);
        badexit("");
      }
      if (check_meme_strand(line, line_num)) {
        free(line);
        badexit("");
      }
      strand_detected = 1;
    } else if (check_line_contains(line, "MOTIF\0")) {
      if (motif_i < -1 && args.w) {
        fprintf(stderr, "%zu)\n", motifs[motif_i]->size);
      }
      motif_i++;
      if (add_motif()) {
        free(line);
        badexit("");
      }
      motifs[motif_i]->file_line_num = line_num;
      parse_meme_name(line, motif_i);
      pos_i = 0;
    } else if (check_line_contains(line, "letter-probability matrix\0")) {
      if (pos_i != 0) {
        free(line);
        fprintf(stderr, "Error: Possible malformed MEME motif (L%zu).",
          line_num);
        badexit("");
      }
      l_p_m_L = line_num;
      live_motif = 1;
    } else if (live_motif) {

      if (!count_nonempty_chars(line) || check_char_is_one_of('-', line) ||
          check_char_is_one_of('*', line)) {
        live_motif = 0;
      } else if (line_num == (l_p_m_L + pos_i + 1)) {

        if (pos_i >= MAX_MOTIF_SIZE / 5 && pos_i < -1) {
          free(line);
          fprintf(stderr, "Error: Motif [%s] is too large (max=%zu)",
            motifs[motif_i]->name, MAX_MOTIF_SIZE / 5);
          badexit("");
        }
        if (add_motif_ppm_column(motifs[motif_i], line, pos_i)) {
          free(line);
          badexit("");
        }
        pos_i++;
        motifs[motif_i]->size = pos_i;

      } else {
        live_motif = 0;
      }

    }
  }
  free(line);
  if (motif_i < -1 && args.w) {
    fprintf(stderr, "%zu)\n", motifs[motif_i]->size);
  }
  if (!motif_info.n) badexit("Error: Failed to detect any motifs in MEME file.");
  if (args.v) {
    fprintf(stderr, "Found %'zu MEME motif(s).\n", motif_info.n);
  }
}

void parse_homer_name(const char *line, const size_t motif_i) {
  size_t name_start = 0, name_end = 0, i = 1, in_between = 0, j = 0;
  while (line[i] != '\0' && line[i] != '\r' && line[i] != '\n') {
    if (line[i] == '\t') {
      if (name_start) {
        name_end = i;
        break;
      } else {
        in_between = 1;
      }
    } else if (in_between) {
      if (!name_start) name_start = i;
    } 
    i++;
  }
  if (!name_start) {
    if (args.w) {
      fprintf(stderr, "Warning: Failed to parse motif name [#%'zu].\n", motif_i + 1);
    }
  } else if (!name_end) {
    if (args.w) {
      fprintf(stderr, "Warning: HOMER motif is missing logodds score [#%'zu].\n", 
        motif_i + 1);
    }
    name_end = i;
  }
  for (size_t k = name_start; k < name_end; k++) {
    motifs[motif_i]->name[j] = line[k];
    j++;
  }
  motifs[motif_i]->name[j] = '\0';
  if (args.w) fprintf(stderr, "    Found motif: %s (size=", motifs[motif_i]->name);
}

void read_homer(void) {
  motif_info.fmt = FMT_HOMER;
  char *line = NULL;
  size_t len = 0;
  ssize_t read;
  size_t line_num = 0, motif_i = -1, pos_i;
  int ready_to_start = 0;
  while ((read = getline(&line, &len, files.m)) != -1) {
    line_num++;
    if (line[0] == '>') {
      ready_to_start = 1;
      if (motif_i < -1 && args.w) {
        fprintf(stderr, "%zu)\n", motifs[motif_i]->size);
      }
      motif_i++;
      if (add_motif()) {
        free(line);
        badexit("");
      }
      motifs[motif_i]->file_line_num = line_num;
      parse_homer_name(line, motif_i);
      pos_i = 0;
    } else if (count_nonempty_chars(line) && ready_to_start) {
      if (pos_i > MAX_MOTIF_SIZE / 5 && pos_i < -1) {
        fprintf(stderr, "Error: Motif [%s] is too large (max=%'zu).\n",
          motifs[motif_i]->name, MAX_MOTIF_SIZE / 5);
      }
      if (add_motif_ppm_column(motifs[motif_i], line, pos_i)) {
        free(line);
        badexit("");
      }
      pos_i++;
      motifs[motif_i]->size = pos_i;
    }
  }
  free(line);
  if (motif_i < -1 && args.w) {
    fprintf(stderr, "%zu)\n", motifs[motif_i]->size);
  }
  if (args.v) {
    fprintf(stderr, "Found %'zu HOMER motif(s).\n", motif_info.n);
  }
}

int get_pwm_max(const motif_t *motif) {
  int max = 0, val;
  for (size_t pos = 0; pos < motif->size; pos++) {
    for (int let = 0; let < 4; let++) {
      val = get_score_i(motif, let, pos);
      if (val > max) max = val;
    }
  }
  return max;
}

int get_pwm_min(const motif_t *motif) {
  int min = 0, val;
  for (size_t pos = 0; pos < motif->size; pos++) {
    for (int let = 0; let < 4; let++) {
      val = get_score_i(motif, let, pos);
      if (val < min) min = val;
    }
  }
  return min;
}

void fill_pwm_rc(motif_t *motif) {
  for (size_t pos = 0; pos < motif->size; pos++) {
    set_score_rc(motif, 'A', motif->size - 1 - pos, get_score(motif, 'T', pos));
    set_score_rc(motif, 'C', motif->size - 1 - pos, get_score(motif, 'G', pos));
    set_score_rc(motif, 'G', motif->size - 1 - pos, get_score(motif, 'C', pos));
    set_score_rc(motif, 'T', motif->size - 1 - pos, get_score(motif, 'A', pos));
  }
}

void trim_motif_name(motif_t *motif) {
  for (size_t i = 0; i < MAX_NAME_SIZE; i++) {
    if (motif->name[i] == ' ' || motif->name[i] == '\t' || motif->name[i] == '\0') {
      motif->name[i] = '\0';
      break;
    }
  }
}

void complete_motifs(void) {
  for (size_t i = 0; i < motif_info.n; i++) {
    motifs[i]->min = get_pwm_min(motifs[i]);
    motifs[i]->max = get_pwm_max(motifs[i]);
    motifs[i]->cdf_offset = motifs[i]->min * motifs[i]->size;
    fill_pwm_rc(motifs[i]);
    motifs[i]->cdf_max = motifs[i]->max - motifs[i]->min;
    motifs[i]->cdf_size = motifs[i]->size * motifs[i]->cdf_max + 1;
    if (args.trim_names) trim_motif_name(motifs[i]);
  }
}

void print_motif(motif_t *motif, const size_t n) {
  fprintf(files.o, "Motif: %s (N%zu L%zu)\n", motif->name, n, motif->file_line_num);
  if (motif->threshold == INT_MAX) {
    fprintf(files.o, "MaxScore=%.2f\tThreshold=%s\n",
      motif->max_score / PWM_INT_MULTIPLIER, "[exceeds max]");
  } else {
    fprintf(files.o, "MaxScore=%.2f\tThreshold=%.2f\n",
      motif->max_score / PWM_INT_MULTIPLIER, motif->threshold / PWM_INT_MULTIPLIER);
  }
  fprintf(files.o, "Motif PWM:\n\tA\tC\tG\tT\n");
  for (size_t i = 0; i < motif->size; i++) {
    fprintf(files.o, "%zu:\t%.2f\t%.2f\t%.2f\t%.2f\n", i + 1,
      get_score(motif, 'A', i) / PWM_INT_MULTIPLIER,
      get_score(motif, 'C', i) / PWM_INT_MULTIPLIER,
      get_score(motif, 'G', i) / PWM_INT_MULTIPLIER,
      get_score(motif, 'T', i) / PWM_INT_MULTIPLIER);
  }
  fprintf(files.o, "Score=%.2f\t-->     p=1\n",
      motif->min_score / PWM_INT_MULTIPLIER);
  fprintf(files.o, "Score=%.2f\t-->     p=%.2g\n",
      (motif->min_score / 2) / PWM_INT_MULTIPLIER,
      score2pval(motif, motif->min_score / 2));
  fprintf(files.o, "Score=0.00\t-->     p=%.2g\n",
      score2pval(motif, 0.0));
  fprintf(files.o, "Score=%.2f\t-->     p=%.2g\n",
      (motif->max_score / 2) / PWM_INT_MULTIPLIER,
      score2pval(motif, motif->max_score / 2));
  fprintf(files.o, "Score=%.2f\t-->     p=%.2g\n",
      motif->max_score / PWM_INT_MULTIPLIER,
      score2pval(motif, motif->max_score));
}

void parse_jaspar_name(const char *line, const size_t motif_i) {
  size_t i = 0, j = 1;
  for (;;) {
    if (line[j] == '\r' || line[j] == '\n' || line[j] == '\0') break;
    motifs[motif_i]->name[i] = line[j];
    i++; j++;
  }
  motifs[motif_i]->name[i] = '\0';
  if (args.w) fprintf(stderr, "    Found motif: %s (size=", motifs[motif_i]->name);
}

int add_jaspar_row(motif_t *motif, const char *line) {
  size_t row_i = -1, left_bracket = -1, right_bracket = -1, i = 0;
  char let;
  for (;;) {
    if (line[i] == '\r' || line[i] == '\n' || line[i] == '\0') break;
    switch (line[i]) {
      case '\t':
      case ' ':
        break;
      case 'a':
      case 'A':
        row_i = 0;
        let = 'A';
        break;
      case 'c':
      case 'C':
        row_i = 1;
        let = 'C';
        break;
      case 'g':
      case 'G':
        row_i = 2;
        let = 'G';
        break;
      case 'u':
      case 'U':
      case 't':
      case 'T':
        row_i = 3;
        let = 'T';
        break;
      case '[':
        left_bracket = i;
        break;
      case ']':
        right_bracket = i;
        break;
    }
    i++;
  }
  if (row_i == -1) {
    fprintf(stderr, "Error: Couldn't find ACGTU in motif [%s] row names.", motif->name);
    return 1;
  }
  if (left_bracket == -1 || right_bracket == -1) {
    fprintf(stderr, "Error: Couldn't find '[]' in motif [%s] row (%zu).",
        motif->name, row_i + 1);
    return 1;
  }
  size_t k = 0, pos_i = -1;
  int prev_line_was_space = 1;
  char prob_c[MOTIF_VALUE_MAX_CHAR];
  ERASE_ARRAY(prob_c, MOTIF_VALUE_MAX_CHAR);
  i = left_bracket + 1;
  for (;;) {
    if (line[i] != ' ' && line[i] != '\t') break;
    i++;
  }
  for (size_t j = i; j < right_bracket; j++) {
    if (line[j] != ' ' && line[j] != '\t') {
      prob_c[k] = line[j];
      k++; prev_line_was_space = 0;
    } else {
      if (!prev_line_was_space) {
        pos_i++;
        if (pos_i + 1 > MAX_MOTIF_SIZE && pos_i < -1) {
          fprintf(stderr, "Error: Motif [%s] has too many columns (need %zu).",
            motif->name, MAX_MOTIF_SIZE); return 1;
        }
        set_score(motif, let, pos_i, atoi(prob_c));
        ERASE_ARRAY(prob_c, MOTIF_VALUE_MAX_CHAR);
        k = 0;
      }
      prev_line_was_space = 1;
    }
  }
  if (!prev_line_was_space) {
    pos_i++;
    if (pos_i > MAX_NAME_SIZE && pos_i < -1) {
      fprintf(stderr, "Error: Motif [%s] has too many columns (need %zu).",
        motif->name, MAX_MOTIF_SIZE); return 1;
    }
    set_score(motif, let, pos_i, atoi(prob_c));
  }
  if (pos_i == -1) {
    fprintf(stderr, "Error: Motif [%s] has an empty row.", motif->name); return 1;
  }
  pos_i++;
  if (motif->size) {
    if (motif->size != pos_i) {
      fprintf(stderr, "Error: Motif [%s] has rows with differing numbers of counts.",
        motif->name); return 1;
    }
  } else {
    motif->size = pos_i;
  }
  return 0;
}

void pcm_to_pwm(motif_t *motif) {
  int nsites = 0, nsites2;
  for (int i = 0; i < 4; i++) {
    nsites += get_score_i(motif, i, 0);
  }
  for (size_t j = 0; j < motif->size; j++) {
    nsites2 = 0;
    for (int i = 0; i < 4; i++) {
      nsites2 += get_score_i(motif, i, j);
    }
    if (abs(nsites2 - nsites) > 1) {
      fprintf(stderr, "Error: Column sums for motif [%s] are not equal.", motif->name);
      badexit("");
    } else if (abs(nsites2 - nsites) == 1 && args.w) {
      fprintf(stderr, "Warning: Found difference of 1 between column sums for motif [%s].",
        motif->name);
    }
  }
  char lets[] = { 'A', 'C', 'G', 'T' };
  for (size_t j = 0; j < motif->size; j++) {
    for (int i = 0; i < 4; i++) {
      set_score(motif, lets[i], j,
          calc_score(
            (args.pseudocount / 4.0 + ((double) get_score_i(motif, i, j))) /
            (args.pseudocount + ((double) nsites)),
            args.bkg[i]));
    }
  }
}

void read_jaspar(void) {
  motif_info.fmt = FMT_JASPAR;
  char *line = NULL;
  size_t len = 0;
  ssize_t read;
  size_t line_num = 0, motif_i = -1, row_i = -1, ready_to_start = 0;
  while ((read = getline(&line, &len, files.m)) != -1) {
    line_num++;
    if (line[0] == '>') {
      ready_to_start = 1;
      if (motif_i < -1 && args.w) {
        fprintf(stderr, "%zu)\n", motifs[motif_i]->size);
      }
      if (motif_i < -1 && row_i != 4) {
        if (row_i < 4) {
          if (args.w) fprintf(stderr, "\n");
          fprintf(stderr, "Error: Motif [%s] has too few rows", motifs[motif_i]->name);
        } else {
          if (args.w) fprintf(stderr, "\n");
          fprintf(stderr, "Error: Motif [%s] has too many rows", motifs[motif_i]->name);
        }
        badexit("");
      }
      motif_i++;
      if (add_motif()) {
        free(line);
        badexit("");
      }
      motifs[motif_i]->file_line_num = line_num;
      parse_jaspar_name(line, motif_i);
      row_i = 0;
    } else if (count_nonempty_chars(line) && ready_to_start) {
      row_i++;
      if (add_jaspar_row(motifs[motif_i], line)) {
        free(line);
        badexit("");
      }
    }
  }
  free(line);
  if (motif_i < -1 && row_i != 4) {
    if (row_i < 4) {
      if (args.w) fprintf(stderr, "\n");
      fprintf(stderr, "Error: Motif [%s] has too few rows", motifs[motif_i]->name);
    } else {
      if (args.w) fprintf(stderr, "\n");
      fprintf(stderr, "Error: Motif [%s] has too many rows", motifs[motif_i]->name);
    }
    badexit("");
  }
  if (motif_i < -1 && args.w) fprintf(stderr, "%zu)\n", motifs[motif_i]->size);
  for (size_t i = 0; i < motif_info.n; i++) {
    pcm_to_pwm(motifs[i]);
  }
  if (args.v) {
    fprintf(stderr, "Found %'zu JASPAR motif(s).\n", motif_info.n);
  }
}

int add_motif_pcm_column(motif_t *motif, const char *line, const size_t pos) {
  double probs[] = {-1.0, -1.0, -1.0, -1.0};
  if (get_line_probs(motif, line, probs, 4)) return 1;
  double pcm_sum = probs[0] + probs[1] + probs[2] + probs[3];
  if (pcm_sum < 0.99) {
    fprintf(stderr, "Error: Motif [%s] PCM row adds up to less than 1", motif->name);
    return 1;
  }
  VEC_ADD(probs, args.pseudocount / 4.0, 4);
  set_score(motif, 'A', pos, calc_score(probs[0] / pcm_sum, args.bkg[0]));
  set_score(motif, 'C', pos, calc_score(probs[1] / pcm_sum, args.bkg[1]));
  set_score(motif, 'G', pos, calc_score(probs[2] / pcm_sum, args.bkg[2]));
  set_score(motif, 'T', pos, calc_score(probs[3] / pcm_sum, args.bkg[3]));
  return 0;
}

void read_hocomoco(void) {
  motif_info.fmt = FMT_HOCOMOCO;
  char *line = NULL;
  size_t len = 0;
  ssize_t read;
  size_t line_num = 0, motif_i = -1, pos_i;
  int ready_to_start = 0;
  while ((read = getline(&line, &len, files.m)) != -1) {
    line_num++;
    if (line[0] == '>') {
      ready_to_start = 1;
      if (motif_i < -1 && args.w) {
        fprintf(stderr, "%zu)\n", motifs[motif_i]->size);
      }
      motif_i++;
      if (add_motif()) {
        free(line);
        badexit("");
      }
      motifs[motif_i]->file_line_num = line_num;
      for (size_t i = 1, j = 0; i < MAX_NAME_SIZE; i++) {
        if (line[i] == '\r' || line[i] == '\n' || line[i] == '\0') {
          motifs[motif_i]->name[j] = '\0';
          break;
        }
        motifs[motif_i]->name[j] = line[i];
        j++;
      }
      if (args.w) fprintf(stderr, "    Found motif: %s (size=", motifs[motif_i]->name);
      pos_i = 0;
    } else if (count_nonempty_chars(line) && ready_to_start) {
      if (pos_i > MAX_MOTIF_SIZE / 5 && pos_i < -1) {
        fprintf(stderr, "Error: Motif [%s] is too large (max=%'zu).\n",
          motifs[motif_i]->name, MAX_MOTIF_SIZE / 5);
      }
      if (add_motif_pcm_column(motifs[motif_i], line, pos_i)) {
        free(line);
        badexit("");
      }
      pos_i++;
      motifs[motif_i]->size = pos_i;
    }
  }
  free(line);
  if (motif_i < -1 && args.w) {
    fprintf(stderr, "%zu)\n", motifs[motif_i]->size);
  }
  if (args.v) {
    fprintf(stderr, "Found %'zu HOCOMOCO motif(s).\n", motif_info.n);
  }
}

void load_motifs(void) {
  switch (detect_motif_fmt()) {
    case FMT_MEME:
      read_meme();
      break;
    case FMT_HOMER:
      read_homer();
      break;
    case FMT_JASPAR:
      read_jaspar();
      break;
    case FMT_HOCOMOCO:
      read_hocomoco();
      break;
    case FMT_UNKNOWN:
      badexit("Error: Failed to detect motif format.");
      break;
  }
  complete_motifs();
  size_t empty_motifs = 0;
  for (size_t i = 0; i < motif_info.n; i++) if (!motifs[i]->size) empty_motifs++;
  if (empty_motifs == motif_info.n) {
    badexit("Error: All parsed motifs are empty.");
  } else if (empty_motifs) {
    fprintf(stderr, "Warning: Found %'zu empty motifs.\n", empty_motifs);
  }
}

void count_bases(void) {
  for (size_t i = 0; i < seq_info.n; i++) {
    for (size_t j = 0; j < seq_sizes[i]; j++) {
      char_counts[seqs[i][j]]++;
    }
  }
}

void count_bases_single(const unsigned char *seq, const size_t len) {
  for (size_t i = 0; i < len; i++) char_counts[seq[i]]++;
}

static inline size_t standard_base_count(void) {
  return
    char_counts['A'] + char_counts['a'] +
    char_counts['C'] + char_counts['c'] +
    char_counts['G'] + char_counts['g'] +
    char_counts['U'] + char_counts['u'] +
    char_counts['T'] + char_counts['t'];
}

double calc_gc(void) {
  double gc = (double) (
    char_counts['G'] + char_counts['C'] + char_counts['g'] + char_counts['c']);
  gc /= standard_base_count();
  return gc;
}

void add_seq_name(char *name, kseq_t *kseq) {
  for (size_t i = 0; i < kseq->name.l; i++) {
    name[i] = kseq->name.s[i];
  }
  if (args.trim_names || !kseq->comment.l) {
    if (kseq->name.l > SEQ_NAME_MAX_CHAR) {
      kseq_destroy(kseq);
      fprintf(stderr, "Error: Sequence name is too large (%zu>%zu).",
        kseq->name.l, SEQ_NAME_MAX_CHAR);
      badexit("");
    }
    name[kseq->name.l] = '\0';
  } else if (kseq->comment.l) {
    if (kseq->name.l + kseq->comment.l + 1 > SEQ_NAME_MAX_CHAR) {
      kseq_destroy(kseq);
      fprintf(stderr, "Error: Sequence name is too large (%zu>%zu).",
        kseq->name.l + kseq->comment.l + 1, SEQ_NAME_MAX_CHAR);
      badexit("");
    }
    name[kseq->name.l] = ' ';
    for (size_t j = 0, i = kseq->name.l + 1; i < kseq->name.l + kseq->comment.l + 1; i++, j++) {
      name[i] = kseq->comment.s[j];
    }
    name[kseq->name.l + kseq->comment.l + 1] = '\0';
  }
}

size_t peak_through_seqs(kseq_t *kseq) {
  size_t name_sizes = 0;
  int ret_val;
  while ((ret_val = kseq_read(kseq)) >= 0) {
    seq_info.n++;
    char **tmp_ptr1 = realloc(seq_names, sizeof(*seq_names) * seq_info.n);
    if (tmp_ptr1 == NULL) {
      kseq_destroy(kseq);
      badexit("Error: Failed to allocate memory for sequence names.");
    } else {
      seq_names = tmp_ptr1;
    }
    size_t *tmp_ptr3 = realloc(seq_sizes, sizeof(*seq_sizes) * seq_info.n);
    if (tmp_ptr3 == NULL) {
      kseq_destroy(kseq);
      badexit("Error: Failed to allocate memory for sequence sizes.");
    } else {
      seq_sizes = tmp_ptr3;
    }
    seq_sizes[seq_info.n - 1] = kseq->seq.l;
    /* TODO: Don't allocate name.l+comment.l if trim_names */
    seq_names[seq_info.n - 1] = malloc(sizeof(char) * kseq->name.l + sizeof(char) * kseq->comment.l + 2);
    name_sizes += kseq->name.l + kseq->comment.l + 2;
    if (seq_names[seq_info.n - 1] == NULL) {
      kseq_destroy(kseq);
      badexit("Error: Failed to allocate memory for sequence name.");
    }
    add_seq_name(seq_names[seq_info.n - 1], kseq);
    const unsigned char *seq_tmp = (unsigned char *) kseq->seq.s;
    for (size_t i = 0; i < kseq->seq.l; i++) {
      char_counts[seq_tmp[i]]++;
    }
  }
  if (ret_val == -2) {
    kseq_destroy(kseq);
    badexit("Error: Failed to parse FASTQ qualities.");
  } else if (ret_val < -2) {
    kseq_destroy(kseq);
    badexit("Error: Failed to read input.");
  } else if (!seq_info.n) {
    kseq_destroy(kseq);
    badexit("Error: Failed to read any sequences from input.");
  }
  gzrewind(files.s);
  kseq_rewind(kseq);
  size_t seq_len_total = 0;
  for (size_t i = 0; i < seq_info.n; i++) seq_len_total += seq_sizes[i];
  if (!seq_len_total) {
    badexit("Error: Only encountered empty sequences.");
  }
  seq_info.total_bases = seq_len_total;
  seq_info.unknowns = seq_len_total - standard_base_count();
  seq_info.gc_pct = calc_gc() * 100.0;
  double unknowns_pct = 100.0 * seq_info.unknowns / seq_len_total;
  if (seq_info.unknowns == seq_len_total) {
    badexit("Error: Failed to read any standard DNA/RNA bases.");
  } else if (unknowns_pct >= 90.0) {
    fprintf(stderr, "!!! Warning: Non-standard base count is extremely high!!! (%.2f%%)\n",
      unknowns_pct);
  } else if (unknowns_pct >= 50.0 && args.v) {
    fprintf(stderr, "Warning: Non-standard base count is very high! (%.2f%%)\n",
      unknowns_pct);
  } else if (unknowns_pct >= 10.0 && args.v) {
    fprintf(stderr, "Warning: Non-standard base count seems high. (%.2f%%)\n",
      unknowns_pct);
  }
  if (char_counts[32] && args.v) {
    fprintf(stderr,
      "Warning: Found spaces (%'zu) in sequences, these will be treated as gaps.\n",
      char_counts[32]);
  }
  size_t max_seq_size = 0;
  for (size_t i = 0; i < seq_info.n; i++) {
    max_seq_size = MAX(max_seq_size, seq_sizes[i]);
  }
  if (args.v) {
    fprintf(stderr, "Found %'zu sequence(s).\n    size=%'zu    GC=%.2f%%\n",
      seq_info.n, seq_len_total, seq_info.gc_pct);
    if (seq_info.unknowns) {
      fprintf(stderr, "Found %'zu (%.2f%%) non-standard bases.\n",
        seq_info.unknowns, unknowns_pct);
    }
    fprintf(stderr, "Approx. memory usage by sequence(s): %'.2f MB\n",
      b2mb(sizeof(unsigned char) * max_seq_size + sizeof(size_t) * seq_info.n +
        sizeof(char) * SEQ_NAME_MAX_CHAR * seq_info.n));
  }
  return max_seq_size;
}

void load_seqs(kseq_t *kseq) {
  size_t name_sizes = 0;
  int ret_val;
  while ((ret_val = kseq_read(kseq)) >= 0) {
    seq_info.n++;
    char **tmp_ptr1 = realloc(seq_names, sizeof(*seq_names) * seq_info.n);
    if (tmp_ptr1 == NULL) {
      kseq_destroy(kseq);
      badexit("Error: Failed to allocate memory for sequence names.");
    } else {
      seq_names = tmp_ptr1;
    }
    unsigned char **tmp_ptr2 = realloc(seqs, sizeof(*seqs) * seq_info.n);
    if (tmp_ptr2 == NULL) {
      kseq_destroy(kseq);
      badexit("Error: Failed to allocate memory for sequences.");
    } else {
      seqs = tmp_ptr2;
    }
    size_t *tmp_ptr3 = realloc(seq_sizes, sizeof(*seq_sizes) * seq_info.n);
    if (tmp_ptr3 == NULL) {
      kseq_destroy(kseq);
      badexit("Error: Failed to allocate memory for sequence sizes.");
    } else {
      seq_sizes = tmp_ptr3;
    }
    seqs[seq_info.n - 1] = (unsigned char *) kseq->seq.s;
    kseq->seq.s = NULL;
    seq_sizes[seq_info.n - 1] = kseq->seq.l;
    seq_names[seq_info.n - 1] = malloc(sizeof(char) * kseq->name.l + sizeof(char) * kseq->comment.l + 2);
    name_sizes += kseq->name.l + kseq->comment.l + 2;
    if (seq_names[seq_info.n - 1] == NULL) {
      kseq_destroy(kseq);
      badexit("Error: Failed to allocate memory for sequence name.");
    }
    add_seq_name(seq_names[seq_info.n - 1], kseq);
  }
  if (ret_val == -2) {
    kseq_destroy(kseq);
    badexit("Error: Failed to parse FASTQ qualities.");
  } else if (ret_val < -2) {
    kseq_destroy(kseq);
    badexit("Error: Failed to read input.");
  } else if (!seq_info.n) {
    kseq_destroy(kseq);
    badexit("Error: Failed to read any sequences from input.");
  }
  kseq_destroy(kseq);
  ERASE_ARRAY(char_counts, 256);
  count_bases();
  size_t seq_len_total = 0;
  for (size_t i = 0; i < seq_info.n; i++) seq_len_total += seq_sizes[i];
  if (!seq_len_total) {
    badexit("Error: Only encountered empty sequences.");
  }
  seq_info.total_bases = seq_len_total;
  seq_info.unknowns = seq_len_total - standard_base_count();
  seq_info.gc_pct = calc_gc() * 100.0;
  double unknowns_pct = 100.0 * seq_info.unknowns / seq_len_total;
  if (seq_info.unknowns == seq_len_total) {
    badexit("Error: Failed to read any standard DNA/RNA bases.");
  } else if (unknowns_pct >= 90.0) {
    fprintf(stderr,
      "!!! Warning: Non-standard base count is extremely high!!! (%.2f%%)\n",
      unknowns_pct);
  } else if (unknowns_pct >= 50.0 && args.v) {
    fprintf(stderr, "Warning: Non-standard base count is very high! (%.2f%%)\n",
      unknowns_pct);
  } else if (unknowns_pct >= 10.0 && args.v) {
    fprintf(stderr, "Warning: Non-standard base count seems high. (%.2f%%)\n",
      unknowns_pct);
  }
  if (char_counts[32] && args.v) {
    fprintf(stderr,
      "Warning: Found spaces (%'zu) in sequences, these will be treated as gaps.\n",
      char_counts[32]);
  }
  if (args.v) {
    fprintf(stderr, "Loaded %'zu sequence(s).\n    size=%'zu    GC=%.2f%%\n",
      seq_info.n, seq_len_total, seq_info.gc_pct);
    if (seq_info.unknowns) {
      fprintf(stderr, "Found %'zu (%.2f%%) non-standard bases.\n",
        seq_info.unknowns, unknowns_pct);
    }
    fprintf(stderr, "Approx. memory usage by sequence(s): %'.2f MB\n",
      b2mb(sizeof(unsigned char) * seq_len_total + sizeof(size_t) * seq_info.n * 2 +
        sizeof(char) * name_sizes));
  }
}

/*
void load_seqs_old(void) {
  size_t line_len = 0, line_num = 0, line_len_no_spaces = 0, seq_i = -1;
  int name_loaded = 0;
  char *line = NULL;
  size_t len = 0;
  ssize_t read;
  while ((read = getline(&line, &len, files.s)) != -1) {
    line_num++;
    if (line[read - 2] == '\r') {
      line[read - 2] = '\0';
      line_len = read - 2;
    } else if (line[read - 1] == '\n') {
      line[read - 1] = '\0';
      line_len = read - 1;
    } else {
      line_len = read;
    }
    if (line[0] == '>') {
      seq_i++;
      char **tmp_ptr1 = realloc(seq_names, sizeof(*seq_names) * (1 + seq_i));
      if (tmp_ptr1 == NULL) {
        free(line);
        badexit("Error: Failed to allocate memory for sequence names.");
      } else {
        seq_names = tmp_ptr1;
      }
      unsigned char **tmp_ptr2 = realloc(seqs, sizeof(*seqs) * (1 + seq_i));
      if (tmp_ptr2 == NULL) {
        free(line);
        badexit("Error: Failed to allocate memory for sequences.");
      } else {
        seqs = tmp_ptr2;
      }
      size_t *tmp_ptr3 = realloc(seq_sizes, sizeof(*seq_sizes) * (1 + seq_i));
      if (tmp_ptr3 == NULL) {
        free(line);
        badexit("Error: Failed to allocate memory for sequence sizes.");
      } else {
        seq_sizes = tmp_ptr3;
      }
      seq_sizes[seq_i] = 0;
      size_t *tmp_ptr6 = realloc(seq_real_sizes, sizeof(*seq_real_sizes) * (1 + seq_i));
      if (tmp_ptr6 == NULL) {
        free(line);
        badexit("Error: Failed to allocate memory for real sequence sizes.");
      } else {
        seq_real_sizes = tmp_ptr6;
      }
      seq_real_sizes[seq_i] = 0;
      size_t *tmp_ptr5 = realloc(seq_line_nums, sizeof(*seq_line_nums) * (1 + seq_i));
      if (tmp_ptr5 == NULL) {
        free(line);
        badexit("Error: Failed to allocate memory for sequence line numbers.");
      } else {
        seq_line_nums = tmp_ptr5;
      }
      seq_line_nums[seq_i] = line_num;
      seq_names[seq_i] = malloc(sizeof(**seq_names) * SEQ_NAME_MAX_CHAR);
      if (seq_names[seq_i] == NULL) {
        free(line);
        badexit("Error: Failed to allocate memory for sequence name.");
      }
      ERASE_ARRAY(seq_names[seq_i], SEQ_NAME_MAX_CHAR);
      for (size_t i = 0; i < SEQ_NAME_MAX_CHAR; i++) {
        if (line[i + 1] == '\0') break;
        seq_names[seq_i][i] = line[i + 1];
      }
      name_loaded = 1;
    } else if (name_loaded) {
      line_len_no_spaces = 0;
      for (size_t i = 0; i < line_len; i++) {
        if (line[i] != ' ') line_len_no_spaces++;
      }
      if (seq_sizes[seq_i] &&
          (seq_sizes[seq_i] + line_len_no_spaces) > seq_real_sizes[seq_i]) {
        size_t rl_size = MAX(SEQ_REALLOC_SIZE, line_len_no_spaces);
        unsigned char *tmp_ptr4 = realloc(seqs[seq_i],
          sizeof(**seqs) * (seq_real_sizes[seq_i] + rl_size));
        if (tmp_ptr4 == NULL) {
          free(line);
          badexit("Error: Failed to allocate memory for sequences.");
        }
        seqs[seq_i] = tmp_ptr4;
        seq_real_sizes[seq_i] += rl_size;
      } else if (!seq_sizes[seq_i]) {
        size_t rl_size = MAX(SEQ_REALLOC_SIZE, line_len_no_spaces);
        seqs[seq_i] = malloc(sizeof(**seqs) * rl_size);
        seq_real_sizes[seq_i] += rl_size;
        if (seqs[seq_i] == NULL) {
          free(line);
          badexit("Error: Failed to allocate memory for sequences.");
        }
      }
      for (size_t j = 0, i = 0; i < line_len; i++) {
        if (line[i] != ' ') {
          seqs[seq_i][j + seq_sizes[seq_i]] = line[i]; j++;
        }
      }
      seq_sizes[seq_i] += line_len_no_spaces;
    }
  }
  seq_info.line_count = line_num;
  free(line);
  if (seq_i == -1) {
    badexit("Error: Sequences don't appear to be fasta-formatted.");
  }
  seq_info.n = seq_i + 1;
  ERASE_ARRAY(char_counts, 256);
  count_bases();
  size_t seq_len_total = 0;
  for (size_t i = 0; i < seq_info.n; i++) seq_len_total += seq_sizes[i];
  if (!seq_len_total) {
    badexit("Error: Only encountered empty sequences.");
  }
  seq_info.total_bases = seq_len_total;
  seq_info.unknowns = seq_len_total - standard_base_count();
  seq_info.gc_pct = calc_gc() * 100.0;
  double unknowns_pct = 100.0 * seq_info.unknowns / seq_len_total;
  if (seq_info.unknowns == seq_len_total) {
    badexit("Error: Failed to read any standard DNA/RNA bases.");
  } else if (unknowns_pct >= 90.0) {
    fprintf(stderr,
      "!!! Warning: Non-standard base count is extremely high!!! (%.2f%%)\n",
      unknowns_pct);
  } else if (unknowns_pct >= 50.0 && args.v) {
    fprintf(stderr, "Warning: Non-standard base count is very high! (%.2f%%)\n",
      unknowns_pct);
  } else if (unknowns_pct >= 10.0 && args.v) {
    fprintf(stderr, "Warning: Non-standard base count seems high. (%.2f%%)\n",
      unknowns_pct);
  }
  if (char_counts[32] && args.v) {
    fprintf(stderr,
      "Internal warning: Found spaces (%'zu) in loaded sequences, alert maintainer.\n",
      char_counts[32]);
  }
  if (args.v) {
    fprintf(stderr, "Loaded %'zu sequence(s).\n    size=%'zu    GC=%.2f%%\n",
      seq_info.n, seq_len_total, seq_info.gc_pct);
    if (seq_info.unknowns) {
      fprintf(stderr, "Found %'zu (%.2f%%) non-standard bases.\n",
        seq_info.unknowns, unknowns_pct);
    }
    fprintf(stderr, "Approx. memory usage by sequence(s): %'.2f MB\n",
      b2mb(sizeof(unsigned char) * seq_len_total + sizeof(size_t) * seq_info.n +
        sizeof(char) * SEQ_NAME_MAX_CHAR * seq_info.n));
  }
}
*/

int char_arrays_are_equal(const char *arr1, const char *arr2, const size_t len) {
  int are_equal = 1;
  for (size_t i = 0; i < len; i++) {
    if (arr1[i] == '\0' && arr2[i] == '\0') {
      break;
    }
    if (arr1[i] != arr2[i]) {
      are_equal = 0;
      break;
    }
  }
  return are_equal;
}

void int_to_char_array(const size_t N, char *arr) {
  ERASE_ARRAY(arr, 128);
  sprintf(arr, "__N%zu", N);
}

int dedup_char_array(char *arr, const size_t arr_max_len, const size_t N) {
  size_t arr_len = 0, dedup_len = 0, success = 0, j = 0;
  char dedup[128];
  ERASE_ARRAY(dedup, 128);
  int_to_char_array(N, dedup);
  for (size_t i = 0; i < arr_max_len; i++) {
    if (arr[i] == '\0') {
      arr_len = i;
      break;
    }
  }
  for (size_t i = 0; i < 128; i++) {
    if (dedup[i] == '\0') {
      dedup_len = i + 1;
      break;
    }
  }
  if (arr_max_len - arr_len >= dedup_len) {
    for (size_t i = arr_len; i < arr_len + dedup_len; i++) {
      arr[i] = dedup[j];
      j++;
    }
    success = 1;
  }
  return success;
}

void find_motif_dupes(void) {
  if (motif_info.n == 1) return;
  size_t *is_dup = malloc(sizeof(size_t) * motif_info.n);
  if (is_dup == NULL) {
    badexit("Error: Failed to allocate memory for motif name duplication check.");
  }
  ERASE_ARRAY(is_dup, motif_info.n);
  for (size_t i = 0; i < motif_info.n - 1; i++) {
    for (size_t j = i + 1; j < motif_info.n; j++) {
      if (is_dup[j]) continue;
      if (char_arrays_are_equal(motifs[i]->name, motifs[j]->name, MAX_NAME_SIZE)) {
        is_dup[i] = 1; is_dup[j] = 1;
      }
    }
  }
  size_t dup_count = 0;
  for (size_t i = 0; i < motif_info.n; i++) dup_count += is_dup[i];
  if (dup_count) {
    if (args.dedup) {
      for (size_t i = 0; i < motif_info.n; i++) {
        if (is_dup[i]) {
          int success = dedup_char_array(motifs[i]->name, MAX_NAME_SIZE, i + 1);
          if (!success) {
            fprintf(stderr,
              "Error: Failed to deduplicate motif #%zu, name is too large.", i + 1);
            free(is_dup);
            badexit("");
          }
        }
      }
    } else {
      fprintf(stderr,
        "Error: Encountered duplicate motif name (use -d to deduplicate).");
      size_t to_print = 5;
      if (to_print > dup_count) to_print = dup_count;
      for (size_t i = 0; i < motif_info.n; i++) {
        if (is_dup[i]) {
          fprintf(stderr, "\n    L%zu #%zu: %s", motifs[i]->file_line_num, i + 1,
            motifs[i]->name);
          to_print--;
          if (!to_print) break;
        }
      }
      if (dup_count > 5) {
        fprintf(stderr, "\n    ...");
        fprintf(stderr, "\n    Found %'zu total non-unique names.", dup_count);
      }
      free(is_dup);
      badexit("");
    }
  }
  free(is_dup);
}

void find_seq_dupes(void) {
  if (seq_info.n == 1) return;
  size_t *is_dup = malloc(sizeof(size_t) * seq_info.n);
  ERASE_ARRAY(is_dup, seq_info.n);
  for (size_t i = 0; i < seq_info.n - 1; i++) {
    for (size_t j = i + 1; j < seq_info.n; j++) {
      if (is_dup[j]) continue;
      if (char_arrays_are_equal(seq_names[i], seq_names[j],
            MAX(strlen(seq_names[i]), strlen(seq_names[j])))) {
        is_dup[i] = 1; is_dup[j] = 1;
      }
    }
  }
  size_t dup_count = 0;
  for (size_t i = 0; i < seq_info.n; i++) dup_count += is_dup[i];
  if (dup_count) {
    if (args.dedup) {
      for (size_t i = 0; i < seq_info.n; i++) {
        if (is_dup[i]) {
          int success = dedup_char_array(seq_names[i], SEQ_NAME_MAX_CHAR, i + 1);
          if (!success) {
            fprintf(stderr,
              "Error: Failed to deduplicate sequence #%zu, name is too large.", i + 1);
            free(is_dup);
            badexit("");
          }
        }
      }
    } else {
      fprintf(stderr,
        "Error: Encountered duplicate sequence name (use -d to deduplicate).");
      size_t to_print = 5;
      if (to_print > dup_count) to_print = dup_count;
      for (size_t i = 0; i < seq_info.n; i++) {
        if (is_dup[i]) {
          fprintf(stderr, "\n    #%zu: %s", i + 1, seq_names[i]);
          to_print--;
          if (!to_print) break;
        }
      }
      if (dup_count > 5) {
        fprintf(stderr, "\n    ...");
        fprintf(stderr, "\n    Found %'zu total non-unique names.", dup_count);
      }
      free(is_dup);
      badexit("");
    }
  }
  free(is_dup);
}

static inline void score_subseq(const motif_t *motif, const unsigned char *seq, const size_t offset, int *score) {
  *score = 0;
  for (size_t i = 0; i < motif->size; i++) {
    *score += get_score(motif, seq[i + offset], i);
  }
}

static inline void score_subseq_rc(const motif_t *motif, const unsigned char *seq, const size_t offset, int *score) {
  *score = 0;
  for (size_t i = 0; i < motif->size; i++) {
    *score += get_score_rc(motif, seq[i + offset], i);
  }
}

void score_seq(const motif_t *motif, const size_t seq_i, const size_t seq_loc) {
  /* Performance notes:
   * - Tried generating unrolled score_subseq functions for each possible motif size,
   *   made it slower
   * - Tried using decreasing loops (i--) in score_subseq, didn't help
   * - Tried instrument profiling, didn't help/made it slower
   * - Slight improvement when using -O3/-Ofast versus -O2
   * - Slight improvements when using local const variables
   * - Slight improvement for RC when stuffing both strands in the same loop
   * - No improvement when calculating pvalue locally with pre-computed offset
   * - Tried providing a vector of reverse cumulative max scores per position to
   *   allow score_subseq to stop early if passing the threshold becomes impossible,
   *   made it slower
   * - Tried allowing for sequences to be the outer loop, and motifs the inner loop;
   *   surprisingly not any faster, even when using multiple threads for multiple
   *   sequences vs a single motif
   * - At this point I think the only thing that might help would be to change the
   *   basic logic of how scanning is performed, maybe see if more cache-friendly
   *   way of doing things are possible (for motif PWM or the actual sequence)
   */
  const unsigned char *seq = seqs[seq_loc];
  const char *seq_name = seq_names[seq_i];
  const size_t seq_size = seq_sizes[seq_i];
  const int mot_size = motif->size;
  if (seq_size < motif->size || motif->threshold == INT_MAX) return;
  int score = INT_MIN;
  if (args.scan_rc) {
    for (size_t i = 0; i <= seq_size - motif->size; i++) {
      score_subseq(motif, seq, i, &score);
      if (score >= motif->threshold) {
        fprintf(files.o, "%s\t%zu\t%zu\t+\t%s\t%.9g\t%.3f\t%.1f\t%.*s\n",
          seq_name,
          i + 1,
          i + motif->size,
          motif->name,
          score2pval(motif, score),
          score / PWM_INT_MULTIPLIER,
          100.0 * score / motif->max_score,
          mot_size,
          seq + i);
      }
      score_subseq_rc(motif, seq, i, &score);
      if (score >= motif->threshold) {
        fprintf(files.o, "%s\t%zu\t%zu\t-\t%s\t%.9g\t%.3f\t%.1f\t%.*s\n",
          seq_name,
          i + 1,
          i + motif->size,
          motif->name,
          score2pval(motif, score),
          score / PWM_INT_MULTIPLIER,
          100.0 * score / motif->max_score,
          mot_size,
          seq + i);
      }
    }
  } else {
    for (size_t i = 0; i <= seq_size - motif->size; i++) {
      score_subseq(motif, seq, i, &score);
      if (score >= motif->threshold) {
        fprintf(files.o, "%s\t%zu\t%zu\t+\t%s\t%.9g\t%.3f\t%.1f\t%.*s\n",
          seq_name,
          i + 1,
          i + motif->size,
          motif->name,
          score2pval(motif, score),
          score / PWM_INT_MULTIPLIER,
          100.0 * score / motif->max_score,
          mot_size,
          seq + i);
      }
    }
  }
}

void print_seq_stats_single(FILE *whereto, const size_t seq_i, const size_t seq_j) {
  ERASE_ARRAY(char_counts, 256);
  count_bases_single(seqs[seq_i], seq_sizes[seq_j]);
  fprintf(whereto, "%zu\t%s\t", seq_j + 1, seq_names[seq_j]);
  fprintf(whereto, "%zu\t", seq_sizes[seq_j]);
  if (!seq_sizes[seq_j]) {
    fprintf(whereto, "nan\t");
  } else {
    fprintf(whereto, "%.2f\t", calc_gc() * 100.0);
  }
  fprintf(whereto, "%zu\n", seq_sizes[seq_j] - standard_base_count());
}

void print_seq_stats(FILE *whereto) {
  for (size_t i = 0; i < seq_info.n; i++) {
    ERASE_ARRAY(char_counts, 256);
    count_bases_single(seqs[i], seq_sizes[i]);
    fprintf(whereto, "%zu\t%s\t", i + 1, seq_names[i]);
    fprintf(whereto, "%zu\t", seq_sizes[i]);
    if (!seq_sizes[i]) {
      fprintf(whereto, "nan\t");
    } else {
      fprintf(whereto, "%.2f\t", calc_gc() * 100.0);
    }
    fprintf(whereto, "%zu\n", seq_sizes[i] - standard_base_count());
  }
}

void add_consensus_motif(const char *consensus) {
  if (add_motif()) badexit("");
  ERASE_ARRAY(motifs[0]->name, MAX_NAME_SIZE);
  size_t i = 0;
  for (;;) {
    motifs[0]->name[i] = consensus[i];
    if (consensus[i] == '\0') {
      motifs[0]->size = i;
      break;
    }
    i++;
  }
  motifs[0]->name[i] = '\0';
  if (motifs[0]->size > MAX_MOTIF_SIZE / 5) {
    fprintf(stderr, "Error: Consensus sequence is too large (%zu>max=%zu).",
      motifs[0]->size, MAX_MOTIF_SIZE / 5);
  }
  size_t let_i;
  for (size_t pos = 0; pos < motifs[0]->size; pos++) {
    let_i = consensus2index[(unsigned char) consensus[pos]];
    if (let_i < 0) {
      fprintf(stderr, "Error: Encountered unknown letter in consensus (%c).",
        consensus[pos]);
      badexit("");
    }
    set_score(motifs[0], 'A', pos,
      calc_score(consensus2probs[let_i * 4 + 0], args.bkg[0]));
    set_score(motifs[0], 'C', pos,
      calc_score(consensus2probs[let_i * 4 + 1], args.bkg[1]));
    set_score(motifs[0], 'G', pos,
      calc_score(consensus2probs[let_i * 4 + 2], args.bkg[2]));
    set_score(motifs[0], 'T', pos,
      calc_score(consensus2probs[let_i * 4 + 3], args.bkg[3]));
  }
  complete_motifs();
}

void print_pb(const double prog) {
  int left = prog * PROGRESS_BAR_WIDTH;
  int right = PROGRESS_BAR_WIDTH - left;
  fprintf(stderr, "\r[%.*s%*s] %3d%%", left, PROGRESS_BAR_STRING, right, "",
      (int) (prog * 100.0));
  fflush(stderr);
}

void *scan_sub_process(void *thread_i) {
  for (size_t i = 0; i < motif_info.n; i++) {
    motif_t *motif = motifs[i];
    if (*((int *) thread_i) == motif->thread) {
      if (args.w && !args.progress) {
        fprintf(stderr, "    Scanning motif: %s\n", motif->name);
      }
      fill_cdf(motif);
      set_threshold(motif);
      for (size_t j = 0; j < seq_info.n; j++) {
        score_seq(motif, j, j);
      }
      if (args.progress) {
        pthread_mutex_lock(&pb_lock);
        pb_counter++;
        print_pb((double) pb_counter / motif_info.n);
        pthread_mutex_unlock(&pb_lock);
      }
    }
  }
  free(thread_i);
  return NULL;
}

int main(int argc, char **argv) {

  if (setlocale(LC_NUMERIC, "en_US") == NULL && args.v) {
    fprintf(stderr, "Warning: setlocale(LC_NUMERIC, \"en_US\") failed.\n");
  }

  motifs = malloc(sizeof(*motifs));
  if (motifs == NULL) {
    badexit("Error: Failed to allocate memory for motifs.");
  }
  seq_names = malloc(sizeof(*seq_names));
  if (seq_names == NULL) {
    badexit("Error: Failed to allocate memory for sequence names.");
  }
  seqs = malloc(sizeof(*seqs));
  if (seqs == NULL) {
    badexit("Error: Failed to allocate memory for sequences.");
  }
  seq_sizes = malloc(sizeof(*seq_sizes));
  if (seq_sizes == NULL) {
    badexit("Error: Failed to allocate memory for sequence sizes.");
  }
  threads = malloc(sizeof(pthread_t));
  if (threads == NULL) {
    badexit("Error: Failed to allocate memory for threads.");
  }

  kseq_t *kseq;
  char *user_bkg, *consensus;
  int has_motifs = 0, has_seqs = 0, has_consensus = 0;
  int use_stdout = 1, use_stdin = 0, use_manual_thresh = 0;
  size_t max_seq_size;

  int opt;

  while ((opt = getopt(argc, argv, "m:1:s:o:b:flt:p:n:j:dgrvwh0")) != -1) {
    switch (opt) {
      case 'm':
        if (has_consensus) {
          badexit("Error: -m and -1 cannot both be used.");
        }
        has_motifs = 1;
        files.m = fopen(optarg, "r");
        if (files.m == NULL) {
          fprintf(stderr, "Error: Failed to open motif file: %s", optarg);
          badexit("");
        }
        files.m_open = 1;
        break;
      case '1':
        if (has_motifs) {
          badexit("Error: -m and -1 cannot both be used.");
        }
        has_consensus = 1;
        consensus = optarg;
        break;
      case 's':
        has_seqs = 1;
        if (optarg[0] == '-' && optarg[1] == '\0') {
          files.s = gzdopen(fileno(stdin), "r");
          use_stdin = 1;
        } else {
          files.s = gzopen(optarg, "r");
          if (files.s == NULL) {
            fprintf(stderr, "Error: Failed to open sequence file: %s", optarg);
            badexit("");
          }
        }
        files.s_open = 1;
        break;
      case 'o':
        use_stdout = 0;
        files.o = fopen(optarg, "w");
        if (files.o == NULL) {
          fprintf(stderr, "Error: Failed to create output file: %s", optarg);
          badexit("");
        }
        files.o_open = 1;
        break;
      case 'b':
        args.use_user_bkg = 1;
        user_bkg = optarg;
        break;
      case 'f':
        args.scan_rc = 0;
        break;
      case 't':
        args.pvalue = atof(optarg);
        use_manual_thresh = 1;
        break;
      case 'p':
        args.pseudocount = atof(optarg);
        if (!args.pseudocount) {
          badexit("Error: -p must be a positive integer.");
        }
        break;
      case 'n':
        args.nsites = atoi(optarg);
        if (!args.nsites) {
          badexit("Error: -n must be a positive integer.");
        }
        break;
      case 'j':
        args.nthreads = atoi(optarg);
        if (!args.nthreads) {
          badexit("Error: -j must be a positive integer.");
        }
      case 'd':
        args.dedup = 1;
        break;
      case 'r':
        args.trim_names = 0;
        break;
      case 'l':
        args.low_mem = 0;
        break;
      case '0':
        args.thresh0 = 1;
        break;
      case 'g':
        args.progress = 1;
        break;
      case 'w':
        args.w = 1;
      case 'v':
        args.v = 1;
        break;
      case 'h':
        usage();
        return EXIT_SUCCESS;
      default:
        return EXIT_FAILURE;
    }
  }

  if (use_manual_thresh && args.thresh0) {
    badexit("Error: Cannot use both -t and -0.");
  } else if (use_manual_thresh && has_consensus) {
    badexit("Error: Cannot use both -1 and -0.");
  }

  if (use_stdout) {
    files.o = stdout;
    files.o_open = 1;
  }

  if (!has_seqs && !has_motifs && !has_consensus) {
    badexit("Error: Missing one of -m, -1, -s args.");
  }

  if (args.use_user_bkg) parse_user_bkg(user_bkg);

  if (has_consensus) {
    args.bkg[0] = 0.25; args.bkg[1] = 0.25; args.bkg[2] = 0.25; args.bkg[3] = 0.25;
    args.pvalue = 1;
    args.nsites = 1000;
    args.pseudocount = 1;
    add_consensus_motif(consensus);
    has_motifs = 1;
    motif_info.is_consensus = 1;
  } else if (has_motifs) {
    load_motifs();
    find_motif_dupes();
  }

  if (has_consensus || !has_seqs || !has_motifs || motif_info.n == 1) {
    if (args.nthreads > 1) {
      fprintf(stderr, "Note: Multi-threading not available for current inputs.\n");
    }
    args.nthreads = 1;
  }

  if (use_stdin || args.nthreads > 1) {
    if (args.low_mem && args.v) {
      fprintf(stderr, "Deactivating low-mem mode.\n");
      args.low_mem = 0;
    }
  }

  if (has_motifs) {
    pthread_t *tmp_threads = realloc(threads, sizeof(pthread_t) * args.nthreads);
    if (tmp_threads == NULL) {
      badexit("Error: Failed to re-allocate memory for threads.");
    }
    threads = tmp_threads;
    for (size_t i = 0; i < motif_info.n; i++) {
      motifs[i]->thread = ((double) i / motif_info.n) * args.nthreads;
    }
  }

  if (has_motifs && !has_seqs) {
    if (args.v) {
      fprintf(stderr,
        "No sequences provided, parsing + printing motifs before exit.\n");
    }
    time_t time1 = time(NULL);
    if (alloc_cdf()) badexit("");
    for (size_t i = 0; i < motif_info.n; i++) {
      fill_cdf(motifs[i]);
      set_threshold(motifs[i]);
      fprintf(files.o, "----------------------------------------\n");
      print_motif(motifs[i], i + 1);
    }
    fprintf(files.o, "----------------------------------------\n");
    free_cdf();
    time_t time2 = time(NULL);
    if (args.v) {
      time_t time3 = difftime(time2, time1);
      if (time3 > 1) {
        fprintf(stderr, "Needed %'zu seconds to print motifs.\n",
          (size_t) time3);
      }
    }
  }

  if (has_seqs) {
    kseq = kseq_init(files.s);
    time_t time1 = time(NULL);
    if (args.v) {
      if (args.low_mem) fprintf(stderr, "Peaking through sequences ...\n");
      else fprintf(stderr, "Reading sequences ...\n");
    }
    if (args.low_mem) {
      max_seq_size = peak_through_seqs(kseq);
    } else {
      load_seqs(kseq);
    }
    find_seq_dupes();
    time_t time2 = time(NULL);
    if (args.v) {
      time_t time3 = difftime(time2, time1);
      if (time3 > 1) {
        if (args.low_mem) {
          fprintf(stderr, "Needed %'zu seconds to peak through sequences.\n",
            (size_t) time3);
        } else {
          fprintf(stderr, "Needed %'zu seconds to load sequences.\n",
            (size_t) time3);
        }
      }
    }
    if (!has_motifs) {
      if (args.v) {
        fprintf(stderr, "No motifs provided, printing sequence stats before exit.\n");
      }
      fprintf(files.o, "##seqnum\tseqname\tsize\tgc_pct\tn_count\n");

      time_t time1 = time(NULL);

      if (args.low_mem) {
        for (size_t j = 0; j < seq_info.n; j++) {
          if (kseq_read(kseq) < 0) {
            badexit("Error: Failed to re-read input file.");
          } else {
            seqs[0] = (unsigned char *) kseq->seq.s;
          }
          print_seq_stats_single(files.o, 0, j);
        }
        kseq_destroy(kseq);
      } else {
        print_seq_stats(files.o);
      }

      time_t time2 = time(NULL);
      if (args.v) {
        time_t time3 = difftime(time2, time1);
        if (time3 > 1) {
          fprintf(stderr, "Needed %'zu seconds to print sequence stats.\n",
            (size_t) time3);
        }
      }

    }
  }

  if (has_seqs && has_motifs) {

    fprintf(files.o, "##minimotif v%s [ ", MINIMOTIF_VERSION);
    for (size_t i = 1; i < argc; i++) {
      fprintf(files.o, "%s ", argv[i]);
    }
    fprintf(files.o, "]\n");
    size_t motif_size = 0;
    for (size_t i = 0; i < motif_info.n; i++) motif_size += motifs[i]->size;
    fprintf(files.o,
      "##MotifCount=%zu MotifSize=%zu SeqCount=%zu SeqSize=%zu GC=%.2f%% Ns=%zu\n",
      motif_info.n, motif_size, seq_info.n, seq_info.total_bases, seq_info.gc_pct,
      seq_info.unknowns);
    fprintf(files.o, 
      "##seqname\tstart\tend\tstrand\tmotif\tpvalue\tscore\tscore_pct\tmatch\n");

    if (args.v) fprintf(stderr, "Scanning ...\n");
    time_t time1 = time(NULL);
    if (alloc_cdf()) badexit("");
    if (args.low_mem) {
      if (args.progress) print_pb(0.0);
      for (size_t i = 0; i < motif_info.n; i++) {
        if (args.w && !args.progress) {
          fprintf(stderr, "    Scanning motif: %s\n", motifs[i]->name);
        }
        fill_cdf(motifs[i]);
        set_threshold(motifs[i]);
        for (size_t j = 0; j < seq_info.n; j++) {
          if (args.w && !args.progress) {
            fprintf(stderr, "        Scanning sequence: %s\n", seq_names[j]);
          }
          if (kseq_read(kseq) < 0) {
            badexit("Error: Failed to re-read input file.");
          } else {
            seqs[0] = (unsigned char *) kseq->seq.s;
          }
          score_seq(motifs[i], j, 0);
        }
        gzrewind(files.s);
        kseq_rewind(kseq);
        if (args.progress) print_pb((i + 1.0) / motif_info.n);
      }
      free(seqs[0]);
      if (args.progress) fprintf(stderr, "\n");
    } else {
      if (args.progress) print_pb(0.0);
      for (size_t t = 0; t < args.nthreads; t++) {
        size_t *thread_i = malloc(sizeof(size_t *));
        if (thread_i == NULL) {
          badexit("Error: Failed to allocate memory for thread index.");
        }
        *thread_i = t;
        pthread_create(&threads[t], NULL, scan_sub_process, thread_i);
      }
      for (size_t t = 0; t < args.nthreads; t++) {
        pthread_join(threads[t], NULL);
      }
      if (args.progress) fprintf(stderr, "\n");
    }
    free_cdf();
    time_t time2 = time(NULL);
    time_t time3 = difftime(time2, time1);
    if (args.v) fprintf(stderr, "Done.\n");
    if (args.v && time3 > 1) fprintf(stderr, "Needed %'zu seconds to scan.\n",
        (size_t) time3);

  }

  close_files();
  free(threads);
  free_motifs();
  free_seqs();

  return EXIT_SUCCESS;

}

