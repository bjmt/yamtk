/*
 *   yamme: De novo motif elicitation
 *   Copyright (C) 2026  Benjamin Jean-Marie Tremblay
 *   (part of yamtk)
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
#include <errno.h>
#include <locale.h>
#include <getopt.h>
#include <math.h>
#include <limits.h>
#include <time.h>
#include <stdint.h>
#include <inttypes.h>
#include <zlib.h>
#include <pthread.h>
#include "kseq.h"
#include "khash.h"
#include "version.h"

KSEQ_INIT(gzFile, gzread)
KHASH_MAP_INIT_INT64(seed_h, uint64_t)
KHASH_SET_INIT_INT64(seed_set_h)

/* ---- Tunable macros ---- */

#ifndef HAMMING_MISMATCH
#define HAMMING_MISMATCH  1
#endif
#ifndef REFINE_PASSES
#define REFINE_PASSES     2
#endif
#ifndef TOP_K_SEEDS
#define TOP_K_SEEDS       4
#endif
#ifndef MAX_MOTIF_WIDTH
#define MAX_MOTIF_WIDTH   30
#endif
#ifndef WIDTH_WARN_AT
#define WIDTH_WARN_AT     20
#endif
#ifndef MIN_REFINE_HITS
#define MIN_REFINE_HITS   20
#endif
#ifndef MIN_IC_BITS
#define MIN_IC_BITS       0.5
#endif
#ifndef CONSENSUS_THRESHOLD
#define CONSENSUS_THRESHOLD 0.75
#endif

/* ---- Constants ---- */

#define MAX_NAME_SIZE           ((uint64_t) 256)
#define MAX_MOTIF_SIZE          ((uint64_t) 250)
#define AMBIGUITY_SCORE                -10000000
#define MIN_BKG_VALUE                      0.001
#define MAX_CDF_SIZE        ((uint64_t) 2097152)
#define PWM_INT_MULTIPLIER                1000.0
#define USER_BKG_MAX_SIZE       ((uint64_t) 256)
#define SEQ_NAME_MAX_CHAR       ((uint64_t) 512)
#define ALLOC_CHUNK_SIZE        ((uint64_t) 256)
#define DEFAULT_PSEUDOCOUNT                    1
#define DEFAULT_SHUFFLE_K                      2
#define DEFAULT_QVALUE_FILTER                1e-3
#define DEFAULT_HIT_PVAL                     1e-4
#define DEFAULT_STOP_PVAL                    1e-3
#define DEFAULT_MIN_W                          6
#define DEFAULT_MAX_W                         15
#define DEFAULT_N_MOTIFS                      10
#define DEFAULT_DEDUP_OVERLAP                0.5
#define MAX_K                                  9
#define FASTA_LINE_LEN                        60
#define PROGRESS_BAR_WIDTH                    60
#define PROGRESS_BAR_STRING \
  "============================================================"

#define VEC_ADD(VEC, X, VEC_LEN) \
  do { for (uint64_t Xi = 0; Xi < VEC_LEN; Xi++) VEC[Xi] += X; } while (0)
#define VEC_DIV(VEC, X, VEC_LEN) \
  do { for (uint64_t Xi = 0; Xi < VEC_LEN; Xi++) VEC[Xi] /= X; } while (0)
#define VEC_MIN(VEC, MIN_RES, VEC_LEN) \
  do { MIN_RES = VEC[0]; for (uint64_t Xi = 1; Xi < VEC_LEN; Xi++) { if (VEC[Xi] < MIN_RES) MIN_RES = VEC[Xi]; } } while (0)
#define ERASE_ARRAY(ARR, LEN) memset(ARR, 0, sizeof(ARR[0]) * (LEN))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define LIKELY(COND)   __builtin_expect(COND, 1)
#define UNLIKELY(COND) __builtin_expect(COND, 0)

#define MALLOC_OR_RET1(OBJ, SIZE) \
  do { OBJ = malloc(SIZE); if (OBJ == NULL) return 1; } while (0)
#define REALLOC_OR_RET1(OBJ, SIZE, TYPE) \
  do { TYPE *TMP = realloc(OBJ, (SIZE) * sizeof(TYPE)); if (TMP == NULL) return 1; OBJ = TMP; } while (0)

/* ---- Motif struct (extended from yamenr.c) ---- */

typedef struct motif_t {
  int         pwm[MAX_MOTIF_SIZE];
  int         pwm_rc[MAX_MOTIF_SIZE];
  double     *cdf;
  int         threshold;
  uint64_t    size;
  uint64_t    cdf_size;
  uint64_t    thread;
  uint64_t    file_line_num;
  int         min;
  int         max;
  int         max_score;
  int         min_score;
  int         cdf_max;
  int         cdf_offset;
  char        name[MAX_NAME_SIZE];
  double     *tmp_pdf;
  double      pwm_probs[50][4];  /* pseudocount-adj probs for output */
  uint64_t    nsites_actual;     /* aligned sites from last refinement */
} motif_t;

/* ---- Args ---- */

typedef struct {
  double   bkg[4];
  double   hit_pval;
  double   stop_pval;
  int      pseudocount;
  int      nthreads;
  int      scan_rc;
  int      mask;
  int      trim_names;
  int      use_user_bkg;
  int      shuffle_k;
  uint64_t seed;
  int      use_seed;
  int      min_w;
  int      max_w;
  int      n_motifs;
  double   dedup_overlap;
  double   qvalue_filter;
  int      progress;
  int      v;
  int      w;
} args_t;

static args_t args = {
  .bkg           = {0.25, 0.25, 0.25, 0.25},
  .hit_pval      = DEFAULT_HIT_PVAL,
  .stop_pval     = DEFAULT_STOP_PVAL,
  .pseudocount   = DEFAULT_PSEUDOCOUNT,
  .nthreads      = 1,
  .scan_rc       = 1,
  .mask          = 0,
  .trim_names    = 1,
  .use_user_bkg  = 0,
  .shuffle_k     = DEFAULT_SHUFFLE_K,
  .seed          = 0,
  .use_seed      = 0,
  .min_w         = DEFAULT_MIN_W,
  .max_w         = DEFAULT_MAX_W,
  .n_motifs      = DEFAULT_N_MOTIFS,
  .dedup_overlap = DEFAULT_DEDUP_OVERLAP,
  .qvalue_filter = DEFAULT_QVALUE_FILTER,
  .progress      = 0,
  .v             = 0,
  .w             = 0
};

static pthread_mutex_t pb_lock = PTHREAD_MUTEX_INITIALIZER;
static uint64_t pb_counter = 0;

static void print_pb(const double prog) {
  const int left = prog * PROGRESS_BAR_WIDTH;
  const int right = PROGRESS_BAR_WIDTH - left;
  fprintf(stderr, "\r[%.*s%*s] %3d%%", left, PROGRESS_BAR_STRING, right, "",
      (int)(prog * 100.0));
  fflush(stderr);
}

/* ---- Files ---- */

typedef struct {
  int    i_open, n_open, tsv_open, meme_open;
  gzFile i;
  gzFile n;
  FILE  *tsv;
  FILE  *meme;
} files_t;

static files_t files = { 0 };

static const char *tsv_path  = "motifs.tsv";
static const char *meme_path = "motifs.meme";

/* ---- Sequence sets ---- */

typedef struct {
  char           **names;
  unsigned char  **seqs;
  uint64_t        *sizes;
  uint64_t         n;
  uint64_t         n_alloc;
  uint64_t         total_bases;
  double           gc_pct;
  uint64_t         unknowns;
} seq_set_t;

static seq_set_t pos_set = {0};
static seq_set_t neg_set = {0};

/* ---- Discovery results ---- */

typedef struct {
  motif_t   motif;
  uint8_t **covered;    /* [si][byte] bitmask over positives */
  double    pvalue;
  double    qvalue;
  uint64_t  sites_pos;
  uint64_t  sites_neg;
  uint64_t  seqs_pos;
  uint64_t  seqs_neg;
  int       dropped;
  uint64_t  discovery_seq;  /* intra-width acceptance order, for deterministic merge */
} disc_result_t;

static disc_result_t *results       = NULL;
static uint64_t       n_results     = 0;
static uint64_t       n_results_alloc = 0;

typedef struct { uint64_t kmer; double pval; } seed_t;

/* ---- Per-thread discovery context ---- */

typedef struct {
  int             tid;            /* 0..nthreads-1, diagnostic only */
  unsigned char **seqs;           /* writable view of positives for this thread */
  /* CDF scratch buffers (point to thread-local arrays; motif->cdf/tmp_pdf
     are aliased here during convert_ppm_to_motif) */
  double         *cdf;
  double         *tmp_pdf;
  uint64_t        cdf_real_size;
  /* Thread-local accepted motifs; merged into global results[] after join */
  disc_result_t  *local_results;
  uint64_t        n_local;
  uint64_t        n_local_alloc;
} thread_ctx_t;

/* ---- Lookup tables (from yamenr.c / yamscan.c) ---- */

static uint64_t char_counts[256];

static const unsigned char char2index[256] = {
  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4, 4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4, 4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
  4,0,4,1,4,4,4,2,4,4,4,4,4,4,4,4, 4,4,4,4,3,3,4,4,4,4,4,4,4,4,4,4,
  4,0,4,1,4,4,4,2,4,4,4,4,4,4,4,4, 4,4,4,4,3,3,4,4,4,4,4,4,4,4,4,4,
  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4, 4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4, 4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4, 4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4, 4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4
};

/* Like char2index but treats lowercase a/c/g/t/u as ambiguity (index 4),
   used when -M masking is enabled to skip lowercase regions during scanning. */
static const unsigned char char2maskindex[256] = {
  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4, 4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4, 4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
  4,0,4,1,4,4,4,2,4,4,4,4,4,4,4,4, 4,4,4,4,3,3,4,4,4,4,4,4,4,4,4,4,
  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4, 4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4, 4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4, 4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4, 4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4, 4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4
};

static const unsigned char shuf_char2index[256] = {
  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4, 4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4, 4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
  4,0,4,1,4,4,4,2,4,4,4,4,4,4,4,4, 4,4,4,4,3,3,4,4,4,4,4,4,4,4,4,4,
  4,0,4,1,4,4,4,2,4,4,4,4,4,4,4,4, 4,4,4,4,3,3,4,4,4,4,4,4,4,4,4,4,
  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4, 4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4, 4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4, 4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4, 4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4
};

static const char index2dna[6] = "ACGTN";

static const uint64_t pow5[16] = {
  1,5,25,125,625,3125,15625,78125,390625,1953125,9765625,48828125,
  244140625,1220703125,6103515625ULL,30517578125ULL
};

/* complement of 2-bit index: A(0)↔T(3), C(1)↔G(2) */
static const int comp4[4] = {3, 2, 1, 0};

/* ---- PRNG (from yamshuf.c via yamenr.c) ---- */

typedef struct { uint64_t s[2]; } xrng_t;
static xrng_t xrng;

/* ---- Peak memory / timing (from yamenr.c) ---- */

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
static long peak_mem(void) { return 0; }
#else
#include <sys/resource.h>
static long peak_mem(void) {
  struct rusage r_mem;
  getrusage(RUSAGE_SELF, &r_mem);
#ifdef __linux__
  return r_mem.ru_maxrss * 1024;
#else
  return r_mem.ru_maxrss;
#endif
}
#endif

static void print_peak_mb(void) {
  long bytes = peak_mem();
  if (bytes > (1 << 30))
    fprintf(stderr, "Approx. peak memory usage: %'.2f GB.\n",
      (((double)bytes/1024.0)/1024.0)/1024.0);
  else if (bytes > (1 << 20))
    fprintf(stderr, "Approx. peak memory usage: %'.2f MB.\n",
      ((double)bytes/1024.0)/1024.0);
  else if (bytes)
    fprintf(stderr, "Approx. peak memory usage: %'.2f KB.\n",
      (double)bytes/1024.0);
}

static void print_time(const uint64_t s, const char *what) {
  if (s > 7200)
    fprintf(stderr, "Needed %'.2f hours to %s.\n", ((double)s/60.0)/60.0, what);
  else if (s > 120)
    fprintf(stderr, "Needed %'.2f minutes to %s.\n", (double)s/60.0, what);
  else if (s > 1)
    fprintf(stderr, "Needed %'" PRIu64 " seconds to %s.\n", s, what);
}

/* ---- Sequence set free (needed by badexit) ---- */

static void free_seq_set(seq_set_t *s) {
  if (s->names) {
    for (uint64_t i = 0; i < s->n; i++) free(s->names[i]);
    free(s->names);
  }
  if (s->seqs) {
    for (uint64_t i = 0; i < s->n; i++) free(s->seqs[i]);
    free(s->seqs);
  }
  free(s->sizes);
  memset(s, 0, sizeof(*s));
}

/* ---- CDF alloc/free (needed by badexit) ---- */

/* Thread context pool (lives in main_me; freed via free_thread_ctxs). */
static thread_ctx_t *thread_ctxs   = NULL;
static int           n_thread_ctxs = 0;

static int init_thread_ctx(thread_ctx_t *ctx, int tid) {
  memset(ctx, 0, sizeof(*ctx));
  ctx->tid = tid;
  ctx->cdf_real_size = 1;
  ctx->cdf     = malloc(sizeof(double));
  ctx->tmp_pdf = malloc(sizeof(double));
  if (!ctx->cdf || !ctx->tmp_pdf) { fprintf(stderr, "Error: init_thread_ctx failed.\n"); return 1; }
  /* seqs / local_results allocated separately once pos_set is loaded */
  return 0;
}

/* Allocate per-thread mutable copy of positives. Caller fills with memcpy. */
static int alloc_ctx_seqs(thread_ctx_t *ctx) {
  ctx->seqs = malloc(sizeof(unsigned char *) * pos_set.n);
  if (!ctx->seqs) return 1;
  for (uint64_t si = 0; si < pos_set.n; si++) {
    ctx->seqs[si] = malloc(pos_set.sizes[si]);
    if (!ctx->seqs[si]) return 1;
  }
  return 0;
}

static void free_thread_ctx(thread_ctx_t *ctx) {
  free(ctx->cdf); ctx->cdf = NULL;
  free(ctx->tmp_pdf); ctx->tmp_pdf = NULL;
  if (ctx->seqs) {
    for (uint64_t si = 0; si < pos_set.n; si++) free(ctx->seqs[si]);
    free(ctx->seqs); ctx->seqs = NULL;
  }
  /* local_results.motif/covered ownership transferred to global results[] on merge */
  free(ctx->local_results); ctx->local_results = NULL;
  ctx->n_local = 0; ctx->n_local_alloc = 0;
}

static void free_thread_ctxs(void) {
  if (!thread_ctxs) return;
  for (int t = 0; t < n_thread_ctxs; t++) free_thread_ctx(&thread_ctxs[t]);
  free(thread_ctxs); thread_ctxs = NULL; n_thread_ctxs = 0;
}

/* ---- Results free (needed by badexit) ---- */

static void free_results(void) {
  if (!results) return;
  for (uint64_t i = 0; i < n_results; i++) {
    if (results[i].covered) {
      for (uint64_t si = 0; si < pos_set.n; si++) free(results[i].covered[si]);
      free(results[i].covered);
    }
  }
  free(results); results = NULL; n_results = 0; n_results_alloc = 0;
}

/* ---- File close ---- */

static void close_files(void) {
  if (files.i_open)    gzclose(files.i);
  if (files.n_open)    gzclose(files.n);
  if (files.tsv_open)  fclose(files.tsv);
  if (files.meme_open) fclose(files.meme);
}

/* ---- badexit ---- */

static void badexit(const char *msg) {
  fprintf(stderr, "%s\nRun yamtk me -h for usage.\n", msg);
  free_results();
  free_seq_set(&pos_set);
  free_seq_set(&neg_set);
  free_thread_ctxs();
  close_files();
  exit(EXIT_FAILURE);
}

/* ---- String helpers (from yamenr.c) ---- */

static inline int str_to_double(char *str, double *res) {
  char *tmp; errno = 0;
  *res = strtod(str, &tmp);
  return (str == tmp || errno != 0 || *tmp != '\0');
}

static inline int str_to_int(char *str, int *res) {
  char *tmp; errno = 0;
  long int r = strtol(str, &tmp, 10);
  if (r > INT_MAX) return 1;
  *res = (int)r;
  return (str == tmp || errno != 0 || *tmp != '\0');
}

static inline int str_to_uint64_t(char *str, uint64_t *res) {
  char *tmp; errno = 0;
  *res = (uint64_t)strtoull(str, &tmp, 10);
  return (str == tmp || errno != 0 || *tmp != '\0');
}

/* ---- Background (from yamenr.c) ---- */

static int check_and_load_bkg(double *bkg) {
  if (bkg[0]==-1.0||bkg[1]==-1.0||bkg[2]==-1.0||bkg[3]==-1.0) {
    fprintf(stderr, "Error: Too few background values found (need 4)."); return 1;
  }
  double mn = 0; VEC_MIN(bkg, mn, 4);
  if (mn < MIN_BKG_VALUE) {
    if (args.v)
      fprintf(stderr, "Warning: Adjusting background values below min (%.2g<%.2g).\n",
        mn, MIN_BKG_VALUE);
    VEC_ADD(bkg, MIN_BKG_VALUE, 4);
  }
  double sum = bkg[0]+bkg[1]+bkg[2]+bkg[3];
  if (fabs(sum-1.0)>0.001 && args.v)
    fprintf(stderr, "Warning: Background doesn't sum to 1.0, adjusting (sum=%.3g).\n", sum);
  VEC_DIV(bkg, sum, 4);
  args.bkg[0]=bkg[0]; args.bkg[1]=bkg[1]; args.bkg[2]=bkg[2]; args.bkg[3]=bkg[3];
  return 0;
}

static void parse_user_bkg(const char *bkg_usr) {
  uint64_t i=0, j=0, bi=0;
  char bc[USER_BKG_MAX_SIZE];
  double b[] = {-1.0,-1.0,-1.0,-1.0};
  ERASE_ARRAY(bc, USER_BKG_MAX_SIZE);
  while (bkg_usr[i] != '\0') {
    if (bkg_usr[i] != ',' && bkg_usr[i] != ' ') { bc[j]=bkg_usr[i]; j++; }
    else if (bkg_usr[i] == ',') {
      if (bi > 2) badexit("Error: Too many background values (need 4).");
      if (str_to_double(bc, &b[bi])) {
        fprintf(stderr, "Error: Failed to parse background value: %s", bc);
        badexit("");
      }
      ERASE_ARRAY(bc, USER_BKG_MAX_SIZE); bi++; j=0;
    }
    i++;
  }
  if (str_to_double(bc, &b[3])) {
    fprintf(stderr, "Error: Failed to parse background value: %s", bc); badexit("");
  }
  if (check_and_load_bkg(b)) badexit("");
  if (args.w)
    fprintf(stderr, "Using background: A=%.3g C=%.3g G=%.3g T=%.3g\n",
      args.bkg[0], args.bkg[1], args.bkg[2], args.bkg[3]);
}

/* ---- Motif PWM math (from yamenr.c) ---- */

static void init_motif(motif_t *motif) {
  ERASE_ARRAY(motif->name, MAX_NAME_SIZE);
  motif->name[0]='m'; motif->name[1]='o'; motif->name[2]='t';
  motif->name[3]='i'; motif->name[4]='f'; motif->name[5]='\0';
  motif->size=0; motif->threshold=0; motif->max_score=0;
  motif->min=0; motif->max=0; motif->cdf_size=0; motif->file_line_num=0;
  motif->min_score=0; motif->cdf_max=0; motif->thread=0;
  motif->nsites_actual=0;
  for (uint64_t i = 0; i < MAX_MOTIF_SIZE; i++) { motif->pwm[i]=0; motif->pwm_rc[i]=0; }
  for (uint64_t i = 4; i < MAX_MOTIF_SIZE; i += 5) {
    motif->pwm[i]=AMBIGUITY_SCORE; motif->pwm_rc[i]=AMBIGUITY_SCORE;
  }
  memset(motif->pwm_probs, 0, sizeof(motif->pwm_probs));
}

static inline void set_score(motif_t *m, const unsigned char let, const uint64_t pos, const int s) {
  m->pwm[char2index[let] + pos*5] = s;
}
static inline void set_score_rc(motif_t *m, const unsigned char let, const uint64_t pos, const int s) {
  m->pwm_rc[char2index[let] + pos*5] = s;
}
static inline int get_score(const motif_t *m, const unsigned char let, const uint64_t pos, const unsigned char *tbl) {
  return m->pwm[tbl[let] + pos*5];
}
static inline int get_score_rc(const motif_t *m, const unsigned char let, const uint64_t pos, const unsigned char *tbl) {
  return m->pwm_rc[tbl[let] + pos*5];
}
static inline int get_score_i(const motif_t *m, const int i, const uint64_t pos) {
  return m->pwm[i + pos*5];
}

/* Parameterized log-odds score (no global nsites/pseudocount) */
static int calc_score_ns(const double prob_i, const double bkg_i, const int ns, const int pc) {
  double x = prob_i * ns;
  x += ((double)pc) / 4.0;
  x /= (double)(ns + pc);
  if (x <= 0.0) x = 1e-12;
  return (int)(log2(x / bkg_i) * PWM_INT_MULTIPLIER);
}

static void fill_pwm_rc(motif_t *m) {
  for (uint64_t pos=0; pos<m->size; pos++) {
    set_score_rc(m,'A',m->size-1-pos,get_score(m,'T',pos,char2index));
    set_score_rc(m,'C',m->size-1-pos,get_score(m,'G',pos,char2index));
    set_score_rc(m,'G',m->size-1-pos,get_score(m,'C',pos,char2index));
    set_score_rc(m,'T',m->size-1-pos,get_score(m,'A',pos,char2index));
  }
}

static int get_pwm_min(const motif_t *m) {
  int mn=0;
  for (uint64_t pos=0; pos<m->size; pos++)
    for (int let=0; let<4; let++) { int v=get_score_i(m,let,pos); if(v<mn) mn=v; }
  return mn;
}

static int get_pwm_max(const motif_t *m) {
  int mx=0;
  for (uint64_t pos=0; pos<m->size; pos++)
    for (int let=0; let<4; let++) { int v=get_score_i(m,let,pos); if(v>mx) mx=v; }
  return mx;
}

/* ---- CDF fill and threshold (from yamenr.c) ---- */

/* Prefix used on -w log lines so multi-thread output is readable.  Empty
   string when n_thread_ctxs <= 1 (no -j parallelism in effect). */
static void w_prefix(const thread_ctx_t *ctx, char *buf, size_t bufsz) {
  if (n_thread_ctxs > 1) snprintf(buf, bufsz, "[t%d] ", ctx->tid);
  else                   buf[0] = '\0';
}

static void fill_cdf(thread_ctx_t *ctx, motif_t *motif, const char *phase) {
  uint64_t max_step, s;
  double pdf_sum = 0.0;
  if (motif->cdf_size > MAX_CDF_SIZE) {
    fprintf(stderr,
      "\nError: CDF size for [%s] too large (%'" PRIu64 ">%'" PRIu64 ").\n",
      motif->name, motif->cdf_size, MAX_CDF_SIZE);
    badexit("");
  }
  if (ctx->cdf_real_size < motif->cdf_size) {
    double *r1 = realloc(ctx->cdf, motif->cdf_size*sizeof(double));
    if (!r1) badexit("Error: Memory re-allocation for CDF failed.");
    ctx->cdf = r1;
    double *r2 = realloc(ctx->tmp_pdf, motif->cdf_size*sizeof(double));
    if (!r2) badexit("Error: Memory re-allocation for PDF failed.");
    ctx->tmp_pdf = r2;
    ctx->cdf_real_size = motif->cdf_size;
  }
  motif->cdf     = ctx->cdf;
  motif->tmp_pdf = ctx->tmp_pdf;
  for (uint64_t i = 0; i < motif->cdf_size; i++) motif->cdf[i] = 1.0;
  for (uint64_t i = 0; i < motif->size; i++) {
    max_step = i * motif->cdf_max;
    for (uint64_t j = 0; j < motif->cdf_size; j++) motif->tmp_pdf[j] = motif->cdf[j];
    ERASE_ARRAY(motif->cdf, max_step + motif->cdf_max + 1);
    for (int j = 0; j < 4; j++) {
      s = (uint64_t)(get_score_i(motif, j, i) - motif->min);
      for (uint64_t k = 0; k <= max_step; k++)
        motif->cdf[k+s] += motif->tmp_pdf[k] * args.bkg[j];
    }
  }
  for (uint64_t i = 0; i < motif->cdf_size; i++) pdf_sum += motif->cdf[i];
  if (fabs(pdf_sum-1.0) > 0.0001) {
    for (uint64_t i = 0; i < motif->cdf_size; i++) motif->cdf[i] /= pdf_sum;
  }
  for (uint64_t i = motif->cdf_size-2; i < (uint64_t)-1; i--)
    motif->cdf[i] += motif->cdf[i+1];
  if (args.w && !args.progress) {
    char pre[16]; w_prefix(ctx, pre, sizeof(pre));
    fprintf(stderr, "        %sCDF [%s | %s] (n=%'" PRIu64 ") done.\n",
      pre, motif->name, phase, motif->cdf_size);
  }
}

static inline double score2pval(const motif_t *motif, const int score) {
  return motif->cdf[score - motif->cdf_offset];
}

static void set_threshold_pval(motif_t *motif, double pval) {
  uint64_t threshold_i = motif->cdf_size;
  for (uint64_t i = 0; i < motif->cdf_size; i++) {
    if (motif->cdf[i] < pval) { threshold_i = i; break; }
  }
  motif->threshold -= motif->min;
  motif->threshold *= (int)motif->size;
  motif->threshold = (int)threshold_i - motif->threshold;
  for (uint64_t i = 0; i < motif->size; i++) {
    int mx = get_score_i(motif, 0, i), mn = mx;
    for (int j = 1; j < 4; j++) {
      int t = get_score_i(motif, j, i);
      if (t > mx) mx = t;
      if (t < mn) mn = t;
    }
    motif->max_score += mx;
    motif->min_score += mn;
  }
  double min_pval = score2pval(motif, motif->max_score);
  if (min_pval / pval > 1.0001) {
    if (args.w)
      fprintf(stderr, "Warning: Min p-value for [%s] exceeds threshold.\n", motif->name);
    motif->threshold = INT_MAX;
  }
}

/* ---- Sequence loading (from yamenr.c) ---- */

static inline uint64_t standard_base_count(void) {
  return char_counts['A']+char_counts['a']+char_counts['C']+char_counts['c']+
         char_counts['G']+char_counts['g']+char_counts['T']+char_counts['t']+
         char_counts['U']+char_counts['u'];
}

static double calc_gc(void) {
  double gc = (double)(char_counts['G']+char_counts['C']+char_counts['g']+char_counts['c']);
  return gc / standard_base_count();
}

static void load_seq_set(kseq_t *kseq, seq_set_t *set, const char *label) {
  int ret;
  while ((ret = kseq_read(kseq)) >= 0) {
    set->n++;
    if (set->n > set->n_alloc) {
      char **t1 = realloc(set->names,
        sizeof(*set->names)*set->n_alloc + sizeof(*set->names)*ALLOC_CHUNK_SIZE);
      if (!t1) { kseq_destroy(kseq); badexit("Error: Failed to alloc seq names."); }
      set->names = t1;
      unsigned char **t2 = realloc(set->seqs,
        sizeof(*set->seqs)*set->n_alloc + sizeof(*set->seqs)*ALLOC_CHUNK_SIZE);
      if (!t2) { kseq_destroy(kseq); badexit("Error: Failed to alloc seqs."); }
      set->seqs = t2;
      uint64_t *t3 = realloc(set->sizes,
        sizeof(*set->sizes)*set->n_alloc + sizeof(*set->sizes)*ALLOC_CHUNK_SIZE);
      if (!t3) { kseq_destroy(kseq); badexit("Error: Failed to alloc seq sizes."); }
      set->sizes = t3;
      set->n_alloc += ALLOC_CHUNK_SIZE;
    }
    set->seqs[set->n-1] = (unsigned char *)kseq->seq.s;
    kseq->seq.s = NULL;
    set->sizes[set->n-1] = kseq->seq.l;
    size_t nl = kseq->name.l, cl = kseq->comment.l;
    if (args.trim_names || !cl) {
      if (nl > SEQ_NAME_MAX_CHAR) { kseq_destroy(kseq); badexit("Error: Seq name too large."); }
    } else {
      if (nl+cl+1 > SEQ_NAME_MAX_CHAR) { kseq_destroy(kseq); badexit("Error: Seq name too large."); }
    }
    set->names[set->n-1] = malloc(nl+cl+2);
    if (!set->names[set->n-1]) { kseq_destroy(kseq); badexit("Error: Failed to alloc seq name."); }
    if (args.trim_names || !cl) {
      memcpy(set->names[set->n-1], kseq->name.s, nl);
      set->names[set->n-1][nl] = '\0';
    } else {
      memcpy(set->names[set->n-1], kseq->name.s, nl);
      set->names[set->n-1][nl] = ' ';
      memcpy(set->names[set->n-1]+nl+1, kseq->comment.s, cl);
      set->names[set->n-1][nl+cl+1] = '\0';
    }
  }
  if (ret == -2) { kseq_destroy(kseq); badexit("Error: Failed to parse FASTQ qualities."); }
  else if (ret < -2) { kseq_destroy(kseq); badexit("Error: Failed to read input."); }
  else if (!set->n) { kseq_destroy(kseq); badexit("Error: No sequences read from input."); }
  kseq_destroy(kseq);
  ERASE_ARRAY(char_counts, 256);
  for (uint64_t i=0; i<set->n; i++)
    for (uint64_t j=0; j<set->sizes[i]; j++) char_counts[set->seqs[i][j]]++;
  uint64_t total=0;
  for (uint64_t i=0; i<set->n; i++) total += set->sizes[i];
  if (!total) badexit("Error: All sequences are empty.");
  set->total_bases = total;
  set->unknowns = total - standard_base_count();
  set->gc_pct = calc_gc() * 100.0;
  if (set->unknowns == total) badexit("Error: No standard DNA/RNA bases found.");
  double unk_pct = 100.0 * set->unknowns / total;
  if (unk_pct >= 90.0)
    fprintf(stderr, "!!! Warning: Non-standard base count extremely high in %s set (%.2f%%).\n",
      label, unk_pct);
  else if (unk_pct >= 50.0 && args.v)
    fprintf(stderr, "Warning: Non-standard base count very high in %s set (%.2f%%).\n",
      label, unk_pct);
  if (args.v)
    fprintf(stderr, "%s set: %'" PRIu64 " seq(s), %'" PRIu64 " bp, GC=%.2f%%\n",
      label, set->n, total, set->gc_pct);
}

/* ---- Compute background from current char_counts ---- */

static void compute_bkg_from_counts(void) {
  uint64_t tot = standard_base_count();
  if (!tot) return;
  double bkg[4];
  bkg[0] = (double)(char_counts['A']+char_counts['a']) / tot;
  bkg[1] = (double)(char_counts['C']+char_counts['c']) / tot;
  bkg[2] = (double)(char_counts['G']+char_counts['g']) / tot;
  bkg[3] = (double)(char_counts['T']+char_counts['t']+
                    char_counts['U']+char_counts['u']) / tot;
  if (check_and_load_bkg(bkg) && args.v)
    fprintf(stderr, "Warning: Could not compute background from sequences; using uniform.\n");
  if (args.v)
    fprintf(stderr, "Computed background: A=%.3g C=%.3g G=%.3g T=%.3g\n",
      args.bkg[0], args.bkg[1], args.bkg[2], args.bkg[3]);
}

/* ---- PRNG (from yamshuf.c via yamenr.c) ---- */

static inline uint64_t splitmix64(uint64_t x) {
  uint64_t z = (x += 0x9E3779B97F4A7C15ULL);
  z = (z ^ (z>>30)) * 0xBF58476D1CE4E5B9ULL;
  z = (z ^ (z>>27)) * 0x94D049BB133111EBULL;
  return z ^ (z>>31);
}

#define ROTL(X,K) (((X)<<(K))|((X)>>(64-(K))))

static inline uint64_t xrand_r(xrng_t *r) {
  const uint64_t s0=r->s[0]; uint64_t s1=r->s[1];
  const uint64_t result=ROTL(s0+s1,17)+s0;
  s1^=s0; r->s[0]=ROTL(s0,49)^s1^(s1<<21); r->s[1]=ROTL(s1,28);
  return result;
}

static inline void sxrand_r(xrng_t *r, uint64_t seed) {
  r->s[0]=splitmix64(seed); r->s[1]=splitmix64(r->s[0]);
}

/* ---- Shuffle functions (from yamshuf.c via yamenr.c) ---- */

static inline void shuf_swap(unsigned char *seq, const uint64_t i, const uint64_t j) {
  const unsigned char t=seq[i]; seq[i]=seq[j]; seq[j]=t;
}

static int shuffle_fisher_yates(unsigned char *seq, const uint64_t len) {
  for (uint64_t i=0, l=len-1; i<l; i++)
    shuf_swap(seq, i, i + xrand_r(&xrng)%(l-i));
  return 0;
}

static inline uint64_t chars2kmer(const unsigned char *seq, const uint64_t k, const uint64_t off) {
  uint64_t kmer=0;
  for (uint64_t j=0, i=k-1; i<(uint64_t)-1; j++, i--)
    kmer += pow5[i] * shuf_char2index[seq[off+j]];
  return kmer;
}

static inline void count_kmers(const unsigned char *seq, const uint64_t size, uint64_t *kt, const uint64_t k) {
  for (uint64_t i=0; i<size-k+1; i++) kt[chars2kmer(seq,k,i)]++;
}

static inline uint64_t cumsum_and_pick(const uint64_t *kmers) {
  const uint64_t k0=kmers[0], k1=k0+kmers[1], k2=k1+kmers[2], k3=k2+kmers[3], k4=k3+kmers[4];
  const uint64_t r=xrand_r(&xrng)%k4;
  return 4-((r<k0)+(r<k1)+(r<k2)+(r<k3));
}

#define COUNT_EDGES(O,T) ((T)[(O)]+(T)[(O)+1]+(T)[(O)+2]+(T)[(O)+3]+(T)[(O)+4])

static int shuffle_euler(unsigned char *seq, const uint64_t size, const uint64_t k,
                          uint64_t *kmer_tab, unsigned char *invalid_vertex,
                          uint64_t *euler_path, uint64_t *next_index) {
  for (uint64_t i=0; i<k-1; i++) seq[i]=index2dna[shuf_char2index[seq[i]]];
  seq[size-1]=index2dna[shuf_char2index[seq[size-1]]];
  uint64_t last_edge=chars2kmer(seq,k,size-k);
  kmer_tab[last_edge]--;
  for (uint64_t i=0,j=0; i<pow5[k-1]; i++,j+=5)
    if (!COUNT_EDGES(j,kmer_tab)) invalid_vertex[i]=1;
  const uint64_t final_vertex=chars2kmer(seq,k-1,size-k+1);
  invalid_vertex[final_vertex]=1;
  if (k>2) {
    for (uint64_t i=0,j=0,jmax=pow5[k-2]; i<pow5[k-1]; i++,j++) {
      if (j==jmax) j=0;
      next_index[i]=j*5;
    }
  }
  for (uint64_t u,i=0; i<pow5[k-1]; i++) {
    u=i;
    while (!invalid_vertex[u]) { euler_path[u]=cumsum_and_pick(kmer_tab+u*5); u=euler_path[u]+next_index[u]; }
    u=i;
    while (!invalid_vertex[u]) { invalid_vertex[u]=1; u=euler_path[u]+next_index[u]; }
  }
  for (uint64_t i=0,edge; i<pow5[k-1]; i++) {
    if (i==final_vertex) continue;
    edge=i*5+euler_path[i];
    if (edge!=last_edge&&kmer_tab[edge]) kmer_tab[edge]--;
  }
  for (uint64_t cv,ne,ki,i=k-2; i<size-2; i++) {
    cv=chars2kmer(seq,k-1,(i+2)-k);
    ki=cv*5;
    if (LIKELY(COUNT_EDGES(ki,kmer_tab))) {
      ne=cumsum_and_pick(kmer_tab+ki);
      kmer_tab[ne+ki]--;
    } else {
      ne=euler_path[cv];
    }
    seq[i+1]=index2dna[ne];
  }
  return 0;
}

static void make_shuffled_negatives(void) {
  if (args.v)
    fprintf(stderr, "Generating shuffled negative set (k=%d) ...\n", args.shuffle_k);
  neg_set.n       = pos_set.n;
  neg_set.n_alloc = pos_set.n;
  neg_set.seqs  = malloc(sizeof(unsigned char *) * pos_set.n);
  if (!neg_set.seqs) badexit("Error: Failed to allocate shuffled negatives.");
  neg_set.sizes = malloc(sizeof(uint64_t) * pos_set.n);
  if (!neg_set.sizes) badexit("Error: Failed to allocate shuffled negative sizes.");
  neg_set.names = NULL;
  uint64_t seed = args.use_seed ? args.seed : (uint64_t)time(NULL);
  sxrand_r(&xrng, seed);
  if (args.v) fprintf(stderr, "RNG seed: %" PRIu64 "\n", seed);
  if (args.shuffle_k == 1) {
    for (uint64_t i=0; i<pos_set.n; i++) {
      neg_set.sizes[i] = pos_set.sizes[i];
      neg_set.seqs[i]  = malloc(pos_set.sizes[i]);
      if (!neg_set.seqs[i]) badexit("Error: Failed to allocate shuffle buffer.");
      memcpy(neg_set.seqs[i], pos_set.seqs[i], pos_set.sizes[i]);
      if (pos_set.sizes[i] > 1) shuffle_fisher_yates(neg_set.seqs[i], pos_set.sizes[i]);
    }
  } else {
    uint64_t ksz = args.shuffle_k;
    uint64_t *kmer_tab    = calloc(pow5[ksz],   sizeof(uint64_t));
    unsigned char *inv_vtx= calloc(pow5[ksz-1], sizeof(unsigned char));
    uint64_t *euler_path  = calloc(pow5[ksz-1], sizeof(uint64_t));
    uint64_t *next_idx    = calloc(pow5[ksz-1], sizeof(uint64_t));
    if (!kmer_tab||!inv_vtx||!euler_path||!next_idx)
      badexit("Error: Failed to allocate Euler shuffle scratch.");
    for (uint64_t i=0; i<pos_set.n; i++) {
      neg_set.sizes[i] = pos_set.sizes[i];
      neg_set.seqs[i]  = malloc(pos_set.sizes[i]);
      if (!neg_set.seqs[i]) badexit("Error: Failed to allocate shuffle buffer.");
      memcpy(neg_set.seqs[i], pos_set.seqs[i], pos_set.sizes[i]);
      if (pos_set.sizes[i] >= ksz) {
        memset(kmer_tab,  0, sizeof(uint64_t)       * pow5[ksz]);
        memset(inv_vtx,   0, sizeof(unsigned char)  * pow5[ksz-1]);
        memset(euler_path,0, sizeof(uint64_t)       * pow5[ksz-1]);
        memset(next_idx,  0, sizeof(uint64_t)       * pow5[ksz-1]);
        count_kmers(neg_set.seqs[i], pos_set.sizes[i], kmer_tab, ksz);
        shuffle_euler(neg_set.seqs[i], pos_set.sizes[i], ksz,
                      kmer_tab, inv_vtx, euler_path, next_idx);
      }
    }
    free(kmer_tab); free(inv_vtx); free(euler_path); free(next_idx);
  }
  neg_set.total_bases = pos_set.total_bases;
  neg_set.gc_pct      = pos_set.gc_pct;
  neg_set.unknowns    = pos_set.unknowns;
  if (args.v)
    fprintf(stderr, "Negative (shuffled) set: %'" PRIu64 " seq(s), %'" PRIu64 " bp\n",
      neg_set.n, neg_set.total_bases);
}

/* ---- Statistical functions (from yamenr.c) ---- */

static double log_choose(double n, double k) {
  return lgamma(n+1) - lgamma(k+1) - lgamma(n-k+1);
}

static double fishers_exact_log_greater(uint64_t a, uint64_t b, uint64_t c, uint64_t d) {
  double r1=a+b, r2=c+d, c1=a+c, N=a+b+c+d;
  double log_denom = log_choose(N, c1);
  uint64_t x_max = (uint64_t)MIN(r1, c1);
  double log_sum = log_choose(r1,a) + log_choose(r2,c1-a) - log_denom;
  for (uint64_t x=a+1; x<=x_max; x++) {
    double term = log_choose(r1,x) + log_choose(r2,c1-x) - log_denom;
    if (term > log_sum) log_sum = term + log1p(exp(log_sum-term));
    else                log_sum = log_sum + log1p(exp(term-log_sum));
  }
  double p = exp(log_sum);
  return (p > 1.0) ? 1.0 : p;
}

typedef struct { double pval; uint64_t idx; } sp_t;
static int cmp_sp(const void *a, const void *b) {
  const sp_t *sa=(const sp_t *)a, *sb=(const sp_t *)b;
  return (sa->pval>sb->pval)-(sa->pval<sb->pval);
}

static void bh_qvalues(const double *p, uint64_t n, double *q) {
  sp_t *sorted = malloc(sizeof(sp_t)*n);
  if (!sorted) { for (uint64_t i=0;i<n;i++) q[i]=1.0; return; }
  for (uint64_t i=0;i<n;i++) { sorted[i].pval=p[i]; sorted[i].idx=i; }
  qsort(sorted, n, sizeof(sp_t), cmp_sp);
  double prev=1.0;
  for (uint64_t i=n; i>0; ) {
    i--;
    double qv = sorted[i].pval * n / (i+1);
    if (qv > prev) qv = prev;
    prev = qv;
    q[sorted[i].idx] = qv > 1.0 ? 1.0 : qv;
  }
  free(sorted);
}

/* ---- Scoring (from yamenr.c) ---- */

static inline void score_subseq(const motif_t *m, const unsigned char *seq,
                                 const uint64_t off, int *s, const unsigned char *tbl) {
  *s=0; for (uint64_t i=0;i<m->size;i++) *s+=get_score(m,seq[i+off],i,tbl);
}
static inline void score_subseq_rc(const motif_t *m, const unsigned char *seq,
                                    const uint64_t off, int *s, int *src,
                                    const unsigned char *tbl) {
  *s=0; *src=0;
  for (uint64_t i=0;i<m->size;i++) {
    *s   += get_score(m,seq[i+off],i,tbl);
    *src += get_score_rc(m,seq[i+off],i,tbl);
  }
}

/* ---- Discovery: k-mer encoding ---- */

/* Encode a w-mer starting at off using 2-bit (4-base) encoding.
   Returns UINT64_MAX if any position is non-ACGT under the given lookup
   table (char2index for normal, char2maskindex when -M is set so lowercase
   positions abort the encoding). */
static uint64_t kmer_encode4(const unsigned char *seq, const uint64_t w,
                              const uint64_t off, const unsigned char *tbl) {
  uint64_t kmer = 0;
  for (uint64_t i = 0; i < w; i++) {
    uint8_t idx = tbl[seq[off+i]];
    if (idx >= 4) return UINT64_MAX;
    kmer = kmer * 4 + idx;
  }
  return kmer;
}

/* Decode a 4-base encoded w-mer into a NUL-terminated ACGT string. */
static void kmer_to_string(uint64_t kmer, const uint64_t w, char *buf) {
  static const char letters[4] = {'A','C','G','T'};
  for (uint64_t i = w; i > 0; i--) {
    buf[i-1] = letters[kmer & 3];
    kmer >>= 2;
  }
  buf[w] = '\0';
}

/* Reverse-complement of a 4-base encoded w-mer */
static uint64_t rc_kmer4(uint64_t kmer, const uint64_t w) {
  uint64_t rc = 0;
  for (uint64_t i = 0; i < w; i++) {
    rc = rc * 4 + (uint64_t)comp4[kmer & 3];
    kmer >>= 2;
  }
  return rc;
}

/* Hamming distance between two 4-base encoded w-mers */
static uint64_t hamming4(uint64_t a, uint64_t b, const uint64_t w) {
  uint64_t mm = 0;
  for (uint64_t i = 0; i < w; i++) {
    if ((a & 3) != (b & 3)) mm++;
    a >>= 2; b >>= 2;
  }
  return mm;
}

/* ---- Discovery: seed enumeration ---- */

static int cmp_seed(const void *a, const void *b) {
  const seed_t *sa=(const seed_t *)a, *sb=(const seed_t *)b;
  return (sa->pval > sb->pval) - (sa->pval < sb->pval);
}

/* Enumerate top-K seeds of width w by per-sequence Fisher's p-value.
   Returns 0 on success. Sets *out_seeds (caller frees) and *out_n. */
static int enumerate_seeds(thread_ctx_t *ctx, uint64_t w, seed_t **out_seeds, uint64_t *out_n) {
  *out_seeds = NULL; *out_n = 0;
  const unsigned char *tbl = args.mask ? char2maskindex : char2index;

  khash_t(seed_h) *pos_h = kh_init(seed_h);
  khash_t(seed_h) *neg_h = kh_init(seed_h);
  if (!pos_h || !neg_h) {
    kh_destroy(seed_h, pos_h); kh_destroy(seed_h, neg_h); return 1;
  }

  /* Count per-sequence canonical k-mer presence in positives */
  for (uint64_t si = 0; si < pos_set.n; si++) {
    const unsigned char *seq = ctx->seqs[si];
    uint64_t seqlen = pos_set.sizes[si];
    if (seqlen < w) continue;

    khash_t(seed_set_h) *seen = kh_init(seed_set_h);
    if (!seen) {
      kh_destroy(seed_h, pos_h); kh_destroy(seed_h, neg_h); return 1;
    }

    for (uint64_t off = 0; off <= seqlen - w; off++) {
      uint64_t kmer = kmer_encode4(seq, w, off, tbl);
      if (kmer == UINT64_MAX) continue;
      uint64_t rc = args.scan_rc ? rc_kmer4(kmer, w) : kmer;
      uint64_t canon = (kmer <= rc) ? kmer : rc;
      int absent;
      kh_put(seed_set_h, seen, canon, &absent);
      if (absent < 0) { kh_destroy(seed_set_h, seen); kh_destroy(seed_h, pos_h); kh_destroy(seed_h, neg_h); return 1; }
    }

    khint_t k;
    for (k = kh_begin(seen); k != kh_end(seen); k++) {
      if (!kh_exist(seen, k)) continue;
      int absent;
      khint_t ph = kh_put(seed_h, pos_h, kh_key(seen, k), &absent);
      if (absent < 0) { kh_destroy(seed_set_h, seen); kh_destroy(seed_h, pos_h); kh_destroy(seed_h, neg_h); return 1; }
      if (absent) kh_val(pos_h, ph) = 0;
      kh_val(pos_h, ph)++;
    }
    kh_destroy(seed_set_h, seen);
  }

  /* Count per-sequence canonical k-mer presence in negatives */
  for (uint64_t si = 0; si < neg_set.n; si++) {
    const unsigned char *seq = neg_set.seqs[si];
    uint64_t seqlen = neg_set.sizes[si];
    if (seqlen < w) continue;

    khash_t(seed_set_h) *seen = kh_init(seed_set_h);
    if (!seen) {
      kh_destroy(seed_h, pos_h); kh_destroy(seed_h, neg_h); return 1;
    }

    for (uint64_t off = 0; off <= seqlen - w; off++) {
      uint64_t kmer = kmer_encode4(seq, w, off, tbl);
      if (kmer == UINT64_MAX) continue;
      uint64_t rc = args.scan_rc ? rc_kmer4(kmer, w) : kmer;
      uint64_t canon = (kmer <= rc) ? kmer : rc;
      int absent;
      kh_put(seed_set_h, seen, canon, &absent);
      if (absent < 0) { kh_destroy(seed_set_h, seen); kh_destroy(seed_h, pos_h); kh_destroy(seed_h, neg_h); return 1; }
    }

    khint_t k;
    for (k = kh_begin(seen); k != kh_end(seen); k++) {
      if (!kh_exist(seen, k)) continue;
      int absent;
      khint_t ph = kh_put(seed_h, neg_h, kh_key(seen, k), &absent);
      if (absent < 0) { kh_destroy(seed_set_h, seen); kh_destroy(seed_h, pos_h); kh_destroy(seed_h, neg_h); return 1; }
      if (absent) kh_val(neg_h, ph) = 0;
      kh_val(neg_h, ph)++;
    }
    kh_destroy(seed_set_h, seen);
  }

  /* Print -w diagnostics before destroying hashes */
  if (args.w) {
    uint64_t pos_uniq = (uint64_t)kh_size(pos_h);
    uint64_t neg_uniq = (uint64_t)kh_size(neg_h);
    uint64_t hbytes = ((uint64_t)kh_n_buckets(pos_h) + kh_n_buckets(neg_h))
                      * (sizeof(uint64_t) + sizeof(uint64_t) + 1);
    long rss = peak_mem();
    char pre[16]; w_prefix(ctx, pre, sizeof(pre));
    fprintf(stderr,
      "    %sEnumerated seeds  w=%" PRIu64 "  pos_unique_kmers=%" PRIu64
      "  neg_unique_kmers=%" PRIu64 "  hash_bytes=%" PRIu64
      "  peak_mem=%.2f MB\n",
      pre, w, pos_uniq, neg_uniq, hbytes, (double)rss / (1024.0*1024.0));
  }

  /* Score and collect candidates */
  uint64_t n_pos = pos_set.n, n_neg = neg_set.n;
  uint64_t n_cands = (uint64_t)kh_size(pos_h);
  if (n_cands == 0) {
    kh_destroy(seed_h, pos_h); kh_destroy(seed_h, neg_h); return 0;
  }

  seed_t *cands = malloc(sizeof(seed_t) * n_cands);
  if (!cands) {
    kh_destroy(seed_h, pos_h); kh_destroy(seed_h, neg_h); return 1;
  }

  uint64_t ci = 0;
  khint_t k;
  for (k = kh_begin(pos_h); k != kh_end(pos_h); k++) {
    if (!kh_exist(pos_h, k)) continue;
    uint64_t kmer = kh_key(pos_h, k);
    uint64_t pc   = kh_val(pos_h, k);
    uint64_t nc   = 0;
    khint_t nh    = kh_get(seed_h, neg_h, kmer);
    if (nh != kh_end(neg_h)) nc = kh_val(neg_h, nh);
    cands[ci].kmer = kmer;
    cands[ci].pval = fishers_exact_log_greater(pc, nc,
                       n_pos > pc ? n_pos - pc : 0,
                       n_neg > nc ? n_neg - nc : 0);
    ci++;
  }

  kh_destroy(seed_h, pos_h);
  kh_destroy(seed_h, neg_h);

  qsort(cands, n_cands, sizeof(seed_t), cmp_seed);

  *out_seeds = cands;
  *out_n     = (n_cands < TOP_K_SEEDS) ? n_cands : TOP_K_SEEDS;
  return 0;
}

/* ---- Discovery: PPM from seed via Hamming alignment ---- */

/* Build int PPM ppm[w][4] by aligning all positive windows within
   HAMMING_MISMATCH of seed_kmer (or its RC).  Returns aligned hit count. */
static int build_ppm_from_seed(thread_ctx_t *ctx, const uint64_t seed_kmer, const uint64_t w,
                                int ppm[][4]) {
  uint64_t rc_seed = args.scan_rc ? rc_kmer4(seed_kmer, w) : UINT64_MAX;
  int nsites = 0;
  const unsigned char *tbl = args.mask ? char2maskindex : char2index;

  for (uint64_t si = 0; si < pos_set.n; si++) {
    const unsigned char *seq = ctx->seqs[si];
    uint64_t seqlen = pos_set.sizes[si];
    if (seqlen < w) continue;

    for (uint64_t off = 0; off <= seqlen - w; off++) {
      uint64_t kmer = kmer_encode4(seq, w, off, tbl);
      if (kmer == UINT64_MAX) continue;

      uint64_t dist_fwd = hamming4(kmer, seed_kmer, w);
      if (dist_fwd <= HAMMING_MISMATCH) {
        for (uint64_t i = 0; i < w; i++) {
          uint8_t idx = tbl[seq[off+i]];
          if (idx < 4) ppm[i][idx]++;
        }
        nsites++;
      } else if (rc_seed != UINT64_MAX && hamming4(kmer, rc_seed, w) <= HAMMING_MISMATCH) {
        /* RC match: add reverse-complement of the window */
        for (uint64_t i = 0; i < w; i++) {
          uint8_t idx = tbl[seq[off + w - 1 - i]];
          if (idx < 4) ppm[i][comp4[idx]]++;
        }
        nsites++;
      }
    }
  }
  return nsites;
}

/* ---- Discovery: PPM → PWM ---- */

/* Convert int count PPM to log-odds PWM in motif, fill probs, CDF, threshold.
   `phase` is a short label used in -w output (e.g. "initial", "refine 1/2"). */
static void convert_ppm_to_motif(thread_ctx_t *ctx, motif_t *motif, int ppm[][4],
                                  const uint64_t w, const int nsites,
                                  const char *phase) {
  motif->nsites_actual = (uint64_t)nsites;
  motif->size          = w;

  for (uint64_t j = 0; j < w; j++) {
    double col = (double)(ppm[j][0]+ppm[j][1]+ppm[j][2]+ppm[j][3]);
    if (col <= 0.0) col = 1.0;
    for (int i = 0; i < 4; i++) {
      double raw = ppm[j][i] / col;
      /* Pseudocount-adjusted probability (stored for MEME/TSV) */
      motif->pwm_probs[j][i] = (raw * nsites + (double)args.pseudocount / 4.0)
                                / (nsites + args.pseudocount);
    }
    set_score(motif,'A',j,calc_score_ns(ppm[j][0]/col,args.bkg[0],nsites,args.pseudocount));
    set_score(motif,'C',j,calc_score_ns(ppm[j][1]/col,args.bkg[1],nsites,args.pseudocount));
    set_score(motif,'G',j,calc_score_ns(ppm[j][2]/col,args.bkg[2],nsites,args.pseudocount));
    set_score(motif,'T',j,calc_score_ns(ppm[j][3]/col,args.bkg[3],nsites,args.pseudocount));
  }

  motif->min        = get_pwm_min(motif);
  motif->max        = get_pwm_max(motif);
  motif->cdf_offset = motif->min * (int)w;
  fill_pwm_rc(motif);
  motif->cdf_max  = motif->max - motif->min;
  motif->cdf_size = w * (uint64_t)motif->cdf_max + 1;

  if (motif->cdf_max == 0 || motif->cdf_size > MAX_CDF_SIZE) {
    motif->threshold = INT_MAX;
    return;
  }

  fill_cdf(ctx, motif, phase);

  motif->threshold  = 0;
  motif->max_score  = 0;
  motif->min_score  = 0;
  set_threshold_pval(motif, args.hit_pval);

  if (motif->threshold == INT_MAX) {
    /* Fallback: relax to top 1% */
    if (args.w && !args.progress) {
      char pre[16]; w_prefix(ctx, pre, sizeof(pre));
      fprintf(stderr, "        %s  (threshold@p=%g unreachable; retrying at p=0.01)\n",
              pre, args.hit_pval);
    }
    motif->threshold = 0;
    motif->max_score = 0;
    motif->min_score = 0;
    set_threshold_pval(motif, 0.01);
  }
}

/* ---- Discovery: refinement ---- */

/* One refinement pass: score positives with current PWM threshold,
   rebuild PPM from hits, re-call convert_ppm_to_motif.
   Returns new nsites (>= MIN_REFINE_HITS) or 0 if bailed. */
static int refine_motif(thread_ctx_t *ctx, motif_t *motif, const uint64_t w,
                         int pass_idx) {
  if (motif->threshold == INT_MAX) return 0;
  int ppm[50][4];
  memset(ppm, 0, sizeof(ppm));
  int nsites = 0;
  const int thr = motif->threshold - 1;
  const unsigned char *tbl = args.mask ? char2maskindex : char2index;

  for (uint64_t si = 0; si < pos_set.n; si++) {
    const unsigned char *seq = ctx->seqs[si];
    uint64_t seqlen = pos_set.sizes[si];
    if (seqlen < w) continue;

    for (uint64_t off = 0; off <= seqlen - w; off++) {
      int score, src;
      if (args.scan_rc) {
        score_subseq_rc(motif, seq, off, &score, &src, tbl);
        if (score > thr) {
          for (uint64_t i = 0; i < w; i++) {
            uint8_t idx = tbl[seq[off+i]];
            if (idx < 4) ppm[i][idx]++;
          }
          nsites++;
        }
        if (src > thr) {
          for (uint64_t i = 0; i < w; i++) {
            uint8_t idx = tbl[seq[off + w - 1 - i]];
            if (idx < 4) ppm[i][comp4[idx]]++;
          }
          nsites++;
        }
      } else {
        score_subseq(motif, seq, off, &score, tbl);
        if (score > thr) {
          for (uint64_t i = 0; i < w; i++) {
            uint8_t idx = tbl[seq[off+i]];
            if (idx < 4) ppm[i][idx]++;
          }
          nsites++;
        }
      }
    }
  }

  if (args.w && !args.progress) {
    char pre[16]; w_prefix(ctx, pre, sizeof(pre));
    fprintf(stderr, "        %sRefinement %d/%d: %d hits\n",
            pre, pass_idx, REFINE_PASSES, nsites);
  }
  if (nsites < MIN_REFINE_HITS) {
    if (args.w && !args.progress) {
      char pre[16]; w_prefix(ctx, pre, sizeof(pre));
      fprintf(stderr, "        %s  (below MIN_REFINE_HITS=%d; bailing)\n",
              pre, MIN_REFINE_HITS);
    }
    return 0;
  }
  char phase[24];
  snprintf(phase, sizeof(phase), "refine %d/%d", pass_idx, REFINE_PASSES);
  convert_ppm_to_motif(ctx, motif, ppm, w, nsites, phase);
  return nsites;
}

/* ---- Discovery: motif evaluation with coverage tracking ---- */

/* Score positives and negatives; count seqs/sites with hits.
   Fills covered[si][byte] bitmask over positives when covered != NULL.
   Returns Fisher's p-value on per-sequence presence. */
static double evaluate_motif(thread_ctx_t *ctx, const motif_t *motif, const uint64_t w,
                              uint64_t *out_sites_pos, uint64_t *out_sites_neg,
                              uint64_t *out_seqs_pos,  uint64_t *out_seqs_neg,
                              uint8_t **covered) {
  const int thr = motif->threshold - 1;
  const unsigned char *tbl = args.mask ? char2maskindex : char2index;
  uint64_t seqs_pos = 0, seqs_neg = 0;
  uint64_t sites_pos = 0, sites_neg = 0;

  for (uint64_t si = 0; si < pos_set.n; si++) {
    const unsigned char *seq = ctx->seqs[si];
    uint64_t seqlen = pos_set.sizes[si];
    int has_hit = 0;
    if (seqlen < w) continue;

    for (uint64_t off = 0; off <= seqlen - w; off++) {
      int score, src;
      if (args.scan_rc) {
        score_subseq_rc(motif, seq, off, &score, &src, tbl);
        if (UNLIKELY(score > thr)) {
          has_hit = 1; sites_pos++;
          if (covered)
            for (uint64_t p = off; p < off+w; p++)
              covered[si][p/8] |= (uint8_t)(1u << (p%8));
        }
        if (UNLIKELY(src > thr)) {
          has_hit = 1; sites_pos++;
          if (covered)
            for (uint64_t p = off; p < off+w; p++)
              covered[si][p/8] |= (uint8_t)(1u << (p%8));
        }
      } else {
        score_subseq(motif, seq, off, &score, tbl);
        if (UNLIKELY(score > thr)) {
          has_hit = 1; sites_pos++;
          if (covered)
            for (uint64_t p = off; p < off+w; p++)
              covered[si][p/8] |= (uint8_t)(1u << (p%8));
        }
      }
    }
    if (has_hit) seqs_pos++;
  }

  for (uint64_t si = 0; si < neg_set.n; si++) {
    const unsigned char *seq = neg_set.seqs[si];
    uint64_t seqlen = neg_set.sizes[si];
    int has_hit = 0;
    if (seqlen < w) continue;

    for (uint64_t off = 0; off <= seqlen - w; off++) {
      int score, src;
      if (args.scan_rc) {
        score_subseq_rc(motif, seq, off, &score, &src, tbl);
        if (UNLIKELY(score > thr)) { has_hit = 1; sites_neg++; }
        if (UNLIKELY(src > thr))   { has_hit = 1; sites_neg++; }
      } else {
        score_subseq(motif, seq, off, &score, tbl);
        if (UNLIKELY(score > thr)) { has_hit = 1; sites_neg++; }
      }
    }
    if (has_hit) seqs_neg++;
  }

  *out_sites_pos = sites_pos;
  *out_sites_neg = sites_neg;
  *out_seqs_pos  = seqs_pos;
  *out_seqs_neg  = seqs_neg;

  uint64_t n_pos = pos_set.n, n_neg = neg_set.n;
  return fishers_exact_log_greater(seqs_pos, seqs_neg,
           n_pos > seqs_pos ? n_pos - seqs_pos : 0,
           n_neg > seqs_neg ? n_neg - seqs_neg : 0);
}

/* ---- Discovery: position masking ---- */

static void mask_positions(thread_ctx_t *ctx, uint8_t **covered) {
  for (uint64_t si = 0; si < pos_set.n; si++) {
    if (!covered[si]) continue;
    uint64_t nbytes = (pos_set.sizes[si] + 7) / 8;
    for (uint64_t b = 0; b < nbytes; b++) {
      if (!covered[si][b]) continue;
      for (int bit = 0; bit < 8; bit++) {
        if (covered[si][b] & (1u << bit)) {
          uint64_t pos = b * 8 + bit;
          if (pos < pos_set.sizes[si]) ctx->seqs[si][pos] = 'N';
        }
      }
    }
  }
}

/* ---- Discovery: add accepted motif to results array ---- */

/* Push an accepted motif onto the thread's local results list.  Results are
   merged (and renamed) by the main thread after all workers complete. */
static int add_result_local(thread_ctx_t *ctx, const motif_t *motif,
                            uint8_t **covered, double pvalue,
                            uint64_t sites_pos, uint64_t sites_neg,
                            uint64_t seqs_pos,  uint64_t seqs_neg,
                            uint64_t discovery_seq) {
  if (ctx->n_local + 1 > ctx->n_local_alloc) {
    uint64_t newcap = ctx->n_local_alloc + ALLOC_CHUNK_SIZE;
    disc_result_t *tmp = realloc(ctx->local_results, newcap * sizeof(disc_result_t));
    if (!tmp) return 1;
    ctx->local_results = tmp;
    ctx->n_local_alloc = newcap;
  }
  uint64_t ri = ctx->n_local++;
  disc_result_t *r = &ctx->local_results[ri];
  r->motif         = *motif;
  r->covered       = covered;
  r->pvalue        = pvalue;
  r->qvalue        = 1.0;
  r->sites_pos     = sites_pos;
  r->sites_neg     = sites_neg;
  r->seqs_pos      = seqs_pos;
  r->seqs_neg      = seqs_neg;
  r->dropped       = 0;
  r->discovery_seq = discovery_seq;
  return 0;
}

/* ---- Discovery: per-width discovery loop ---- */

static void discover_for_width(thread_ctx_t *ctx, const uint64_t w) {
  /* Refresh ctx->seqs from the pristine pos_set.seqs (each width pass starts
     from unmasked positives; mask_positions then mutates ctx->seqs only) */
  for (uint64_t si = 0; si < pos_set.n; si++)
    memcpy(ctx->seqs[si], pos_set.seqs[si], pos_set.sizes[si]);

  if (args.v && !args.progress) {
    char pre[16]; w_prefix(ctx, pre, sizeof(pre));
    fprintf(stderr, "  %sWidth %'" PRIu64 " ...\n", pre, w);
  }

  uint64_t intra = 0;  /* intra-width acceptance counter for deterministic merge */
  for (int ki = 0; ki < args.n_motifs; ki++) {
    seed_t *seeds = NULL; uint64_t n_seeds = 0;
    if (enumerate_seeds(ctx, w, &seeds, &n_seeds)) badexit("Error: enumerate_seeds failed.");
    if (n_seeds == 0) { free(seeds); break; }
    if (args.w && !args.progress) {
      char pre[16]; w_prefix(ctx, pre, sizeof(pre));
      fprintf(stderr, "    %sIteration %d/%d (w=%" PRIu64 "): trying %"
              PRIu64 " seed(s)\n", pre, ki + 1, args.n_motifs, w, n_seeds);
    }

    int accepted = 0;
    for (uint64_t si = 0; si < n_seeds && !accepted; si++) {
      char kbuf[64]; kmer_to_string(seeds[si].kmer, w, kbuf);
      if (args.w && !args.progress) {
        char pre[16]; w_prefix(ctx, pre, sizeof(pre));
        fprintf(stderr, "      %sSeed %" PRIu64 "/%" PRIu64
                " kmer=%s fisher_log_p=%.3g\n",
                pre, si + 1, n_seeds, kbuf, seeds[si].pval);
      }
      int ppm[50][4]; memset(ppm, 0, sizeof(ppm));
      int nsites = build_ppm_from_seed(ctx, seeds[si].kmer, w, ppm);
      if (nsites < MIN_REFINE_HITS) {
        if (args.w && !args.progress) {
          char pre[16]; w_prefix(ctx, pre, sizeof(pre));
          fprintf(stderr, "        %sPPM sites=%d (< MIN_REFINE_HITS=%d); skipping seed\n",
                  pre, nsites, MIN_REFINE_HITS);
        }
        continue;
      }
      if (args.w && !args.progress) {
        char pre[16]; w_prefix(ctx, pre, sizeof(pre));
        fprintf(stderr, "        %sInitial PPM sites=%d\n", pre, nsites);
      }

      motif_t motif; init_motif(&motif);
      /* Temporary name; final names assigned post-merge in main thread */
      snprintf(motif.name, MAX_NAME_SIZE, "w%" PRIu64 "_%" PRIu64, w, intra);
      convert_ppm_to_motif(ctx, &motif, ppm, w, nsites, "initial");
      if (motif.threshold == INT_MAX) {
        if (args.w && !args.progress) {
          char pre[16]; w_prefix(ctx, pre, sizeof(pre));
          fprintf(stderr, "        %sThreshold unreachable; skipping seed\n", pre);
        }
        continue;
      }

      for (int ri = 0; ri < REFINE_PASSES; ri++) {
        if (!refine_motif(ctx, &motif, w, ri + 1)) break;
      }
      if (motif.threshold == INT_MAX) continue;

      /* Allocate coverage bitmask */
      uint8_t **covered = calloc(pos_set.n, sizeof(uint8_t *));
      if (!covered) badexit("Error: Failed to alloc covered array.");
      int cov_ok = 1;
      for (uint64_t sj = 0; sj < pos_set.n; sj++) {
        uint64_t nbytes = (pos_set.sizes[sj] + 7) / 8;
        covered[sj] = calloc(nbytes, 1);
        if (!covered[sj]) { cov_ok = 0; break; }
      }
      if (!cov_ok) {
        for (uint64_t sj = 0; sj < pos_set.n; sj++) free(covered[sj]);
        free(covered);
        badexit("Error: Failed to alloc coverage bitmask.");
      }

      uint64_t sp, sn, ep, en;
      double pval = evaluate_motif(ctx, &motif, w, &sp, &sn, &ep, &en, covered);

      if (pval > args.stop_pval) {
        if (args.w && !args.progress) {
          char pre[16]; w_prefix(ctx, pre, sizeof(pre));
          fprintf(stderr, "        %sRejected: p=%.3g > stop_pval=%.3g\n",
                  pre, pval, args.stop_pval);
        }
        for (uint64_t sj = 0; sj < pos_set.n; sj++) free(covered[sj]);
        free(covered);
        continue;
      }

      if (add_result_local(ctx, &motif, covered, pval, sp, sn, ep, en, intra)) {
        for (uint64_t sj = 0; sj < pos_set.n; sj++) free(covered[sj]);
        free(covered);
        badexit("Error: Failed to add result.");
      }

      mask_positions(ctx, covered);
      accepted = 1;
      intra++;

      if (args.v && !args.progress) {
        char pre[16]; w_prefix(ctx, pre, sizeof(pre));
        fprintf(stderr, "    %sAccepted %s  w=%" PRIu64 "  p=%.3g  sites=%" PRIu64 "  seqs_pos=%" PRIu64 "/%" PRIu64 "\n",
          pre, motif.name, w, pval, sp, ep, pos_set.n);
      }
    }
    free(seeds);
    if (!accepted) {
      if (args.w && !args.progress) {
        char pre[16]; w_prefix(ctx, pre, sizeof(pre));
        fprintf(stderr, "    %sNo seed accepted; ending iterations for width %" PRIu64 "\n",
                pre, w);
      }
      break;
    }
  }
}

/* ---- Discovery: width sweep ---- */

/* Merge per-thread local_results into the global results[] in deterministic
   width-then-discovery_seq order, and assign final motif_<N> names. */
static int cmp_local_result(const void *a, const void *b) {
  const disc_result_t *ra = (const disc_result_t *)a;
  const disc_result_t *rb = (const disc_result_t *)b;
  if (ra->motif.size != rb->motif.size)
    return (ra->motif.size > rb->motif.size) - (ra->motif.size < rb->motif.size);
  return (ra->discovery_seq > rb->discovery_seq) - (ra->discovery_seq < rb->discovery_seq);
}

static void merge_thread_results(void) {
  uint64_t total = 0;
  for (int t = 0; t < n_thread_ctxs; t++) total += thread_ctxs[t].n_local;
  if (total == 0) return;
  if (total > n_results_alloc) {
    disc_result_t *tmp = realloc(results, total * sizeof(disc_result_t));
    if (!tmp) badexit("Error: Failed to alloc merged results.");
    results = tmp;
    n_results_alloc = total;
  }
  for (int t = 0; t < n_thread_ctxs; t++) {
    for (uint64_t i = 0; i < thread_ctxs[t].n_local; i++) {
      results[n_results++] = thread_ctxs[t].local_results[i];
    }
  }
  qsort(results, n_results, sizeof(disc_result_t), cmp_local_result);
  for (uint64_t i = 0; i < n_results; i++)
    snprintf(results[i].motif.name, MAX_NAME_SIZE, "motif_%" PRIu64, i + 1);
}

/* Tiny thread-safe FIFO of widths. Pre-populated before workers start. */
typedef struct {
  uint64_t       *widths;
  uint64_t        head;
  uint64_t        tail;
  pthread_mutex_t mu;
} width_queue_t;

typedef struct {
  thread_ctx_t   *ctx;
  width_queue_t  *q;
} worker_arg_t;

static void *worker(void *arg) {
  worker_arg_t *wa = (worker_arg_t *)arg;
  for (;;) {
    uint64_t w;
    pthread_mutex_lock(&wa->q->mu);
    if (wa->q->head >= wa->q->tail) {
      pthread_mutex_unlock(&wa->q->mu);
      break;
    }
    w = wa->q->widths[wa->q->head++];
    pthread_mutex_unlock(&wa->q->mu);
    discover_for_width(wa->ctx, w);
    if (args.progress) {
      pthread_mutex_lock(&pb_lock);
      pb_counter++;
      print_pb((double)pb_counter / (double)(args.max_w - args.min_w + 1));
      pthread_mutex_unlock(&pb_lock);
    }
  }
  return NULL;
}

static void discover_all(void) {
  /* Build width queue */
  uint64_t n_widths = (uint64_t)(args.max_w - args.min_w + 1);
  if (args.progress) print_pb(0.0);
  width_queue_t q;
  q.widths = malloc(sizeof(uint64_t) * n_widths);
  if (!q.widths) badexit("Error: Failed to alloc width queue.");
  for (uint64_t i = 0; i < n_widths; i++) q.widths[i] = (uint64_t)args.min_w + i;
  q.head = 0; q.tail = n_widths;
  if (pthread_mutex_init(&q.mu, NULL) != 0) badexit("Error: mutex init failed.");

  if (n_thread_ctxs <= 1) {
    /* Serial path: same code, single context, no pthread overhead */
    worker_arg_t wa = { .ctx = &thread_ctxs[0], .q = &q };
    worker(&wa);
  } else {
    pthread_t *threads = malloc(sizeof(pthread_t) * n_thread_ctxs);
    worker_arg_t *wargs = malloc(sizeof(worker_arg_t) * n_thread_ctxs);
    if (!threads || !wargs) badexit("Error: Failed to alloc thread arrays.");
    for (int t = 0; t < n_thread_ctxs; t++) {
      wargs[t].ctx = &thread_ctxs[t];
      wargs[t].q   = &q;
      if (pthread_create(&threads[t], NULL, worker, &wargs[t]) != 0)
        badexit("Error: pthread_create failed.");
    }
    for (int t = 0; t < n_thread_ctxs; t++) pthread_join(threads[t], NULL);
    free(threads); free(wargs);
  }

  pthread_mutex_destroy(&q.mu);
  free(q.widths);
  if (args.progress) fprintf(stderr, "\n");

  merge_thread_results();
}

/* ---- Rebuild coverage on unmasked sequences ---- */

/* After discover_all() the positives are fully restored.  Re-scan every
   accepted motif on those unmasked sequences so that dedup compares
   coverage bitmasks built under identical, unmasked conditions. */
static void rebuild_coverage(void) {
  /* Borrow thread 0's ctx but force it to read from the pristine pos_set.seqs
     (rebuild happens in the main thread after all discovery is complete). */
  thread_ctx_t *ctx = &thread_ctxs[0];
  unsigned char **saved = ctx->seqs;
  ctx->seqs = pos_set.seqs;
  for (uint64_t ri = 0; ri < n_results; ri++) {
    disc_result_t *r = &results[ri];
    for (uint64_t si = 0; si < pos_set.n; si++) {
      uint64_t nbytes = (pos_set.sizes[si] + 7) / 8;
      memset(r->covered[si], 0, nbytes);
    }
    uint64_t sites_pos, sites_neg, seqs_pos, seqs_neg;
    double pval = evaluate_motif(ctx, &r->motif, r->motif.size,
                                 &sites_pos, &sites_neg,
                                 &seqs_pos,  &seqs_neg, r->covered);
    r->pvalue    = pval;
    r->sites_pos = sites_pos;
    r->sites_neg = sites_neg;
    r->seqs_pos  = seqs_pos;
    r->seqs_neg  = seqs_neg;
  }
  ctx->seqs = saved;
}

/* ---- Cross-width dedup ---- */

static int cmp_result_pval(const void *a, const void *b) {
  const disc_result_t *ra=(const disc_result_t *)a, *rb=(const disc_result_t *)b;
  return (ra->pvalue > rb->pvalue) - (ra->pvalue < rb->pvalue);
}

static void dedup_results(void) {
  if (n_results < 2) return;
  qsort(results, n_results, sizeof(disc_result_t), cmp_result_pval);

  for (uint64_t i = 0; i < n_results; i++) {
    if (results[i].dropped) continue;
    for (uint64_t j = i+1; j < n_results; j++) {
      if (results[j].dropped) continue;
      uint64_t inter = 0, sz_i = 0, sz_j = 0;
      for (uint64_t si = 0; si < pos_set.n; si++) {
        uint64_t nbytes = (pos_set.sizes[si] + 7) / 8;
        for (uint64_t b = 0; b < nbytes; b++) {
          uint8_t ci = results[i].covered[si][b];
          uint8_t cj = results[j].covered[si][b];
          inter += (uint64_t)__builtin_popcount(ci & cj);
          sz_i  += (uint64_t)__builtin_popcount(ci);
          sz_j  += (uint64_t)__builtin_popcount(cj);
        }
      }
      uint64_t sz_min = sz_i < sz_j ? sz_i : sz_j;
      if (sz_min > 0 && (double)inter / sz_min > args.dedup_overlap)
        results[j].dropped = 1;
    }
  }
}

/* ---- Output: consensus ---- */

/* Map IUPAC base-set mask to letter.  Bit positions: 0=A, 1=C, 2=G, 3=T. */
static char iupac_from_mask(int mask) {
  switch (mask) {
    case 0x1: return 'A';
    case 0x2: return 'C';
    case 0x4: return 'G';
    case 0x8: return 'T';
    case 0x3: return 'M';  /* A,C */
    case 0x5: return 'R';  /* A,G */
    case 0x9: return 'W';  /* A,T */
    case 0x6: return 'S';  /* C,G */
    case 0xA: return 'Y';  /* C,T */
    case 0xC: return 'K';  /* G,T */
    case 0x7: return 'V';  /* A,C,G */
    case 0xB: return 'H';  /* A,C,T */
    case 0xD: return 'D';  /* A,G,T */
    case 0xE: return 'B';  /* C,G,T */
    default:  return 'N';
  }
}

/* Render a column as the smallest IUPAC code whose constituent bases together
   account for at least CONSENSUS_THRESHOLD of the probability mass. */
static char consensus_char(const double probs[4]) {
  int order[4] = {0, 1, 2, 3};
  /* Insertion sort descending by probability (4 elements; trivial). */
  for (int i = 1; i < 4; i++) {
    int j = i;
    while (j > 0 && probs[order[j]] > probs[order[j-1]]) {
      int t = order[j]; order[j] = order[j-1]; order[j-1] = t;
      j--;
    }
  }
  double cum = 0.0;
  int mask = 0;
  for (int k = 0; k < 4; k++) {
    cum += probs[order[k]];
    mask |= 1 << order[k];
    if (cum >= CONSENSUS_THRESHOLD) break;
  }
  return iupac_from_mask(mask);
}

static void build_consensus(const motif_t *m, char *buf) {
  for (uint64_t j = 0; j < m->size; j++) buf[j] = consensus_char(m->pwm_probs[j]);
  buf[m->size] = '\0';
}

/* ---- Output: IC trimming of low-information flanking columns ---- */

static double column_ic(const double probs[4]) {
  double ic = 2.0;  /* log2(|alphabet|) */
  for (int i = 0; i < 4; i++) {
    if (probs[i] > 0) ic += probs[i] * log2(probs[i]);
  }
  return ic < 0 ? 0 : ic;
}

/* Walk inward from each end of pwm_probs; strip flanking columns whose IC is
   below MIN_IC_BITS.  Updates motif.size and shifts the kept columns down to
   index 0.  No-ops if trimming would leave fewer than 3 informative columns. */
static void ic_trim_motif(motif_t *m) {
  uint64_t left = 0;
  while (left < m->size && column_ic(m->pwm_probs[left]) < MIN_IC_BITS) left++;
  uint64_t right = m->size;
  while (right > left && column_ic(m->pwm_probs[right - 1]) < MIN_IC_BITS) right--;
  uint64_t new_size = right - left;
  if (new_size < 3) return;
  if (left == 0 && right == m->size) return;
  for (uint64_t i = 0; i < new_size; i++) {
    for (int j = 0; j < 4; j++) m->pwm_probs[i][j] = m->pwm_probs[i + left][j];
  }
  m->size = new_size;
}

static void trim_all_results(void) {
  for (uint64_t i = 0; i < n_results; i++) {
    if (results[i].dropped) continue;
    ic_trim_motif(&results[i].motif);
  }
}

/* ---- Output: TSV ---- */

static void write_tsv(int argc, char **argv) {
  FILE *f = files.tsv;

  fprintf(f, "##yamme v%s [ ", YAMTK_VERSION);
  for (int i = 1; i < argc; i++) fprintf(f, "%s ", argv[i]);
  fprintf(f, "]\n");
  fprintf(f,
    "##motif\trank\twidth\tconsensus\tnsites"
    "\tseqs_pos\tseqs_neg\tsites_pos\tsites_neg"
    "\tn_pos\tn_neg\tpvalue\tqvalue\n");

  char cons[MAX_MOTIF_WIDTH + 2];
  uint64_t rank = 0;
  for (uint64_t ri = 0; ri < n_results; ri++) {
    if (results[ri].dropped) continue;
    rank++;
    build_consensus(&results[ri].motif, cons);
    fprintf(f,
      "%s\t%" PRIu64 "\t%" PRIu64 "\t%s\t%" PRIu64
      "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64
      "\t%" PRIu64 "\t%" PRIu64 "\t%.9g\t%.9g\n",
      results[ri].motif.name,
      rank,
      results[ri].motif.size,
      cons,
      results[ri].motif.nsites_actual,
      results[ri].seqs_pos,
      results[ri].seqs_neg,
      results[ri].sites_pos,
      results[ri].sites_neg,
      pos_set.n,
      neg_set.n,
      results[ri].pvalue,
      results[ri].qvalue);
  }
}

/* ---- Output: MEME ---- */

static void write_meme(void) {
  FILE *f = files.meme;
  char cons[MAX_MOTIF_WIDTH + 2];

  fprintf(f, "MEME version 4\n\n");
  fprintf(f, "ALPHABET= ACGT\n\n");
  fprintf(f, "strands: %s\n\n", args.scan_rc ? "+ -" : "+");
  fprintf(f, "Background letter frequencies:\n");
  fprintf(f, "A %.6f C %.6f G %.6f T %.6f\n\n",
    args.bkg[0], args.bkg[1], args.bkg[2], args.bkg[3]);

  for (uint64_t ri = 0; ri < n_results; ri++) {
    if (results[ri].dropped) continue;
    const motif_t *m = &results[ri].motif;
    build_consensus(m, cons);
    fprintf(f, "MOTIF %s %s\n\n", m->name, cons);
    fprintf(f, "letter-probability matrix: alength= 4 w= %" PRIu64
               " nsites= %" PRIu64 " E= %.4e\n",
      m->size, m->nsites_actual, results[ri].qvalue);
    for (uint64_t j = 0; j < m->size; j++) {
      fprintf(f, " %.6f %.6f %.6f %.6f\n",
        m->pwm_probs[j][0], m->pwm_probs[j][1],
        m->pwm_probs[j][2], m->pwm_probs[j][3]);
    }
    fprintf(f, "\n");
  }
}

/* ---- Usage ---- */

static void usage(void) {
  printf(
    "yamtk v%s  Copyright (C) %s  Benjamin Jean-Marie Tremblay\n"
    "Usage:  yamtk me [options] -i positives.fa[.gz]\n"
    "\n"
    " -i <str>   Positives FASTA/FASTQ ('-' = stdin, requires -n).\n"
    " -n <str>   Negatives FASTA (default: shuffle of positives).\n"
    " -o <str>   Output TSV file (default: motifs.tsv; '-' = stdout).\n"
    " -O <str>   Output MEME motif file (default: motifs.meme; '-' = stdout).\n"
    " -k <int>   Min motif width (default: %d, min 3).\n"
    " -K <int>   Max motif width (default: %d, max %d).\n"
    " -N <int>   Max motifs to discover (default: %d).\n"
    " -t <dbl>   Per-motif stopping p-value at discovery (default: %g).\n"
    " -P <dbl>   Hit-scoring p-value threshold (default: %g).\n"
    " -D <dbl>   Cross-width dedup overlap threshold (default: %g).\n"
    " -q <dbl>   BH q-value filter on surviving motifs (default: %g).\n"
    " -S <int>   Shuffle k-mer size for default negatives (default: %d, max %d).\n"
    " -b A,C,G,T Background (default: computed from positives).\n"
    " -p <int>   Pseudocount for PWM generation (default: %d).\n"
    " -R         Disable reverse-strand scoring.\n"
    " -M         Mask lower-case bases (skip scanning at those positions).\n"
    " -s <uint>  RNG seed (default: time-seeded).\n"
    " -j <int>   Threads (default: 1).\n"
    " -g         Show progress bar.\n"
    " -v / -w / -h   Verbose / very-verbose / help.\n"
    , YAMTK_VERSION, YAMTK_YEAR,
    DEFAULT_MIN_W, DEFAULT_MAX_W, MAX_MOTIF_WIDTH, DEFAULT_N_MOTIFS,
    DEFAULT_STOP_PVAL, DEFAULT_HIT_PVAL, DEFAULT_DEDUP_OVERLAP,
    DEFAULT_QVALUE_FILTER, DEFAULT_SHUFFLE_K, MAX_K, DEFAULT_PSEUDOCOUNT
  );
}

/* ---- main_me ---- */

int main_me(int argc, char **argv) {
  int opt;
  int has_pos = 0, has_neg = 0;
  int use_stdin_pos = 0;
  char *user_bkg  = NULL;
  char *seed_str  = NULL;

  struct timespec ts_program;
  clock_gettime(CLOCK_MONOTONIC, &ts_program);

  while ((opt = getopt(argc, argv, "i:n:o:O:k:K:N:t:P:D:S:b:p:q:RMs:j:gvwh")) != -1) {
    switch (opt) {
      case 'i':
        if (files.i_open) badexit("Error: -i specified more than once.");
        has_pos = 1;
        if (optarg[0]=='-' && optarg[1]=='\0') {
          files.i = gzdopen(fileno(stdin), "r");
          use_stdin_pos = 1;
        } else {
          files.i = gzopen(optarg, "r");
          if (!files.i) {
            fprintf(stderr, "Error: Cannot open -i file: %s\n", optarg);
            badexit("");
          }
        }
        files.i_open = 1;
        break;
      case 'n':
        if (files.n_open) badexit("Error: -n specified more than once.");
        has_neg = 1;
        files.n = gzopen(optarg, "r");
        if (!files.n) {
          fprintf(stderr, "Error: Cannot open -n file: %s\n", optarg);
          badexit("");
        }
        files.n_open = 1;
        break;
      case 'o':
        tsv_path = optarg;
        break;
      case 'O':
        if (optarg[0] == '\0') meme_path = NULL;
        else meme_path = optarg;
        break;
      case 'k':
        if (str_to_int(optarg, &args.min_w)) badexit("Error: Failed to parse -k value.");
        break;
      case 'K':
        if (str_to_int(optarg, &args.max_w)) badexit("Error: Failed to parse -K value.");
        break;
      case 'N':
        if (str_to_int(optarg, &args.n_motifs)) badexit("Error: Failed to parse -N value.");
        if (args.n_motifs < 0) badexit("Error: -N must be >= 0.");
        break;
      case 't':
        if (str_to_double(optarg, &args.stop_pval)) badexit("Error: Failed to parse -t value.");
        if (args.stop_pval <= 0.0 || args.stop_pval > 1.0) badexit("Error: -t must be in (0,1].");
        break;
      case 'P':
        if (str_to_double(optarg, &args.hit_pval)) badexit("Error: Failed to parse -P value.");
        if (args.hit_pval <= 0.0 || args.hit_pval > 1.0) badexit("Error: -P must be in (0,1].");
        break;
      case 'D':
        if (str_to_double(optarg, &args.dedup_overlap)) badexit("Error: Failed to parse -D value.");
        if (args.dedup_overlap < 0.0 || args.dedup_overlap > 1.0) badexit("Error: -D must be in [0,1].");
        break;
      case 'q':
        if (str_to_double(optarg, &args.qvalue_filter)) badexit("Error: Failed to parse -q value.");
        if (args.qvalue_filter < 0.0 || args.qvalue_filter > 1.0) badexit("Error: -q must be in [0,1].");
        break;
      case 'S':
        if (str_to_int(optarg, &args.shuffle_k)) badexit("Error: Failed to parse -S value.");
        if (args.shuffle_k < 1 || args.shuffle_k > MAX_K)
          badexit("Error: -S must be between 1 and 9.");
        break;
      case 'b':
        user_bkg = optarg;
        args.use_user_bkg = 1;
        break;
      case 'p':
        if (str_to_int(optarg, &args.pseudocount)) badexit("Error: Failed to parse -p value.");
        if (args.pseudocount < 0) badexit("Error: -p must be >= 0.");
        break;
      case 'R':
        args.scan_rc = 0;
        break;
      case 'M':
        args.mask = 1;
        break;
      case 's':
        seed_str = optarg;
        args.use_seed = 1;
        break;
      case 'j':
        if (str_to_int(optarg, &args.nthreads)) badexit("Error: Failed to parse -j value.");
        if (args.nthreads < 1) badexit("Error: -j must be >= 1.");
        break;
      case 'g':
        args.progress = 1;
        break;
      case 'w':
        args.w = 1;
        /* fall through */
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

  if (setlocale(LC_NUMERIC, "en_US") == NULL && args.v)
    fprintf(stderr, "Warning: setlocale(LC_NUMERIC, \"en_US\") failed.\n");

  if (!has_pos) badexit("Error: Missing -i (positives FASTA).");
  if (!has_neg && use_stdin_pos) badexit("Error: Cannot shuffle stdin input; provide -n.");

  if (args.min_w < 3) badexit("Error: -k must be >= 3.");
  if (args.max_w > MAX_MOTIF_WIDTH) {
    fprintf(stderr, "Error: -K must be <= %d.\n", MAX_MOTIF_WIDTH); badexit("");
  }
  if (args.min_w > args.max_w) badexit("Error: -k must be <= -K.");
  if (args.v && args.max_w > WIDTH_WARN_AT)
    fprintf(stderr, "Warning: -K > %d may use significant memory for k-mer enumeration.\n",
      WIDTH_WARN_AT);

  if (args.n_motifs == 0) {
    /* Short-circuit: no motifs requested */
    if (tsv_path && tsv_path[0]=='-' && tsv_path[1]=='\0') {
      files.tsv = stdout;
    } else if (tsv_path) {
      files.tsv = fopen(tsv_path, "w");
      if (!files.tsv) {
        fprintf(stderr, "Error: Cannot open output file: %s\n", tsv_path); badexit("");
      }
      files.tsv_open = 1;
    }
    if (files.tsv) {
      fprintf(files.tsv, "##yamme v%s [ ", YAMTK_VERSION);
      for (int i = 1; i < argc; i++) fprintf(files.tsv, "%s ", argv[i]);
      fprintf(files.tsv, "]\n");
      fprintf(files.tsv,
        "##motif\trank\twidth\tconsensus\tnsites"
        "\tseqs_pos\tseqs_neg\tsites_pos\tsites_neg"
        "\tn_pos\tn_neg\tpvalue\tqvalue\n");
    }
    close_files();
    return EXIT_SUCCESS;
  }

  if (args.use_seed && seed_str) {
    if (str_to_uint64_t(seed_str, &args.seed))
      badexit("Error: Failed to parse -s value.");
  }

  if (user_bkg) {
    parse_user_bkg(user_bkg);
  }

  /* Load positives */
  if (args.v) fprintf(stderr, "Loading positives ...\n");
  time_t t0 = time(NULL);
  pos_set.n_alloc = ALLOC_CHUNK_SIZE;
  pos_set.names   = malloc(sizeof(*pos_set.names) * ALLOC_CHUNK_SIZE);
  pos_set.seqs    = malloc(sizeof(*pos_set.seqs)  * ALLOC_CHUNK_SIZE);
  pos_set.sizes   = malloc(sizeof(*pos_set.sizes) * ALLOC_CHUNK_SIZE);
  if (!pos_set.names||!pos_set.seqs||!pos_set.sizes) badexit("Error: alloc pos_set failed.");
  kseq_t *kseq_pos = kseq_init(files.i);
  load_seq_set(kseq_pos, &pos_set, "positive");
  if (args.v) print_time((uint64_t)difftime(time(NULL), t0), "load positives");

  /* Compute background from positives unless user supplied */
  if (!args.use_user_bkg) compute_bkg_from_counts();

  /* Load or shuffle negatives */
  if (has_neg) {
    if (args.v) fprintf(stderr, "Loading negatives ...\n");
    t0 = time(NULL);
    neg_set.n_alloc = ALLOC_CHUNK_SIZE;
    neg_set.names   = malloc(sizeof(*neg_set.names) * ALLOC_CHUNK_SIZE);
    neg_set.seqs    = malloc(sizeof(*neg_set.seqs)  * ALLOC_CHUNK_SIZE);
    neg_set.sizes   = malloc(sizeof(*neg_set.sizes) * ALLOC_CHUNK_SIZE);
    if (!neg_set.names||!neg_set.seqs||!neg_set.sizes) badexit("Error: alloc neg_set failed.");
    kseq_t *kseq_neg = kseq_init(files.n);
    load_seq_set(kseq_neg, &neg_set, "negative");
    if (args.v) print_time((uint64_t)difftime(time(NULL), t0), "load negatives");
  } else {
    make_shuffled_negatives();
  }

  /* Open output files */
  if (tsv_path && tsv_path[0]=='-' && tsv_path[1]=='\0') {
    files.tsv = stdout;
  } else if (tsv_path) {
    files.tsv = fopen(tsv_path, "w");
    if (!files.tsv) {
      fprintf(stderr, "Error: Cannot open TSV output: %s\n", tsv_path); badexit("");
    }
    files.tsv_open = 1;
  }

  if (meme_path && meme_path[0]=='-' && meme_path[1]=='\0') {
    files.meme = stdout;
  } else if (meme_path) {
    files.meme = fopen(meme_path, "w");
    if (!files.meme) {
      fprintf(stderr, "Error: Cannot open MEME output: %s\n", meme_path); badexit("");
    }
    files.meme_open = 1;
  }

  /* Allocate thread contexts (one per worker).  Cap by the number of widths
     to discover so we never create idle workers. */
  int n_widths = args.max_w - args.min_w + 1;
  n_thread_ctxs = args.nthreads;
  if (n_thread_ctxs < 1) n_thread_ctxs = 1;
  if (n_thread_ctxs > n_widths) n_thread_ctxs = n_widths;
  if (n_thread_ctxs < 1) n_thread_ctxs = 1;
  thread_ctxs = malloc(sizeof(thread_ctx_t) * n_thread_ctxs);
  if (!thread_ctxs) badexit("Error: Failed to alloc thread contexts.");
  for (int t = 0; t < n_thread_ctxs; t++) {
    if (init_thread_ctx(&thread_ctxs[t], t)) badexit("Error: init_thread_ctx failed.");
    if (alloc_ctx_seqs(&thread_ctxs[t])) badexit("Error: alloc_ctx_seqs failed.");
  }
  if (args.v && n_thread_ctxs > 1)
    fprintf(stderr, "Discovery: %d worker thread(s) (width-level parallelism).\n",
            n_thread_ctxs);

  /* Discover motifs */
  if (args.v) fprintf(stderr, "Discovering motifs (widths %d-%d, max %d each) ...\n",
    args.min_w, args.max_w, args.n_motifs);
  t0 = time(NULL);
  discover_all();
  if (args.v) print_time((uint64_t)difftime(time(NULL), t0), "discover motifs");
  if (args.v) fprintf(stderr, "Found %'" PRIu64 " candidate motif(s) before dedup.\n", n_results);

  /* Rebuild coverage bitmasks on unmasked positives so dedup is accurate */
  if (n_results > 0) rebuild_coverage();

  /* Cross-width dedup */
  if (n_results > 1) {
    if (args.v) fprintf(stderr, "Deduplicating ...\n");
    dedup_results();
  }

  /* BH correction on survivors */
  uint64_t n_surv = 0;
  for (uint64_t ri = 0; ri < n_results; ri++) if (!results[ri].dropped) n_surv++;
  if (args.v) fprintf(stderr, "%'" PRIu64 " motif(s) after dedup.\n", n_surv);

  if (n_surv > 0) {
    double *pv = malloc(sizeof(double) * n_surv);
    double *qv = malloc(sizeof(double) * n_surv);
    uint64_t *idx = malloc(sizeof(uint64_t) * n_surv);
    if (!pv||!qv||!idx) badexit("Error: Failed to alloc BH arrays.");
    uint64_t ci = 0;
    for (uint64_t ri = 0; ri < n_results; ri++) {
      if (!results[ri].dropped) { pv[ci]=results[ri].pvalue; idx[ci]=ri; ci++; }
    }
    bh_qvalues(pv, n_surv, qv);
    for (uint64_t i = 0; i < n_surv; i++) results[idx[i]].qvalue = qv[i];
    free(pv); free(qv); free(idx);
  }

  /* Drop motifs whose BH q-value exceeds the user-set filter (default 1.0, no filter) */
  if (args.qvalue_filter < 1.0) {
    uint64_t n_kept = 0;
    for (uint64_t ri = 0; ri < n_results; ri++) {
      if (results[ri].dropped) continue;
      if (results[ri].qvalue > args.qvalue_filter) results[ri].dropped = 1;
      else n_kept++;
    }
    if (args.v) fprintf(stderr,
      "%'" PRIu64 " motif(s) survive q-value filter (-q %.3g).\n",
      n_kept, args.qvalue_filter);
  }

  /* Renumber survivors so motif names are sequential (motif_1, motif_2, ...)
     instead of preserving gaps left by dedup / q-value drops. */
  {
    uint64_t out_i = 0;
    for (uint64_t ri = 0; ri < n_results; ri++) {
      if (results[ri].dropped) continue;
      snprintf(results[ri].motif.name, MAX_NAME_SIZE, "motif_%" PRIu64, ++out_i);
    }
  }

  /* IC-trim low-information flanking columns before writing output */
  trim_all_results();

  /* Write output */
  if (files.tsv) write_tsv(argc, argv);
  if (files.meme) write_meme();

  /* Clean up */
  free_results();
  free_seq_set(&pos_set);
  free_seq_set(&neg_set);
  free_thread_ctxs();
  close_files();

  if (args.v) {
    struct timespec ts_end;
    clock_gettime(CLOCK_MONOTONIC, &ts_end);
    double elapsed = (double)(ts_end.tv_sec - ts_program.tv_sec)
                   + (double)(ts_end.tv_nsec - ts_program.tv_nsec) / 1e9;
    fprintf(stderr, "Done.\n");
    fprintf(stderr, "Total runtime: %.3fs\n", elapsed);
    print_peak_mb();
  }

  return EXIT_SUCCESS;
}
