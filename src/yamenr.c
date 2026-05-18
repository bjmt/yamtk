/*
 *   yamenr: Motif enrichment between positive and negative sequence sets
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
#include <pthread.h>
#include <stdint.h>
#include <inttypes.h>
#include <zlib.h>
#include "kseq.h"
#include "khash.h"
#include "version.h"

KSEQ_INIT(gzFile, gzread)
KHASH_MAP_INIT_STR(seq_str_h, uint64_t);
KHASH_SET_INIT_STR(motif_str_h);

#define MAX_NAME_SIZE           ((uint64_t) 256)
#define MAX_MOTIF_SIZE          ((uint64_t) 250)
#define AMBIGUITY_SCORE                -10000000
#define MIN_BKG_VALUE                      0.001
#define MAX_CDF_SIZE        ((uint64_t) 2097152)
#define PWM_INT_MULTIPLIER                1000.0
#define USER_BKG_MAX_SIZE       ((uint64_t) 256)
#define MEME_BKG_MAX_SIZE       ((uint64_t) 256)
#define MOTIF_VALUE_MAX_CHAR    ((uint64_t) 256)
#define SEQ_NAME_MAX_CHAR       ((uint64_t) 512)
#define ALLOC_CHUNK_SIZE        ((uint64_t) 256)
#define SEQ_REALLOC_SIZE                  524288
#define DEFAULT_NSITES                      1000
#define DEFAULT_PVALUE                    0.0001
#define DEFAULT_PSEUDOCOUNT                    1
#define DEFAULT_SHUFFLE_K                      2
#define DEFAULT_QVALUE_FILTER                0.1
#define MAX_K                                  9
#define FASTA_LINE_LEN                        60
#define BED_FIELD_MAX_CHAR      ((uint64_t) 256)
#define BED_ALLOC_CHUNK_SIZE    ((uint64_t) 256)
#define PROGRESS_BAR_WIDTH                    60
#define PROGRESS_BAR_STRING \
  "============================================================"

#define VEC_ADD(VEC, X, VEC_LEN) \
  do { for (uint64_t Xi = 0; Xi < VEC_LEN; Xi++) VEC[Xi] += X; } while (0)
#define VEC_DIV(VEC, X, VEC_LEN) \
  do { for (uint64_t Xi = 0; Xi < VEC_LEN; Xi++) VEC[Xi] /= X; } while (0)
#define VEC_SUM(VEC, SUM_RES, VEC_LEN) \
  do { SUM_RES = 0; for (uint64_t Xi = 0; Xi < VEC_LEN; Xi++) SUM_RES += VEC[Xi]; } while (0)
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

enum MOTIF_FMT {
  FMT_MEME = 1, FMT_HOMER = 2, FMT_JASPAR = 3, FMT_HOCOMOCO = 4, FMT_UNKNOWN = 5
};

enum TEST_MODE { TEST_SEQS = 0, TEST_SITES = 1, TEST_RANKSUM = 2 };

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

/* ---- Args ---- */

typedef struct args_t {
  double   bkg[4];
  double   pvalue;
  double   qvalue_filter;
  int      nsites;
  int      pseudocount;
  int      nthreads;
  int      scan_rc;
  int      mask;
  int      trim_names;
  int      use_user_bkg;
  int      shuffle_k;
  uint64_t seed;
  int      use_seed;
  int      dedup;
  int      test_mode;
  int      progress;
  int      v;
  int      w;
  int      low_mem;       /* -l: streaming, no all-seqs-in-memory */
  int      use_bed_pos;   /* -x */
  int      use_bed_neg;   /* -X (requires -n and -x) */
  const char *pos_path;   /* retained for -l shuffled re-open */
} args_t;

static args_t args = {
  .bkg          = {0.25, 0.25, 0.25, 0.25},
  .pvalue       = DEFAULT_PVALUE,
  .qvalue_filter = DEFAULT_QVALUE_FILTER,
  .nsites       = DEFAULT_NSITES,
  .pseudocount  = DEFAULT_PSEUDOCOUNT,
  .nthreads     = 1,
  .scan_rc      = 1,
  .mask         = 0,
  .trim_names   = 1,
  .use_user_bkg = 0,
  .shuffle_k    = DEFAULT_SHUFFLE_K,
  .seed         = 0,
  .use_seed     = 0,
  .dedup        = 0,
  .test_mode    = TEST_SEQS,
  .progress     = 0,
  .v            = 0,
  .w            = 0,
  .low_mem      = 0,
  .use_bed_pos  = 0,
  .use_bed_neg  = 0,
  .pos_path     = NULL
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

typedef struct files_t {
  int    m_open, i_open, n_open, o_open, x_pos_open, x_neg_open;
  FILE  *m;
  gzFile i;
  gzFile n;
  FILE  *o;
  gzFile x_pos;
  gzFile x_neg;
} files_t;

static files_t files = { .m_open=0, .i_open=0, .n_open=0, .o_open=0,
                         .x_pos_open=0, .x_neg_open=0 };

static khash_t(seq_str_h) *seq_hash_tab;

static void close_files(void) {
  if (files.m_open) fclose(files.m);
  if (files.i_open) gzclose(files.i);
  if (files.n_open) gzclose(files.n);
  if (files.o_open) fclose(files.o);
  if (files.x_pos_open) gzclose(files.x_pos);
  if (files.x_neg_open) gzclose(files.x_neg);
}

/* ---- Sequence sets ---- */

typedef struct seq_set_t {
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

/* Forward declarations for helpers defined later in the file. */
static void badexit(const char *msg);
static inline int str_to_uint64_t(char *str, uint64_t *res);

/* ---- BED reader ----
   Parametric on enr_bed_t so we can hold two independent BED sets (pos / neg).
   Closely modeled on yamseq.c's BED code; range_names are dropped (yamenr
   doesn't output region names) and seq_idx is added to cache per-region
   sequence-set indices once binding has run. */

typedef struct enr_bed_t {
  uint64_t  *starts;       /* 0-based                          */
  uint64_t  *ends;         /* 0-based, exclusive               */
  char      *strands;      /* '.', '+', '-'                    */
  char     **seq_names;
  uint64_t   n_regions;
  uint64_t   n_alloc;
  /* Filled by enr_bed_bind_seqs(): per-region index into seq_set,
     or UINT64_MAX if the BED seq is not present (fatal in fully-loaded mode). */
  uint64_t  *seq_idx;
  /* Filled by enr_bed_build_hash(): name -> head region index;
     seq_next[i] -> next region for the same seq (UINT64_MAX terminates). */
  khash_t(seq_str_h) *seq_hash;
  uint64_t  *seq_next;
} enr_bed_t;

static enr_bed_t bed_pos = {0};
static enr_bed_t bed_neg = {0};

static void enr_bed_free(enr_bed_t *b) {
  if (b->n_alloc) {
    for (uint64_t i = 0; i < b->n_regions; i++) free(b->seq_names[i]);
    free(b->seq_names);
    free(b->starts);
    free(b->ends);
    free(b->strands);
  }
  free(b->seq_idx);
  free(b->seq_next);
  if (b->seq_hash) kh_destroy(seq_str_h, b->seq_hash);
  memset(b, 0, sizeof(*b));
}

static void enr_bed_grow(enr_bed_t *b) {
  if (b->n_regions + 1 <= b->n_alloc) return;
  const uint64_t na = b->n_alloc + BED_ALLOC_CHUNK_SIZE;
  char     **p1 = realloc(b->seq_names, sizeof(*b->seq_names) * na);
  uint64_t  *p2 = realloc(b->starts,    sizeof(*b->starts)    * na);
  uint64_t  *p3 = realloc(b->ends,      sizeof(*b->ends)      * na);
  char      *p4 = realloc(b->strands,   sizeof(*b->strands)   * na);
  if (!p1 || !p2 || !p3 || !p4) badexit("Error: Failed to grow BED storage.");
  b->seq_names = p1; b->starts = p2; b->ends = p3; b->strands = p4;
  b->n_alloc = na;
}

static uint64_t bed_count_nonempty_chars(const char *line) {
  uint64_t n = 0;
  for (uint64_t i = 0; line[i] != '\0'; i++)
    if (line[i] != ' ' && line[i] != '\t' && line[i] != '\r' && line[i] != '\n') n++;
  return n;
}

static uint64_t bed_count_fields(const char *line) {
  uint64_t r = 1, i = 0;
  while (line[i] != '\0') { if (line[i] == '\t') r++; i++; }
  return r;
}

static uint64_t bed_field_size(const char *line, const uint64_t k) {
  uint64_t r = 0, i = 0, n = 0;
  while (line[i] != '\0') {
    if (line[i] == '\t') n++;
    else if (n + 1 == k) r++;
    else if (n + 1 > k) break;
    i++;
  }
  return r;
}

static uint64_t bed_field_start(const char *line, const uint64_t k) {
  uint64_t i = 0, n = 0;
  while (line[i] != '\0') {
    if (line[i] == '\t') n++;
    else if (n + 1 == k) break;
    i++;
  }
  return i;
}

static uint64_t bed_field_end(const char *line, const uint64_t k) {
  uint64_t i = 0, n = 0;
  while (line[i] != '\0') {
    if (line[i] == '\t') n++;
    if (n == k) { i--; break; }
    i++;
  }
  return i;
}

static inline uint64_t bed_parse_field(const char *line, const uint64_t k,
                                       char *field, const int no_spaces) {
  uint64_t s = bed_field_start(line, k);
  uint64_t e = bed_field_end(line, k);
  uint64_t sz = bed_field_size(line, k);
  uint64_t fi = 0;
  ERASE_ARRAY(field, BED_FIELD_MAX_CHAR);
  if (sz > 0 && sz < BED_FIELD_MAX_CHAR) {
    for (uint64_t i = s; i <= e; i++) {
      if (no_spaces && (line[i]==' '||line[i]=='\t'||line[i]=='\r'||line[i]=='\n')) sz--;
      else field[fi++] = line[i];
    }
  }
  return sz;
}

static void enr_bed_read(enr_bed_t *b, gzFile gz, const char *label) {
  kstream_t *kbed = ks_init(gz);
  if (!kbed) badexit("Error: Failed to init BED stream.");
  kstring_t line = { 0, 0, 0 };
  int rv;
  uint64_t ln = 0;
  char fld[BED_FIELD_MAX_CHAR];
  while ((rv = ks_getuntil(kbed, '\n', &line, 0)) >= 0) {
    ln++;
    if (bed_count_nonempty_chars(line.s) == 0) continue;
    if (line.s[0] == '#') continue;
    if (line.l >= 7 && memcmp(line.s, "browser", 7) == 0) continue;
    if (line.l >= 5 && memcmp(line.s, "track",   5) == 0) continue;
    const uint64_t nf = bed_count_fields(line.s);
    if (nf < 3) {
      fprintf(stderr, "Error: %s BED line %" PRIu64 " has %" PRIu64
        " fields; need >=3.\n", label, ln, nf);
      ks_destroy(kbed); free(line.s); badexit("");
    }
    enr_bed_grow(b);

    /* strand (col 6, optional) */
    if (nf >= 6) {
      uint64_t sz = bed_parse_field(line.s, 6, fld, 1);
      if (sz != 1 || (fld[0] != '+' && fld[0] != '-' && fld[0] != '.')) {
        fprintf(stderr, "Error: %s BED line %" PRIu64 " strand must be +/-/. (got '%s').\n",
          label, ln, fld);
        ks_destroy(kbed); free(line.s); badexit("");
      }
      b->strands[b->n_regions] = fld[0];
    } else {
      b->strands[b->n_regions] = '.';
    }

    /* start (col 2) */
    if (bed_parse_field(line.s, 2, fld, 1) == 0) {
      fprintf(stderr, "Error: %s BED line %" PRIu64 " empty start.\n", label, ln);
      ks_destroy(kbed); free(line.s); badexit("");
    }
    uint64_t tmp;
    if (str_to_uint64_t(fld, &tmp)) {
      fprintf(stderr, "Error: %s BED line %" PRIu64 " bad start '%s'.\n", label, ln, fld);
      ks_destroy(kbed); free(line.s); badexit("");
    }
    b->starts[b->n_regions] = tmp;

    /* end (col 3) */
    if (bed_parse_field(line.s, 3, fld, 1) == 0) {
      fprintf(stderr, "Error: %s BED line %" PRIu64 " empty end.\n", label, ln);
      ks_destroy(kbed); free(line.s); badexit("");
    }
    if (str_to_uint64_t(fld, &tmp)) {
      fprintf(stderr, "Error: %s BED line %" PRIu64 " bad end '%s'.\n", label, ln, fld);
      ks_destroy(kbed); free(line.s); badexit("");
    }
    b->ends[b->n_regions] = tmp;
    if (b->starts[b->n_regions] >= b->ends[b->n_regions]) {
      fprintf(stderr, "Error: %s BED line %" PRIu64 " has start >= end.\n", label, ln);
      ks_destroy(kbed); free(line.s); badexit("");
    }

    /* seq name (col 1) */
    const uint64_t sz1 = bed_parse_field(line.s, 1, fld, 0);
    if (sz1 == 0 || sz1 >= BED_FIELD_MAX_CHAR) {
      fprintf(stderr, "Error: %s BED line %" PRIu64 " empty/too-long seq name.\n", label, ln);
      ks_destroy(kbed); free(line.s); badexit("");
    }
    b->seq_names[b->n_regions] = malloc(sz1 + 1);
    if (!b->seq_names[b->n_regions]) {
      ks_destroy(kbed); free(line.s); badexit("Error: BED alloc failed.");
    }
    memcpy(b->seq_names[b->n_regions], fld, sz1);
    b->seq_names[b->n_regions][sz1] = '\0';
    /* Trim at first space (kseq drops everything after the first space in
       FASTA names; BED seq names must match the same canonical form). */
    for (uint64_t i = 0; i < sz1; i++) {
      if (b->seq_names[b->n_regions][i] == ' ') {
        b->seq_names[b->n_regions][i] = '\0'; break;
      }
    }
    b->n_regions++;
  }
  if (rv == -3) { ks_destroy(kbed); free(line.s); badexit("Error: Failed to read BED stream."); }
  ks_destroy(kbed);
  free(line.s);
  if (!b->n_regions) badexit("Error: BED has no usable records.");
  if (args.v) fprintf(stderr, "Read %" PRIu64 " %s BED region(s).\n", b->n_regions, label);
}

/* Resolve each BED region's seq name against a seq_set; verify bounds.
   Fatal on unknown seq or out-of-bounds end. Sets b->seq_idx[i]. */
static void enr_bed_bind_seqs(enr_bed_t *b, const seq_set_t *s, const char *label) {
  /* Build a name->index map for the seq_set. */
  khash_t(seq_str_h) *h = kh_init(seq_str_h);
  if (!h) badexit("Error: hash init failed.");
  int absent;
  for (uint64_t i = 0; i < s->n; i++) {
    khint_t k = kh_put(seq_str_h, h, s->names[i], &absent);
    if (absent == -1) { kh_destroy(seq_str_h, h); badexit("Error: hash put failed."); }
    if (absent == 0) {
      fprintf(stderr, "Error: Duplicate sequence name '%s' in %s set.\n", s->names[i], label);
      kh_destroy(seq_str_h, h); badexit("");
    }
    kh_val(h, k) = i;
  }
  b->seq_idx = malloc(sizeof(*b->seq_idx) * b->n_regions);
  if (!b->seq_idx) { kh_destroy(seq_str_h, h); badexit("Error: alloc failed."); }
  for (uint64_t i = 0; i < b->n_regions; i++) {
    khint_t k = kh_get(seq_str_h, h, b->seq_names[i]);
    if (k == kh_end(h)) {
      fprintf(stderr, "Error: %s BED references unknown sequence '%s'.\n", label, b->seq_names[i]);
      kh_destroy(seq_str_h, h); badexit("");
    }
    uint64_t si = kh_val(h, k);
    b->seq_idx[i] = si;
    if (b->ends[i] > s->sizes[si]) {
      fprintf(stderr,
        "Error: %s BED region '%s:%" PRIu64 "-%" PRIu64 "' exceeds seq length %" PRIu64 ".\n",
        label, b->seq_names[i], b->starts[i], b->ends[i], s->sizes[si]);
      kh_destroy(seq_str_h, h); badexit("");
    }
  }
  kh_destroy(seq_str_h, h);
}

/* Streaming-mode helper: build name -> region-list linked structure. */
static void enr_bed_build_hash(enr_bed_t *b) {
  b->seq_hash = kh_init(seq_str_h);
  if (!b->seq_hash) badexit("Error: BED hash init failed.");
  b->seq_next = malloc(sizeof(*b->seq_next) * b->n_regions);
  if (!b->seq_next) badexit("Error: BED next alloc failed.");
  for (uint64_t i = 0; i < b->n_regions; i++) b->seq_next[i] = UINT64_MAX;
  int absent;
  for (uint64_t i = 0; i < b->n_regions; i++) {
    khint_t k = kh_put(seq_str_h, b->seq_hash, b->seq_names[i], &absent);
    if (absent == -1) badexit("Error: BED hash put failed.");
    if (absent == 0) {
      b->seq_next[i] = kh_val(b->seq_hash, k);
      kh_val(b->seq_hash, k) = i;
    } else {
      kh_val(b->seq_hash, k) = i;
    }
  }
}

static uint64_t enr_bed_collect(const enr_bed_t *b, const char *seq_name, uint64_t *out) {
  if (!b->seq_hash) return 0;
  khint_t k = kh_get(seq_str_h, b->seq_hash, seq_name);
  if (k == kh_end(b->seq_hash)) return 0;
  uint64_t n = 0;
  uint64_t idx = kh_val(b->seq_hash, k);
  while (idx != UINT64_MAX) {
    out[n++] = idx;
    idx = b->seq_next[idx];
  }
  return n;
}

/* ---- Motif structs ---- */

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
  char        consensus[MAX_MOTIF_SIZE/5 + 1];
  double     *tmp_pdf;
} motif_t;

static motif_t **motifs = NULL;

typedef struct motif_info_t {
  int      fmt;
  uint64_t n;
  uint64_t n_alloc;
} motif_info_t;

static motif_info_t motif_info = { .fmt=0, .n=0, .n_alloc=0 };

static void free_motifs(void) {
  if (!motifs) return;
  for (uint64_t i = 0; i < motif_info.n; i++) free(motifs[i]);
  free(motifs);
  motifs = NULL;
}

/* ---- CDF ---- */

static uint64_t  *cdf_real_size = NULL;
static double   **cdf           = NULL;
static double   **tmp_pdf       = NULL;

static int alloc_cdf(void) {
  cdf_real_size = malloc(sizeof(uint64_t) * args.nthreads);
  if (!cdf_real_size) { fprintf(stderr, "Error: alloc_cdf failed."); return 1; }
  for (uint64_t i = 0; i < (uint64_t)args.nthreads; i++) cdf_real_size[i] = 1;
  cdf = malloc(sizeof(double *) * args.nthreads);
  if (!cdf) { fprintf(stderr, "Error: alloc_cdf failed."); return 1; }
  tmp_pdf = malloc(sizeof(double *) * args.nthreads);
  if (!tmp_pdf) { fprintf(stderr, "Error: alloc_cdf failed."); return 1; }
  for (uint64_t i = 0; i < (uint64_t)args.nthreads; i++) {
    cdf[i] = malloc(sizeof(double));
    if (!cdf[i]) { fprintf(stderr, "Error: alloc_cdf failed."); return 1; }
    tmp_pdf[i] = malloc(sizeof(double));
    if (!tmp_pdf[i]) { fprintf(stderr, "Error: alloc_cdf failed."); return 1; }
  }
  return 0;
}

static void free_cdf(void) {
  if (!cdf) return;
  for (uint64_t i = 0; i < (uint64_t)args.nthreads; i++) {
    free(cdf[i]); free(tmp_pdf[i]);
  }
  free(cdf); cdf = NULL;
  free(tmp_pdf); tmp_pdf = NULL;
  free(cdf_real_size); cdf_real_size = NULL;
}

/* ---- Result arrays ---- */

static uint64_t *seq_hits_pos  = NULL;
static uint64_t *seq_hits_neg  = NULL;
static uint64_t *site_hits_pos = NULL;
static uint64_t *site_hits_neg = NULL;
static int      *max_scores_arr = NULL;  /* [n_motifs * (n_pos + n_neg)], ranksum only */

/* Unit counts: number of independent observations for the test. Equals
   pos_set.n / neg_set.n when no BED is in use, else the count of BED
   regions in bed_pos / bed_neg. Set once after sequences are loaded. */
static uint64_t pos_unit_n = 0;
static uint64_t neg_unit_n = 0;

static pthread_t *threads = NULL;

/* ---- Lookup tables (from yamscan.c) ---- */

static uint64_t char_counts[256];

static const unsigned char char2index[256] = {
  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4, 4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4, 4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
  4,0,4,1,4,4,4,2,4,4,4,4,4,4,4,4, 4,4,4,4,3,3,4,4,4,4,4,4,4,4,4,4,  /* A,C,G,T,U */
  4,0,4,1,4,4,4,2,4,4,4,4,4,4,4,4, 4,4,4,4,3,3,4,4,4,4,4,4,4,4,4,4,  /* a,c,g,t,u */
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

/* ---- Helpers ---- */

static void badexit(const char *msg) {
  fprintf(stderr, "%s\nRun yamtk enr -h for usage.\n", msg);
  free_motifs();
  free_seq_set(&pos_set);
  free_seq_set(&neg_set);
  enr_bed_free(&bed_pos);
  enr_bed_free(&bed_neg);
  free_cdf();
  free(threads);
  free(seq_hits_pos); free(seq_hits_neg);
  free(site_hits_pos); free(site_hits_neg);
  free(max_scores_arr);
  if (seq_hash_tab) kh_destroy(seq_str_h, seq_hash_tab);
  close_files();
  exit(EXIT_FAILURE);
}

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

static int check_line_contains(const char *line, const char *sub) {
  const uint64_t n = strlen(sub);
  if (strlen(line) < n) return 0;
  for (uint64_t i = 0; i < n; i++) if (line[i] != sub[i]) return 0;
  return 1;
}

static uint64_t count_nonempty_chars(const char *line) {
  uint64_t n = 0, i = 0;
  for (;;) {
    switch (line[i]) {
      case ' ': case '\t': case '\r': case '\v': case '\f': case '\n': break;
      case '\0': return n;
      default: n++;
    }
    i++;
  }
  return n;
}

static int check_char_is_one_of(const char c, const char *list) {
  for (uint64_t i = 0; i < strlen(list); i++) if (list[i] == c) return 1;
  return 0;
}

/* True iff the first non-whitespace character looks like the start of a
   numeric PPM value (digit or decimal point). Used to detect end of a
   MEME letter-probability matrix block when a non-data line (e.g. JASPAR's
   `URL ...`) follows. */
static int is_ppm_data_line(const char *line) {
  uint64_t i = 0;
  while (line[i] == ' ' || line[i] == '\t') i++;
  return (line[i] >= '0' && line[i] <= '9') || line[i] == '.';
}

/* ---- Background ---- */

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
  if (args.w) {
    fprintf(stderr, "Using background: A=%.3g C=%.3g G=%.3g T=%.3g\n",
      args.bkg[0], args.bkg[1], args.bkg[2], args.bkg[3]);
  }
}

/* ---- Motif PWM helpers ---- */

static void init_motif(motif_t *motif) {
  ERASE_ARRAY(motif->name, MAX_NAME_SIZE);
  motif->name[0]='m'; motif->name[1]='o'; motif->name[2]='t';
  motif->name[3]='i'; motif->name[4]='f'; motif->name[5]='\0';
  for (uint64_t i = 0; i < sizeof(motif->consensus); i++) motif->consensus[i] = '\0';
  motif->size=0; motif->threshold=0; motif->max_score=0;
  motif->min=0; motif->max=0; motif->cdf_size=0; motif->file_line_num=0;
  motif->min_score=0; motif->cdf_max=0; motif->thread=0;
  for (uint64_t i = 0; i < MAX_MOTIF_SIZE; i++) { motif->pwm[i]=0; motif->pwm_rc[i]=0; }
  for (uint64_t i = 4; i < MAX_MOTIF_SIZE; i += 5) {
    motif->pwm[i]=AMBIGUITY_SCORE; motif->pwm_rc[i]=AMBIGUITY_SCORE;
  }
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

static int calc_score(const double prob_i, const double bkg_i) {
  double x = prob_i * args.nsites;
  x += ((double)args.pseudocount) / 4.0;
  x /= (double)(args.nsites + args.pseudocount);
  return (int)(log2(x / bkg_i) * PWM_INT_MULTIPLIER);
}

static int normalize_probs(double *probs, const char *name) {
  double sum = probs[0]+probs[1]+probs[2]+probs[3];
  if (fabs(sum-1.0) > 0.1) {
    fprintf(stderr, "Error: Position for [%s] does not add up to 1 (sum=%.3g)", name, sum);
    return 1;
  }
  if (fabs(sum-1.0) > 0.02) VEC_DIV(probs, sum, 4);
  return 0;
}

static int get_line_probs(const motif_t *motif, const char *line, double *probs, const uint64_t n) {
  uint64_t i=0, j=0, which_i=(uint64_t)-1;
  int prev_space=1;
  char pos_i[MOTIF_VALUE_MAX_CHAR];
  ERASE_ARRAY(pos_i, MOTIF_VALUE_MAX_CHAR);
  for (; line[i]==' '||line[i]=='\t'; i++) {}
  while (line[i]!='\0'&&line[i]!='\r'&&line[i]!='\n') {
    if (line[i]!=' '&&line[i]!='\t') {
      pos_i[j++]=line[i]; prev_space=0;
    } else {
      if (!prev_space) {
        which_i++;
        if (which_i > n-1) {
          fprintf(stderr, "Error: Motif [%s] has too many columns.", motif->name); return 1;
        }
        if (str_to_double(pos_i, &probs[which_i])) {
          fprintf(stderr, "Error: Failed to parse probability for [%s]: %s", motif->name, pos_i);
          return 1;
        }
        ERASE_ARRAY(pos_i, MOTIF_VALUE_MAX_CHAR); j=0;
      }
      prev_space=1;
    }
    i++;
  }
  if (!prev_space) {
    which_i++;
    if (which_i > n-1) {
      fprintf(stderr, "Error: Motif [%s] has too many columns.", motif->name); return 1;
    }
    if (str_to_double(pos_i, &probs[which_i])) {
      fprintf(stderr, "Error: Failed to parse probability for [%s]: %s", motif->name, pos_i);
      return 1;
    }
  }
  if (which_i == (uint64_t)-1) {
    fprintf(stderr, "Error: Motif [%s] has an empty row.", motif->name); return 1;
  }
  if (which_i < n-1) {
    fprintf(stderr, "Error: Motif [%s] has too few columns.", motif->name); return 1;
  }
  return 0;
}

/* IUPAC consensus from a single column's ACGT probability vector.
   Sort bases by descending probability and emit the smallest IUPAC code
   whose cumulative probability reaches CONSENSUS_THRESHOLD. */
#define CONSENSUS_THRESHOLD 0.75
static char iupac_from_mask(int mask) {
  switch (mask) {
    case 0x1: return 'A'; case 0x2: return 'C';
    case 0x4: return 'G'; case 0x8: return 'T';
    case 0x3: return 'M'; case 0x5: return 'R';
    case 0x9: return 'W'; case 0x6: return 'S';
    case 0xA: return 'Y'; case 0xC: return 'K';
    case 0x7: return 'V'; case 0xB: return 'H';
    case 0xD: return 'D'; case 0xE: return 'B';
    default:  return 'N';
  }
}
static char consensus_char(const double probs[4]) {
  int order[4] = {0, 1, 2, 3};
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

static int add_motif_ppm_column(motif_t *motif, const char *line, const uint64_t pos) {
  double probs[4] = {-1.0,-1.0,-1.0,-1.0};
  if (get_line_probs(motif, line, probs, 4)) return 1;
  if (normalize_probs(probs, motif->name)) return 1;
  if (pos < sizeof(motif->consensus) - 1) motif->consensus[pos] = consensus_char(probs);
  set_score(motif, 'A', pos, calc_score(probs[0], args.bkg[0]));
  set_score(motif, 'C', pos, calc_score(probs[1], args.bkg[1]));
  set_score(motif, 'G', pos, calc_score(probs[2], args.bkg[2]));
  set_score(motif, 'T', pos, calc_score(probs[3], args.bkg[3]));
  return 0;
}

static int add_motif(void) {
  motif_info.n++;
  const uint64_t last_i = motif_info.n - 1;
  if (motif_info.n > motif_info.n_alloc) {
    motif_t **tmp = realloc(motifs,
      sizeof(*motifs)*motif_info.n_alloc + sizeof(*motifs)*ALLOC_CHUNK_SIZE);
    if (!tmp) { fprintf(stderr, "Error: Failed to allocate memory for motifs."); return 1; }
    motifs = tmp;
    motif_info.n_alloc += ALLOC_CHUNK_SIZE;
  }
  motifs[last_i] = malloc(sizeof(motif_t));
  if (!motifs[last_i]) { fprintf(stderr, "Error: Failed to allocate motif."); return 1; }
  init_motif(motifs[last_i]);
  return 0;
}

/* ---- CDF computation (from yamscan.c) ---- */

static void fill_cdf(motif_t *motif) {
  uint64_t max_step, s;
  double pdf_sum = 0.0;
  if (args.w && args.nthreads == 1)
    fprintf(stderr, "        Generating CDF for [%s] (n=%'" PRIu64 ") ... ",
      motif->name, motif->cdf_size);
  if (motif->cdf_size > MAX_CDF_SIZE) {
    fprintf(stderr,
      "\nError: CDF size for [%s] too large (%'" PRIu64 ">%'" PRIu64 ").\n",
      motif->name, motif->cdf_size, MAX_CDF_SIZE);
    badexit("");
  }
  if (cdf_real_size[motif->thread] < motif->cdf_size) {
    double *r1 = realloc(cdf[motif->thread], motif->cdf_size*sizeof(double));
    if (!r1) badexit("Error: Memory re-allocation for CDF failed.");
    cdf[motif->thread] = r1;
    double *r2 = realloc(tmp_pdf[motif->thread], motif->cdf_size*sizeof(double));
    if (!r2) badexit("Error: Memory re-allocation for PDF failed.");
    tmp_pdf[motif->thread] = r2;
    cdf_real_size[motif->thread] = motif->cdf_size;
  }
  motif->cdf = cdf[motif->thread];
  motif->tmp_pdf = tmp_pdf[motif->thread];
  for (uint64_t i = 0; i < motif->cdf_size; i++) motif->cdf[i] = 1.0;
  for (uint64_t i = 0; i < motif->size; i++) {
    max_step = i * motif->cdf_max;
    for (uint64_t j = 0; j < motif->cdf_size; j++) motif->tmp_pdf[j] = motif->cdf[j];
    ERASE_ARRAY(motif->cdf, max_step + motif->cdf_max + 1);
    for (int j = 0; j < 4; j++) {
      s = get_score_i(motif, j, i) - motif->min;
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
  if (args.w && args.nthreads == 1) fprintf(stderr, "done.\n");
}

static inline double score2pval(const motif_t *motif, const int score) {
  return motif->cdf[score - motif->cdf_offset];
}

static void set_threshold(motif_t *motif) {
  uint64_t threshold_i = motif->cdf_size;
  for (uint64_t i = 0; i < motif->cdf_size; i++) {
    if (motif->cdf[i] < args.pvalue) { threshold_i = i; break; }
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
  if (min_pval / args.pvalue > 1.0001) {
    if (args.w)
      fprintf(stderr, "Warning: Min p-value for [%s] exceeds threshold, motif won't be scored.\n",
        motif->name);
    motif->threshold = INT_MAX;
  }
}

/* ---- Motif format detection ---- */

static int detect_motif_fmt(void) {
  int jh=0, fmt=0, tabs=0;
  char *line = NULL; size_t len=0; ssize_t read;
  while ((read = getline(&line, &len, files.m)) != -1) {
    if (!count_nonempty_chars(line)) continue;
    if (check_line_contains(line, "MEME version \0")) {
      if (args.w) fprintf(stderr, "Detected MEME format.\n");
      fmt = FMT_MEME; break;
    }
    if (jh) {
      if (line[0]=='A' && check_char_is_one_of('[',line) && check_char_is_one_of(']',line)) {
        fmt = FMT_JASPAR;
        if (args.w) fprintf(stderr, "Detected JASPAR format.\n");
        break;
      } else if (tabs) {
        fmt = FMT_HOMER;
        if (args.w) fprintf(stderr, "Detected HOMER format.\n");
        break;
      } else {
        if (check_char_is_one_of('-', line)) badexit("Error: yamenr cannot read HOCOMOCO PWMs.");
        fmt = FMT_HOCOMOCO;
        if (args.w) fprintf(stderr, "Detected HOCOMOCO format.\n");
        break;
      }
    } else if (line[0] == '>') {
      if (check_char_is_one_of('\t', line)) tabs = 1;
      jh = 1;
    }
  }
  rewind(files.m);
  free(line);
  return fmt ? fmt : FMT_UNKNOWN;
}

/* ---- MEME reader ---- */

static int check_meme_alph(const char *line, const uint64_t ln) {
  if (check_line_contains(line, "ALPHABET= ACDEFGHIKLMNPQRSTVWY\0")) {
    fprintf(stderr, "Error: Detected protein alphabet (L%" PRIu64 ").", ln); return 1;
  }
  return 0;
}

static int check_meme_strand(const char *line, const uint64_t ln) {
  uint64_t fwd=0, rev=0, i=0;
  for (; line[i]; i++) { if (line[i]=='+') fwd++; if (line[i]=='-') rev++; }
  if (args.scan_rc && fwd && !rev && args.v)
    fprintf(stderr, "Warning: MEME motifs are forward-only (L%" PRIu64 ").\n", ln);
  if (!args.scan_rc && fwd && rev && args.v)
    fprintf(stderr, "Warning: MEME motifs are both-strand (L%" PRIu64 ").\n", ln);
  return 0;
}

static int get_meme_bkg(const char *line, const uint64_t ln) {
  if (args.use_user_bkg) return 0;
  double bp[] = {-1.0,-1.0,-1.0,-1.0};
  uint64_t i=1, li=0, j=0, empty=0;
  char bc[MEME_BKG_MAX_SIZE]; ERASE_ARRAY(bc, MEME_BKG_MAX_SIZE);
  if (line[0] != 'A') {
    fprintf(stderr, "Error: Expected 'A' at start of background line (L%" PRIu64 ").", ln);
    return 1;
  }
  while (line[i]!='\0'&&line[i]!='\n'&&line[i]!='\r') {
    if (line[i]!=' '&&line[i]!='\t') {
      if      (line[i]=='C') { if (!empty||li!=0){fprintf(stderr,"Error: MEME bkg parse error (L%" PRIu64 ").",ln);return 1;} if(str_to_double(bc,&bp[li])){fprintf(stderr,"Error: Bad bkg value.");return 1;} ERASE_ARRAY(bc,MEME_BKG_MAX_SIZE);li=1;j=0; }
      else if (line[i]=='G') { if (!empty||li!=1){fprintf(stderr,"Error: MEME bkg parse error (L%" PRIu64 ").",ln);return 1;} if(str_to_double(bc,&bp[li])){fprintf(stderr,"Error: Bad bkg value.");return 1;} ERASE_ARRAY(bc,MEME_BKG_MAX_SIZE);li=2;j=0; }
      else if (line[i]=='T'||line[i]=='U') { if (!empty||li!=2){fprintf(stderr,"Error: MEME bkg parse error (L%" PRIu64 ").",ln);return 1;} if(str_to_double(bc,&bp[li])){fprintf(stderr,"Error: Bad bkg value.");return 1;} ERASE_ARRAY(bc,MEME_BKG_MAX_SIZE);li=3;j=0; }
      else if (check_char_is_one_of(line[i],"0123456789.\0")) { bc[j++]=line[i]; }
      else { fprintf(stderr,"Error: Unexpected char in MEME background (L%" PRIu64 ").",ln);return 1; }
      empty=0;
    } else { empty=1; }
    i++;
  }
  if (bc[0]!='\0') { if(str_to_double(bc,&bp[li])){fprintf(stderr,"Error: Bad bkg value.");return 1;} }
  if (check_and_load_bkg(bp)) return 1;
  if (args.w)
    fprintf(stderr, "MEME background: A=%.3g C=%.3g G=%.3g T=%.3g\n",
      args.bkg[0],args.bkg[1],args.bkg[2],args.bkg[3]);
  return 0;
}

static void parse_meme_name(const char *line, const uint64_t mi) {
  /* Capture identifier AND any altname; trim_motif_name (run when args.trim_names)
     truncates at the first space later, so default behavior is unchanged. */
  uint64_t i=5, j=0;
  int prev_space=1;
  while (line[i]&&line[i]!='\r'&&line[i]!='\n'&&j<MAX_NAME_SIZE-1) {
    if (line[i]==' '||line[i]=='\t') {
      if (!prev_space && j>0 && j<MAX_NAME_SIZE-1) motifs[mi]->name[j++]=' ';
      prev_space=1;
    } else {
      motifs[mi]->name[j++]=line[i];
      prev_space=0;
    }
    i++;
  }
  if (j>0 && motifs[mi]->name[j-1]==' ') j--;
  motifs[mi]->name[j]='\0';
  if (args.w) fprintf(stderr, "    Found motif: %s (size=", motifs[mi]->name);
}

static void read_meme(void) {
  motif_info.fmt = FMT_MEME;
  char *line=NULL; size_t len=0; ssize_t read;
  uint64_t ln=0, lpm=0, bkg_L=0, mi=(uint64_t)-1, pi=(uint64_t)-1;
  int alph=0, strand=0, live=0;
  int warned_bkg=0, warned_alph=0, warned_strand=0;
  while ((read = getline(&line,&len,files.m)) != -1) {
    ln++;
    if (check_line_contains(line,"Background letter frequencies\0")) {
      /* Only honor the first; concatenated MEME files repeat the header per chunk. */
      if (!bkg_L && mi==(uint64_t)-1) bkg_L=ln;
      else if (!warned_bkg) {
        fprintf(stderr,
          "Warning: Multiple 'Background letter frequencies' lines in MEME file "
          "(L%" PRIu64 "); using the first.\n", ln);
        warned_bkg = 1;
      }
    } else if (bkg_L && bkg_L==ln-1) {
      if (get_meme_bkg(line,ln)) { free(line); badexit(""); }
    } else if (check_line_contains(line,"ALPHABET\0")) {
      if (!alph && mi==(uint64_t)-1) { if(check_meme_alph(line,ln)){free(line);badexit("");} alph=1; }
      else if (!warned_alph) {
        fprintf(stderr,
          "Warning: Multiple ALPHABET lines in MEME file (L%" PRIu64 "); using the first.\n", ln);
        warned_alph = 1;
      }
    } else if (check_line_contains(line,"strands:\0")) {
      if (!strand && mi==(uint64_t)-1) { check_meme_strand(line,ln); strand=1; }
      else if (!warned_strand) {
        fprintf(stderr,
          "Warning: Multiple 'strands:' lines in MEME file (L%" PRIu64 "); using the first.\n", ln);
        warned_strand = 1;
      }
    } else if (check_line_contains(line,"MOTIF\0")) {
      if (mi!=(uint64_t)-1&&args.w) fprintf(stderr,"%" PRIu64 ")\n",motifs[mi]->size);
      mi++;
      if (add_motif()) { free(line); badexit(""); }
      motifs[mi]->file_line_num=ln;
      parse_meme_name(line,mi);
      pi=0;
    } else if (check_line_contains(line,"letter-probability matrix\0")) {
      if (pi==0) { lpm=ln; live=1; }
    } else if (live) {
      if (!count_nonempty_chars(line)||!is_ppm_data_line(line)) {
        live=0;
      } else if (ln==(lpm+pi+1)) {
        if (pi>=MAX_MOTIF_SIZE/5&&pi<(uint64_t)-1) {
          fprintf(stderr,"Error: Motif [%s] too large.",motifs[mi]->name);
          free(line); badexit("");
        }
        if (add_motif_ppm_column(motifs[mi],line,pi)) { free(line); badexit(""); }
        pi++;
        motifs[mi]->size=pi;
      } else { live=0; }
    }
  }
  free(line);
  if (mi!=(uint64_t)-1&&args.w) fprintf(stderr,"%" PRIu64 ")\n",motifs[mi]->size);
  if (!motif_info.n) badexit("Error: No motifs found in MEME file.");
  if (args.v) fprintf(stderr,"Found %'" PRIu64 " MEME motif(s).\n",motif_info.n);
}

/* ---- HOMER reader ---- */

static void parse_homer_name(const char *line, const uint64_t mi) {
  uint64_t ns=0,ne=0,i=1,in_b=0,j=0;
  while (line[i]&&line[i]!='\r'&&line[i]!='\n') {
    if (line[i]=='\t') { if(ns){ne=i;break;}else{in_b=1;} }
    else if (in_b && !ns) ns=i;
    i++;
  }
  if (!ns) { if(args.w) fprintf(stderr,"Warning: Failed to parse HOMER motif name.\n"); }
  else if (!ne) ne=i;
  for (uint64_t k=ns; k<ne&&j<MAX_NAME_SIZE-1; k++) motifs[mi]->name[j++]=line[k];
  motifs[mi]->name[j]='\0';
  if (args.w) fprintf(stderr,"    Found motif: %s (size=",motifs[mi]->name);
}

static void read_homer(void) {
  motif_info.fmt = FMT_HOMER;
  char *line=NULL; size_t len=0; ssize_t read;
  uint64_t ln=0, mi=(uint64_t)-1, pi; int ready=0;
  while ((read=getline(&line,&len,files.m))!=-1) {
    ln++;
    if (line[0]=='>') {
      ready=1;
      if (mi!=(uint64_t)-1&&args.w) fprintf(stderr,"%" PRIu64 ")\n",motifs[mi]->size);
      mi++;
      if (add_motif()) { free(line); badexit(""); }
      motifs[mi]->file_line_num=ln;
      parse_homer_name(line,mi); pi=0;
    } else if (count_nonempty_chars(line)&&ready) {
      if (pi>MAX_MOTIF_SIZE/5&&pi<(uint64_t)-1) {
        fprintf(stderr,"Error: Motif [%s] too large.",motifs[mi]->name);
        free(line); badexit("");
      }
      if (add_motif_ppm_column(motifs[mi],line,pi)) { free(line); badexit(""); }
      pi++; motifs[mi]->size=pi;
    }
  }
  free(line);
  if (mi!=(uint64_t)-1&&args.w) fprintf(stderr,"%" PRIu64 ")\n",motifs[mi]->size);
  if (args.v) fprintf(stderr,"Found %'" PRIu64 " HOMER motif(s).\n",motif_info.n);
}

/* ---- JASPAR reader ---- */

static void parse_jaspar_name(const char *line, const uint64_t mi) {
  uint64_t i=0, j=1;
  for (;line[j]&&line[j]!='\r'&&line[j]!='\n'&&i<MAX_NAME_SIZE-1; j++) motifs[mi]->name[i++]=line[j];
  motifs[mi]->name[i]='\0';
  if (args.w) fprintf(stderr,"    Found motif: %s (size=",motifs[mi]->name);
}

static int add_jaspar_row(motif_t *motif, const char *line) {
  uint64_t ri=(uint64_t)-1, lb=(uint64_t)-1, rb=(uint64_t)-1, i=0;
  char let='N';
  for (; line[i]&&line[i]!='\r'&&line[i]!='\n'; i++) {
    switch(line[i]) {
      case 'A': case 'a': ri=0; let='A'; break;
      case 'C': case 'c': ri=1; let='C'; break;
      case 'G': case 'g': ri=2; let='G'; break;
      case 'T': case 't': case 'U': case 'u': ri=3; let='T'; break;
      case '[': lb=i; break;
      case ']': rb=i; break;
    }
  }
  if (ri==(uint64_t)-1) { fprintf(stderr,"Error: No ACGTU in JASPAR row for [%s].",motif->name); return 1; }
  if (lb==(uint64_t)-1||rb==(uint64_t)-1) { fprintf(stderr,"Error: No [] in JASPAR row for [%s].",motif->name); return 1; }
  uint64_t k=0, pi=(uint64_t)-1; int prev_sp=1, tmp=0;
  char pc[MOTIF_VALUE_MAX_CHAR]; ERASE_ARRAY(pc, MOTIF_VALUE_MAX_CHAR);
  for (i=lb+1; line[i]==' '||line[i]=='\t'; i++) {}
  for (uint64_t j2=i; j2<rb; j2++) {
    if (line[j2]!=' '&&line[j2]!='\t') { pc[k++]=line[j2]; prev_sp=0; }
    else {
      if (!prev_sp) {
        pi++;
        if (pi+1>MAX_MOTIF_SIZE/5&&pi<(uint64_t)-1) {
          fprintf(stderr,"Error: Motif [%s] too many columns.",motif->name); return 1;
        }
        if (str_to_int(pc,&tmp)) {
          fprintf(stderr,"Error: Bad count in JASPAR motif [%s]: %s",motif->name,pc); return 1;
        }
        set_score(motif,let,pi,tmp);
        ERASE_ARRAY(pc,MOTIF_VALUE_MAX_CHAR); k=0;
      }
      prev_sp=1;
    }
  }
  if (!prev_sp) {
    pi++;
    if (str_to_int(pc,&tmp)) {
      fprintf(stderr,"Error: Bad count in JASPAR motif [%s]: %s",motif->name,pc); return 1;
    }
    set_score(motif,let,pi,tmp);
  }
  if (pi==(uint64_t)-1) { fprintf(stderr,"Error: Empty JASPAR row for [%s].",motif->name); return 1; }
  pi++;
  if (motif->size) { if (motif->size!=pi) { fprintf(stderr,"Error: JASPAR rows differ in [%s].",motif->name); return 1; } }
  else { motif->size=pi; }
  return 0;
}

static void pcm_to_pwm(motif_t *motif) {
  int nsites=0, nsites2;
  for (int i=0;i<4;i++) nsites+=get_score_i(motif,i,0);
  for (uint64_t j=0;j<motif->size;j++) {
    nsites2=0;
    for (int i=0;i<4;i++) nsites2+=get_score_i(motif,i,j);
    if (abs(nsites2-nsites)>1) {
      fprintf(stderr,"Error: Column sums differ for motif [%s].",motif->name); badexit("");
    }
  }
  /* Compute consensus from raw counts before they're overwritten with scores. */
  for (uint64_t j = 0; j < motif->size && j < sizeof(motif->consensus) - 1; j++) {
    double probs[4];
    double col = 0;
    for (int i = 0; i < 4; i++) col += (double)get_score_i(motif, i, j);
    if (col <= 0) col = 1;
    for (int i = 0; i < 4; i++) probs[i] = (double)get_score_i(motif, i, j) / col;
    motif->consensus[j] = consensus_char(probs);
  }
  const char lets[4]={'A','C','G','T'};
  for (uint64_t j=0;j<motif->size;j++)
    for (int i=0;i<4;i++)
      set_score(motif,lets[i],j,
        calc_score((args.pseudocount/4.0+((double)get_score_i(motif,i,j)))/(args.pseudocount+((double)nsites)),args.bkg[i]));
}

static void read_jaspar(void) {
  motif_info.fmt = FMT_JASPAR;
  char *line=NULL; size_t len=0; ssize_t read;
  uint64_t ln=0, mi=(uint64_t)-1, ri=(uint64_t)-1; int ready=0;
  while ((read=getline(&line,&len,files.m))!=-1) {
    ln++;
    if (line[0]=='>') {
      ready=1;
      if (mi!=(uint64_t)-1) {
        if (args.w) fprintf(stderr,"%" PRIu64 ")\n",motifs[mi]->size);
        if (ri!=4) { fprintf(stderr,"Error: JASPAR motif [%s] has wrong row count.",motifs[mi]->name); free(line); badexit(""); }
      }
      mi++;
      if (add_motif()) { free(line); badexit(""); }
      motifs[mi]->file_line_num=ln;
      parse_jaspar_name(line,mi); ri=0;
    } else if (count_nonempty_chars(line)&&ready) {
      ri++;
      if (add_jaspar_row(motifs[mi],line)) { free(line); badexit(""); }
    }
  }
  free(line);
  if (mi!=(uint64_t)-1 && ri!=4) {
    fprintf(stderr,"Error: JASPAR motif [%s] has wrong row count.",motifs[mi]->name); badexit("");
  }
  if (mi!=(uint64_t)-1&&args.w) fprintf(stderr,"%" PRIu64 ")\n",motifs[mi]->size);
  for (uint64_t i=0;i<motif_info.n;i++) pcm_to_pwm(motifs[i]);
  if (args.v) fprintf(stderr,"Found %'" PRIu64 " JASPAR motif(s).\n",motif_info.n);
}

/* ---- HOCOMOCO reader ---- */

static int add_motif_pcm_column(motif_t *motif, const char *line, const uint64_t pos) {
  double probs[4]={-1.0,-1.0,-1.0,-1.0};
  if (get_line_probs(motif,line,probs,4)) return 1;
  double pcm_sum=probs[0]+probs[1]+probs[2]+probs[3];
  if (pcm_sum<0.99) { fprintf(stderr,"Error: PCM row sums < 1 for [%s].",motif->name); return 1; }
  if (pos < sizeof(motif->consensus) - 1) {
    double cprobs[4] = { probs[0]/pcm_sum, probs[1]/pcm_sum, probs[2]/pcm_sum, probs[3]/pcm_sum };
    motif->consensus[pos] = consensus_char(cprobs);
  }
  VEC_ADD(probs, args.pseudocount/4.0, 4);
  set_score(motif,'A',pos,calc_score(probs[0]/pcm_sum,args.bkg[0]));
  set_score(motif,'C',pos,calc_score(probs[1]/pcm_sum,args.bkg[1]));
  set_score(motif,'G',pos,calc_score(probs[2]/pcm_sum,args.bkg[2]));
  set_score(motif,'T',pos,calc_score(probs[3]/pcm_sum,args.bkg[3]));
  return 0;
}

static void read_hocomoco(void) {
  motif_info.fmt = FMT_HOCOMOCO;
  char *line=NULL; size_t len=0; ssize_t read;
  uint64_t ln=0, mi=(uint64_t)-1, pi; int ready=0;
  while ((read=getline(&line,&len,files.m))!=-1) {
    ln++;
    if (line[0]=='>') {
      ready=1;
      if (mi!=(uint64_t)-1&&args.w) fprintf(stderr,"%" PRIu64 ")\n",motifs[mi]->size);
      mi++;
      if (add_motif()) { free(line); badexit(""); }
      motifs[mi]->file_line_num=ln;
      for (uint64_t i=1, j=0; i<MAX_NAME_SIZE; i++) {
        if (line[i]=='\r'||line[i]=='\n'||line[i]=='\0') { motifs[mi]->name[j]='\0'; break; }
        motifs[mi]->name[j++]=line[i];
        if (j==MAX_NAME_SIZE-1) { motifs[mi]->name[j]='\0'; break; }
      }
      if (args.w) fprintf(stderr,"    Found motif: %s (size=",motifs[mi]->name);
      pi=0;
    } else if (count_nonempty_chars(line)&&ready) {
      if (pi>MAX_MOTIF_SIZE/5&&pi<(uint64_t)-1) {
        fprintf(stderr,"Error: Motif [%s] too large.",motifs[mi]->name);
        free(line); badexit("");
      }
      if (add_motif_pcm_column(motifs[mi],line,pi)) { free(line); badexit(""); }
      pi++; motifs[mi]->size=pi;
    }
  }
  free(line);
  if (mi!=(uint64_t)-1&&args.w) fprintf(stderr,"%" PRIu64 ")\n",motifs[mi]->size);
  if (args.v) fprintf(stderr,"Found %'" PRIu64 " HOCOMOCO motif(s).\n",motif_info.n);
}

/* ---- complete_motifs + load_motifs ---- */

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

static void trim_motif_name(motif_t *m) {
  for (uint64_t i=0; i<MAX_NAME_SIZE; i++) {
    if (m->name[i]==' '||m->name[i]=='\t'||m->name[i]=='\0') { m->name[i]='\0'; break; }
  }
}

static void complete_motifs(void) {
  for (uint64_t i=0; i<motif_info.n; i++) {
    motifs[i]->min = get_pwm_min(motifs[i]);
    motifs[i]->max = get_pwm_max(motifs[i]);
    motifs[i]->cdf_offset = motifs[i]->min * (int)motifs[i]->size;
    fill_pwm_rc(motifs[i]);
    motifs[i]->cdf_max = motifs[i]->max - motifs[i]->min;
    motifs[i]->cdf_size = motifs[i]->size * motifs[i]->cdf_max + 1;
    if (args.trim_names) trim_motif_name(motifs[i]);
  }
}

static void load_motifs(void) {
  switch (detect_motif_fmt()) {
    case FMT_MEME:     read_meme();     break;
    case FMT_HOMER:    read_homer();    break;
    case FMT_JASPAR:   read_jaspar();   break;
    case FMT_HOCOMOCO: read_hocomoco(); break;
    default: badexit("Error: Failed to detect motif format.");
  }
  if (motif_info.n > 100000)
    fprintf(stderr, "Warning: Very large motif count may be slow!\n");
  complete_motifs();
  uint64_t empty=0;
  for (uint64_t i=0; i<motif_info.n; i++) if (!motifs[i]->size) empty++;
  if (empty == motif_info.n) badexit("Error: All motifs are empty.");
  else if (empty) fprintf(stderr, "Warning: %'" PRIu64 " empty motif(s).\n", empty);
}

/* ---- Low-memory streaming state ---- */

static uint64_t *stream_pos_positions = NULL;
static uint64_t *stream_neg_positions = NULL;

/* ---- Sequence loading ---- */

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
  /* compute GC and unknowns */
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

/* ---- PRNG (from yamshuf.c) ---- */

typedef struct { uint64_t s[2]; } xrng_t;

static xrng_t xrng;

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

/* ---- Shuffle functions (from yamshuf.c) ---- */

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
  244140625,1220703125,6103515625,30517578125
};

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

/* ---- Build shuffled negative set ---- */

static void make_shuffled_negatives(void) {
  /* When -x is in use, the null set is built from BED region slices so the
     k-mer composition of the null matches what's actually being scored. */
  const int use_bed = args.use_bed_pos;
  const uint64_t n_src = use_bed ? bed_pos.n_regions : pos_set.n;
  if (args.v)
    fprintf(stderr, "Generating shuffled negative set (k=%d, %s) ...\n",
      args.shuffle_k, use_bed ? "BED-restricted positives" : "full positives");
  neg_set.n = n_src;
  neg_set.n_alloc = n_src;
  neg_set.seqs = malloc(sizeof(unsigned char *) * n_src);
  if (!neg_set.seqs) badexit("Error: Failed to allocate shuffled negatives.");
  neg_set.sizes = malloc(sizeof(uint64_t) * n_src);
  if (!neg_set.sizes) badexit("Error: Failed to allocate shuffled negative sizes.");
  neg_set.names = NULL;  /* not needed for scoring */
  /* Seed the PRNG once */
  uint64_t seed = args.use_seed ? args.seed : (uint64_t)time(NULL);
  sxrand_r(&xrng, seed);
  if (args.v)
    fprintf(stderr, "RNG seed: %" PRIu64 "\n", seed);
  const uint64_t ksz = (uint64_t) args.shuffle_k;
  uint64_t      *kmer_tab   = NULL;
  unsigned char *inv_vtx    = NULL;
  uint64_t      *euler_path = NULL;
  uint64_t      *next_idx   = NULL;
  if (ksz > 1) {
    kmer_tab   = calloc(pow5[ksz],   sizeof(uint64_t));
    inv_vtx    = calloc(pow5[ksz-1], sizeof(unsigned char));
    euler_path = calloc(pow5[ksz-1], sizeof(uint64_t));
    next_idx   = calloc(pow5[ksz-1], sizeof(uint64_t));
    if (!kmer_tab||!inv_vtx||!euler_path||!next_idx)
      badexit("Error: Failed to allocate Euler shuffle scratch.");
  }
  uint64_t total_bp = 0;
  for (uint64_t i = 0; i < n_src; i++) {
    /* Resolve the (seq, slice) for this region or the whole sequence. */
    const unsigned char *src;
    uint64_t L;
    if (use_bed) {
      const uint64_t si = bed_pos.seq_idx[i];
      src = pos_set.seqs[si] + bed_pos.starts[i];
      L   = bed_pos.ends[i] - bed_pos.starts[i];
    } else {
      src = pos_set.seqs[i];
      L   = pos_set.sizes[i];
    }
    neg_set.sizes[i] = L;
    neg_set.seqs[i]  = malloc(L ? L : 1);
    if (!neg_set.seqs[i]) badexit("Error: Failed to allocate shuffle buffer.");
    if (L) memcpy(neg_set.seqs[i], src, L);
    if (ksz == 1) {
      if (L > 1) shuffle_fisher_yates(neg_set.seqs[i], L);
    } else if (L >= ksz) {
      memset(kmer_tab,   0, sizeof(uint64_t)      * pow5[ksz]);
      memset(inv_vtx,    0, sizeof(unsigned char) * pow5[ksz-1]);
      memset(euler_path, 0, sizeof(uint64_t)      * pow5[ksz-1]);
      memset(next_idx,   0, sizeof(uint64_t)      * pow5[ksz-1]);
      count_kmers(neg_set.seqs[i], L, kmer_tab, ksz);
      shuffle_euler(neg_set.seqs[i], L, ksz, kmer_tab, inv_vtx, euler_path, next_idx);
    }
    total_bp += L;
  }
  free(kmer_tab); free(inv_vtx); free(euler_path); free(next_idx);
  neg_set.total_bases = total_bp;
  neg_set.gc_pct      = pos_set.gc_pct;     /* approximate; only used in -v messages */
  neg_set.unknowns    = pos_set.unknowns;   /* approximate likewise */
  if (args.v)
    fprintf(stderr, "Negative (shuffled) set: %'" PRIu64 " seq(s), %'" PRIu64 " bp\n",
      neg_set.n, neg_set.total_bases);
}

/* ---- Statistical functions ---- */

static double log_choose(double n, double k) {
  return lgamma(n+1) - lgamma(k+1) - lgamma(n-k+1);
}

/* One-tailed Fisher's exact P(X >= a) for 2x2 [a b / c d]. */
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

typedef struct { int score; int is_pos; } ranked_t;

static int cmp_ranked(const void *a, const void *b) {
  const ranked_t *ra=(const ranked_t *)a, *rb=(const ranked_t *)b;
  return (ra->score>rb->score)-(ra->score<rb->score);
}

/* One-sided Mann-Whitney U (pos > neg); scores layout: [pos | neg]. */
static void mann_whitney_u_greater(const int *scores, int n_pos, int n_neg,
                                    double *out_pval, double *out_auc) {
  int N = n_pos + n_neg;
  ranked_t *arr = malloc(sizeof(ranked_t) * N);
  if (!arr) { *out_pval=1.0; *out_auc=0.5; return; }
  for (int i=0; i<n_pos; i++) { arr[i].score=scores[i];       arr[i].is_pos=1; }
  for (int i=0; i<n_neg; i++) { arr[n_pos+i].score=scores[n_pos+i]; arr[n_pos+i].is_pos=0; }
  qsort(arr, N, sizeof(ranked_t), cmp_ranked);
  /* Assign average ranks to ties; accumulate tie correction. */
  double rank_sum_pos=0.0, tie_correction=0.0;
  int i=0;
  while (i < N) {
    int j=i;
    while (j<N && arr[j].score==arr[i].score) j++;
    double avg_rank = (double)(i+1+j) / 2.0;  /* 1-indexed average */
    uint64_t tie_size = j-i;
    if (tie_size > 1)
      tie_correction += (double)(tie_size*tie_size*tie_size - tie_size);
    for (int k=i; k<j; k++) if (arr[k].is_pos) rank_sum_pos += avg_rank;
    i=j;
  }
  free(arr);
  double U_pos = rank_sum_pos - (double)n_pos*(n_pos+1)/2.0;
  *out_auc = U_pos / ((double)n_pos * n_neg);
  double mu = (double)n_pos*n_neg/2.0;
  double var = ((double)n_pos*n_neg/12.0) *
               ((double)(N+1) - tie_correction/((double)N*((double)N-1)));
  if (var <= 0.0) { *out_pval = (*out_auc >= 0.5) ? 0.0 : 1.0; return; }
  double z = (U_pos - mu) / sqrt(var);
  /* One-sided p: P(Z > z) = 0.5 * erfc(z/sqrt(2)) */
  double p = 0.5 * erfc(z / M_SQRT2);
  *out_pval = (p < 0.0) ? 0.0 : (p > 1.0) ? 1.0 : p;
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

/* ---- Duplicate-name detection (from yamscan.c) ---- */

static void int_to_char_array(const uint64_t N, char *arr) {
  ERASE_ARRAY(arr, 128);
  snprintf(arr, 128, "__N%" PRIu64 "", N);
}

static int dedup_char_array(char *arr, const uint64_t arr_max_len, const uint64_t N) {
  uint64_t arr_len = 0, dedup_len = 0, success = 0, j = 0;
  char dedup[128];
  ERASE_ARRAY(dedup, 128);
  int_to_char_array(N, dedup);
  for (uint64_t i = 0; i < arr_max_len; i++) {
    if (arr[i] == '\0') { arr_len = i; break; }
  }
  for (uint64_t i = 0; i < 128; i++) {
    if (dedup[i] == '\0') { dedup_len = i + 1; break; }
  }
  if (arr_max_len - arr_len >= dedup_len) {
    for (uint64_t i = arr_len; i < arr_len + dedup_len; i++) {
      arr[i] = dedup[j++];
    }
    success = 1;
  }
  return success;
}

static void find_motif_dupes(void) {
  if (motif_info.n == 1) return;
  uint64_t *is_dup = malloc(sizeof(uint64_t) * motif_info.n);
  if (!is_dup) badexit("Error: Failed to allocate memory for motif name duplication check.");
  ERASE_ARRAY(is_dup, motif_info.n);
  khash_t(motif_str_h) *motif_hash_tab = kh_init(motif_str_h);
  khint64_t k;
  int absent;
  for (uint64_t i = 0; i < motif_info.n; i++) {
    k = kh_put(motif_str_h, motif_hash_tab, motifs[i]->name, &absent);
    if (absent == -1) { free(is_dup); badexit("Error: Failed to hash motif names."); }
    else if (absent == 0) is_dup[i] = 1;
  }
  kh_destroy(motif_str_h, motif_hash_tab);
  uint64_t dup_count = 0;
  for (uint64_t i = 0; i < motif_info.n; i++) dup_count += is_dup[i];
  if (dup_count) {
    if (args.dedup) {
      for (uint64_t i = 0; i < motif_info.n; i++) {
        if (is_dup[i]) {
          if (!dedup_char_array(motifs[i]->name, MAX_NAME_SIZE, i + 1)) {
            fprintf(stderr,
              "Error: Failed to deduplicate motif #%" PRIu64 ", name is too large.", i + 1);
            free(is_dup); badexit("");
          }
        }
      }
    } else {
      fprintf(stderr, "Error: Encountered duplicate motif name (use -d to deduplicate).");
      uint64_t to_print = (dup_count < 5) ? dup_count : 5;
      for (uint64_t i = 0; i < motif_info.n; i++) {
        if (is_dup[i]) {
          fprintf(stderr, "\n    L%" PRIu64 " #%" PRIu64 ": %s",
            motifs[i]->file_line_num, i + 1, motifs[i]->name);
          if (!--to_print) break;
        }
      }
      if (dup_count > 5) {
        fprintf(stderr, "\n    ...");
        fprintf(stderr, "\n    Found %'" PRIu64 " total non-unique names.", dup_count);
      }
      free(is_dup); badexit("");
    }
  }
  free(is_dup);
}

/* ---- Scoring scan ---- */

static inline void score_subseq(const motif_t *m, const unsigned char *seq, const uint64_t off, int *s, const unsigned char *tbl) {
  *s=0; for (uint64_t i=0;i<m->size;i++) *s+=get_score(m,seq[i+off],i,tbl);
}
static inline void score_subseq_rc(const motif_t *m, const unsigned char *seq, const uint64_t off, int *s, int *src, const unsigned char *tbl) {
  *s=0; *src=0;
  for (uint64_t i=0;i<m->size;i++) {
    *s   += get_score(m,seq[i+off],i,tbl);
    *src += get_score_rc(m,seq[i+off],i,tbl);
  }
}

/* Count hits and track max PWM score for [start, end) of seq, with explicit
   strand mask. do_fwd / do_rc allow per-region strand control under -x/-X. */
static void score_range_scan(const motif_t *motif, const unsigned char *seq,
                              const uint64_t start, const uint64_t end,
                              const int do_fwd, const int do_rc,
                              uint64_t *out_hits, int *out_max) {
  const uint64_t mot_size = motif->size;
  if (end <= start || (end - start) < mot_size) return;
  /* Use INT_MAX-1 as effective threshold when threshold==INT_MAX (no hits, but max computed). */
  const int threshold = motif->threshold - 1;
  const unsigned char *tbl = args.mask ? char2maskindex : char2index;
  const uint64_t last = end - mot_size;
  int score, score_rc;
  if (do_fwd && do_rc) {
    for (uint64_t i = start; i <= last; i++) {
      score_subseq_rc(motif, seq, i, &score, &score_rc, tbl);
      if (UNLIKELY(score    > threshold)) (*out_hits)++;
      if (UNLIKELY(score_rc > threshold)) (*out_hits)++;
      if (score    > *out_max) *out_max = score;
      if (score_rc > *out_max) *out_max = score_rc;
    }
  } else if (do_fwd) {
    for (uint64_t i = start; i <= last; i++) {
      score_subseq(motif, seq, i, &score, tbl);
      if (UNLIKELY(score > threshold)) (*out_hits)++;
      if (score > *out_max) *out_max = score;
    }
  } else if (do_rc) {
    /* Reverse-only: compute both but ignore forward score. */
    for (uint64_t i = start; i <= last; i++) {
      score_subseq_rc(motif, seq, i, &score, &score_rc, tbl);
      if (UNLIKELY(score_rc > threshold)) (*out_hits)++;
      if (score_rc > *out_max) *out_max = score_rc;
    }
  }
}

/* Count hits and track max PWM score for one (motif, full sequence) pair.
   Wrapper around score_range_scan using args.scan_rc as the strand mask. */
static void score_seq_scan(const motif_t *motif, const unsigned char *seq,
                            const uint64_t seq_size, uint64_t *out_hits, int *out_max) {
  score_range_scan(motif, seq, 0, seq_size, 1, args.scan_rc, out_hits, out_max);
}

/* Decide do_fwd / do_rc for a BED region given its strand char and -R setting. */
static inline void strand_mask_for(const char strand, int *do_fwd, int *do_rc) {
  if (strand == '+')      { *do_fwd = 1;            *do_rc = 0; }
  else if (strand == '-') { *do_fwd = 0;            *do_rc = args.scan_rc ? 1 : 0; }
  else                    { *do_fwd = 1;            *do_rc = args.scan_rc ? 1 : 0; }
  /* If user passed -R and the region is strand '-', we have no work to do
     (no rc, no fwd). Make that explicit: skip such regions in callers via
     do_fwd|do_rc == 0. */
}

/* Per-region scoring position count: (width - mot_size + 1) * strand_factor.
   strand_factor = 2 when both strands active, 1 otherwise. Returns 0 when
   region is shorter than the motif or no strand is active. */
static inline uint64_t region_positions(const uint64_t start, const uint64_t end,
                                        const uint64_t mot_size, const char strand) {
  if (end <= start || (end - start) < mot_size) return 0;
  int do_fwd, do_rc;
  strand_mask_for(strand, &do_fwd, &do_rc);
  const uint64_t s = do_fwd + do_rc;
  if (s == 0) return 0;
  return (end - start - mot_size + 1) * s;
}

/* ---- Low-memory streaming scoring (single-threaded) ---- */

/* Set up CDF + threshold for every motif on thread 0. Called once before
   any streaming scoring. */
static void setup_motifs_low_mem(void) {
  for (uint64_t mi = 0; mi < motif_info.n; mi++) {
    motifs[mi]->thread = 0;
    fill_cdf(motifs[mi]);
    set_threshold(motifs[mi]);
  }
}

/* Score one sequence against every loaded motif, updating per-motif counters.
   `is_pos` selects which accumulators get updated. */
static void score_one_seq_all_motifs(const unsigned char *seq, const uint64_t L,
                                     const int is_pos) {
  const int strands = args.scan_rc ? 2 : 1;
  uint64_t *seq_hits  = is_pos ? seq_hits_pos  : seq_hits_neg;
  uint64_t *site_hits = is_pos ? site_hits_pos : site_hits_neg;
  uint64_t *positions = is_pos ? stream_pos_positions : stream_neg_positions;
  for (uint64_t mi = 0; mi < motif_info.n; mi++) {
    motif_t *m = motifs[mi];
    uint64_t hits = 0; int mx = INT_MIN;
    score_seq_scan(m, seq, L, &hits, &mx);
    site_hits[mi] += hits;
    if (hits > 0) seq_hits[mi]++;
    if (positions && L >= m->size)
      positions[mi] += (L - m->size + 1) * (uint64_t) strands;
  }
}

/* BED-restricted streaming variant: score every region for this sequence
   against every motif, treating each region as one independent unit. */
static void score_one_seq_all_motifs_bed(const unsigned char *seq, const uint64_t L,
                                          const enr_bed_t *bed,
                                          const uint64_t *regions, const uint64_t nr,
                                          const int is_pos) {
  uint64_t *seq_hits  = is_pos ? seq_hits_pos  : seq_hits_neg;
  uint64_t *site_hits = is_pos ? site_hits_pos : site_hits_neg;
  uint64_t *positions = is_pos ? stream_pos_positions : stream_neg_positions;
  (void) L;  /* bounds were validated at bed bind / read time when possible;
                in streaming we accept the BED start/end as-is — out-of-bounds
                is handled by score_range_scan returning early. */
  for (uint64_t ri = 0; ri < nr; ri++) {
    const uint64_t rid = regions[ri];
    int do_fwd, do_rc;
    strand_mask_for(bed->strands[rid], &do_fwd, &do_rc);
    for (uint64_t mi = 0; mi < motif_info.n; mi++) {
      motif_t *m = motifs[mi];
      uint64_t hits = 0; int mx = INT_MIN;
      score_range_scan(m, seq, bed->starts[rid], bed->ends[rid],
                       do_fwd, do_rc, &hits, &mx);
      site_hits[mi] += hits;
      if (hits > 0) seq_hits[mi]++;
      if (positions)
        positions[mi] += region_positions(bed->starts[rid], bed->ends[rid],
                                          m->size, bed->strands[rid]);
    }
  }
}

/* Walk pos / neg FASTA via kseq, calling score_one_seq_all_motifs per record
   and accumulating set-level GC / total_bases / unknowns. When `bed` is
   non-NULL, scoring is restricted to that BED's regions for the streamed
   sequence (one unit per region). */
static void stream_score_set(kseq_t *kseq, seq_set_t *set, const int is_pos,
                             const char *label, const enr_bed_t *bed) {
  ERASE_ARRAY(char_counts, 256);
  set->n = 0;
  set->total_bases = 0;
  set->unknowns = 0;
  set->gc_pct = 0.0;
  uint64_t *region_buf = NULL;
  if (bed) {
    region_buf = malloc(sizeof(*region_buf) * bed->n_regions);
    if (!region_buf) badexit("Error: Failed to alloc BED region scratch.");
  }
  int ret;
  while ((ret = kseq_read(kseq)) >= 0) {
    if (kseq->seq.l == 0) continue;
    set->n++;
    const unsigned char *seq = (const unsigned char *) kseq->seq.s;
    const uint64_t L = kseq->seq.l;
    for (uint64_t j = 0; j < L; j++) char_counts[seq[j]]++;
    if (bed) {
      const uint64_t nr = enr_bed_collect(bed, kseq->name.s ? kseq->name.s : "", region_buf);
      if (nr) score_one_seq_all_motifs_bed(seq, L, bed, region_buf, nr, is_pos);
    } else {
      score_one_seq_all_motifs(seq, L, is_pos);
    }
    if (args.progress && (set->n % 100) == 0) {
      fprintf(stderr, "\r%s seqs streamed: %" PRIu64, label, set->n);
      fflush(stderr);
    }
  }
  free(region_buf);
  if (args.progress) fprintf(stderr, "\r%s seqs streamed: %" PRIu64 "\n", label, set->n);
  if (ret == -2) badexit("Error: Failed to parse FASTQ qualities.");
  if (ret < -2)  badexit("Error: Failed to read input.");
  if (!set->n)   { fprintf(stderr, "Error: No sequences read from %s input.\n", label); badexit(""); }

  uint64_t total = 0;
  for (uint64_t c = 0; c < 256; c++) total += char_counts[c];
  if (!total) { fprintf(stderr, "Error: All %s sequences are empty.\n", label); badexit(""); }
  set->total_bases = total;
  set->unknowns = total - standard_base_count();
  set->gc_pct = calc_gc() * 100.0;
  if (set->unknowns == total) {
    fprintf(stderr, "Error: No standard DNA/RNA bases found in %s set.\n", label);
    badexit("");
  }
  double unk_pct = 100.0 * (double) set->unknowns / (double) total;
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

/* ---- Thread worker ---- */

static void *enr_sub_process(void *arg) {
  uint64_t thread_i = *(uint64_t *)arg;
  free(arg);
  const uint64_t n_total = pos_unit_n + neg_unit_n;
  for (uint64_t mi=0; mi<motif_info.n; mi++) {
    motif_t *motif = motifs[mi];
    if (motif->thread != thread_i) continue;
    if (args.w && !args.progress) fprintf(stderr, "    Scoring motif: %s\n", motif->name);
    fill_cdf(motif);
    set_threshold(motif);
    uint64_t sp=0, sn=0, hp=0, hn=0;
    /* Positives: iterate BED regions when -x, else full sequences. */
    if (args.use_bed_pos) {
      for (uint64_t ri=0; ri<bed_pos.n_regions; ri++) {
        const uint64_t si = bed_pos.seq_idx[ri];
        int do_fwd, do_rc;
        strand_mask_for(bed_pos.strands[ri], &do_fwd, &do_rc);
        uint64_t hits=0; int mx=INT_MIN;
        score_range_scan(motif, pos_set.seqs[si], bed_pos.starts[ri], bed_pos.ends[ri],
                         do_fwd, do_rc, &hits, &mx);
        sp += hits;
        if (hits > 0) hp++;
        if (args.test_mode == TEST_RANKSUM && max_scores_arr)
          max_scores_arr[mi*n_total + ri] = mx;
      }
    } else {
      for (uint64_t si=0; si<pos_set.n; si++) {
        uint64_t hits=0; int mx=INT_MIN;
        score_seq_scan(motif, pos_set.seqs[si], pos_set.sizes[si], &hits, &mx);
        sp += hits;
        if (hits > 0) hp++;
        if (args.test_mode == TEST_RANKSUM && max_scores_arr)
          max_scores_arr[mi*n_total + si] = mx;
      }
    }
    /* Negatives: iterate BED regions when -X, else full sequences. */
    if (args.use_bed_neg) {
      for (uint64_t ri=0; ri<bed_neg.n_regions; ri++) {
        const uint64_t si = bed_neg.seq_idx[ri];
        int do_fwd, do_rc;
        strand_mask_for(bed_neg.strands[ri], &do_fwd, &do_rc);
        uint64_t hits=0; int mx=INT_MIN;
        score_range_scan(motif, neg_set.seqs[si], bed_neg.starts[ri], bed_neg.ends[ri],
                         do_fwd, do_rc, &hits, &mx);
        sn += hits;
        if (hits > 0) hn++;
        if (args.test_mode == TEST_RANKSUM && max_scores_arr)
          max_scores_arr[mi*n_total + pos_unit_n + ri] = mx;
      }
    } else {
      for (uint64_t si=0; si<neg_set.n; si++) {
        uint64_t hits=0; int mx=INT_MIN;
        score_seq_scan(motif, neg_set.seqs[si], neg_set.sizes[si], &hits, &mx);
        sn += hits;
        if (hits > 0) hn++;
        if (args.test_mode == TEST_RANKSUM && max_scores_arr)
          max_scores_arr[mi*n_total + pos_unit_n + si] = mx;
      }
    }
    seq_hits_pos[mi]  = hp;
    seq_hits_neg[mi]  = hn;
    site_hits_pos[mi] = sp;
    site_hits_neg[mi] = sn;
    if (args.progress) {
      pthread_mutex_lock(&pb_lock);
      pb_counter++;
      print_pb((double)pb_counter / (double)motif_info.n);
      pthread_mutex_unlock(&pb_lock);
    }
  }
  return NULL;
}

/* ---- Usage ---- */

static void usage(void) {
  printf(
    "yamtk v%s  Copyright (C) %s  Benjamin Jean-Marie Tremblay\n"
    "Usage:  yamtk enr [options] -i positives.fa[.gz] -m motifs.txt\n"
    "\n"
    " -i <str>   Positives FASTA/FASTQ ('-' = stdin, requires -n).\n"
    " -n <str>   Negatives FASTA. If omitted, positives are shuffled (see -k, -s).\n"
    " -m <str>   Motif file (MEME/JASPAR/HOMER/HOCOMOCO).\n"
    " -o <str>   Output TSV file (default: stdout).\n"
    " -b A,C,G,T Background (default: from MEME bkg or uniform).\n"
    " -p <int>   Pseudocount for PWM generation (default: %d).\n"
    " -N <int>   Motif sites for PPM->PCM conversion (default: %d).\n"
    " -d         Deduplicate motif names (default: abort on duplicates).\n"
    " -r         Do not trim motif names to the first word.\n"
    " -t <dbl>   P-value threshold for a hit (default: %g).\n"
    " -T <str>   Test mode: 'seqs' (default), 'sites', or 'ranksum'.\n"
    "            seqs    = Fisher's exact on per-sequence hit presence.\n"
    "            sites   = Fisher's exact on per-position hit rate.\n"
    "            ranksum = Threshold-free Mann-Whitney U on max PWM score.\n"
    " -R         Disable reverse-strand scoring.\n"
    " -M         Mask lower-case bases (skip scanning at those positions).\n"
    " -k <int>   Shuffle k-mer size when -n is absent (default: %d).\n"
    " -s <uint>  RNG seed for shuffling (default: time-seeded).\n"
    " -q <dbl>   Only report rows with q-value <= this (default: %g).\n"
    " -j <int>   Threads (default: 1).\n"
    " -l         Low-memory streaming mode (one seq at a time). Incompatible\n"
    "            with -T ranksum, stdin (-i -), and forces -j 1.\n"
    " -x <bed>   BED file restricting scoring to ranges in positives. Each\n"
    "            BED region is treated as one independent unit for the test.\n"
    "            BED col-6 strand: '.'=both (resp. -R), '+'=fwd, '-'=rev.\n"
    " -X <bed>   BED file restricting scoring to ranges in negatives. Requires\n"
    "            -n and -x; if absent, full -n FASTA is used. If -x is set and\n"
    "            -n is absent, the shuffled null is built from BED slices.\n"
    " -g         Show progress bar.\n"
    " -v / -w / -h   Verbose / very-verbose / help.\n"
    , YAMTK_VERSION, YAMTK_YEAR,
    DEFAULT_PSEUDOCOUNT, DEFAULT_NSITES, DEFAULT_PVALUE, DEFAULT_SHUFFLE_K,
    DEFAULT_QVALUE_FILTER
  );
}

/* ---- Result sort ---- */

typedef struct {
  uint64_t motif_i;
  double   pvalue;
  double   qvalue;
  double   effect;
  double   log2_effect;
} result_t;

static int cmp_result(const void *a, const void *b) {
  const result_t *ra=(const result_t *)a, *rb=(const result_t *)b;
  if (ra->qvalue != rb->qvalue) return (ra->qvalue > rb->qvalue) - (ra->qvalue < rb->qvalue);
  if (ra->pvalue != rb->pvalue) return (ra->pvalue > rb->pvalue) - (ra->pvalue < rb->pvalue);
  return (ra->motif_i > rb->motif_i) - (ra->motif_i < rb->motif_i);
}

/* ---- main_enr ---- */

int main_enr(int argc, char **argv) {

  struct timespec ts_program;
  clock_gettime(CLOCK_MONOTONIC, &ts_program);

  /* Initial allocations */
  motifs = malloc(sizeof(*motifs) * ALLOC_CHUNK_SIZE);
  if (!motifs) badexit("Error: Failed to allocate motifs.");
  motif_info.n_alloc = ALLOC_CHUNK_SIZE;

  threads = malloc(sizeof(pthread_t));
  if (!threads) badexit("Error: Failed to allocate threads.");

  int opt;
  int has_pos=0, has_neg=0, has_motifs=0;
  int use_stdout=1, use_stdin_pos=0;
  char *user_bkg = NULL;
  char *test_mode_str = NULL;
  char *seed_str = NULL;
  char neg_source_str[512]; neg_source_str[0]='\0';

  while ((opt = getopt(argc, argv, "i:n:m:o:b:p:N:t:T:RMdrk:s:q:j:lgx:X:vwh")) != -1) {
    switch (opt) {
      case 'i':
        if (files.i_open) badexit("Error: -i specified more than once.");
        has_pos=1;
        if (optarg[0]=='-'&&optarg[1]=='\0') {
          files.i = gzdopen(fileno(stdin), "r");
          use_stdin_pos=1;
        } else {
          files.i = gzopen(optarg, "r");
          if (!files.i) {
            fprintf(stderr, "Error: Failed to open positives file \"%s\" [%s]",
              optarg, strerror(errno));
            badexit("");
          }
          args.pos_path = optarg;
        }
        files.i_open=1;
        break;
      case 'n':
        if (files.n_open) badexit("Error: -n specified more than once.");
        has_neg=1;
        files.n = gzopen(optarg, "r");
        if (!files.n) {
          fprintf(stderr, "Error: Failed to open negatives file \"%s\" [%s]",
            optarg, strerror(errno));
          badexit("");
        }
        files.n_open=1;
        snprintf(neg_source_str, sizeof(neg_source_str), "%s", "fasta");
        break;
      case 'm':
        if (files.m_open) badexit("Error: -m specified more than once.");
        has_motifs=1;
        files.m = fopen(optarg, "r");
        if (!files.m) {
          fprintf(stderr, "Error: Failed to open motif file \"%s\" [%s]",
            optarg, strerror(errno));
          badexit("");
        }
        files.m_open=1;
        break;
      case 'o':
        if (files.o_open) badexit("Error: -o specified more than once.");
        use_stdout=0;
        files.o = fopen(optarg, "w");
        if (!files.o) {
          fprintf(stderr, "Error: Failed to create output file \"%s\" [%s]",
            optarg, strerror(errno));
          badexit("");
        }
        files.o_open=1;
        break;
      case 'b':
        if (user_bkg) badexit("Error: -b specified more than once.");
        user_bkg = optarg;
        args.use_user_bkg = 1;
        break;
      case 'p':
        if (str_to_int(optarg, &args.pseudocount)) badexit("Error: Failed to parse -p value.");
        if (!args.pseudocount) badexit("Error: -p must be a positive integer.");
        break;
      case 'N':
        if (str_to_int(optarg, &args.nsites)) badexit("Error: Failed to parse -N value.");
        if (!args.nsites) badexit("Error: -N must be a positive integer.");
        break;
      case 'd':
        if (args.dedup) badexit("Error: -d specified more than once.");
        args.dedup = 1;
        break;
      case 'r':
        if (!args.trim_names) badexit("Error: -r specified more than once.");
        args.trim_names = 0;
        break;
      case 't':
        if (str_to_double(optarg, &args.pvalue)) badexit("Error: Failed to parse -t value.");
        if (args.pvalue <= 0.0 || args.pvalue > 1.0) badexit("Error: -t must be in (0,1].");
        break;
      case 'T':
        test_mode_str = optarg;
        break;
      case 'R':
        args.scan_rc = 0;
        break;
      case 'M':
        args.mask = 1;
        break;
      case 'k':
        if (str_to_int(optarg, &args.shuffle_k)) badexit("Error: Failed to parse -k value.");
        if (args.shuffle_k < 1 || args.shuffle_k > MAX_K)
          badexit("Error: -k must be between 1 and 9.");
        break;
      case 's':
        seed_str = optarg;
        args.use_seed = 1;
        break;
      case 'q':
        if (str_to_double(optarg, &args.qvalue_filter)) badexit("Error: Failed to parse -q value.");
        if (args.qvalue_filter < 0.0 || args.qvalue_filter > 1.0) badexit("Error: -q must be in [0,1].");
        break;
      case 'j':
        if (str_to_int(optarg, &args.nthreads)) badexit("Error: Failed to parse -j value.");
        if (args.nthreads < 1) badexit("Error: -j must be >= 1.");
        break;
      case 'g':
        args.progress = 1;
        break;
      case 'l':
        args.low_mem = 1;
        break;
      case 'x':
        if (files.x_pos_open) badexit("Error: -x specified more than once.");
        files.x_pos = gzopen(optarg, "r");
        if (!files.x_pos) {
          fprintf(stderr, "Error: Failed to open positive BED file \"%s\" [%s]",
            optarg, strerror(errno));
          badexit("");
        }
        files.x_pos_open = 1;
        args.use_bed_pos = 1;
        break;
      case 'X':
        if (files.x_neg_open) badexit("Error: -X specified more than once.");
        files.x_neg = gzopen(optarg, "r");
        if (!files.x_neg) {
          fprintf(stderr, "Error: Failed to open negative BED file \"%s\" [%s]",
            optarg, strerror(errno));
          badexit("");
        }
        files.x_neg_open = 1;
        args.use_bed_neg = 1;
        break;
      case 'w':
        args.w=1;
        /* fall through */
      case 'v':
        args.v=1;
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
  if (!has_motifs) badexit("Error: Missing -m (motif file).");
  if (!has_neg && use_stdin_pos) badexit("Error: Cannot shuffle stdin input; provide -n.");
  if (args.use_bed_neg && !has_neg)
    badexit("Error: -X requires -n (external negative FASTA).");
  if (args.use_bed_neg && !args.use_bed_pos)
    badexit("Error: -X must be paired with -x.");
  if (args.use_bed_pos && use_stdin_pos)
    badexit("Error: -x is incompatible with stdin input (cannot index sequence names).");
  if (args.low_mem && args.use_bed_pos && !has_neg)
    badexit("Error: -l with -x requires -n (shuffled-from-BED is not streamable).");

  if (test_mode_str) {
    if      (strcmp(test_mode_str,"seqs")==0)    args.test_mode=TEST_SEQS;
    else if (strcmp(test_mode_str,"sites")==0)   args.test_mode=TEST_SITES;
    else if (strcmp(test_mode_str,"ranksum")==0) args.test_mode=TEST_RANKSUM;
    else badexit("Error: -T must be one of: seqs, sites, ranksum.");
  }

  if (args.low_mem) {
    if (args.test_mode == TEST_RANKSUM)
      badexit("Error: -l is incompatible with -T ranksum. Use -T seqs or -T sites.");
    if (use_stdin_pos)
      badexit("Error: -l is incompatible with stdin input (cannot rewind for shuffle).");
    if (args.nthreads > 1) {
      if (args.v)
        fprintf(stderr, "Warning: -l forces -j 1 (ignoring -j %d).\n", args.nthreads);
      args.nthreads = 1;
    }
  }

  if (args.use_seed && seed_str) {
    if (str_to_uint64_t(seed_str, &args.seed))
      badexit("Error: Failed to parse -s value.");
  }

  if (use_stdout) { files.o=stdout; files.o_open=1; }

  if (user_bkg) parse_user_bkg(user_bkg);

  seq_hash_tab = kh_init(seq_str_h);
  if (!seq_hash_tab) badexit("Error: Failed to init sequence name hash.");

  /* Load motifs */
  if (args.v) fprintf(stderr, "Loading motifs ...\n");
  time_t t0 = time(NULL);
  load_motifs();
  find_motif_dupes();
  if (args.v) print_time((uint64_t)difftime(time(NULL),t0), "load motifs");

  /* Cap threads to motif count */
  if ((uint64_t)args.nthreads > motif_info.n) {
    if (args.v && args.nthreads > 1)
      fprintf(stderr, "Note: Reducing threads to motif count (%'" PRIu64 ").\n", motif_info.n);
    args.nthreads = (int)motif_info.n;
  }

  /* Assign motifs to threads */
  pthread_t *tmp_threads = realloc(threads, sizeof(pthread_t)*args.nthreads);
  if (!tmp_threads) badexit("Error: Failed to alloc threads.");
  threads = tmp_threads;
  for (uint64_t i=0; i<motif_info.n; i++)
    motifs[i]->thread = (uint64_t)(((double)i/motif_info.n)*args.nthreads);

  if (alloc_cdf()) badexit("Error: alloc_cdf failed.");

  if (args.low_mem) {
    /* Allocate per-motif accumulators */
    seq_hits_pos  = calloc(motif_info.n, sizeof(uint64_t));
    seq_hits_neg  = calloc(motif_info.n, sizeof(uint64_t));
    site_hits_pos = calloc(motif_info.n, sizeof(uint64_t));
    site_hits_neg = calloc(motif_info.n, sizeof(uint64_t));
    if (!seq_hits_pos||!seq_hits_neg||!site_hits_pos||!site_hits_neg)
      badexit("Error: Failed to alloc result arrays.");
    if (args.test_mode == TEST_SITES) {
      stream_pos_positions = calloc(motif_info.n, sizeof(uint64_t));
      stream_neg_positions = calloc(motif_info.n, sizeof(uint64_t));
      if (!stream_pos_positions||!stream_neg_positions)
        badexit("Error: Failed to alloc streaming position counts.");
    }
    setup_motifs_low_mem();

    /* Read BED(s) before streaming so per-seq lookups are O(1) inside the loop. */
    if (args.use_bed_pos) {
      if (args.v) fprintf(stderr, "Reading positive BED ...\n");
      enr_bed_read(&bed_pos, files.x_pos, "positive");
      enr_bed_build_hash(&bed_pos);
    }
    if (args.use_bed_neg) {
      if (args.v) fprintf(stderr, "Reading negative BED ...\n");
      enr_bed_read(&bed_neg, files.x_neg, "negative");
      enr_bed_build_hash(&bed_neg);
    }

    /* Positives: stream + score */
    if (args.v) fprintf(stderr, "Streaming positives ...\n");
    t0 = time(NULL);
    kseq_t *kseq_pos = kseq_init(files.i);
    stream_score_set(kseq_pos, &pos_set, /*is_pos=*/1, "Positive",
                     args.use_bed_pos ? &bed_pos : NULL);
    kseq_destroy(kseq_pos);
    if (args.v) print_time((uint64_t)difftime(time(NULL),t0), "stream+score positives");

    /* Negatives: stream + score */
    if (has_neg) {
      if (args.v) fprintf(stderr, "Streaming negatives ...\n");
      t0 = time(NULL);
      kseq_t *kseq_neg = kseq_init(files.n);
      stream_score_set(kseq_neg, &neg_set, /*is_pos=*/0, "Negative",
                       args.use_bed_neg ? &bed_neg : NULL);
      kseq_destroy(kseq_neg);
      if (args.v) print_time((uint64_t)difftime(time(NULL),t0), "stream+score negatives");
      snprintf(neg_source_str, sizeof(neg_source_str), "%s", "fasta");
    } else {
      /* Shuffled negatives: re-open positives, shuffle each seq in place into
         a reusable scratch buffer, score the shuffled bases. */
      if (!args.pos_path) badexit("Error: -l with shuffled negatives requires a file -i (no stdin).");
      uint64_t seed = args.use_seed ? args.seed : (uint64_t) time(NULL);
      sxrand_r(&xrng, seed);
      if (args.v) {
        fprintf(stderr, "RNG seed: %" PRIu64 "\n", seed);
        fprintf(stderr, "Streaming + shuffling negatives (k=%d) ...\n", args.shuffle_k);
      }
      gzFile fi2 = gzopen(args.pos_path, "r");
      if (!fi2) {
        fprintf(stderr, "Error: Failed to re-open positives file \"%s\" [%s]\n",
          args.pos_path, strerror(errno));
        badexit("");
      }
      kseq_t *kseq_neg = kseq_init(fi2);

      /* Euler-only scratch (k > 1) */
      const uint64_t ksz = (uint64_t) args.shuffle_k;
      uint64_t      *kmer_tab   = NULL;
      unsigned char *inv_vtx    = NULL;
      uint64_t      *euler_path = NULL;
      uint64_t      *next_idx   = NULL;
      if (ksz > 1) {
        kmer_tab   = calloc(pow5[ksz],   sizeof(uint64_t));
        inv_vtx    = calloc(pow5[ksz-1], sizeof(unsigned char));
        euler_path = calloc(pow5[ksz-1], sizeof(uint64_t));
        next_idx   = calloc(pow5[ksz-1], sizeof(uint64_t));
        if (!kmer_tab||!inv_vtx||!euler_path||!next_idx)
          badexit("Error: Failed to allocate Euler shuffle scratch.");
      }

      /* Growable per-seq scratch */
      uint64_t       scratch_cap = 64 * 1024;
      unsigned char *scratch     = malloc(scratch_cap);
      if (!scratch) badexit("Error: Failed to allocate streaming shuffle scratch.");

      t0 = time(NULL);
      ERASE_ARRAY(char_counts, 256);
      neg_set.n = 0;
      neg_set.total_bases = 0;
      neg_set.unknowns = 0;
      neg_set.gc_pct = 0.0;
      int ret_neg;
      while ((ret_neg = kseq_read(kseq_neg)) >= 0) {
        const uint64_t L = kseq_neg->seq.l;
        if (L == 0) continue;
        if (L > scratch_cap) {
          uint64_t new_cap = scratch_cap;
          while (new_cap < L) new_cap *= 2;
          unsigned char *tmp = realloc(scratch, new_cap);
          if (!tmp) badexit("Error: Failed to grow streaming shuffle scratch.");
          scratch = tmp;
          scratch_cap = new_cap;
        }
        memcpy(scratch, kseq_neg->seq.s, L);
        if (ksz == 1) {
          if (L > 1) shuffle_fisher_yates(scratch, L);
        } else if (L >= ksz) {
          memset(kmer_tab,   0, sizeof(uint64_t)      * pow5[ksz]);
          memset(inv_vtx,    0, sizeof(unsigned char) * pow5[ksz-1]);
          memset(euler_path, 0, sizeof(uint64_t)      * pow5[ksz-1]);
          memset(next_idx,   0, sizeof(uint64_t)      * pow5[ksz-1]);
          count_kmers(scratch, L, kmer_tab, ksz);
          shuffle_euler(scratch, L, ksz, kmer_tab, inv_vtx, euler_path, next_idx);
        }
        neg_set.n++;
        for (uint64_t j = 0; j < L; j++) char_counts[scratch[j]]++;
        score_one_seq_all_motifs(scratch, L, /*is_pos=*/0);
      }
      if (ret_neg == -2) badexit("Error: Failed to parse FASTQ qualities (shuffled negatives).");
      if (ret_neg < -2)  badexit("Error: Failed to read positives for shuffling.");

      kseq_destroy(kseq_neg);
      gzclose(fi2);
      free(scratch);
      free(kmer_tab); free(inv_vtx); free(euler_path); free(next_idx);

      uint64_t total = 0;
      for (uint64_t c = 0; c < 256; c++) total += char_counts[c];
      neg_set.total_bases = total;
      if (total) {
        neg_set.unknowns = total - standard_base_count();
        neg_set.gc_pct   = calc_gc() * 100.0;
      }
      if (args.v)
        fprintf(stderr, "Negative (shuffled) set: %'" PRIu64 " seq(s), %'" PRIu64 " bp\n",
          neg_set.n, neg_set.total_bases);
      if (args.v) print_time((uint64_t)difftime(time(NULL),t0), "stream+shuffle+score negatives");
      snprintf(neg_source_str, sizeof(neg_source_str),
        "shuffle k=%d seed=%" PRIu64,
        args.shuffle_k, args.use_seed ? args.seed : (uint64_t) 0);
    }
    free_cdf();
    /* In streaming mode pos_set.n / neg_set.n are only known now. Lock in
       unit counts before falling through to compute_stats. */
    pos_unit_n = args.use_bed_pos ? bed_pos.n_regions : pos_set.n;
    neg_unit_n = args.use_bed_neg ? bed_neg.n_regions : neg_set.n;
    /* Fall through to the existing stats / output block. */
    goto compute_stats;
  }

  /* Load positives */
  if (args.v) fprintf(stderr, "Loading positives ...\n");
  t0 = time(NULL);
  kseq_t *kseq_pos = kseq_init(files.i);
  pos_set.n_alloc = ALLOC_CHUNK_SIZE;
  pos_set.names = malloc(sizeof(*pos_set.names)*ALLOC_CHUNK_SIZE);
  pos_set.seqs  = malloc(sizeof(*pos_set.seqs) *ALLOC_CHUNK_SIZE);
  pos_set.sizes = malloc(sizeof(*pos_set.sizes)*ALLOC_CHUNK_SIZE);
  if (!pos_set.names||!pos_set.seqs||!pos_set.sizes) badexit("Error: alloc pos_set failed.");
  load_seq_set(kseq_pos, &pos_set, "positive");
  if (args.v) print_time((uint64_t)difftime(time(NULL),t0), "load positives");

  /* Read + bind the positive BED, if any. Validation is fatal on unknown
     seq names or out-of-bounds ranges; this is statistically load-bearing. */
  if (args.use_bed_pos) {
    if (args.v) fprintf(stderr, "Reading positive BED ...\n");
    enr_bed_read(&bed_pos, files.x_pos, "positive");
    enr_bed_bind_seqs(&bed_pos, &pos_set, "positive");
  }

  /* Load or shuffle negatives */
  if (has_neg) {
    if (args.v) fprintf(stderr, "Loading negatives ...\n");
    t0 = time(NULL);
    kseq_t *kseq_neg = kseq_init(files.n);
    neg_set.n_alloc = ALLOC_CHUNK_SIZE;
    neg_set.names = malloc(sizeof(*neg_set.names)*ALLOC_CHUNK_SIZE);
    neg_set.seqs  = malloc(sizeof(*neg_set.seqs) *ALLOC_CHUNK_SIZE);
    neg_set.sizes = malloc(sizeof(*neg_set.sizes)*ALLOC_CHUNK_SIZE);
    if (!neg_set.names||!neg_set.seqs||!neg_set.sizes) badexit("Error: alloc neg_set failed.");
    load_seq_set(kseq_neg, &neg_set, "negative");
    if (args.v) print_time((uint64_t)difftime(time(NULL),t0), "load negatives");
    if (args.use_bed_neg) {
      if (args.v) fprintf(stderr, "Reading negative BED ...\n");
      enr_bed_read(&bed_neg, files.x_neg, "negative");
      enr_bed_bind_seqs(&bed_neg, &neg_set, "negative");
    }
    /* Warn if length distributions differ greatly */
    if (args.v) {
      double ratio = (pos_set.total_bases > neg_set.total_bases)
        ? (double)pos_set.total_bases/neg_set.total_bases
        : (double)neg_set.total_bases/pos_set.total_bases;
      if (ratio > 2.0)
        fprintf(stderr,
          "Warning: Median sequence lengths differ between sets (ratio=%.2fx); "
          "enrichment may be biased for -T sites.\n", ratio);
    }
  } else {
    make_shuffled_negatives();
    uint64_t seed = args.use_seed ? args.seed : 0;
    snprintf(neg_source_str, sizeof(neg_source_str),
      "shuffle k=%d seed=%" PRIu64, args.shuffle_k, seed);
  }

  /* Now that both sets and any BEDs are loaded, lock in the unit counts.
     Each BED region (or whole sequence, sans BED) is one independent
     observation for the enrichment test. */
  pos_unit_n = args.use_bed_pos ? bed_pos.n_regions : pos_set.n;
  neg_unit_n = args.use_bed_neg ? bed_neg.n_regions : neg_set.n;

  /* Warn for small ranksum samples */
  if (args.test_mode == TEST_RANKSUM && args.v) {
    if (pos_unit_n < 10 || neg_unit_n < 10)
      fprintf(stderr,
        "Warning: Small sample sizes (pos=%'" PRIu64 ", neg=%'" PRIu64 ") may make\n"
        "  normal approximation for ranksum unreliable.\n",
        pos_unit_n, neg_unit_n);
  }

  /* Allocate result matrices */
  seq_hits_pos  = calloc(motif_info.n, sizeof(uint64_t));
  seq_hits_neg  = calloc(motif_info.n, sizeof(uint64_t));
  site_hits_pos = calloc(motif_info.n, sizeof(uint64_t));
  site_hits_neg = calloc(motif_info.n, sizeof(uint64_t));
  if (!seq_hits_pos||!seq_hits_neg||!site_hits_pos||!site_hits_neg)
    badexit("Error: Failed to alloc result arrays.");

  if (args.test_mode == TEST_RANKSUM) {
    uint64_t n_total = pos_unit_n + neg_unit_n;
    uint64_t ranksum_bytes = motif_info.n * n_total * sizeof(int);
    /* Hard cap: refuse rather than OOM-crash. 4 GB matches the practical
       single-allocation ceiling on most 64-bit systems and is plenty for
       any realistic enrichment run. */
    static const uint64_t RANKSUM_MAX_BYTES = (uint64_t) 4 * 1024 * 1024 * 1024;
    if (ranksum_bytes > RANKSUM_MAX_BYTES) {
      fprintf(stderr,
        "Error: ranksum max-score matrix would need %.2f GB (%" PRIu64 " motifs x %" PRIu64
        " sequences); cap is %.0f GB.\n"
        "Use -T seqs or -T sites for inputs this large.\n",
        (double) ranksum_bytes / (1024.0 * 1024.0 * 1024.0),
        motif_info.n, n_total,
        (double) RANKSUM_MAX_BYTES / (1024.0 * 1024.0 * 1024.0));
      badexit("");
    }
    if (args.v)
      fprintf(stderr, "Allocating %.2f MB for ranksum max-score matrix.\n",
        (double)ranksum_bytes / (1024.0*1024.0));
    max_scores_arr = malloc(ranksum_bytes);
    if (!max_scores_arr) badexit("Error: Failed to alloc ranksum matrix.");
    /* Initialize to INT_MIN (sequence too short contributes lowest rank). */
    for (uint64_t i=0; i<motif_info.n*n_total; i++) max_scores_arr[i]=INT_MIN;
  }

  /* Scan */
  if (args.v) fprintf(stderr, "Scanning ...\n");
  t0 = time(NULL);
  if (args.progress) print_pb(0.0);
  for (uint64_t t=0; t<(uint64_t)args.nthreads; t++) {
    uint64_t *ti = malloc(sizeof(*ti));
    if (!ti) badexit("Error: Failed to alloc thread index.");
    *ti = t;
    pthread_create(&threads[t], NULL, enr_sub_process, ti);
  }
  for (uint64_t t=0; t<(uint64_t)args.nthreads; t++)
    pthread_join(threads[t], NULL);
  if (args.progress) fprintf(stderr, "\n");
  free_cdf();
  if (args.v) {
    print_time((uint64_t)difftime(time(NULL),t0), "scan");
  }

compute_stats:
  /* Compute statistics */
  if (args.v) fprintf(stderr, "Computing statistics ...\n");

  /* Precompute per-motif total scoring positions for -T sites */
  uint64_t *pos_positions = NULL, *neg_positions = NULL;
  int positions_borrowed = 0;
  if (args.test_mode == TEST_SITES) {
    if (args.low_mem) {
      pos_positions = stream_pos_positions;
      neg_positions = stream_neg_positions;
      positions_borrowed = 1;
    } else {
      pos_positions = malloc(sizeof(uint64_t)*motif_info.n);
      neg_positions = malloc(sizeof(uint64_t)*motif_info.n);
      if (!pos_positions||!neg_positions) badexit("Error: Failed to alloc positions.");
      const int strands_full = args.scan_rc ? 2 : 1;
      for (uint64_t mi=0; mi<motif_info.n; mi++) {
        const uint64_t sz = motifs[mi]->size;
        pos_positions[mi]=0; neg_positions[mi]=0;
        if (args.use_bed_pos) {
          for (uint64_t ri=0; ri<bed_pos.n_regions; ri++)
            pos_positions[mi] += region_positions(bed_pos.starts[ri], bed_pos.ends[ri],
                                                  sz, bed_pos.strands[ri]);
        } else {
          for (uint64_t si=0; si<pos_set.n; si++)
            if (pos_set.sizes[si] >= sz)
              pos_positions[mi] += (pos_set.sizes[si]-sz+1) * strands_full;
        }
        if (args.use_bed_neg) {
          for (uint64_t ri=0; ri<bed_neg.n_regions; ri++)
            neg_positions[mi] += region_positions(bed_neg.starts[ri], bed_neg.ends[ri],
                                                  sz, bed_neg.strands[ri]);
        } else {
          for (uint64_t si=0; si<neg_set.n; si++)
            if (neg_set.sizes[si] >= sz)
              neg_positions[mi] += (neg_set.sizes[si]-sz+1) * strands_full;
        }
      }
    }
  }

  double *pvalues = malloc(sizeof(double)*motif_info.n);
  double *qvalues = malloc(sizeof(double)*motif_info.n);
  double *effects = malloc(sizeof(double)*motif_info.n);
  double *log2eff = malloc(sizeof(double)*motif_info.n);
  if (!pvalues||!qvalues||!effects||!log2eff) badexit("Error: Failed to alloc stat arrays.");

  const double eps = 1e-9;

  for (uint64_t mi=0; mi<motif_info.n; mi++) {
    uint64_t a, b, c, d;
    double eff, l2e, pval;
    if (args.test_mode == TEST_SEQS) {
      a = seq_hits_pos[mi];
      b = seq_hits_neg[mi];
      c = pos_unit_n - a;
      d = neg_unit_n - b;
      pval = fishers_exact_log_greater(a, b, c, d);
      double pos_rate = (double)a / (double)pos_unit_n;
      double neg_rate = (double)b / (double)neg_unit_n;
      double neg_denom = (neg_rate > eps) ? neg_rate : (0.5 / (double)neg_unit_n);
      eff = pos_rate / neg_denom;
      l2e = log2(eff);
    } else if (args.test_mode == TEST_SITES) {
      a = site_hits_pos[mi];
      b = site_hits_neg[mi];
      uint64_t pp = pos_positions[mi], np = neg_positions[mi];
      c = (a <= pp) ? pp - a : 0;
      d = (b <= np) ? np - b : 0;
      pval = fishers_exact_log_greater(a, b, c, d);
      double pos_rate = (pp > 0) ? (double)a/(double)pp : 0.0;
      double neg_denom = (np > 0 && b > 0) ? (double)b/(double)np : (0.5 / ((double)(np>0?np:1)));
      eff = (pos_rate > 0.0) ? pos_rate / neg_denom : 0.0;
      l2e = (eff > 0.0) ? log2(eff) : -INFINITY;
    } else { /* ranksum */
      const uint64_t n_total = pos_unit_n + neg_unit_n;
      double auc;
      mann_whitney_u_greater(&max_scores_arr[mi*n_total], (int)pos_unit_n, (int)neg_unit_n,
                              &pval, &auc);
      eff = auc;
      double auc_denom = (1.0-auc > eps) ? (1.0-auc) : eps;
      l2e = log2(auc / auc_denom);
    }
    pvalues[mi] = pval;
    effects[mi] = eff;
    log2eff[mi] = l2e;
  }

  bh_qvalues(pvalues, motif_info.n, qvalues);

  if (!positions_borrowed) { free(pos_positions); free(neg_positions); }
  free(stream_pos_positions); stream_pos_positions = NULL;
  free(stream_neg_positions); stream_neg_positions = NULL;

  /* Sort results */
  result_t *results = malloc(sizeof(result_t)*motif_info.n);
  if (!results) badexit("Error: Failed to alloc results.");
  for (uint64_t i=0; i<motif_info.n; i++) {
    results[i].motif_i   = i;
    results[i].pvalue    = pvalues[i];
    results[i].qvalue    = qvalues[i];
    results[i].effect    = effects[i];
    results[i].log2_effect = log2eff[i];
  }
  qsort(results, motif_info.n, sizeof(result_t), cmp_result);

  free(pvalues); free(qvalues); free(effects); free(log2eff);

  /* Write output */
  const char *test_label  = (args.test_mode==TEST_SEQS)?"seqs":(args.test_mode==TEST_SITES)?"sites":"ranksum";
  const char *effect_label = (args.test_mode==TEST_SEQS)?"seq_fold":(args.test_mode==TEST_SITES)?"site_fold":"auc";

  fprintf(files.o, "##yamenr v%s [ ", YAMTK_VERSION);
  for (int i=1; i<argc; i++) fprintf(files.o, "%s ", argv[i]);
  fprintf(files.o, "]\n");
  fprintf(files.o,
    "##MotifCount=%'" PRIu64 " PosSeqs=%'" PRIu64 " NegSeqs=%'" PRIu64
    " PosUnits=%'" PRIu64 " NegUnits=%'" PRIu64
    " NegSource=%s TestMode=%s Effect=%s\n",
    motif_info.n, pos_set.n, neg_set.n, pos_unit_n, neg_unit_n,
    neg_source_str, test_label, effect_label);
  fprintf(files.o,
    "##motif\tmotif_id\tconsensus\tpos_n\tpos_seq_hits\tpos_site_hits"
    "\tneg_n\tneg_seq_hits\tneg_site_hits\teffect\tlog2_effect\tpvalue\tqvalue\n");

  for (uint64_t ri=0; ri<motif_info.n; ri++) {
    uint64_t mi = results[ri].motif_i;
    if (results[ri].qvalue > args.qvalue_filter) continue;
    fprintf(files.o,
      "%s\t%" PRIu64 "\t%s\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64
      "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64
      "\t%.6g\t%.6g\t%.9g\t%.9g\n",
      motifs[mi]->name, mi+1, motifs[mi]->consensus,
      pos_unit_n, seq_hits_pos[mi], site_hits_pos[mi],
      neg_unit_n, seq_hits_neg[mi], site_hits_neg[mi],
      results[ri].effect, results[ri].log2_effect,
      results[ri].pvalue, results[ri].qvalue);
  }

  /* Clean up */
  free(results);
  free(max_scores_arr); max_scores_arr=NULL;
  free(seq_hits_pos); free(seq_hits_neg);
  free(site_hits_pos); free(site_hits_neg);
  seq_hits_pos=seq_hits_neg=site_hits_pos=site_hits_neg=NULL;
  free(threads); threads=NULL;
  free_motifs();
  free_seq_set(&pos_set);
  free_seq_set(&neg_set);
  kh_destroy(seq_str_h, seq_hash_tab); seq_hash_tab = NULL;
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
