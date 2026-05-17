/*
 *   yamref: Motif refinement against positive sequences
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
#include "kseq.h"
#include "version.h"

KSEQ_INIT(gzFile, gzread)

/* ---- Constants ---- */

#define MAX_NAME_SIZE           ((uint64_t) 256)
#define MAX_MOTIF_SIZE          ((uint64_t) 250)
#define MAX_MOTIF_WIDTH                       50  /* practical width cap (matches yamme) */
#define MIN_REFINE_HITS                       20
#define MIN_IC_BITS                          0.5
#define CONSENSUS_THRESHOLD                 0.75
#define AMBIGUITY_SCORE                -10000000
#define MIN_BKG_VALUE                      0.001
#define MAX_CDF_SIZE        ((uint64_t) 2097152)
#define PWM_INT_MULTIPLIER                1000.0
#define USER_BKG_MAX_SIZE       ((uint64_t) 256)
#define MEME_BKG_MAX_SIZE       ((uint64_t) 256)
#define MOTIF_VALUE_MAX_CHAR    ((uint64_t) 256)
#define SEQ_NAME_MAX_CHAR       ((uint64_t) 512)
#define ALLOC_CHUNK_SIZE        ((uint64_t) 256)
#define DEFAULT_PSEUDOCOUNT                    1
#define DEFAULT_NSITES                      1000
#define DEFAULT_HIT_PVAL                     1e-3
#define DEFAULT_R_PASSES                       2
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

/* ---- Motif struct (mirrors yamme) ---- */

typedef struct motif_t {
  int         pwm[MAX_MOTIF_SIZE];
  int         pwm_rc[MAX_MOTIF_SIZE];
  double     *cdf;
  int         threshold;
  uint64_t    size;
  uint64_t    cdf_size;
  uint64_t    file_line_num;
  int         min;
  int         max;
  int         max_score;
  int         min_score;
  int         cdf_max;
  int         cdf_offset;
  char        name[MAX_NAME_SIZE];
  double     *tmp_pdf;
  double      pwm_probs[MAX_MOTIF_WIDTH][4];  /* pseudocount-adj probs for output */
  uint64_t    nsites_actual;
  int         dropped;
} motif_t;

enum MOTIF_FMT {
  FMT_MEME     = 1,
  FMT_HOMER    = 2,
  FMT_JASPAR   = 3,
  FMT_HOCOMOCO = 4,
  FMT_UNKNOWN  = 5
};

typedef struct {
  int       fmt;
  uint64_t  n;
  uint64_t  n_alloc;
} motif_info_t;

static motif_info_t motif_info = { .fmt = 0, .n = 0, .n_alloc = 0 };
static motif_t **motifs = NULL;

/* ---- Args ---- */

typedef struct {
  double   bkg[4];
  double   hit_pval;
  int      pseudocount;
  int      nsites;
  int      r_passes;
  int      extend;
  int      auto_extend;
  int      ic_trim;
  int      quality_gate;
  double   ic_min;
  int      scan_rc;
  int      mask;
  int      use_user_bkg;
  int      progress;
  int      v;
  int      w;
} args_t;

static args_t args = {
  .bkg          = {0.25, 0.25, 0.25, 0.25},
  .hit_pval     = DEFAULT_HIT_PVAL,
  .pseudocount  = DEFAULT_PSEUDOCOUNT,
  .nsites       = DEFAULT_NSITES,
  .r_passes     = DEFAULT_R_PASSES,
  .extend       = 0,
  .auto_extend  = 0,
  .ic_trim      = 0,
  .quality_gate = 0,
  .ic_min       = MIN_IC_BITS,
  .scan_rc      = 1,
  .mask         = 0,
  .use_user_bkg = 0,
  .progress     = 0,
  .v            = 0,
  .w            = 0
};

static void print_pb(const double prog) {
  const int left = prog * PROGRESS_BAR_WIDTH;
  const int right = PROGRESS_BAR_WIDTH - left;
  fprintf(stderr, "\r[%.*s%*s] %3d%%", left, PROGRESS_BAR_STRING, right, "",
      (int)(prog * 100.0));
  fflush(stderr);
}

/* ---- Files ---- */

typedef struct {
  int    m_open, i_open, o_open;
  FILE  *m;
  gzFile i;
  FILE  *o;
} files_t;

static files_t files = { 0 };

static const char *out_path = NULL;  /* default stdout */

/* ---- Sequence set ---- */

typedef struct {
  char           **names;
  unsigned char  **seqs;
  uint64_t        *sizes;
  uint64_t         n;
  uint64_t         n_alloc;
  uint64_t         total_bases;
  uint64_t         unknowns;
} seq_set_t;

static seq_set_t pos_set = {0};

/* ---- Lookup tables ---- */

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

/* complement of 2-bit index: A(0)<->T(3), C(1)<->G(2) */
static const int comp4[4] = {3, 2, 1, 0};

/* IUPAC consensus letter -> uniform probability over its constituent bases. */
static const double consensus2probs[] = {
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

static const int consensus2index[256] = {
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

/* ---- Single shared CDF scratch (single-threaded) ---- */

static double  *cdf_scratch     = NULL;
static double  *tmp_pdf_scratch = NULL;
static uint64_t cdf_scratch_size = 0;

/* ---- Memory & timing helpers ---- */

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

/* ---- Cleanup ---- */

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

static void free_motifs(void) {
  if (!motifs) return;
  for (uint64_t i = 0; i < motif_info.n; i++) free(motifs[i]);
  free(motifs); motifs = NULL; motif_info.n = 0; motif_info.n_alloc = 0;
}

static void free_cdf(void) {
  free(cdf_scratch); cdf_scratch = NULL;
  free(tmp_pdf_scratch); tmp_pdf_scratch = NULL;
  cdf_scratch_size = 0;
}

static void close_files(void) {
  if (files.m_open) fclose(files.m);
  if (files.i_open) gzclose(files.i);
  if (files.o_open) fclose(files.o);
}

static void badexit(const char *msg) {
  fprintf(stderr, "%s\nRun yamtk ref -h for usage.\n", msg);
  free_motifs();
  free_seq_set(&pos_set);
  free_cdf();
  close_files();
  exit(EXIT_FAILURE);
}

/* ---- String helpers ---- */

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
  if (args.w)
    fprintf(stderr, "Using background: A=%.3g C=%.3g G=%.3g T=%.3g\n",
      args.bkg[0], args.bkg[1], args.bkg[2], args.bkg[3]);
}

/* ---- Motif init/add ---- */

static void init_motif(motif_t *motif) {
  ERASE_ARRAY(motif->name, MAX_NAME_SIZE);
  motif->name[0]='m'; motif->name[1]='o'; motif->name[2]='t';
  motif->name[3]='i'; motif->name[4]='f'; motif->name[5]='\0';
  motif->size=0; motif->threshold=0; motif->max_score=0;
  motif->min=0; motif->max=0; motif->cdf_size=0; motif->file_line_num=0;
  motif->min_score=0; motif->cdf_max=0;
  motif->nsites_actual=0; motif->dropped=0;
  motif->cdf = NULL; motif->tmp_pdf = NULL;
  for (uint64_t i = 0; i < MAX_MOTIF_SIZE; i++) { motif->pwm[i]=0; motif->pwm_rc[i]=0; }
  for (uint64_t i = 4; i < MAX_MOTIF_SIZE; i += 5) {
    motif->pwm[i]=AMBIGUITY_SCORE; motif->pwm_rc[i]=AMBIGUITY_SCORE;
  }
  memset(motif->pwm_probs, 0, sizeof(motif->pwm_probs));
}

static int add_motif(void) {
  motif_info.n++;
  const uint64_t last_i = motif_info.n - 1;
  if (motif_info.n > motif_info.n_alloc) {
    motif_t **tmp_ptr = realloc(motifs,
      sizeof(*motifs) * motif_info.n_alloc + sizeof(*motifs) * ALLOC_CHUNK_SIZE);
    if (!tmp_ptr) {
      fprintf(stderr, "Error: Failed to allocate memory for motifs."); return 1;
    }
    motifs = tmp_ptr;
    motif_info.n_alloc += ALLOC_CHUNK_SIZE;
  }
  motifs[last_i] = malloc(sizeof(motif_t));
  if (!motifs[last_i]) {
    fprintf(stderr, "Error: Failed to allocate memory for motif."); return 1;
  }
  init_motif(motifs[last_i]);
  return 0;
}

/* ---- PWM math ---- */

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
  if (x <= 0.0) x = 1e-12;
  return (int)(log2(x / bkg_i) * PWM_INT_MULTIPLIER);
}

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

/* ---- CDF and threshold (single-threaded) ---- */

static void fill_cdf(motif_t *motif) {
  uint64_t max_step, s;
  double pdf_sum = 0.0;
  if (motif->cdf_size > MAX_CDF_SIZE) {
    fprintf(stderr,
      "\nError: CDF size for [%s] too large (%'" PRIu64 ">%'" PRIu64 ").\n",
      motif->name, motif->cdf_size, MAX_CDF_SIZE);
    badexit("");
  }
  if (cdf_scratch_size < motif->cdf_size) {
    double *r1 = realloc(cdf_scratch, motif->cdf_size*sizeof(double));
    if (!r1) badexit("Error: Memory re-allocation for CDF failed.");
    cdf_scratch = r1;
    double *r2 = realloc(tmp_pdf_scratch, motif->cdf_size*sizeof(double));
    if (!r2) badexit("Error: Memory re-allocation for PDF failed.");
    tmp_pdf_scratch = r2;
    cdf_scratch_size = motif->cdf_size;
  }
  motif->cdf     = cdf_scratch;
  motif->tmp_pdf = tmp_pdf_scratch;
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
  motif->max_score = 0; motif->min_score = 0;
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

/* ---- Motif parser helpers ---- */

static int check_line_contains(const char *line, const char *substring) {
  uint64_t i = 0, j = 0;
  while (substring[j] != '\0') {
    if (line[i] == '\0') return 0;
    if (line[i] == substring[j]) { i++; j++; }
    else { i++; j = 0; }
  }
  return 1;
}

static uint64_t count_nonempty_chars(const char *line) {
  uint64_t i = 0, n = 0;
  while (line[i] != '\0') {
    if (line[i] != ' ' && line[i] != '\t' && line[i] != '\r' && line[i] != '\n') n++;
    i++;
  }
  return n;
}

static int check_char_is_one_of(const char c, const char *list) {
  uint64_t i = 0;
  while (list[i] != '\0') { if (c == list[i]) return 1; i++; }
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

static int detect_motif_fmt(void) {
  int jaspar_or_hocomoco = 0, file_fmt = 0, has_tabs = 0;
  char *line = NULL;
  size_t len = 0;
  ssize_t r;
  while ((r = getline(&line, &len, files.m)) != -1) {
    if (!count_nonempty_chars(line)) continue;
    if (check_line_contains(line, "MEME version \0")) {
      if (args.w) fprintf(stderr, "Detected MEME format (version %c).\n", line[13]);
      file_fmt = FMT_MEME; break;
    }
    if (jaspar_or_hocomoco) {
      if (line[0] == 'A' &&
          check_char_is_one_of('[', line) &&
          check_char_is_one_of(']', line)) {
        file_fmt = FMT_JASPAR;
        if (args.w) fprintf(stderr, "Detected JASPAR format.\n");
        break;
      } else {
        if (line[0] == 'A' ||
            check_char_is_one_of('[', line) ||
            check_char_is_one_of(']', line)) {
          free(line); badexit("Error: Detected malformed JASPAR format.");
        }
        if (has_tabs) {
          file_fmt = FMT_HOMER;
          if (args.w) fprintf(stderr, "Detected HOMER format.\n");
          break;
        } else {
          if (check_char_is_one_of('-', line)) {
            free(line); badexit("Error: yamref cannot read HOCOMOCO PWMs.");
          }
          file_fmt = FMT_HOCOMOCO;
          if (args.w) fprintf(stderr, "Detected HOCOMOCO format.\n");
          break;
        }
      }
    } else if (line[0] == '>') {
      if (check_char_is_one_of('\t', line)) has_tabs = 1;
      jaspar_or_hocomoco = 1;
    }
  }
  rewind(files.m);
  free(line);
  if (!file_fmt) file_fmt = FMT_UNKNOWN;
  return file_fmt;
}

static int normalize_probs(double *probs, const char *name) {
  double sum = probs[0] + probs[1] + probs[2] + probs[3];
  if (fabs(sum - 1.0) > 0.1) {
    fprintf(stderr,
      "Error: Position for [%s] does not add up to 1 (sum=%.3g)",
      name, sum);
    return 1;
  }
  if (fabs(sum - 1.0) > 0.02) {
    if (args.w)
      fprintf(stderr,
        "Warning: Position for [%s] does not add up to 1, adjusting (sum=%.3g)\n",
        name, sum);
    VEC_DIV(probs, sum, 4);
  }
  return 0;
}

static int get_line_probs(const motif_t *motif, const char *line, double *probs, const uint64_t n) {
  uint64_t i = 0, j = 0, which_i = -1;
  int prev_was_space = 1;
  char pos_i[MOTIF_VALUE_MAX_CHAR];
  ERASE_ARRAY(pos_i, MOTIF_VALUE_MAX_CHAR);
  while (line[i] == ' ' || line[i] == '\t') i++;
  while (line[i] != '\0' && line[i] != '\r' && line[i] != '\n') {
    if (line[i] != ' ' && line[i] != '\t') {
      pos_i[j] = line[i]; j++; prev_was_space = 0;
    } else {
      if (!prev_was_space) {
        which_i++;
        if (which_i > n - 1) {
          fprintf(stderr, "Error: Motif [%s] has too many columns (need %" PRIu64 ").",
            motif->name, n); return 1;
        }
        if (str_to_double(pos_i, &probs[which_i])) {
          fprintf(stderr, "Error: Failed to parse value for motif: %s.\n", motif->name);
          fprintf(stderr, "  Line: %s  Bad value: %s", line, pos_i);
          return 1;
        }
        ERASE_ARRAY(pos_i, MOTIF_VALUE_MAX_CHAR);
        j = 0;
      }
      prev_was_space = 1;
    }
    i++;
  }
  if (!prev_was_space) {
    which_i++;
    if (which_i > n - 1) {
      fprintf(stderr, "Error: Motif [%s] has too many columns (need %" PRIu64 ").",
        motif->name, n); return 1;
    }
    if (str_to_double(pos_i, &probs[which_i])) {
      fprintf(stderr, "Error: Failed to parse value for motif: %s.\n", motif->name);
      fprintf(stderr, "  Line: %s  Bad value: %s", line, pos_i);
      return 1;
    }
  }
  if (which_i == -1) {
    fprintf(stderr, "Error: Motif [%s] has an empty row.", motif->name); return 1;
  }
  if (which_i < n - 1) {
    fprintf(stderr, "Error: Motif [%s] has too few columns (need %" PRIu64 ").",
      motif->name, n); return 1;
  }
  return 0;
}

static int add_motif_ppm_column(motif_t *motif, const char *line, const uint64_t pos) {
  double probs[4] = {-1.0, -1.0, -1.0, -1.0};
  if (get_line_probs(motif, line, probs, 4)) return 1;
  if (normalize_probs(probs, motif->name)) return 1;
  if (pos >= MAX_MOTIF_WIDTH) {
    fprintf(stderr, "Error: Motif [%s] exceeds max width (%d).", motif->name, MAX_MOTIF_WIDTH);
    return 1;
  }
  for (int i = 0; i < 4; i++) motif->pwm_probs[pos][i] = probs[i];
  set_score(motif, 'A', pos, calc_score(probs[0], args.bkg[0]));
  set_score(motif, 'C', pos, calc_score(probs[1], args.bkg[1]));
  set_score(motif, 'G', pos, calc_score(probs[2], args.bkg[2]));
  set_score(motif, 'T', pos, calc_score(probs[3], args.bkg[3]));
  return 0;
}

/* ---- MEME reader ---- */

static int check_meme_alph(const char *line, const uint64_t line_num) {
  if (check_line_contains(line, "ALPHABET= ACDEFGHIKLMNPQRSTVWY\0")) {
    fprintf(stderr, "Error: Detected protein alphabet (L%" PRIu64 ").", line_num);
    return 1;
  }
  return 0;
}

static int get_meme_bkg(const char *line, const uint64_t line_num) {
  if (args.use_user_bkg) return 0;
  double bkg_probs[] = {-1.0, -1.0, -1.0, -1.0};
  uint64_t i = 1, let_i = 0, j = 0, empty = 0;
  char bkg_char[MEME_BKG_MAX_SIZE];
  ERASE_ARRAY(bkg_char, MEME_BKG_MAX_SIZE);
  if (line[0] != 'A') {
    fprintf(stderr, "Error: Expected first character of background line to be 'A' (L%" PRIu64 ").",
      line_num); return 1;
  }
  while (line[i] != '\0' && line[i] != '\n' && line[i] != '\r') {
    if (let_i > 3) {
      fprintf(stderr, "Error: Parsed too many background values in MEME file (L%" PRIu64 ").",
        line_num); return 1;
    }
    if (line[i] != ' ' && line[i] != '\t') {
      if (line[i] == 'C' || line[i] == 'G' || line[i] == 'T' || line[i] == 'U') {
        int expect = (line[i]=='C') ? 0 : (line[i]=='G') ? 1 : 2;
        if (!empty) {
          fprintf(stderr, "Error: Expected whitespace before letter (L%" PRIu64 ").", line_num); return 1;
        }
        if (let_i != (uint64_t)expect) {
          fprintf(stderr, "Error: Wrong letter order in MEME background (L%" PRIu64 ").", line_num); return 1;
        }
        if (str_to_double(bkg_char, &bkg_probs[let_i])) {
          fprintf(stderr, "Error: Failed to parse background value: %s", bkg_char); return 1;
        }
        ERASE_ARRAY(bkg_char, MEME_BKG_MAX_SIZE);
        let_i = expect + 1; j = 0;
      } else if (check_char_is_one_of(line[i], "0123456789.\0")) {
        bkg_char[j] = line[i]; j++;
      } else {
        fprintf(stderr, "Error: Unexpected char (%c) in MEME background (L%" PRIu64 ").",
          line[i], line_num); return 1;
      }
      empty = 0;
    } else {
      empty = 1;
    }
    i++;
  }
  if (bkg_char[0] != '\0') {
    if (str_to_double(bkg_char, &bkg_probs[let_i])) {
      fprintf(stderr, "Error: Failed to parse background value: %s", bkg_char); return 1;
    }
  }
  if (check_and_load_bkg(bkg_probs)) return 1;
  if (args.w)
    fprintf(stderr, "Found MEME bkg: A=%.3g C=%.3g G=%.3g T=%.3g\n",
      args.bkg[0], args.bkg[1], args.bkg[2], args.bkg[3]);
  return 0;
}

static void parse_meme_name(const char *line, const uint64_t motif_i) {
  /* Capture identifier AND any altname; trim_motif_name (run when args.trim_names)
     truncates at the first space later, so default behavior is unchanged. */
  uint64_t i = 5, j = 0;
  int prev_space = 1;
  while (line[i] != '\0' && line[i] != '\r' && line[i] != '\n' && j < MAX_NAME_SIZE - 1) {
    if (line[i] == ' ' || line[i] == '\t') {
      if (!prev_space && j > 0 && j < MAX_NAME_SIZE - 1) motifs[motif_i]->name[j++] = ' ';
      prev_space = 1;
    } else {
      motifs[motif_i]->name[j++] = line[i];
      prev_space = 0;
    }
    i++;
  }
  if (j > 0 && motifs[motif_i]->name[j - 1] == ' ') j--;
  motifs[motif_i]->name[j] = '\0';
}

static void read_meme(void) {
  motif_info.fmt = FMT_MEME;
  char *line = NULL;
  size_t len = 0;
  ssize_t r;
  uint64_t line_num = 0, l_p_m_L = 0, bkg_L = 0, motif_i = -1, pos_i = 0;
  int alph_detected = 0, live_motif = 0;
  int warned_bkg = 0, warned_alph = 0;
  while ((r = getline(&line, &len, files.m)) != -1) {
    line_num++;
    if (check_line_contains(line, "Background letter frequencies\0")) {
      /* Only honor the first occurrence; concatenated MEME files (e.g.
         JASPAR2026_CORE_*_meme.txt) repeat the header per chunk. */
      if (!bkg_L) bkg_L = line_num;
      else if (!warned_bkg) {
        fprintf(stderr,
          "Warning: Multiple 'Background letter frequencies' lines in MEME file "
          "(L%" PRIu64 "); using the first.\n", line_num);
        warned_bkg = 1;
      }
    } else if (bkg_L && bkg_L == line_num - 1) {
      if (get_meme_bkg(line, line_num)) { free(line); badexit(""); }
    } else if (check_line_contains(line, "ALPHABET\0")) {
      if (!alph_detected) {
        if (check_meme_alph(line, line_num)) { free(line); badexit(""); }
        alph_detected = 1;
      } else if (!warned_alph) {
        fprintf(stderr,
          "Warning: Multiple ALPHABET lines in MEME file (L%" PRIu64 "); using the first.\n",
          line_num);
        warned_alph = 1;
      }
    } else if (check_line_contains(line, "MOTIF\0")) {
      motif_i++;
      if (add_motif()) { free(line); badexit(""); }
      motifs[motif_i]->file_line_num = line_num;
      parse_meme_name(line, motif_i);
      pos_i = 0;
    } else if (check_line_contains(line, "letter-probability matrix\0")) {
      if (pos_i != 0) { free(line); badexit("Error: Malformed MEME motif."); }
      l_p_m_L = line_num;
      live_motif = 1;
    } else if (live_motif) {
      if (!count_nonempty_chars(line) || !is_ppm_data_line(line)) {
        live_motif = 0;
      } else if (line_num == (l_p_m_L + pos_i + 1)) {
        if (pos_i >= MAX_MOTIF_WIDTH) {
          fprintf(stderr, "Error: Motif [%s] is too large (max=%d)",
            motifs[motif_i]->name, MAX_MOTIF_WIDTH);
          free(line); badexit("");
        }
        if (add_motif_ppm_column(motifs[motif_i], line, pos_i)) { free(line); badexit(""); }
        pos_i++;
        motifs[motif_i]->size = pos_i;
      } else {
        live_motif = 0;
      }
    }
  }
  free(line);
  if (!motif_info.n) badexit("Error: Failed to detect any motifs in MEME file.");
  if (args.v) fprintf(stderr, "Found %'" PRIu64 " MEME motif(s).\n", motif_info.n);
}

/* ---- HOMER reader ---- */

static void parse_homer_name(const char *line, const uint64_t motif_i) {
  uint64_t name_start = 0, name_end = 0, i = 1, in_between = 0, j = 0;
  while (line[i] != '\0' && line[i] != '\r' && line[i] != '\n') {
    if (line[i] == '\t') {
      if (name_start) { name_end = i; break; }
      else in_between = 1;
    } else if (in_between) {
      if (!name_start) name_start = i;
    }
    i++;
  }
  if (!name_end) name_end = i;
  if (!name_start) { motifs[motif_i]->name[0] = '\0'; return; }
  for (uint64_t k = name_start; k < name_end && j < MAX_NAME_SIZE - 1; k++) {
    motifs[motif_i]->name[j] = line[k]; j++;
  }
  motifs[motif_i]->name[j] = '\0';
}

static void read_homer(void) {
  motif_info.fmt = FMT_HOMER;
  char *line = NULL;
  size_t len = 0;
  ssize_t r;
  uint64_t line_num = 0, motif_i = -1, pos_i = 0;
  int ready_to_start = 0;
  while ((r = getline(&line, &len, files.m)) != -1) {
    line_num++;
    if (line[0] == '>') {
      ready_to_start = 1;
      motif_i++;
      if (add_motif()) { free(line); badexit(""); }
      motifs[motif_i]->file_line_num = line_num;
      parse_homer_name(line, motif_i);
      pos_i = 0;
    } else if (count_nonempty_chars(line) && ready_to_start) {
      if (pos_i >= MAX_MOTIF_WIDTH) {
        fprintf(stderr, "Error: Motif [%s] is too large (max=%d).",
          motifs[motif_i]->name, MAX_MOTIF_WIDTH);
        free(line); badexit("");
      }
      if (add_motif_ppm_column(motifs[motif_i], line, pos_i)) { free(line); badexit(""); }
      pos_i++;
      motifs[motif_i]->size = pos_i;
    }
  }
  free(line);
  if (args.v) fprintf(stderr, "Found %'" PRIu64 " HOMER motif(s).\n", motif_info.n);
}

/* ---- JASPAR reader (count matrix) ---- */

static void parse_jaspar_name(const char *line, const uint64_t motif_i) {
  uint64_t i = 0, j = 1;
  while (line[j] != '\r' && line[j] != '\n' && line[j] != '\0' && i < MAX_NAME_SIZE - 1) {
    motifs[motif_i]->name[i] = line[j];
    i++; j++;
  }
  motifs[motif_i]->name[i] = '\0';
}

static int add_jaspar_row(motif_t *motif, const char *line, double counts[MAX_MOTIF_WIDTH][4],
                          uint64_t *width_out) {
  uint64_t left_bracket = -1, right_bracket = -1, i = 0;
  int row_i = -1;
  for (;;) {
    if (line[i] == '\r' || line[i] == '\n' || line[i] == '\0') break;
    switch (line[i]) {
      case 'a': case 'A': row_i = 0; break;
      case 'c': case 'C': row_i = 1; break;
      case 'g': case 'G': row_i = 2; break;
      case 'u': case 'U': case 't': case 'T': row_i = 3; break;
      case '[': left_bracket = i; break;
      case ']': right_bracket = i; break;
    }
    i++;
  }
  if (row_i == -1) {
    fprintf(stderr, "Error: Couldn't find ACGTU in motif [%s] row names.", motif->name);
    return 1;
  }
  if (left_bracket == (uint64_t)-1 || right_bracket == (uint64_t)-1) {
    fprintf(stderr, "Error: Couldn't find '[]' in motif [%s] row.", motif->name);
    return 1;
  }
  uint64_t k = 0, pos_i = -1;
  int prev_was_space = 1;
  char prob_c[MOTIF_VALUE_MAX_CHAR];
  ERASE_ARRAY(prob_c, MOTIF_VALUE_MAX_CHAR);
  i = left_bracket + 1;
  while (i < right_bracket && (line[i] == ' ' || line[i] == '\t')) i++;
  for (uint64_t j = i; j < right_bracket; j++) {
    if (line[j] != ' ' && line[j] != '\t') {
      prob_c[k] = line[j]; k++; prev_was_space = 0;
    } else {
      if (!prev_was_space) {
        pos_i++;
        if (pos_i >= MAX_MOTIF_WIDTH) {
          fprintf(stderr, "Error: Motif [%s] has too many columns (max %d).",
            motif->name, MAX_MOTIF_WIDTH); return 1;
        }
        double v;
        if (str_to_double(prob_c, &v)) {
          fprintf(stderr, "Error: Failed to parse count for motif: %s. Bad value: %s",
            motif->name, prob_c); return 1;
        }
        counts[pos_i][row_i] = v;
        ERASE_ARRAY(prob_c, MOTIF_VALUE_MAX_CHAR);
        k = 0;
      }
      prev_was_space = 1;
    }
  }
  if (!prev_was_space) {
    pos_i++;
    if (pos_i >= MAX_MOTIF_WIDTH) {
      fprintf(stderr, "Error: Motif [%s] has too many columns (max %d).",
        motif->name, MAX_MOTIF_WIDTH); return 1;
    }
    double v;
    if (str_to_double(prob_c, &v)) {
      fprintf(stderr, "Error: Failed to parse count for motif: %s. Bad value: %s",
        motif->name, prob_c); return 1;
    }
    counts[pos_i][row_i] = v;
  }
  if (pos_i == (uint64_t)-1) {
    fprintf(stderr, "Error: Motif [%s] has an empty row.", motif->name); return 1;
  }
  pos_i++;
  if (*width_out && *width_out != pos_i) {
    fprintf(stderr, "Error: Motif [%s] has rows with differing widths.", motif->name);
    return 1;
  }
  *width_out = pos_i;
  return 0;
}

static void pcm_to_pwm(motif_t *motif, double counts[MAX_MOTIF_WIDTH][4]) {
  for (uint64_t j = 0; j < motif->size; j++) {
    double col = counts[j][0] + counts[j][1] + counts[j][2] + counts[j][3];
    if (col <= 0.0) col = 1.0;
    for (int i = 0; i < 4; i++) {
      double raw = counts[j][i] / col;
      double adj = (raw * args.nsites + (double)args.pseudocount / 4.0) /
                   (args.nsites + args.pseudocount);
      motif->pwm_probs[j][i] = adj;
    }
    const char lets[4] = {'A','C','G','T'};
    for (int i = 0; i < 4; i++) {
      set_score(motif, lets[i], j, calc_score(counts[j][i] / col, args.bkg[i]));
    }
  }
}

static void read_jaspar(void) {
  motif_info.fmt = FMT_JASPAR;
  char *line = NULL;
  size_t len = 0;
  ssize_t r;
  uint64_t line_num = 0, motif_i = -1, row_i = 0, ready_to_start = 0;
  static double counts[MAX_MOTIF_WIDTH][4];
  uint64_t cur_width = 0;
  while ((r = getline(&line, &len, files.m)) != -1) {
    line_num++;
    if (line[0] == '>') {
      ready_to_start = 1;
      if (motif_i < -1 && row_i != 4) {
        fprintf(stderr, "Error: Motif [%s] has %s rows", motifs[motif_i]->name,
          row_i < 4 ? "too few" : "too many");
        free(line); badexit("");
      }
      if (motif_i < -1) {
        motifs[motif_i]->size = cur_width;
        pcm_to_pwm(motifs[motif_i], counts);
      }
      motif_i++;
      if (add_motif()) { free(line); badexit(""); }
      motifs[motif_i]->file_line_num = line_num;
      parse_jaspar_name(line, motif_i);
      row_i = 0;
      cur_width = 0;
      memset(counts, 0, sizeof(counts));
    } else if (count_nonempty_chars(line) && ready_to_start) {
      row_i++;
      if (add_jaspar_row(motifs[motif_i], line, counts, &cur_width)) {
        free(line); badexit("");
      }
    }
  }
  free(line);
  if (motif_i < -1 && row_i != 4) {
    fprintf(stderr, "Error: Motif [%s] has %s rows", motifs[motif_i]->name,
      row_i < 4 ? "too few" : "too many");
    badexit("");
  }
  if (motif_info.n) {
    motifs[motif_info.n - 1]->size = cur_width;
    pcm_to_pwm(motifs[motif_info.n - 1], counts);
  }
  if (args.v) fprintf(stderr, "Found %'" PRIu64 " JASPAR motif(s).\n", motif_info.n);
}

/* ---- HOCOMOCO reader (PCM rows: rows = positions, cols = ACGT counts) ---- */

static int add_motif_pcm_column(motif_t *motif, const char *line, const uint64_t pos) {
  double probs[4] = {-1.0, -1.0, -1.0, -1.0};
  if (get_line_probs(motif, line, probs, 4)) return 1;
  double pcm_sum = probs[0] + probs[1] + probs[2] + probs[3];
  if (pcm_sum < 0.99) {
    fprintf(stderr, "Error: Motif [%s] PCM row adds up to less than 1", motif->name);
    return 1;
  }
  if (pos >= MAX_MOTIF_WIDTH) {
    fprintf(stderr, "Error: Motif [%s] exceeds max width (%d).", motif->name, MAX_MOTIF_WIDTH);
    return 1;
  }
  for (int i = 0; i < 4; i++) {
    double raw = probs[i] / pcm_sum;
    double adj = (raw * args.nsites + (double)args.pseudocount / 4.0) /
                 (args.nsites + args.pseudocount);
    motif->pwm_probs[pos][i] = adj;
  }
  /* HOCOMOCO uses pseudocount addition on the raw counts (preserves yamscan behavior). */
  double adj_probs[4];
  for (int i = 0; i < 4; i++) adj_probs[i] = probs[i] + args.pseudocount / 4.0;
  set_score(motif, 'A', pos, calc_score(adj_probs[0] / pcm_sum, args.bkg[0]));
  set_score(motif, 'C', pos, calc_score(adj_probs[1] / pcm_sum, args.bkg[1]));
  set_score(motif, 'G', pos, calc_score(adj_probs[2] / pcm_sum, args.bkg[2]));
  set_score(motif, 'T', pos, calc_score(adj_probs[3] / pcm_sum, args.bkg[3]));
  return 0;
}

static void read_hocomoco(void) {
  motif_info.fmt = FMT_HOCOMOCO;
  char *line = NULL;
  size_t len = 0;
  ssize_t r;
  uint64_t line_num = 0, motif_i = -1, pos_i = 0;
  int ready_to_start = 0;
  while ((r = getline(&line, &len, files.m)) != -1) {
    line_num++;
    if (line[0] == '>') {
      ready_to_start = 1;
      motif_i++;
      if (add_motif()) { free(line); badexit(""); }
      motifs[motif_i]->file_line_num = line_num;
      for (uint64_t i = 1, j = 0; i < MAX_NAME_SIZE; i++) {
        if (line[i] == '\r' || line[i] == '\n' || line[i] == '\0') {
          motifs[motif_i]->name[j] = '\0'; break;
        }
        motifs[motif_i]->name[j] = line[i]; j++;
        if (j == MAX_NAME_SIZE - 1) { motifs[motif_i]->name[j] = '\0'; break; }
      }
      pos_i = 0;
    } else if (count_nonempty_chars(line) && ready_to_start) {
      if (pos_i >= MAX_MOTIF_WIDTH) {
        fprintf(stderr, "Error: Motif [%s] is too large (max=%d).",
          motifs[motif_i]->name, MAX_MOTIF_WIDTH);
        free(line); badexit("");
      }
      if (add_motif_pcm_column(motifs[motif_i], line, pos_i)) { free(line); badexit(""); }
      pos_i++;
      motifs[motif_i]->size = pos_i;
    }
  }
  free(line);
  if (args.v) fprintf(stderr, "Found %'" PRIu64 " HOCOMOCO motif(s).\n", motif_info.n);
}

/* ---- Motif loading entrypoint ---- */

static void load_motifs(void) {
  switch (detect_motif_fmt()) {
    case FMT_MEME:     read_meme();     break;
    case FMT_HOMER:    read_homer();    break;
    case FMT_JASPAR:   read_jaspar();   break;
    case FMT_HOCOMOCO: read_hocomoco(); break;
    case FMT_UNKNOWN:
      badexit("Error: Failed to detect motif format.");
  }
  uint64_t empty_motifs = 0;
  for (uint64_t i = 0; i < motif_info.n; i++) if (!motifs[i]->size) empty_motifs++;
  if (empty_motifs == motif_info.n) badexit("Error: All parsed motifs are empty.");
  if (empty_motifs && args.v)
    fprintf(stderr, "Warning: %'" PRIu64 " empty motif(s) ignored.\n", empty_motifs);
}

/* ---- Consensus motif (single seed from a literal IUPAC string) ---- */

static void add_consensus_motif(const char *consensus) {
  if (add_motif()) badexit("");
  motif_t *m = motifs[0];
  ERASE_ARRAY(m->name, MAX_NAME_SIZE);
  uint64_t i = 0;
  for (; consensus[i] != '\0' && i < MAX_NAME_SIZE - 1; i++) m->name[i] = consensus[i];
  m->name[i] = '\0';
  m->size = i;
  if (m->size == 0) badexit("Error: Consensus is empty.");
  if (m->size > MAX_MOTIF_WIDTH) {
    fprintf(stderr, "Error: Consensus is too large (%" PRIu64 " > max=%d).",
      m->size, MAX_MOTIF_WIDTH);
    badexit("");
  }
  for (uint64_t pos = 0; pos < m->size; pos++) {
    int let_i = consensus2index[(unsigned char)consensus[pos]];
    if (let_i == -1) {
      fprintf(stderr, "Error: Unknown letter in consensus (%c).", consensus[pos]);
      badexit("");
    }
    for (int j = 0; j < 4; j++) m->pwm_probs[pos][j] = consensus2probs[let_i * 4 + j];
  }
  if (args.v) fprintf(stderr, "Built consensus seed: %s (w=%" PRIu64 ")\n", m->name, m->size);
}

/* ---- Sequence loading ---- */

static inline uint64_t standard_base_count(void) {
  return char_counts['A']+char_counts['a']+char_counts['C']+char_counts['c']+
         char_counts['G']+char_counts['g']+char_counts['T']+char_counts['t']+
         char_counts['U']+char_counts['u'];
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
    size_t nl = kseq->name.l;
    if (nl > SEQ_NAME_MAX_CHAR) { kseq_destroy(kseq); badexit("Error: Seq name too large."); }
    set->names[set->n-1] = malloc(nl+1);
    if (!set->names[set->n-1]) { kseq_destroy(kseq); badexit("Error: Failed to alloc seq name."); }
    memcpy(set->names[set->n-1], kseq->name.s, nl);
    set->names[set->n-1][nl] = '\0';
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
  if (set->unknowns == total) badexit("Error: No standard DNA/RNA bases found.");
  if (args.v)
    fprintf(stderr, "%s set: %'" PRIu64 " seq(s), %'" PRIu64 " bp\n",
      label, set->n, total);
}

/* ---- Scoring ---- */

static inline void score_subseq(const motif_t *m, const unsigned char *seq,
                                 const uint64_t off, int *s, const unsigned char *tbl) {
  *s=0;
  for (uint64_t i=0; i<m->size; i++) *s += get_score(m, seq[i+off], i, tbl);
}
static inline void score_subseq_rc(const motif_t *m, const unsigned char *seq,
                                    const uint64_t off, int *s, int *src,
                                    const unsigned char *tbl) {
  *s=0; *src=0;
  for (uint64_t i=0; i<m->size; i++) {
    *s   += get_score(m, seq[i+off], i, tbl);
    *src += get_score_rc(m, seq[i+off], i, tbl);
  }
}

/* ---- Build motif from PPM ---- */

static void convert_ppm_to_motif(motif_t *motif, int ppm[][4], const uint64_t w,
                                  const int nsites) {
  motif->nsites_actual = (uint64_t)nsites;
  motif->size          = w;

  for (uint64_t j = 0; j < w; j++) {
    double col = (double)(ppm[j][0]+ppm[j][1]+ppm[j][2]+ppm[j][3]);
    if (col <= 0.0) col = 1.0;
    for (int i = 0; i < 4; i++) {
      double raw = ppm[j][i] / col;
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

  fill_cdf(motif);

  motif->threshold  = 0;
  motif->max_score  = 0;
  motif->min_score  = 0;
  set_threshold_pval(motif, args.hit_pval);

  if (motif->threshold == INT_MAX) {
    if (args.w)
      fprintf(stderr, "  threshold@p=%g unreachable for [%s]; retrying at p=0.01\n",
        args.hit_pval, motif->name);
    motif->threshold = 0;
    motif->max_score = 0;
    motif->min_score = 0;
    set_threshold_pval(motif, 0.01);
  }
}

/* ---- Refinement pass (with optional flank extension) ---- */

/* Score positives with current PWM threshold (over the seed's current `size`),
   gather count matrix over the extended window [off-flank, off+size+flank).
   Hits where the extended window walks off either end are skipped entirely.
   Returns new nsites (>= MIN_REFINE_HITS) or 0 if bailed. */
static int refine_motif_ext(motif_t *motif, const int flank) {
  if (motif->threshold == INT_MAX) return 0;
  const uint64_t w_in  = motif->size;
  const uint64_t w_out = w_in + 2 * (uint64_t)flank;
  if (w_out > MAX_MOTIF_WIDTH) {
    fprintf(stderr, "Error: Extended width (%" PRIu64 ") exceeds max (%d) for motif [%s].\n",
      w_out, MAX_MOTIF_WIDTH, motif->name);
    return 0;
  }
  int ppm[MAX_MOTIF_WIDTH][4];
  memset(ppm, 0, sizeof(ppm));
  int nsites = 0;
  const int thr = motif->threshold - 1;
  const unsigned char *tbl = args.mask ? char2maskindex : char2index;

  for (uint64_t si = 0; si < pos_set.n; si++) {
    const unsigned char *seq = pos_set.seqs[si];
    uint64_t seqlen = pos_set.sizes[si];
    if (seqlen < w_in) continue;

    for (uint64_t off = 0; off <= seqlen - w_in; off++) {
      int score, src;
      if (args.scan_rc) {
        score_subseq_rc(motif, seq, off, &score, &src, tbl);
        if (score > thr) {
          /* Boundary policy: skip whole hit if extended window runs off the seq */
          if ((int64_t)off >= flank &&
              off + w_in + (uint64_t)flank <= seqlen) {
            for (uint64_t i = 0; i < w_out; i++) {
              uint8_t idx = tbl[seq[off - flank + i]];
              if (idx < 4) ppm[i][idx]++;
            }
            nsites++;
          }
        }
        if (src > thr) {
          if ((int64_t)off >= flank &&
              off + w_in + (uint64_t)flank <= seqlen) {
            for (uint64_t i = 0; i < w_out; i++) {
              uint8_t idx = tbl[seq[off - flank + w_out - 1 - i]];
              if (idx < 4) ppm[i][comp4[idx]]++;
            }
            nsites++;
          }
        }
      } else {
        score_subseq(motif, seq, off, &score, tbl);
        if (score > thr) {
          if ((int64_t)off >= flank &&
              off + w_in + (uint64_t)flank <= seqlen) {
            for (uint64_t i = 0; i < w_out; i++) {
              uint8_t idx = tbl[seq[off - flank + i]];
              if (idx < 4) ppm[i][idx]++;
            }
            nsites++;
          }
        }
      }
    }
  }

  if (nsites < MIN_REFINE_HITS) {
    if (args.v)
      fprintf(stderr, "  (below MIN_REFINE_HITS=%d; bailing)\n", MIN_REFINE_HITS);
    return 0;
  }
  convert_ppm_to_motif(motif, ppm, w_out, nsites);
  return nsites;
}

/* ---- Build CDF/threshold from already-populated pwm_probs ---- */

/* Build integer PWM scores + CDF + threshold from the motif's pwm_probs.
   We round-trip through a synthetic count matrix at scale SCALE to reuse
   convert_ppm_to_motif. That call would otherwise overwrite pwm_probs with
   pseudocount-adjusted probabilities computed against the synthetic nsites
   (which mismatches the real refinement nsites and pollutes seed IC reporting
   plus the -r 0 output motif). Save and restore pwm_probs around the call. */
static void rebuild_motif_from_probs(motif_t *motif) {
  int ppm[MAX_MOTIF_WIDTH][4];
  memset(ppm, 0, sizeof(ppm));
  const int SCALE = 1000;
  double saved_probs[MAX_MOTIF_WIDTH][4];
  memcpy(saved_probs, motif->pwm_probs, sizeof(saved_probs));
  for (uint64_t j = 0; j < motif->size; j++) {
    for (int i = 0; i < 4; i++) {
      double v = motif->pwm_probs[j][i] * SCALE;
      ppm[j][i] = (int)(v + 0.5);
    }
  }
  convert_ppm_to_motif(motif, ppm, motif->size, SCALE);
  memcpy(motif->pwm_probs, saved_probs, sizeof(saved_probs));
}

/* ---- Output: consensus + IC trim ---- */

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

static void build_consensus(const motif_t *m, char *buf) {
  for (uint64_t j = 0; j < m->size; j++) buf[j] = consensus_char(m->pwm_probs[j]);
  buf[m->size] = '\0';
}

static double column_ic(const double probs[4]) {
  double ic = 2.0;
  for (int i = 0; i < 4; i++) {
    if (probs[i] > 0) ic += probs[i] * log2(probs[i]);
  }
  return ic < 0 ? 0 : ic;
}

static double motif_total_ic(const motif_t *m) {
  double sum = 0;
  for (uint64_t j = 0; j < m->size; j++) sum += column_ic(m->pwm_probs[j]);
  return sum;
}

/* Count hits without mutating the motif. Returns total + and - strand hits
   above the motif's current threshold. */
static uint64_t count_hits(const motif_t *m) {
  if (m->threshold == INT_MAX) return 0;
  uint64_t n = 0;
  const int thr = m->threshold - 1;
  const unsigned char *tbl = args.mask ? char2maskindex : char2index;
  for (uint64_t si = 0; si < pos_set.n; si++) {
    const unsigned char *seq = pos_set.seqs[si];
    uint64_t seqlen = pos_set.sizes[si];
    if (seqlen < m->size) continue;
    for (uint64_t off = 0; off <= seqlen - m->size; off++) {
      int score, src;
      if (args.scan_rc) {
        score_subseq_rc(m, seq, off, &score, &src, tbl);
        if (score > thr) n++;
        if (src   > thr) n++;
      } else {
        score_subseq(m, seq, off, &score, tbl);
        if (score > thr) n++;
      }
    }
  }
  return n;
}

static void ic_trim_motif(motif_t *m, double min_ic) {
  uint64_t left = 0;
  while (left < m->size && column_ic(m->pwm_probs[left]) < min_ic) left++;
  uint64_t right = m->size;
  while (right > left && column_ic(m->pwm_probs[right - 1]) < min_ic) right--;
  uint64_t new_size = right - left;
  if (new_size < 3) return;
  if (left == 0 && right == m->size) return;
  for (uint64_t i = 0; i < new_size; i++) {
    for (int j = 0; j < 4; j++) m->pwm_probs[i][j] = m->pwm_probs[i + left][j];
  }
  m->size = new_size;
}

/* ---- Output: MEME v4 ---- */

static void write_meme(int argc, char **argv) {
  FILE *f = files.o;
  char cons[MAX_MOTIF_WIDTH + 2];

  /* fprintf(f, "##yamref v%s [ ", YAMTK_VERSION); */
  /* for (int i = 1; i < argc; i++) fprintf(f, "%s ", argv[i]); */
  /* fprintf(f, "]\n"); */

  fprintf(f, "MEME version 4\n\n");
  fprintf(f, "ALPHABET= ACGT\n\n");
  fprintf(f, "strands: %s\n\n", args.scan_rc ? "+ -" : "+");
  fprintf(f, "Background letter frequencies:\n");
  fprintf(f, "A %.6f C %.6f G %.6f T %.6f\n\n",
    args.bkg[0], args.bkg[1], args.bkg[2], args.bkg[3]);

  for (uint64_t ri = 0; ri < motif_info.n; ri++) {
    motif_t *m = motifs[ri];
    if (m->dropped) continue;
    build_consensus(m, cons);
    fprintf(f, "MOTIF %s %s\n\n", m->name, cons);
    fprintf(f, "letter-probability matrix: alength= 4 w= %" PRIu64 " nsites= %" PRIu64 "\n",
      m->size, m->nsites_actual ? m->nsites_actual : (uint64_t)args.nsites);
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
    "Usage:  yamtk ref [options] [ -m motifs.txt | -1 CONSENSUS ] -i positives.fa[.gz]\n"
    "\n"
    " -m <str>   Seed motif file (MEME/JASPAR/HOMER/HOCOMOCO).\n"
    " -1 <str>   Seed from a single IUPAC consensus string (e.g. CACGTG).\n"
    " -i <str>   Positives FASTA/FASTQ ('-' = stdin).\n"
    " -o <str>   MEME output file (default: stdout).\n"
    " -t <dbl>   Hit-scoring p-value (default: %g).\n"
    " -r <int>   Refinement passes (default: %d; 0 = trim/extend only).\n"
    " -e <int>   Extend by N flanking positions per side on pass 1 (default: 0).\n"
    " -E         Auto-extend: grow flanks 1 column at a time until both sides'\n"
    "            new column IC falls below -I. Implies -T.\n"
    " -T         IC-trim flanks after refinement.\n"
    " -I <dbl>   IC threshold for -E stopping and -T trimming (default: %.2g).\n"
    " -Q         Drop motifs whose refined total IC < seed total IC.\n"
    " -b A,C,G,T Background (default: from MEME bkg or uniform).\n"
    " -p <int>   Pseudocount for PWM generation (default: %d).\n"
    " -R         Disable reverse-strand scoring.\n"
    " -M         Mask lower-case bases (skip scanning at those positions).\n"
    " -g         Show progress bar.\n"
    " -v / -w / -h   Verbose / very-verbose / help.\n"
    , YAMTK_VERSION, YAMTK_YEAR,
    DEFAULT_HIT_PVAL, DEFAULT_R_PASSES, MIN_IC_BITS, DEFAULT_PSEUDOCOUNT
  );
}

/* ---- main_ref ---- */

int main_ref(int argc, char **argv) {
  int opt;
  int has_m = 0, has_i = 0;
  char *user_bkg = NULL;
  char *consensus = NULL;

  struct timespec ts_program;
  clock_gettime(CLOCK_MONOTONIC, &ts_program);

  while ((opt = getopt(argc, argv, "m:1:i:o:t:r:e:ETI:Qb:p:RMgvwh")) != -1) {
    switch (opt) {
      case 'm':
        if (files.m_open) badexit("Error: -m specified more than once.");
        if (consensus) badexit("Error: -m and -1 cannot both be used.");
        has_m = 1;
        files.m = fopen(optarg, "r");
        if (!files.m) {
          fprintf(stderr, "Error: Cannot open -m file: %s\n", optarg); badexit("");
        }
        files.m_open = 1;
        break;
      case '1':
        if (has_m) badexit("Error: -m and -1 cannot both be used.");
        if (consensus) badexit("Error: -1 specified more than once.");
        consensus = optarg;
        break;
      case 'i':
        if (files.i_open) badexit("Error: -i specified more than once.");
        has_i = 1;
        if (optarg[0]=='-' && optarg[1]=='\0') {
          files.i = gzdopen(fileno(stdin), "r");
        } else {
          files.i = gzopen(optarg, "r");
          if (!files.i) {
            fprintf(stderr, "Error: Cannot open -i file: %s\n", optarg); badexit("");
          }
        }
        files.i_open = 1;
        break;
      case 'o':
        out_path = optarg;
        break;
      case 't':
        if (str_to_double(optarg, &args.hit_pval)) badexit("Error: Failed to parse -t value.");
        if (args.hit_pval <= 0.0 || args.hit_pval > 1.0) badexit("Error: -t must be in (0,1].");
        break;
      case 'r':
        if (str_to_int(optarg, &args.r_passes)) badexit("Error: Failed to parse -r value.");
        if (args.r_passes < 0) badexit("Error: -r must be >= 0.");
        break;
      case 'e':
        if (str_to_int(optarg, &args.extend)) badexit("Error: Failed to parse -e value.");
        if (args.extend < 0) badexit("Error: -e must be >= 0.");
        break;
      case 'E':
        args.auto_extend = 1; break;
      case 'T':
        args.ic_trim = 1; break;
      case 'I':
        if (str_to_double(optarg, &args.ic_min)) badexit("Error: Failed to parse -I value.");
        if (args.ic_min < 0.0 || args.ic_min > 2.0)
          badexit("Error: -I must be in [0, 2].");
        break;
      case 'Q':
        args.quality_gate = 1; break;
      case 'b':
        user_bkg = optarg; args.use_user_bkg = 1; break;
      case 'p':
        if (str_to_int(optarg, &args.pseudocount)) badexit("Error: Failed to parse -p value.");
        if (args.pseudocount < 0) badexit("Error: -p must be >= 0.");
        break;
      case 'R':
        args.scan_rc = 0; break;
      case 'M':
        args.mask = 1; break;
      case 'g':
        args.progress = 1; break;
      case 'v':
        args.v = 1; break;
      case 'w':
        args.v = 1; args.w = 1; break;
      case 'h':
        usage();
        free_motifs(); free_seq_set(&pos_set); free_cdf(); close_files();
        return EXIT_SUCCESS;
      default:
        usage();
        free_motifs(); free_seq_set(&pos_set); free_cdf(); close_files();
        return EXIT_FAILURE;
    }
  }

  setlocale(LC_NUMERIC, "");

  if (!has_m && !consensus) badexit("Error: Missing required argument -m or -1.");
  if (!has_i) badexit("Error: Missing required argument -i.");
  if (args.auto_extend && args.extend > 0)
    badexit("Error: Cannot specify both -e and -E.");

  if (user_bkg) parse_user_bkg(user_bkg);

  /* Open output */
  if (out_path && !(out_path[0]=='-' && out_path[1]=='\0')) {
    files.o = fopen(out_path, "w");
    if (!files.o) {
      fprintf(stderr, "Error: Cannot open -o file: %s\n", out_path); badexit("");
    }
    files.o_open = 1;
  } else {
    files.o = stdout;
  }

  /* Build seed motif(s) */
  time_t t0 = time(NULL);
  if (consensus) {
    if (args.v) fprintf(stderr, "Building consensus seed motif ...\n");
    add_consensus_motif(consensus);
  } else {
    if (args.v) fprintf(stderr, "Loading motifs from %s ...\n", "input");
    load_motifs();
  }
  if (args.v) print_time((uint64_t)difftime(time(NULL), t0), "load motifs");

  /* Validate extension fits */
  if (args.extend > 0) {
    for (uint64_t i = 0; i < motif_info.n; i++) {
      if (motifs[i]->size + 2*(uint64_t)args.extend > MAX_MOTIF_WIDTH) {
        fprintf(stderr, "Error: Motif [%s] width %" PRIu64 " + 2*%d > max (%d).\n",
          motifs[i]->name, motifs[i]->size, args.extend, MAX_MOTIF_WIDTH);
        badexit("");
      }
    }
  }

  /* Load positives */
  t0 = time(NULL);
  if (args.v) fprintf(stderr, "Loading positives ...\n");
  kseq_t *kseq = kseq_init(files.i);
  load_seq_set(kseq, &pos_set, "Positive");
  if (args.v) print_time((uint64_t)difftime(time(NULL), t0), "load positives");

  /* Refine each motif */
  uint64_t n_kept = 0;
  if (args.progress) print_pb(0.0);
  for (uint64_t i = 0; i < motif_info.n; i++) {
    motif_t *m = motifs[i];
    if (args.v && !args.progress)
      fprintf(stderr, "[%s] width=%" PRIu64 ", refining ...\n", m->name, m->size);

    /* Build initial CDF/threshold from seed probs. This is required before any
       scan (refine_motif_ext expects motif->threshold set). */
    rebuild_motif_from_probs(m);

    if (m->threshold == INT_MAX) {
      fprintf(stderr, "Warning: [%s] threshold unreachable; dropping.\n", m->name);
      m->dropped = 1;
      continue;
    }

    /* Snapshot seed quality before refinement */
    uint64_t seed_w   = m->size;
    double   seed_ic  = motif_total_ic(m);
    uint64_t seed_hit = count_hits(m);

    int hits = 1;  /* sentinel for r=0 path */
    int dropped = 0;

    /* Convergence-tracking buffers: pwm_probs of the previous pass (or seed
       before pass 1). Convergence triggers only when width didn't change. */
    double prev_probs[MAX_MOTIF_WIDTH][4];
    int    prev_size;
    memcpy(prev_probs, m->pwm_probs, sizeof(prev_probs));
    prev_size = (int)m->size;

    if (args.auto_extend) {
      /* Auto-extend: grow flanks 1 column per step until both new flanks fall
         below MIN_IC_BITS. Each step does a full refine at width+2. */
      int steps = 0;
      while (m->size + 2 <= MAX_MOTIF_WIDTH) {
        motif_t saved = *m;
        hits = refine_motif_ext(m, 1);
        if (hits == 0) { *m = saved; break; }
        double left_ic  = column_ic(m->pwm_probs[0]);
        double right_ic = column_ic(m->pwm_probs[m->size - 1]);
        if (args.v)
          fprintf(stderr, "  auto-extend step %d: w=%" PRIu64
                          " hits=%d leftIC=%.2f rightIC=%.2f\n",
                  ++steps, m->size, hits, left_ic, right_ic);
        if (left_ic < args.ic_min && right_ic < args.ic_min) {
          *m = saved;
          if (args.v) fprintf(stderr, "  (both flanks < %.2g bits; reverting last step)\n",
                              args.ic_min);
          break;
        }
      }
      /* Refresh convergence baseline after auto-extend (width may have grown). */
      memcpy(prev_probs, m->pwm_probs, sizeof(prev_probs));
      prev_size = (int)m->size;
      /* Additional refinement passes at the final width */
      for (int pass = 1; pass <= args.r_passes; pass++) {
        hits = refine_motif_ext(m, 0);
        if (args.w)
          fprintf(stderr, "  refine pass %d/%d [%s] flank=0: %d hits\n",
                  pass, args.r_passes, m->name, hits);
        if (hits == 0) {
          fprintf(stderr, "Warning: [%s] dropped (insufficient hits at refine pass %d).\n",
            m->name, pass);
          dropped = 1; break;
        }
        if (prev_size == (int)m->size &&
            memcmp(prev_probs, m->pwm_probs, sizeof(prev_probs)) == 0) {
          if (args.v) fprintf(stderr, "  (PWM unchanged; stopping early)\n");
          break;
        }
        memcpy(prev_probs, m->pwm_probs, sizeof(prev_probs));
        prev_size = (int)m->size;
      }
      if (dropped) { m->dropped = 1; continue; }
      /* -E implies trim: asymmetric IC drop on one flank should not leave a
         dead column behind. */
      ic_trim_motif(m, args.ic_min);
    } else {
      for (int pass = 1; pass <= args.r_passes; pass++) {
        int flank = (pass == 1) ? args.extend : 0;
        hits = refine_motif_ext(m, flank);
        if (args.v)
          fprintf(stderr, "  refine pass %d/%d [%s] flank=%d: %d hits\n",
                  pass, args.r_passes, m->name, flank, hits);
        if (hits == 0) {
          fprintf(stderr, "Warning: [%s] dropped (insufficient hits at pass %d).\n",
            m->name, pass);
          dropped = 1;
          break;
        }
        /* Convergence: only meaningful when this pass didn't change the width
           (i.e., flank=0). After pass 1 with -e>0 the width changed, so reset
           the baseline and continue. */
        if (flank == 0 && prev_size == (int)m->size &&
            memcmp(prev_probs, m->pwm_probs, sizeof(prev_probs)) == 0) {
          if (args.v) fprintf(stderr, "  (PWM unchanged; stopping early)\n");
          break;
        }
        memcpy(prev_probs, m->pwm_probs, sizeof(prev_probs));
        prev_size = (int)m->size;
      }
      if (dropped) { m->dropped = 1; continue; }

      /* If -r 0 -e N>0: refresh PPM by gathering counts once without scoring. */
      if (args.r_passes == 0 && args.extend > 0) {
        hits = refine_motif_ext(m, args.extend);
        if (args.v)
          fprintf(stderr, "  extend-only [%s] flank=%d: %d hits\n",
                  m->name, args.extend, hits);
        if (hits == 0) {
          fprintf(stderr, "Warning: [%s] dropped (insufficient hits during extend-only).\n",
            m->name);
          m->dropped = 1;
          continue;
        }
      }

      if (args.ic_trim) ic_trim_motif(m, args.ic_min);
    }

    /* Refinement summary */
    double   refined_ic  = motif_total_ic(m);
    uint64_t refined_hit = m->nsites_actual;
    if (!args.progress)
      fprintf(stderr,
        "[%s] w: %" PRIu64 "->%" PRIu64 " | IC: %.2f->%.2f bits | hits: %"
        PRIu64 "->%" PRIu64 "\n",
        m->name, seed_w, m->size, seed_ic, refined_ic, seed_hit, refined_hit);

    if (args.quality_gate && refined_ic < seed_ic) {
      if (!args.progress)
        fprintf(stderr, "  (refined IC < seed IC; dropping under -Q)\n");
      m->dropped = 1;
      if (args.progress) print_pb((double)(i + 1) / (double)motif_info.n);
      continue;
    }

    n_kept++;
    if (args.progress) print_pb((double)(i + 1) / (double)motif_info.n);
  }
  if (args.progress) fprintf(stderr, "\n");

  if (args.v) fprintf(stderr, "Refined %" PRIu64 "/%" PRIu64 " motif(s).\n",
                       n_kept, motif_info.n);

  if (n_kept == 0) {
    fprintf(stderr, "Error: No motifs survived refinement.\n");
    free_motifs(); free_seq_set(&pos_set); free_cdf(); close_files();
    return EXIT_FAILURE;
  }

  write_meme(argc, argv);

  free_motifs();
  free_seq_set(&pos_set);
  free_cdf();
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
