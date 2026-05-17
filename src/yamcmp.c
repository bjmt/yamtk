/*
 *   yamcmp: Compare query motifs against a target motif database (TOMTOM-style)
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
#include <pthread.h>
#include "version.h"

/* ---- Constants ---- */

#define MAX_NAME_SIZE         ((uint64_t) 256)
#define MAX_MOTIF_WIDTH                       50
#define USER_BKG_MAX_SIZE     ((uint64_t) 256)
#define MEME_BKG_MAX_SIZE     ((uint64_t) 256)
#define MOTIF_VALUE_MAX_CHAR  ((uint64_t) 256)
#define ALLOC_CHUNK_SIZE      ((uint64_t) 256)
#define MIN_BKG_VALUE                      0.001
#define CONSENSUS_THRESHOLD                 0.75
#define DEFAULT_QVALUE_FILTER                0.1
#define DEFAULT_MIN_OVERLAP                    5
#define DEFAULT_PSEUDOCOUNT                  0.0
#define DEFAULT_NSITES                ((uint64_t) 20)  /* TOMTOM default when MEME has no nsites= */
#define PCC_BINS                              51   /* PCC discretization, step = 0.04 */
#define PROGRESS_BAR_WIDTH                    60
#define PROGRESS_BAR_STRING \
  "============================================================"

/* Null mode is runtime-selectable via -e:
     default = empirical (TOMTOM-style; pools actual target-database columns)
     -e      = parametric enumeration of the simplex grid at resolution
               ENUM_GRID_K, weighted by Dirichlet-Multinomial under
               ENUM_DIRICHLET_N. Use when the target db is small (< ~500
               columns) and the empirical histogram becomes too sparse. */
#define ENUM_DIRICHLET_N                       4   /* Dirichlet concentration: alpha = N * bkg.
                                                      N=4, uniform bkg => Dirichlet(1,1,1,1)
                                                      = uniform on the 4-base simplex. */
#define ENUM_GRID_K                            5   /* Simplex grid resolution; #cols = C(K+3,3) */
#define ENUM_N_COLS \
  (((ENUM_GRID_K + 1) * (ENUM_GRID_K + 2) * (ENUM_GRID_K + 3)) / 6)

#define VEC_ADD(VEC, X, VEC_LEN) \
  do { for (uint64_t Xi = 0; Xi < VEC_LEN; Xi++) VEC[Xi] += X; } while (0)
#define VEC_DIV(VEC, X, VEC_LEN) \
  do { for (uint64_t Xi = 0; Xi < VEC_LEN; Xi++) VEC[Xi] /= X; } while (0)
#define VEC_MIN(VEC, MIN_RES, VEC_LEN) \
  do { MIN_RES = VEC[0]; for (uint64_t Xi = 1; Xi < VEC_LEN; Xi++) { if (VEC[Xi] < MIN_RES) MIN_RES = VEC[Xi]; } } while (0)
#define ERASE_ARRAY(ARR, LEN) memset(ARR, 0, sizeof(ARR[0]) * (LEN))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define MAX(x, y) (((x) > (y)) ? (x) : (y))

/* ---- Motif struct ---- */

typedef struct motif_t {
  uint64_t  size;
  uint64_t  file_line_num;
  uint64_t  nsites;                      /* from MEME 'nsites=' / JASPAR PCM total;
                                            DEFAULT_NSITES if unknown */
  char      name[MAX_NAME_SIZE];
  double    ppm[MAX_MOTIF_WIDTH][4];     /* query/target PPM (rows sum to 1) */
  double    ppm_rc[MAX_MOTIF_WIDTH][4];  /* reverse complement PPM */
  /* per-query precomputed null PMFs (set only for query motifs).
     pmf_fwd[start][L] is the convolved PMF of forward-query columns
     [start, start+L), length L*(PCC_BINS-1)+1. Same for pmf_rc with q->ppm_rc.
     L ranges 1..(size-start); pmf_*[start][0] is unused. */
  double ***pmf_fwd;
  double ***pmf_rc;
  uint64_t  pmf_wq;                      /* size at time of cache build (for cleanup) */
  uint64_t  thread;
} motif_t;

enum MOTIF_FMT {
  FMT_MEME     = 1,
  FMT_HOMER    = 2,
  FMT_JASPAR   = 3,
  FMT_HOCOMOCO = 4,
  FMT_UNKNOWN  = 5
};

typedef struct {
  motif_t **motifs;
  uint64_t  n;
  uint64_t  n_alloc;
  int       fmt;
} motif_set_t;

static motif_set_t query_set  = { .motifs = NULL, .n = 0, .n_alloc = 0, .fmt = 0 };
static motif_set_t target_set = { .motifs = NULL, .n = 0, .n_alloc = 0, .fmt = 0 };

/* Parser writes into *cur_set (swapped between queries and targets). */
static motif_set_t *cur_set = NULL;

/* ---- Args ---- */

typedef struct args_t {
  double   bkg[4];
  double   qvalue_filter;
  double   pseudocount;
  uint64_t nsites_override;  /* 0 = no override (use parsed-or-DEFAULT_NSITES per motif) */
  int      nthreads;
  int      min_overlap;
  int      scan_rc;
  int      trim_names;
  int      use_user_bkg;
  int      enum_null;        /* 0 = empirical (default), 1 = enum */
  int      progress;         /* show progress bar (-g) */
  int      v;
  int      w;
  int      metric_idx;
} args_t;

static args_t args = {
  .bkg            = {0.25, 0.25, 0.25, 0.25},
  .qvalue_filter  = DEFAULT_QVALUE_FILTER,
  .pseudocount    = DEFAULT_PSEUDOCOUNT,
  .nsites_override = 0,
  .nthreads       = 1,
  .min_overlap    = DEFAULT_MIN_OVERLAP,
  .scan_rc        = 1,
  .trim_names     = 1,
  .use_user_bkg   = 0,
  .enum_null      = 0,
  .progress       = 0,
  .v              = 0,
  .w              = 0,
  .metric_idx     = 0
};

/* ---- Files ---- */

typedef struct {
  int   q_open, t_open, o_open;
  FILE *q;   /* query motif file */
  FILE *t;   /* target motif file */
  FILE *o;
} files_t;

static files_t files = { 0 };

/* ---- Memory & timing ---- */

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

static void free_pmf_cache(double ***cache, uint64_t wq) {
  if (!cache) return;
  for (uint64_t s = 0; s < wq; s++) {
    if (!cache[s]) continue;
    for (uint64_t L = 1; L <= wq - s; L++) free(cache[s][L]);
    free(cache[s]);
  }
  free(cache);
}

static void free_motif_set(motif_set_t *set) {
  if (!set->motifs) return;
  for (uint64_t i = 0; i < set->n; i++) {
    motif_t *m = set->motifs[i];
    if (!m) continue;
    free_pmf_cache(m->pmf_fwd, m->pmf_wq);
    free_pmf_cache(m->pmf_rc,  m->pmf_wq);
    free(m);
  }
  free(set->motifs);
  set->motifs = NULL; set->n = 0; set->n_alloc = 0;
}

static void close_files(void) {
  if (files.q_open) fclose(files.q);
  if (files.t_open) fclose(files.t);
  if (files.o_open) fclose(files.o);
}

static void badexit(const char *msg) {
  fprintf(stderr, "%s\nRun yamtk cmp -h for usage.\n", msg);
  free_motif_set(&query_set);
  free_motif_set(&target_set);
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

/* ---- Motif management ---- */

static void init_motif(motif_t *motif) {
  ERASE_ARRAY(motif->name, MAX_NAME_SIZE);
  motif->name[0]='m'; motif->name[1]='o'; motif->name[2]='t';
  motif->name[3]='i'; motif->name[4]='f'; motif->name[5]='\0';
  motif->size = 0;
  motif->file_line_num = 0;
  motif->nsites = DEFAULT_NSITES;
  memset(motif->ppm, 0, sizeof(motif->ppm));
  memset(motif->ppm_rc, 0, sizeof(motif->ppm_rc));
  motif->pmf_fwd = NULL;
  motif->pmf_rc  = NULL;
  motif->pmf_wq  = 0;
  motif->thread = 0;
}

static int add_motif(void) {
  cur_set->n++;
  const uint64_t last_i = cur_set->n - 1;
  if (cur_set->n > cur_set->n_alloc) {
    motif_t **tmp_ptr = realloc(cur_set->motifs,
      sizeof(*cur_set->motifs) * cur_set->n_alloc + sizeof(*cur_set->motifs) * ALLOC_CHUNK_SIZE);
    if (!tmp_ptr) {
      fprintf(stderr, "Error: Failed to allocate memory for motifs."); return 1;
    }
    cur_set->motifs = tmp_ptr;
    cur_set->n_alloc += ALLOC_CHUNK_SIZE;
  }
  cur_set->motifs[last_i] = malloc(sizeof(motif_t));
  if (!cur_set->motifs[last_i]) {
    fprintf(stderr, "Error: Failed to allocate memory for motif."); return 1;
  }
  init_motif(cur_set->motifs[last_i]);
  return 0;
}

/* Fill ppm_rc from ppm: ppm_rc[size-1-i][b] = ppm[i][comp(b)]; comp A<->T, C<->G */
static void fill_ppm_rc(motif_t *m) {
  for (uint64_t i = 0; i < m->size; i++) {
    m->ppm_rc[m->size - 1 - i][0] = m->ppm[i][3];  /* A <- T */
    m->ppm_rc[m->size - 1 - i][1] = m->ppm[i][2];  /* C <- G */
    m->ppm_rc[m->size - 1 - i][2] = m->ppm[i][1];  /* G <- C */
    m->ppm_rc[m->size - 1 - i][3] = m->ppm[i][0];  /* T <- A */
  }
}

/* ---- Parser helpers ---- */

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

static int detect_motif_fmt(FILE *f) {
  int jaspar_or_hocomoco = 0, file_fmt = 0, has_tabs = 0;
  char *line = NULL;
  size_t len = 0;
  ssize_t r;
  while ((r = getline(&line, &len, f)) != -1) {
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
            free(line); badexit("Error: yamcmp cannot read HOCOMOCO PWMs.");
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
  rewind(f);
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
  for (int i = 0; i < 4; i++) motif->ppm[pos][i] = probs[i];
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
  motif_t *m = cur_set->motifs[motif_i];
  while (line[i] != '\0' && line[i] != '\r' && line[i] != '\n' && j < MAX_NAME_SIZE - 1) {
    if (line[i] == ' ' || line[i] == '\t') {
      if (!prev_space && j > 0 && j < MAX_NAME_SIZE - 1) m->name[j++] = ' ';
      prev_space = 1;
    } else {
      m->name[j++] = line[i];
      prev_space = 0;
    }
    i++;
  }
  if (j > 0 && m->name[j - 1] == ' ') j--;
  m->name[j] = '\0';
}

static void read_meme(FILE *f, int read_bkg) {
  cur_set->fmt = FMT_MEME;
  char *line = NULL;
  size_t len = 0;
  ssize_t r;
  uint64_t line_num = 0, l_p_m_L = 0, bkg_L = 0, motif_i = -1, pos_i = 0;
  int alph_detected = 0, live_motif = 0;
  int warned_bkg = 0, warned_alph = 0;
  while ((r = getline(&line, &len, f)) != -1) {
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
    } else if (bkg_L && bkg_L == line_num - 1 && read_bkg) {
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
      cur_set->motifs[motif_i]->file_line_num = line_num;
      parse_meme_name(line, motif_i);
      pos_i = 0;
    } else if (check_line_contains(line, "letter-probability matrix\0")) {
      if (pos_i != 0) { free(line); badexit("Error: Malformed MEME motif."); }
      l_p_m_L = line_num;
      live_motif = 1;
      /* Parse nsites= from this header line if present (for pseudocount math). */
      const char *p = strstr(line, "nsites=");
      if (p) {
        p += 7;
        while (*p == ' ' || *p == '\t') p++;
        char *end = NULL;
        unsigned long ns = strtoul(p, &end, 10);
        if (end != p && ns > 0) cur_set->motifs[motif_i]->nsites = (uint64_t)ns;
      }
    } else if (live_motif) {
      if (!count_nonempty_chars(line) || !is_ppm_data_line(line)) {
        live_motif = 0;
      } else if (line_num == (l_p_m_L + pos_i + 1)) {
        if (pos_i >= MAX_MOTIF_WIDTH) {
          fprintf(stderr, "Error: Motif [%s] is too large (max=%d)",
            cur_set->motifs[motif_i]->name, MAX_MOTIF_WIDTH);
          free(line); badexit("");
        }
        if (add_motif_ppm_column(cur_set->motifs[motif_i], line, pos_i)) { free(line); badexit(""); }
        pos_i++;
        cur_set->motifs[motif_i]->size = pos_i;
      } else {
        live_motif = 0;
      }
    }
  }
  free(line);
  if (!cur_set->n) badexit("Error: Failed to detect any motifs in MEME file.");
  if (args.v) fprintf(stderr, "Found %'" PRIu64 " MEME motif(s).\n", cur_set->n);
}

/* ---- HOMER reader ---- */

static void parse_homer_name(const char *line, const uint64_t motif_i) {
  uint64_t name_start = 0, name_end = 0, i = 1, in_between = 0, j = 0;
  motif_t *m = cur_set->motifs[motif_i];
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
  if (!name_start) { m->name[0] = '\0'; return; }
  for (uint64_t k = name_start; k < name_end && j < MAX_NAME_SIZE - 1; k++) {
    m->name[j] = line[k]; j++;
  }
  m->name[j] = '\0';
}

static void read_homer(FILE *f) {
  cur_set->fmt = FMT_HOMER;
  char *line = NULL;
  size_t len = 0;
  ssize_t r;
  uint64_t line_num = 0, motif_i = -1, pos_i = 0;
  int ready_to_start = 0;
  while ((r = getline(&line, &len, f)) != -1) {
    line_num++;
    if (line[0] == '>') {
      ready_to_start = 1;
      motif_i++;
      if (add_motif()) { free(line); badexit(""); }
      cur_set->motifs[motif_i]->file_line_num = line_num;
      parse_homer_name(line, motif_i);
      pos_i = 0;
    } else if (count_nonempty_chars(line) && ready_to_start) {
      if (pos_i >= MAX_MOTIF_WIDTH) {
        fprintf(stderr, "Error: Motif [%s] is too large (max=%d).",
          cur_set->motifs[motif_i]->name, MAX_MOTIF_WIDTH);
        free(line); badexit("");
      }
      if (add_motif_ppm_column(cur_set->motifs[motif_i], line, pos_i)) { free(line); badexit(""); }
      pos_i++;
      cur_set->motifs[motif_i]->size = pos_i;
    }
  }
  free(line);
  if (args.v) fprintf(stderr, "Found %'" PRIu64 " HOMER motif(s).\n", cur_set->n);
}

/* ---- JASPAR reader (count matrix) ---- */

static void parse_jaspar_name(const char *line, const uint64_t motif_i) {
  uint64_t i = 0, j = 1;
  motif_t *m = cur_set->motifs[motif_i];
  while (line[j] != '\r' && line[j] != '\n' && line[j] != '\0' && i < MAX_NAME_SIZE - 1) {
    m->name[i] = line[j];
    i++; j++;
  }
  m->name[i] = '\0';
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

static void counts_to_ppm(motif_t *motif, double counts[MAX_MOTIF_WIDTH][4]) {
  /* Use the first column's count total as nsites (PCM columns should sum
     consistently). */
  if (motif->size > 0) {
    double total = counts[0][0] + counts[0][1] + counts[0][2] + counts[0][3];
    if (total > 0.5) motif->nsites = (uint64_t)(total + 0.5);
  }
  for (uint64_t j = 0; j < motif->size; j++) {
    double col = counts[j][0] + counts[j][1] + counts[j][2] + counts[j][3];
    if (col <= 0.0) col = 1.0;
    for (int i = 0; i < 4; i++) motif->ppm[j][i] = counts[j][i] / col;
  }
}

static void read_jaspar(FILE *f) {
  cur_set->fmt = FMT_JASPAR;
  char *line = NULL;
  size_t len = 0;
  ssize_t r;
  uint64_t line_num = 0, motif_i = -1, row_i = 0, ready_to_start = 0;
  static double counts[MAX_MOTIF_WIDTH][4];
  uint64_t cur_width = 0;
  while ((r = getline(&line, &len, f)) != -1) {
    line_num++;
    if (line[0] == '>') {
      ready_to_start = 1;
      if (motif_i < -1 && row_i != 4) {
        fprintf(stderr, "Error: Motif [%s] has %s rows", cur_set->motifs[motif_i]->name,
          row_i < 4 ? "too few" : "too many");
        free(line); badexit("");
      }
      if (motif_i < -1) {
        cur_set->motifs[motif_i]->size = cur_width;
        counts_to_ppm(cur_set->motifs[motif_i], counts);
      }
      motif_i++;
      if (add_motif()) { free(line); badexit(""); }
      cur_set->motifs[motif_i]->file_line_num = line_num;
      parse_jaspar_name(line, motif_i);
      row_i = 0;
      cur_width = 0;
      memset(counts, 0, sizeof(counts));
    } else if (count_nonempty_chars(line) && ready_to_start) {
      row_i++;
      if (add_jaspar_row(cur_set->motifs[motif_i], line, counts, &cur_width)) {
        free(line); badexit("");
      }
    }
  }
  free(line);
  if (motif_i < -1 && row_i != 4) {
    fprintf(stderr, "Error: Motif [%s] has %s rows", cur_set->motifs[motif_i]->name,
      row_i < 4 ? "too few" : "too many");
    badexit("");
  }
  if (cur_set->n) {
    cur_set->motifs[cur_set->n - 1]->size = cur_width;
    counts_to_ppm(cur_set->motifs[cur_set->n - 1], counts);
  }
  if (args.v) fprintf(stderr, "Found %'" PRIu64 " JASPAR motif(s).\n", cur_set->n);
}

/* ---- HOCOMOCO reader ---- */

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
  for (int i = 0; i < 4; i++) motif->ppm[pos][i] = probs[i] / pcm_sum;
  return 0;
}

static void read_hocomoco(FILE *f) {
  cur_set->fmt = FMT_HOCOMOCO;
  char *line = NULL;
  size_t len = 0;
  ssize_t r;
  uint64_t line_num = 0, motif_i = -1, pos_i = 0;
  int ready_to_start = 0;
  while ((r = getline(&line, &len, f)) != -1) {
    line_num++;
    if (line[0] == '>') {
      ready_to_start = 1;
      motif_i++;
      if (add_motif()) { free(line); badexit(""); }
      cur_set->motifs[motif_i]->file_line_num = line_num;
      for (uint64_t i = 1, j = 0; i < MAX_NAME_SIZE; i++) {
        if (line[i] == '\r' || line[i] == '\n' || line[i] == '\0') {
          cur_set->motifs[motif_i]->name[j] = '\0'; break;
        }
        cur_set->motifs[motif_i]->name[j] = line[i]; j++;
        if (j == MAX_NAME_SIZE - 1) { cur_set->motifs[motif_i]->name[j] = '\0'; break; }
      }
      pos_i = 0;
    } else if (count_nonempty_chars(line) && ready_to_start) {
      if (pos_i >= MAX_MOTIF_WIDTH) {
        fprintf(stderr, "Error: Motif [%s] is too large (max=%d).",
          cur_set->motifs[motif_i]->name, MAX_MOTIF_WIDTH);
        free(line); badexit("");
      }
      if (add_motif_pcm_column(cur_set->motifs[motif_i], line, pos_i)) { free(line); badexit(""); }
      pos_i++;
      cur_set->motifs[motif_i]->size = pos_i;
    }
  }
  free(line);
  if (args.v) fprintf(stderr, "Found %'" PRIu64 " HOCOMOCO motif(s).\n", cur_set->n);
}

static void trim_motif_name(motif_t *motif) {
  for (uint64_t i = 0; i < MAX_NAME_SIZE; i++) {
    if (motif->name[i] == ' ' || motif->name[i] == '\t' || motif->name[i] == '\0') {
      motif->name[i] = '\0';
      break;
    }
  }
}

static void load_motifs_from(FILE *f, motif_set_t *into, const char *which, int read_bkg) {
  cur_set = into;
  switch (detect_motif_fmt(f)) {
    case FMT_MEME:     read_meme(f, read_bkg); break;
    case FMT_HOMER:    read_homer(f);          break;
    case FMT_JASPAR:   read_jaspar(f);         break;
    case FMT_HOCOMOCO: read_hocomoco(f);       break;
    case FMT_UNKNOWN:
      fprintf(stderr, "Error: Failed to detect %s motif format.", which);
      badexit("");
  }
  if (args.trim_names)
    for (uint64_t i = 0; i < into->n; i++) trim_motif_name(into->motifs[i]);
  uint64_t empty_motifs = 0;
  for (uint64_t i = 0; i < into->n; i++) if (!into->motifs[i]->size) empty_motifs++;
  if (empty_motifs == into->n) {
    fprintf(stderr, "Error: All parsed %s motifs are empty.", which); badexit("");
  }
  if (empty_motifs && args.v)
    fprintf(stderr, "Warning: %'" PRIu64 " empty %s motif(s) ignored.\n", empty_motifs, which);
  /* If user supplied -N, force every motif's nsites to that value
     (overrides any value parsed from MEME / JASPAR). */
  if (args.nsites_override > 0) {
    for (uint64_t i = 0; i < into->n; i++) into->motifs[i]->nsites = args.nsites_override;
  }
  /* Apply TOMTOM-style pseudocount smoothing (no-op if pseudocount == 0):
       new_p[j] = (old_p[j] * nsites + pseudo * bkg[j]) / (nsites + pseudo) */
  if (args.pseudocount > 0.0) {
    for (uint64_t i = 0; i < into->n; i++) {
      motif_t *m = into->motifs[i];
      double denom = (double)m->nsites + args.pseudocount;
      for (uint64_t p = 0; p < m->size; p++) {
        for (int j = 0; j < 4; j++) {
          m->ppm[p][j] = (m->ppm[p][j] * (double)m->nsites
                          + args.pseudocount * args.bkg[j]) / denom;
        }
      }
    }
  }
  for (uint64_t i = 0; i < into->n; i++) fill_ppm_rc(into->motifs[i]);
  cur_set = NULL;
}

/* ---- Consensus ---- */

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

static void build_consensus_range(const double ppm[MAX_MOTIF_WIDTH][4],
                                   const uint64_t start, const uint64_t len, char *buf) {
  for (uint64_t j = 0; j < len; j++) buf[j] = consensus_char(ppm[start + j]);
  buf[len] = '\0';
}

/* ---- Column-similarity metrics ---- */

/* Pearson correlation of (q - 0.25) and (t - 0.25). Range [-1, 1].
   Falls back to 0 if either vector has zero variance. */
static double pcc_column(const double q[4], const double t[4]) {
  double qd[4], td[4];
  double sqq = 0, stt = 0, sqt = 0;
  for (int i = 0; i < 4; i++) { qd[i] = q[i] - 0.25; td[i] = t[i] - 0.25; }
  for (int i = 0; i < 4; i++) {
    sqq += qd[i] * qd[i];
    stt += td[i] * td[i];
    sqt += qd[i] * td[i];
  }
  if (sqq < 1e-12 || stt < 1e-12) return 0.0;
  return sqt / sqrt(sqq * stt);
}

typedef struct {
  const char *name;
  double (*column)(const double q[4], const double t[4]);
  double score_min;
  double score_max;
  int    higher_is_better;
} metric_t;

static const metric_t METRICS[] = {
  { "pcc", pcc_column, -1.0, 1.0, 1 },
  /* future: { "ed", ... }, { "sw", ... }, { "allr", ... } */
};
#define N_METRICS ((int)(sizeof(METRICS)/sizeof(METRICS[0])))

static int find_metric(const char *name) {
  for (int i = 0; i < N_METRICS; i++)
    if (strcmp(METRICS[i].name, name) == 0) return i;
  return -1;
}

/* ---- Discretization ---- */

static inline int score_to_bin(double s, const metric_t *m) {
  double t = (s - m->score_min) / (m->score_max - m->score_min);  /* [0, 1] */
  if (t < 0) t = 0; else if (t > 1) t = 1;
  int bin = (int)floor(t * (PCC_BINS - 1) + 0.5);
  if (bin < 0) bin = 0; else if (bin > PCC_BINS - 1) bin = PCC_BINS - 1;
  return bin;
}

/* Sum-of-L column scores → bin index in [0, L*(PCC_BINS-1)]. */
static inline int sum_score_to_bin(double s, int L, const metric_t *m) {
  double per_col = (s / L - m->score_min) / (m->score_max - m->score_min);  /* avg in [0, 1] */
  if (per_col < 0) per_col = 0; else if (per_col > 1) per_col = 1;
  int bin = (int)floor(per_col * (PCC_BINS - 1) * L + 0.5);
  int hi = L * (PCC_BINS - 1);
  if (bin < 0) bin = 0; else if (bin > hi) bin = hi;
  return bin;
}

/* ---- Per-query column null PMFs and L-fold convolution cache ---- */

/* --- Enum-mode tables (built once if args.enum_null) --- */

typedef struct {
  double col[4];
  double weight;   /* normalized so sum over all enum cols == 1 */
} enum_col_t;

static enum_col_t enum_cols[ENUM_N_COLS];

static void build_enum_cols(void) {
  double alpha[4];
  for (int i = 0; i < 4; i++) alpha[i] = (double)ENUM_DIRICHLET_N * args.bkg[i];
  if (alpha[0] < 0.1) alpha[0] = 0.1;
  if (alpha[1] < 0.1) alpha[1] = 0.1;
  if (alpha[2] < 0.1) alpha[2] = 0.1;
  if (alpha[3] < 0.1) alpha[3] = 0.1;
  double alpha_sum = alpha[0] + alpha[1] + alpha[2] + alpha[3];
  /* Dirichlet-Multinomial log PMF:
       log P(c | alpha, K) = log(K!) - sum log(c_i!)
                           + sum [log Gamma(c_i + alpha_i) - log Gamma(alpha_i)]
                           - [log Gamma(K + sum alpha) - log Gamma(sum alpha)] */
  int idx = 0;
  for (int a = 0; a <= ENUM_GRID_K; a++)
    for (int b = 0; b <= ENUM_GRID_K - a; b++)
      for (int c = 0; c <= ENUM_GRID_K - a - b; c++) {
        int d = ENUM_GRID_K - a - b - c;
        int counts[4] = {a, b, c, d};
        enum_cols[idx].col[0] = (double)a / ENUM_GRID_K;
        enum_cols[idx].col[1] = (double)b / ENUM_GRID_K;
        enum_cols[idx].col[2] = (double)c / ENUM_GRID_K;
        enum_cols[idx].col[3] = (double)d / ENUM_GRID_K;
        double lw = lgamma((double)ENUM_GRID_K + 1.0);
        for (int i = 0; i < 4; i++) {
          lw -= lgamma((double)counts[i] + 1.0);
          lw += lgamma((double)counts[i] + alpha[i]) - lgamma(alpha[i]);
        }
        lw -= lgamma((double)ENUM_GRID_K + alpha_sum) - lgamma(alpha_sum);
        enum_cols[idx].weight = lw;  /* store log temporarily */
        idx++;
      }
  /* Normalize: exp-then-sum-then-divide (numerically stable via max-shift) */
  double max_lw = enum_cols[0].weight;
  for (int i = 1; i < ENUM_N_COLS; i++)
    if (enum_cols[i].weight > max_lw) max_lw = enum_cols[i].weight;
  double zsum = 0.0;
  for (int i = 0; i < ENUM_N_COLS; i++) {
    enum_cols[i].weight = exp(enum_cols[i].weight - max_lw);
    zsum += enum_cols[i].weight;
  }
  for (int i = 0; i < ENUM_N_COLS; i++) enum_cols[i].weight /= zsum;
}

/* --- Per-column null PMF: empirical (default) or enum (-e) --- */

static void build_col_pmf_empirical(const double q[4], double pmf_out[PCC_BINS]) {
  for (int i = 0; i < PCC_BINS; i++) pmf_out[i] = 0.0;
  const metric_t *m = &METRICS[args.metric_idx];
  uint64_t total = 0;
  for (uint64_t ti = 0; ti < target_set.n; ti++) {
    motif_t *t = target_set.motifs[ti];
    for (uint64_t pos = 0; pos < t->size; pos++) {
      double v = m->column(q, t->ppm[pos]);
      pmf_out[score_to_bin(v, m)] += 1.0;
      total++;
    }
  }
  if (total == 0) return;
  double inv = 1.0 / (double)total;
  for (int i = 0; i < PCC_BINS; i++) pmf_out[i] *= inv;
}

static void build_col_pmf_enum(const double q[4], double pmf_out[PCC_BINS]) {
  for (int i = 0; i < PCC_BINS; i++) pmf_out[i] = 0.0;
  const metric_t *m = &METRICS[args.metric_idx];
  for (int i = 0; i < ENUM_N_COLS; i++) {
    double v = m->column(q, enum_cols[i].col);
    pmf_out[score_to_bin(v, m)] += enum_cols[i].weight;
  }
}

static void build_col_pmf(const double q[4], double pmf_out[PCC_BINS]) {
  if (args.enum_null) build_col_pmf_enum(q, pmf_out);
  else                build_col_pmf_empirical(q, pmf_out);
}

/* In-place convolution: result = pmf_a * pmf_b, where pmf_a has length len_a
   and pmf_b has length PCC_BINS. result is written into out (length len_a+PCC_BINS-1). */
static void convolve_with_col(const double *pmf_a, int len_a,
                               const double pmf_b[PCC_BINS], double *out) {
  int len_out = len_a + PCC_BINS - 1;
  for (int i = 0; i < len_out; i++) out[i] = 0.0;
  for (int i = 0; i < len_a; i++) {
    double a = pmf_a[i];
    if (a == 0.0) continue;
    for (int j = 0; j < PCC_BINS; j++) out[i + j] += a * pmf_b[j];
  }
}

/* Build cache[start][L] = PMF of convolved query columns [start, start+L)
   for one orientation (qppm = q->ppm or q->ppm_rc). PMF lengths are
   L*(PCC_BINS-1)+1. Indices: start in [0, wq), L in [1, wq-start].
   Built incrementally per start: cache[start][L] = cache[start][L-1] ⊛ col_pmf[start+L-1]. */
static int build_pmf_cache(const motif_t *q, const double (*qppm)[4],
                            double ****cache_out) {
  uint64_t wq = q->size;
  /* Per-column PMFs for this orientation */
  double col_pmf[MAX_MOTIF_WIDTH][PCC_BINS];
  for (uint64_t i = 0; i < wq; i++) build_col_pmf(qppm[i], col_pmf[i]);

  double ***cache = calloc(wq, sizeof(double **));
  if (!cache) return 1;
  for (uint64_t start = 0; start < wq; start++) {
    uint64_t max_L = wq - start;
    cache[start] = calloc(max_L + 1, sizeof(double *));  /* L in [1, max_L]; [0] unused */
    if (!cache[start]) goto fail;
    /* L = 1: copy the column PMF directly */
    cache[start][1] = malloc(PCC_BINS * sizeof(double));
    if (!cache[start][1]) goto fail;
    memcpy(cache[start][1], col_pmf[start], PCC_BINS * sizeof(double));
    /* L = 2..max_L: convolve previous window with next column */
    for (uint64_t L = 2; L <= max_L; L++) {
      int prev_len = (int)(L - 1) * (PCC_BINS - 1) + 1;
      int new_len  = (int)L       * (PCC_BINS - 1) + 1;
      cache[start][L] = malloc(new_len * sizeof(double));
      if (!cache[start][L]) goto fail;
      convolve_with_col(cache[start][L - 1], prev_len, col_pmf[start + L - 1], cache[start][L]);
    }
  }
  *cache_out = cache;
  return 0;
fail:
  free_pmf_cache(cache, wq);
  return 1;
}

/* For one query: build per-(start, L) PMF caches for both orientations
   (only forward if !args.scan_rc). Stored on the motif. */
static int build_query_pmfs(motif_t *q) {
  q->pmf_wq = q->size;
  if (build_pmf_cache(q, q->ppm, &q->pmf_fwd)) return 1;
  if (args.scan_rc) {
    if (build_pmf_cache(q, q->ppm_rc, &q->pmf_rc)) return 1;
  }
  return 0;
}

/* Tail probability P(S' >= observed_bin) from a pmf of length L*(PCC_BINS-1)+1. */
static double tail_prob(const double *pmf, int L, int observed_bin) {
  int len = L * (PCC_BINS - 1) + 1;
  if (observed_bin < 0) observed_bin = 0;
  if (observed_bin >= len) return 0.0;
  /* For "higher is better" metrics, tail is upper. Sum from observed to end. */
  double s = 0.0;
  for (int i = observed_bin; i < len; i++) s += pmf[i];
  if (s < 0.0) s = 0.0;
  if (s > 1.0) s = 1.0;
  return s;
}

/* ---- Alignment scoring ---- */

/* Score a single alignment: query (already with chosen orientation, given by q_ppm)
   slid against target (always +) at offset d_in_q (the query column corresponding
   to target column 0). Returns the sum of column metric scores and L (overlap). */
static double align_score(const double q_ppm[MAX_MOTIF_WIDTH][4], int wq,
                           const double t_ppm[MAX_MOTIF_WIDTH][4], int wt,
                           int d_in_q, int *out_L, int *out_q_start, int *out_t_start) {
  /* Target column j corresponds to query column (d_in_q + j).
     Overlap: max(0, d_in_q) <= j_q (in q-coords) < min(wq, wt + d_in_q)
     Equivalently in target indexing: j_t in [max(0, -d_in_q), min(wt, wq - d_in_q)) */
  int t_start = (d_in_q >= 0) ? 0 : -d_in_q;
  int t_end   = (wq - d_in_q < wt) ? wq - d_in_q : wt;
  int L = t_end - t_start;
  if (L <= 0) { *out_L = 0; return 0.0; }
  int q_start = d_in_q + t_start;
  double s = 0.0;
  const metric_t *m = &METRICS[args.metric_idx];
  for (int j = 0; j < L; j++) {
    s += m->column(q_ppm[q_start + j], t_ppm[t_start + j]);
  }
  *out_L = L;
  *out_q_start = q_start;
  *out_t_start = t_start;
  return s;
}

/* Iterate offsets and orientations; return best alignment for (q, t).
   Bonferroni denominator = number of (offset × orientation) evaluations
   with overlap >= min_overlap. */
typedef struct {
  double score;
  int    L;
  int    offset;  /* d_in_q: query column corresponding to target column 0 */
  int    orientation;  /* 0 = +, 1 = - (RC query) */
  int    q_start_oriented; /* start col index in the oriented query (for PMF lookup) */
  int    q_start_in_orig;  /* offset back to the ORIGINAL (forward) query coords (for output) */
  int    t_start;
  int    n_tested;
} align_result_t;

static void best_alignment(motif_t *q, motif_t *t, align_result_t *out) {
  const metric_t *m = &METRICS[args.metric_idx];
  int wq = (int)q->size;
  int wt = (int)t->size;
  int mo = args.min_overlap;
  if (mo > wq) mo = wq;
  if (mo > wt) mo = wt;
  /* d_in_q ranges so that overlap >= mo. d in [mo - wt, wq - mo]. */
  int d_lo = mo - wt;
  int d_hi = wq - mo;
  out->score = -INFINITY;
  out->L = 0;
  out->offset = 0;
  out->orientation = 0;
  out->q_start_oriented = 0;
  out->q_start_in_orig = 0;
  out->t_start = 0;
  out->n_tested = 0;
  int orientations = args.scan_rc ? 2 : 1;
  for (int orient = 0; orient < orientations; orient++) {
    const double (*qppm)[4] = (orient == 0) ? q->ppm : q->ppm_rc;
    for (int d = d_lo; d <= d_hi; d++) {
      int L, q_start, t_start;
      double s = align_score(qppm, wq, t->ppm, wt, d, &L, &q_start, &t_start);
      if (L < mo) continue;
      out->n_tested++;
      int better;
      if (m->higher_is_better) better = s > out->score;
      else                      better = s < out->score || out->n_tested == 1;
      if (better || out->n_tested == 1) {
        out->score = s;
        out->L = L;
        out->offset = d;
        out->orientation = orient;
        out->q_start_oriented = q_start;
        /* Translate q_start (in oriented query) back to original-query coords.
           For RC, oriented column i corresponds to original column (wq - 1 - i). */
        if (orient == 0) {
          out->q_start_in_orig = q_start;
        } else {
          out->q_start_in_orig = wq - 1 - (q_start + L - 1);
        }
        out->t_start = t_start;
      }
    }
  }
}

/* ---- Result types and sort ---- */

typedef struct {
  uint64_t target_i;
  double   score;
  double   pvalue;
  double   qvalue;
  int      L;
  int      offset;        /* d_in_q: target column 0 maps to query column `offset` */
  int      orientation;
  int      q_start_in_orig;
  int      t_start;
} pair_result_t;

static int cmp_pair(const void *a, const void *b) {
  const pair_result_t *ra = (const pair_result_t *)a, *rb = (const pair_result_t *)b;
  if (ra->qvalue != rb->qvalue) return (ra->qvalue > rb->qvalue) - (ra->qvalue < rb->qvalue);
  if (ra->pvalue != rb->pvalue) return (ra->pvalue > rb->pvalue) - (ra->pvalue < rb->pvalue);
  return (ra->target_i > rb->target_i) - (ra->target_i < rb->target_i);
}

/* ---- BH q-values (per query) ---- */

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

/* ---- Per-query worker output storage ---- */

static pair_result_t **query_results = NULL;   /* [q_i] -> array of length target_set.n */
static pthread_mutex_t pb_lock = PTHREAD_MUTEX_INITIALIZER;
static uint64_t pb_counter = 0;                /* progress: completed queries */

static void print_pb(const double prog) {
  const int left = prog * PROGRESS_BAR_WIDTH;
  const int right = PROGRESS_BAR_WIDTH - left;
  fprintf(stderr, "\r[%.*s%*s] %3d%%", left, PROGRESS_BAR_STRING, right, "",
      (int) (prog * 100.0));
  fflush(stderr);
}

/* ---- Thread worker ---- */

static void *cmp_sub_process(void *arg) {
  uint64_t thread_i = *(uint64_t *)arg;
  free(arg);
  const metric_t *m = &METRICS[args.metric_idx];
  for (uint64_t qi = 0; qi < query_set.n; qi++) {
    motif_t *q = query_set.motifs[qi];
    if (q->thread != thread_i) continue;
    if (args.w && !args.progress)
      fprintf(stderr, "    Building null PMFs for query: %s\n", q->name);
    if (build_query_pmfs(q)) {
      fprintf(stderr, "Error: Failed to build PMFs for query %s.", q->name);
      exit(EXIT_FAILURE);
    }
    pair_result_t *rows = query_results[qi];
    double *pvals = malloc(sizeof(double) * target_set.n);
    double *qvals = malloc(sizeof(double) * target_set.n);
    if (!pvals || !qvals) {
      fprintf(stderr, "Error: Failed to alloc pvals for query %s.", q->name);
      exit(EXIT_FAILURE);
    }
    for (uint64_t ti = 0; ti < target_set.n; ti++) {
      motif_t *t = target_set.motifs[ti];
      align_result_t a;
      best_alignment(q, t, &a);
      double pval = 1.0;
      if (a.n_tested > 0 && a.L > 0 && (uint64_t)a.L <= q->pmf_wq) {
        double ***cache = (a.orientation == 0) ? q->pmf_fwd : q->pmf_rc;
        if (cache && a.q_start_oriented >= 0 &&
            (uint64_t)a.q_start_oriented < q->pmf_wq &&
            cache[a.q_start_oriented] && cache[a.q_start_oriented][a.L]) {
          int observed_bin = sum_score_to_bin(a.score, a.L, m);
          double p_align = tail_prob(cache[a.q_start_oriented][a.L], a.L, observed_bin);
          /* Bonferroni-correct for number of offsets tested */
          pval = p_align * (double)a.n_tested;
          if (pval > 1.0) pval = 1.0;
          if (pval < 0.0) pval = 0.0;
        }
      }
      rows[ti].target_i        = ti;
      rows[ti].score           = a.score;
      rows[ti].pvalue          = pval;
      rows[ti].qvalue          = 1.0;
      rows[ti].L               = a.L;
      rows[ti].offset          = a.offset;
      rows[ti].orientation     = a.orientation;
      rows[ti].q_start_in_orig = a.q_start_in_orig;
      rows[ti].t_start         = a.t_start;
      pvals[ti] = pval;
    }
    bh_qvalues(pvals, target_set.n, qvals);
    for (uint64_t ti = 0; ti < target_set.n; ti++) rows[ti].qvalue = qvals[ti];
    free(pvals); free(qvals);
    /* Free the per-query PMF caches after we're done with this query */
    free_pmf_cache(q->pmf_fwd, q->pmf_wq); q->pmf_fwd = NULL;
    free_pmf_cache(q->pmf_rc,  q->pmf_wq); q->pmf_rc  = NULL;
    q->pmf_wq = 0;
    if (args.progress) {
      pthread_mutex_lock(&pb_lock);
      pb_counter++;
      print_pb((double)pb_counter / (double)query_set.n);
      pthread_mutex_unlock(&pb_lock);
    }
  }
  return NULL;
}

/* ---- Output ---- */

static void write_results(FILE *o, int argc, char **argv) {
  fprintf(o, "##yamcmp v%s [", YAMTK_VERSION);
  for (int i = 0; i < argc; i++) fprintf(o, " %s", argv[i]);
  fprintf(o, " ]\n");
  fprintf(o, "##query_id\ttarget_id\toffset\tp_value\tq_value\tscore\toverlap\t"
             "query_consensus\ttarget_consensus\torientation\n");

  char qcons[MAX_MOTIF_WIDTH + 1];
  char tcons[MAX_MOTIF_WIDTH + 1];
  for (uint64_t qi = 0; qi < query_set.n; qi++) {
    motif_t *q = query_set.motifs[qi];
    pair_result_t *rows = query_results[qi];
    qsort(rows, target_set.n, sizeof(pair_result_t), cmp_pair);
    for (uint64_t ri = 0; ri < target_set.n; ri++) {
      pair_result_t *r = &rows[ri];
      if (r->qvalue > args.qvalue_filter) continue;
      motif_t *t = target_set.motifs[r->target_i];
      /* Build the aligned-consensus substrings for the overlap region */
      const double (*qppm)[4] = (r->orientation == 0) ? q->ppm : q->ppm_rc;
      /* In oriented query coords, the overlap starts at:
           q_start_oriented = (orient == 0) ? r->q_start_in_orig
                                            : (q->size - 1 - (r->q_start_in_orig + r->L - 1)) */
      uint64_t q_start_oriented;
      if (r->orientation == 0) q_start_oriented = (uint64_t)r->q_start_in_orig;
      else q_start_oriented = q->size - r->q_start_in_orig - (uint64_t)r->L;
      build_consensus_range(qppm, q_start_oriented, (uint64_t)r->L, qcons);
      build_consensus_range(t->ppm, (uint64_t)r->t_start, (uint64_t)r->L, tcons);
      fprintf(o, "%s\t%s\t%d\t%.6g\t%.6g\t%.6g\t%d\t%s\t%s\t%c\n",
        q->name, t->name, r->offset, r->pvalue, r->qvalue, r->score, r->L,
        qcons, tcons, r->orientation == 0 ? '+' : '-');
    }
  }
}

/* ---- Usage ---- */

static void usage(void) {
  printf(
    "yamtk v%s  Copyright (C) %s  Benjamin Jean-Marie Tremblay\n"
    "Usage:  yamtk cmp [options] -m queries.meme -t targets.meme\n"
    "\n"
    " -m <str>   Query motif file (MEME/JASPAR/HOMER/HOCOMOCO).\n"
    " -t <str>   Target motif database (same formats).\n"
    " -o <str>   Output TSV file (default: stdout).\n"
    " -d <str>   Column-similarity metric: pcc (default: pcc).\n"
    " -n <int>   Minimum overlap columns (default: %d).\n"
    " -R         Disable reverse-strand scoring.\n"
    " -q <dbl>   Only report rows with q-value <= this (default: %g).\n"
    " -b A,C,G,T Background (default: from query MEME or uniform).\n"
    " -p <dbl>   Pseudocount added to PPMs before scoring (default: %g).\n"
    " -N <int>   Override every motif's nsites (default: parsed nsites,\n"
    "            falling back to %" PRIu64 ").\n"
    " -r         Do not trim motif names to the first word.\n"
    " -e         Use parametric enumeration null (Dirichlet-Multinomial over\n"
    "            a 56-column grid) instead of the empirical null. Use when\n"
    "            the target db is small (< ~500 columns).\n"
    " -j <int>   Threads (default: 1).\n"
    " -g         Show progress bar.\n"
    " -v / -w / -h   Verbose / very-verbose / help.\n"
    , YAMTK_VERSION, YAMTK_YEAR,
    DEFAULT_MIN_OVERLAP, DEFAULT_QVALUE_FILTER, DEFAULT_PSEUDOCOUNT, DEFAULT_NSITES
  );
}

/* ---- main_cmp ---- */

int main_cmp(int argc, char **argv) {
  int opt;
  int has_q = 0, has_t = 0;
  char *user_bkg = NULL;
  char *metric_name = NULL;
  const char *output_path = NULL;
  pthread_t *threads = NULL;

  setlocale(LC_NUMERIC, "");

  struct timespec ts_program;
  clock_gettime(CLOCK_MONOTONIC, &ts_program);

  while ((opt = getopt(argc, argv, "m:t:o:d:n:Rq:b:p:N:rej:gvwh")) != -1) {
    switch (opt) {
      case 'm':
        if (files.q_open) badexit("Error: -m specified more than once.");
        has_q = 1;
        files.q = fopen(optarg, "r");
        if (!files.q) {
          fprintf(stderr, "Error: Cannot open -m file: %s [%s]", optarg, strerror(errno));
          badexit("");
        }
        files.q_open = 1;
        break;
      case 't':
        if (files.t_open) badexit("Error: -t specified more than once.");
        has_t = 1;
        files.t = fopen(optarg, "r");
        if (!files.t) {
          fprintf(stderr, "Error: Cannot open -t file: %s [%s]", optarg, strerror(errno));
          badexit("");
        }
        files.t_open = 1;
        break;
      case 'o':
        if (files.o_open) badexit("Error: -o specified more than once.");
        output_path = optarg;
        files.o = fopen(optarg, "w");
        if (!files.o) {
          fprintf(stderr, "Error: Cannot create output file: %s [%s]", optarg, strerror(errno));
          badexit("");
        }
        files.o_open = 1;
        break;
      case 'd':
        metric_name = optarg;
        break;
      case 'n':
        if (str_to_int(optarg, &args.min_overlap)) badexit("Error: Failed to parse -n value.");
        if (args.min_overlap < 1) badexit("Error: -n must be >= 1.");
        break;
      case 'R':
        args.scan_rc = 0;
        break;
      case 'q':
        if (str_to_double(optarg, &args.qvalue_filter)) badexit("Error: Failed to parse -q value.");
        if (args.qvalue_filter < 0.0 || args.qvalue_filter > 1.0)
          badexit("Error: -q must be in [0,1].");
        break;
      case 'b':
        if (user_bkg) badexit("Error: -b specified more than once.");
        user_bkg = optarg;
        args.use_user_bkg = 1;
        break;
      case 'p':
        if (str_to_double(optarg, &args.pseudocount)) badexit("Error: Failed to parse -p value.");
        if (args.pseudocount < 0.0) badexit("Error: -p must be >= 0.");
        break;
      case 'N': {
        int n;
        if (str_to_int(optarg, &n)) badexit("Error: Failed to parse -N value.");
        if (n < 1) badexit("Error: -N must be >= 1.");
        args.nsites_override = (uint64_t)n;
        break;
      }
      case 'r':
        args.trim_names = 0;
        break;
      case 'e':
        args.enum_null = 1;
        break;
      case 'g':
        args.progress = 1;
        break;
      case 'j':
        if (str_to_int(optarg, &args.nthreads)) badexit("Error: Failed to parse -j value.");
        if (args.nthreads < 1) badexit("Error: -j must be >= 1.");
        break;
      case 'w':
        args.w = 1; /* fall through */
      case 'v':
        args.v = 1;
        break;
      case 'h':
        usage();
        close_files();
        return EXIT_SUCCESS;
      default:
        usage();
        close_files();
        return EXIT_FAILURE;
    }
  }

  if (!has_q) badexit("Error: Missing -m (query motifs).");
  if (!has_t) badexit("Error: Missing -t (target motifs).");

  if (metric_name) {
    args.metric_idx = find_metric(metric_name);
    if (args.metric_idx < 0) {
      fprintf(stderr, "Error: Unknown metric '%s'. Available: ", metric_name);
      for (int i = 0; i < N_METRICS; i++) fprintf(stderr, "%s%s", i ? ", " : "", METRICS[i].name);
      badexit("");
    }
  }

  if (user_bkg) parse_user_bkg(user_bkg);

  time_t t0 = time(NULL);

  /* Load queries (read MEME bkg only if user didn't override) */
  if (args.v) fprintf(stderr, "Loading query motifs ...\n");
  load_motifs_from(files.q, &query_set, "query", 1);
  if (!query_set.n) badexit("Error: No query motifs loaded.");

  /* Load targets (don't override args.bkg from target file) */
  if (args.v) fprintf(stderr, "Loading target motifs ...\n");
  load_motifs_from(files.t, &target_set, "target", 0);
  if (!target_set.n) badexit("Error: No target motifs loaded.");

  if (args.v) print_time((uint64_t)difftime(time(NULL), t0), "load motifs");

  /* args.bkg is now final (user override or query-MEME value). */
  if (args.enum_null) {
    build_enum_cols();
    if (args.w)
      fprintf(stderr, "Using enumeration null over %d grid columns.\n", ENUM_N_COLS);
  } else if (args.w) {
    uint64_t total_cols = 0;
    for (uint64_t ti = 0; ti < target_set.n; ti++) total_cols += target_set.motifs[ti]->size;
    fprintf(stderr, "Using empirical null over %" PRIu64 " target columns from %" PRIu64
      " target motifs.\n", total_cols, target_set.n);
  }

  /* Cap threads to #queries */
  if ((uint64_t)args.nthreads > query_set.n) {
    if (args.v && args.nthreads > 1)
      fprintf(stderr, "Reducing -j from %d to %" PRIu64 " (one thread per query max).\n",
        args.nthreads, query_set.n);
    args.nthreads = (int)query_set.n;
  }
  threads = malloc(sizeof(pthread_t) * args.nthreads);
  if (!threads) badexit("Error: Failed to alloc threads.");
  for (uint64_t i = 0; i < query_set.n; i++)
    query_set.motifs[i]->thread = (uint64_t)(((double)i / query_set.n) * args.nthreads);

  /* Allocate per-query result rows */
  query_results = malloc(sizeof(*query_results) * query_set.n);
  if (!query_results) { free(threads); badexit("Error: Failed to alloc results array."); }
  for (uint64_t qi = 0; qi < query_set.n; qi++) {
    query_results[qi] = malloc(sizeof(pair_result_t) * target_set.n);
    if (!query_results[qi]) {
      for (uint64_t k = 0; k < qi; k++) free(query_results[k]);
      free(query_results); query_results = NULL;
      free(threads);
      badexit("Error: Failed to alloc result rows.");
    }
  }

  if (args.v) fprintf(stderr, "Comparing %" PRIu64 " queries x %" PRIu64 " targets ...\n",
    query_set.n, target_set.n);
  t0 = time(NULL);

  if (args.progress) print_pb(0.0);
  for (uint64_t i = 0; i < (uint64_t)args.nthreads; i++) {
    uint64_t *ti = malloc(sizeof(*ti));
    if (!ti) badexit("Error: Failed to alloc thread index.");
    *ti = i;
    pthread_create(&threads[i], NULL, cmp_sub_process, ti);
  }
  for (uint64_t i = 0; i < (uint64_t)args.nthreads; i++) pthread_join(threads[i], NULL);
  if (args.progress) fprintf(stderr, "\n");

  if (args.v) print_time((uint64_t)difftime(time(NULL), t0), "compare");
  free(threads);

  /* Write output */
  FILE *out = files.o_open ? files.o : stdout;
  (void)output_path;
  write_results(out, argc, argv);

  /* Free per-query result rows */
  for (uint64_t qi = 0; qi < query_set.n; qi++) free(query_results[qi]);
  free(query_results); query_results = NULL;

  free_motif_set(&query_set);
  free_motif_set(&target_set);
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
