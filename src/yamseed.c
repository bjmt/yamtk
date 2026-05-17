/*
 *   yamseed: Seed input sequences with samples from input motifs
 *   Copyright (C) 2026  Benjamin Jean-Marie Tremblay
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
#include <getopt.h>
#include <math.h>
#include <limits.h>
#include <time.h>
#include <stdint.h>
#include <inttypes.h>
#include <zlib.h>
#include <ctype.h>
#include "kseq.h"
#include "khash.h"
#include "version.h"

KSEQ_INIT(gzFile, gzread)
KHASH_MAP_INIT_STR(seq_str_h, uint64_t)
KHASH_MAP_INIT_STR(motif_str_h, uint64_t)

/* ---- Constants ---- */

#define MAX_NAME_SIZE           ((uint64_t) 256)
#define MAX_MOTIF_WIDTH                       50
#define MOTIF_VALUE_MAX_CHAR    ((uint64_t) 256)
#define ALLOC_CHUNK_SIZE        ((uint64_t) 256)
#define FASTA_LINE_LEN                        60
#define DEFAULT_SEED                           4
#define DEFAULT_NSITES                      1000
#define DEFAULT_PSEUDOCOUNT                    1
#define MAX_RETRIES                          100
#define BED_FIELD_MAX_CHAR    ((uint64_t) 256)
#define BED_NAME_MAX_CHAR     ((uint64_t) 512)
#define BED_ALLOC_CHUNK_SIZE  ((uint64_t) 256)

#define ERASE_ARRAY(ARR, LEN) memset(ARR, 0, sizeof(ARR[0]) * (LEN))
#define VEC_DIV(VEC, X, VEC_LEN) \
  do { for (uint64_t Xi = 0; Xi < VEC_LEN; Xi++) VEC[Xi] /= X; } while (0)

#define ROTL(x, k) (((x) << (k)) | ((x) >> (64 - (k))))

static const char index2dna[6] = "ACGTN";

/* ---- PRNG (xoroshiro128++, seeded via splitmix64) ---- */

typedef struct { uint64_t s[2]; } xrng_t;

static inline uint64_t splitmix64(uint64_t x) {
  uint64_t z = (x += 0x9e3779b97f4a7c15);
  z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9;
  z = (z ^ (z >> 27)) * 0x94d049bb133111eb;
  return z ^ (z >> 31);
}

static inline uint64_t xrand_r(xrng_t *r) {
  const uint64_t s0 = r->s[0];
  uint64_t s1 = r->s[1];
  const uint64_t result = ROTL(s0 + s1, 17) + s0;
  s1 ^= s0;
  r->s[0] = ROTL(s0, 49) ^ s1 ^ (s1 << 21);
  r->s[1] = ROTL(s1, 28);
  return result;
}

static inline void sxrand_r(xrng_t *r, uint64_t seed) {
  r->s[0] = splitmix64(seed);
  r->s[1] = splitmix64(r->s[0]);
}

/* Uniform double in [0, 1). 53-bit mantissa precision. */
static inline double xrand_double(xrng_t *r) {
  return (xrand_r(r) >> 11) * (1.0 / (double)(1ULL << 53));
}

static xrng_t xrng;

/* ---- Insertion counters (for -v / -w summary) ---- */
static uint64_t total_insertions = 0;
static uint64_t *per_motif_count = NULL;

/* ---- Motif struct (slim — no PWM scoring; PPM-only) ---- */

typedef struct motif_t {
  char     name[MAX_NAME_SIZE];
  uint64_t size;
  uint64_t file_line_num;
  double   pwm_probs[MAX_MOTIF_WIDTH][4];  /* 0=A 1=C 2=G 3=T */
  int      dropped;
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

typedef struct args_t {
  double   lambda;        /* -f, per-bp Poisson rate (random mode) */
  uint64_t min_spacing;   /* -M */
  int      seed;          /* -s */
  int      trim_names;    /* 1 = trim (default); -r sets to 0 */
  int      use_random;    /* set by -f */
  int      use_bed;       /* set by -x */
  int      no_rc;         /* -R */
  int      progress;      /* -g */
  int      v;
  int      w;
} args_t;

static args_t args = {
  .lambda      = 0.0,
  .min_spacing = 0,
  .seed        = DEFAULT_SEED,
  .trim_names  = 1,
  .use_random  = 0,
  .use_bed     = 0,
  .no_rc       = 0,
  .progress    = 0,
  .v           = 0,
  .w           = 0
};

/* ---- Files ---- */

typedef struct files_t {
  int     m_open;
  int     i_open;
  int     o_open;
  int     x_open;
  int     O_open;
  FILE   *m;
  gzFile  i;
  FILE   *o;
  gzFile  x;
  FILE   *O;
} files_t;

static files_t files = {
  .m_open = 0,
  .i_open = 0,
  .o_open = 0,
  .x_open = 0,
  .O_open = 0
};

/* ---- Cleanup ---- */

static void close_files(void) {
  if (files.m_open) fclose(files.m);
  if (files.i_open) gzclose(files.i);
  if (files.o_open) fclose(files.o);
  if (files.x_open) gzclose(files.x);
  if (files.O_open) fclose(files.O);
}

static void free_motifs(void) {
  if (!motifs) return;
  for (uint64_t i = 0; i < motif_info.n; i++) free(motifs[i]);
  free(motifs); motifs = NULL;
  motif_info.n = 0; motif_info.n_alloc = 0;
}

static void free_bed(void);

static void badexit(const char *msg) {
  if (msg && msg[0]) fprintf(stderr, "%s\n", msg);
  fprintf(stderr, "Run yamtk seed -h to see usage.\n");
  free(per_motif_count); per_motif_count = NULL;
  free_motifs();
  free_bed();
  close_files();
  exit(EXIT_FAILURE);
}

/* ---- String helpers ---- */

static inline int str_to_int(char *str, int *res) {
  char *tmp; errno = 0;
  long int res_long = strtol(str, &tmp, 10);
  if (res_long > INT_MAX) return 1;
  *res = (int) res_long;
  if (str == tmp || errno != 0 || *tmp != '\0') return 1;
  return 0;
}

static inline int str_to_uint64_t(char *str, uint64_t *res) {
  char *tmp; errno = 0;
  *res = (uint64_t) strtoull(str, &tmp, 10);
  if (str == tmp || errno != 0 || *tmp != '\0') return 1;
  return 0;
}

static inline int str_to_double(char *str, double *res) {
  char *tmp; errno = 0;
  *res = strtod(str, &tmp);
  if (str == tmp || errno != 0 || *tmp != '\0') return 1;
  return 0;
}

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

static int is_ppm_data_line(const char *line) {
  uint64_t i = 0;
  while (line[i] == ' ' || line[i] == '\t') i++;
  return (line[i] >= '0' && line[i] <= '9') || line[i] == '.';
}

/* ---- Motif init/add ---- */

static void init_motif(motif_t *motif) {
  ERASE_ARRAY(motif->name, MAX_NAME_SIZE);
  motif->name[0]='m'; motif->name[1]='o'; motif->name[2]='t';
  motif->name[3]='i'; motif->name[4]='f'; motif->name[5]='\0';
  motif->size = 0;
  motif->file_line_num = 0;
  motif->dropped = 0;
  memset(motif->pwm_probs, 0, sizeof(motif->pwm_probs));
}

static int add_motif(void) {
  motif_info.n++;
  const uint64_t last_i = motif_info.n - 1;
  if (motif_info.n > motif_info.n_alloc) {
    motif_t **tmp_ptr = realloc(motifs,
      sizeof(*motifs) * motif_info.n_alloc + sizeof(*motifs) * ALLOC_CHUNK_SIZE);
    if (!tmp_ptr) {
      fprintf(stderr, "Error: Failed to allocate memory for motifs.\n");
      return 1;
    }
    motifs = tmp_ptr;
    motif_info.n_alloc += ALLOC_CHUNK_SIZE;
  }
  motifs[last_i] = malloc(sizeof(motif_t));
  if (!motifs[last_i]) {
    fprintf(stderr, "Error: Failed to allocate memory for motif.\n");
    return 1;
  }
  init_motif(motifs[last_i]);
  return 0;
}

static void trim_motif_name(motif_t *m) {
  for (uint64_t i = 0; i < MAX_NAME_SIZE; i++) {
    if (m->name[i] == ' ' || m->name[i] == '\t' || m->name[i] == '\0') {
      m->name[i] = '\0';
      break;
    }
  }
}

/* ---- PPM column ingest (no PWM scoring) ---- */

static int normalize_probs(double *probs, const char *name) {
  double sum = probs[0] + probs[1] + probs[2] + probs[3];
  if (fabs(sum - 1.0) > 0.1) {
    fprintf(stderr,
      "Error: Position for [%s] does not add up to 1 (sum=%.3g)\n", name, sum);
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
          fprintf(stderr, "Error: Motif [%s] has too many columns (need %" PRIu64 ").\n",
            motif->name, n);
          return 1;
        }
        if (str_to_double(pos_i, &probs[which_i])) {
          fprintf(stderr, "Error: Failed to parse value for motif: %s.\n", motif->name);
          fprintf(stderr, "  Line: %s  Bad value: %s\n", line, pos_i);
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
      fprintf(stderr, "Error: Motif [%s] has too many columns (need %" PRIu64 ").\n",
        motif->name, n);
      return 1;
    }
    if (str_to_double(pos_i, &probs[which_i])) {
      fprintf(stderr, "Error: Failed to parse value for motif: %s.\n", motif->name);
      fprintf(stderr, "  Line: %s  Bad value: %s\n", line, pos_i);
      return 1;
    }
  }
  if (which_i == (uint64_t)-1) {
    fprintf(stderr, "Error: Motif [%s] has an empty row.\n", motif->name);
    return 1;
  }
  if (which_i < n - 1) {
    fprintf(stderr, "Error: Motif [%s] has too few columns (need %" PRIu64 ").\n",
      motif->name, n);
    return 1;
  }
  return 0;
}

/* MEME/HOMER PPM column: probs are already proportions; store directly. */
static int add_motif_ppm_column(motif_t *motif, const char *line, const uint64_t pos) {
  double probs[4] = {-1.0, -1.0, -1.0, -1.0};
  if (get_line_probs(motif, line, probs, 4)) return 1;
  if (normalize_probs(probs, motif->name)) return 1;
  if (pos >= MAX_MOTIF_WIDTH) {
    fprintf(stderr, "Error: Motif [%s] exceeds max width (%d).\n", motif->name, MAX_MOTIF_WIDTH);
    return 1;
  }
  for (int i = 0; i < 4; i++) motif->pwm_probs[pos][i] = probs[i];
  return 0;
}

/* JASPAR/HOCOMOCO: counts → pseudocounted probs (so no zero-probability slots). */
static int add_motif_pcm_column(motif_t *motif, const double counts[4], const uint64_t pos) {
  double col = counts[0] + counts[1] + counts[2] + counts[3];
  if (col <= 0.0) col = 1.0;
  if (pos >= MAX_MOTIF_WIDTH) {
    fprintf(stderr, "Error: Motif [%s] exceeds max width (%d).\n", motif->name, MAX_MOTIF_WIDTH);
    return 1;
  }
  for (int i = 0; i < 4; i++) {
    double raw = counts[i] / col;
    double adj = (raw * DEFAULT_NSITES + (double)DEFAULT_PSEUDOCOUNT / 4.0) /
                 (DEFAULT_NSITES + DEFAULT_PSEUDOCOUNT);
    motif->pwm_probs[pos][i] = adj;
  }
  return 0;
}

/* ---- Format detection ---- */

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
            free(line); badexit("Error: yamseed cannot read HOCOMOCO PWMs.");
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

/* ---- MEME reader ---- */

static int check_meme_alph(const char *line, const uint64_t line_num) {
  if (check_line_contains(line, "ALPHABET= ACDEFGHIKLMNPQRSTVWY\0")) {
    fprintf(stderr, "Error: Detected protein alphabet (L%" PRIu64 ").\n", line_num);
    return 1;
  }
  return 0;
}

static void parse_meme_name(const char *line, const uint64_t motif_i) {
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
  uint64_t line_num = 0, l_p_m_L = 0, motif_i = -1, pos_i = 0;
  int alph_detected = 0, live_motif = 0, warned_alph = 0;
  while ((r = getline(&line, &len, files.m)) != -1) {
    line_num++;
    if (check_line_contains(line, "ALPHABET\0")) {
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
          fprintf(stderr, "Error: Motif [%s] is too large (max=%d).\n",
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
  if (args.v) fprintf(stderr, "Found %" PRIu64 " MEME motif(s).\n", motif_info.n);
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
        fprintf(stderr, "Error: Motif [%s] is too large (max=%d).\n",
          motifs[motif_i]->name, MAX_MOTIF_WIDTH);
        free(line); badexit("");
      }
      if (add_motif_ppm_column(motifs[motif_i], line, pos_i)) { free(line); badexit(""); }
      pos_i++;
      motifs[motif_i]->size = pos_i;
    }
  }
  free(line);
  if (args.v) fprintf(stderr, "Found %" PRIu64 " HOMER motif(s).\n", motif_info.n);
}

/* ---- JASPAR reader ---- */

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
    fprintf(stderr, "Error: Couldn't find ACGTU in motif [%s] row names.\n", motif->name);
    return 1;
  }
  if (left_bracket == (uint64_t)-1 || right_bracket == (uint64_t)-1) {
    fprintf(stderr, "Error: Couldn't find '[]' in motif [%s] row.\n", motif->name);
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
          fprintf(stderr, "Error: Motif [%s] has too many columns (max %d).\n",
            motif->name, MAX_MOTIF_WIDTH);
          return 1;
        }
        double v;
        if (str_to_double(prob_c, &v)) {
          fprintf(stderr, "Error: Failed to parse count for motif: %s. Bad value: %s\n",
            motif->name, prob_c);
          return 1;
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
      fprintf(stderr, "Error: Motif [%s] has too many columns (max %d).\n",
        motif->name, MAX_MOTIF_WIDTH);
      return 1;
    }
    double v;
    if (str_to_double(prob_c, &v)) {
      fprintf(stderr, "Error: Failed to parse count for motif: %s. Bad value: %s\n",
        motif->name, prob_c);
      return 1;
    }
    counts[pos_i][row_i] = v;
  }
  if (pos_i == (uint64_t)-1) {
    fprintf(stderr, "Error: Motif [%s] has an empty row.\n", motif->name);
    return 1;
  }
  pos_i++;
  if (*width_out && *width_out != pos_i) {
    fprintf(stderr, "Error: Motif [%s] has rows with differing widths.\n", motif->name);
    return 1;
  }
  *width_out = pos_i;
  return 0;
}

static int finalize_jaspar_motif(motif_t *motif, double counts[MAX_MOTIF_WIDTH][4]) {
  for (uint64_t j = 0; j < motif->size; j++) {
    if (add_motif_pcm_column(motif, counts[j], j)) return 1;
  }
  return 0;
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
      if (motif_i < (uint64_t)-1 && row_i != 4) {
        fprintf(stderr, "Error: Motif [%s] has %s rows\n", motifs[motif_i]->name,
          row_i < 4 ? "too few" : "too many");
        free(line); badexit("");
      }
      if (motif_i < (uint64_t)-1) {
        motifs[motif_i]->size = cur_width;
        if (finalize_jaspar_motif(motifs[motif_i], counts)) { free(line); badexit(""); }
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
  if (motif_i < (uint64_t)-1 && row_i != 4) {
    fprintf(stderr, "Error: Motif [%s] has %s rows\n", motifs[motif_i]->name,
      row_i < 4 ? "too few" : "too many");
    badexit("");
  }
  if (motif_info.n) {
    motifs[motif_info.n - 1]->size = cur_width;
    if (finalize_jaspar_motif(motifs[motif_info.n - 1], counts)) badexit("");
  }
  if (args.v) fprintf(stderr, "Found %" PRIu64 " JASPAR motif(s).\n", motif_info.n);
}

/* ---- HOCOMOCO reader (PCM rows: rows = positions, cols = ACGT counts) ---- */

static int add_motif_hocomoco_row(motif_t *motif, const char *line, const uint64_t pos) {
  double probs[4] = {-1.0, -1.0, -1.0, -1.0};
  if (get_line_probs(motif, line, probs, 4)) return 1;
  double pcm_sum = probs[0] + probs[1] + probs[2] + probs[3];
  if (pcm_sum < 0.99) {
    fprintf(stderr, "Error: Motif [%s] PCM row adds up to less than 1\n", motif->name);
    return 1;
  }
  return add_motif_pcm_column(motif, probs, pos);
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
        fprintf(stderr, "Error: Motif [%s] is too large (max=%d).\n",
          motifs[motif_i]->name, MAX_MOTIF_WIDTH);
        free(line); badexit("");
      }
      if (add_motif_hocomoco_row(motifs[motif_i], line, pos_i)) { free(line); badexit(""); }
      pos_i++;
      motifs[motif_i]->size = pos_i;
    }
  }
  free(line);
  if (args.v) fprintf(stderr, "Found %" PRIu64 " HOCOMOCO motif(s).\n", motif_info.n);
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
    fprintf(stderr, "Warning: %" PRIu64 " empty motif(s) ignored.\n", empty_motifs);
  if (args.trim_names) {
    for (uint64_t i = 0; i < motif_info.n; i++) trim_motif_name(motifs[i]);
  }
}

/* ---- FASTA output ---- */

static void write_seq(const unsigned char *seq, const uint64_t size, const char *name) {
  fprintf(files.o, ">%s\n", name);
  for (uint64_t i = 0; i < size; i += FASTA_LINE_LEN) {
    fprintf(files.o, "%.*s\n", FASTA_LINE_LEN, seq + i);
  }
}

/* ---- Sampling + random-mode seeding ---- */

static unsigned char sample_base(const double probs[4], xrng_t *r) {
  double u = xrand_double(r);
  double c = probs[0]; if (u < c) return 0;
  c += probs[1];       if (u < c) return 1;
  c += probs[2];       if (u < c) return 2;
  return 3;
}

static void overwrite_motif(unsigned char *seq, const uint64_t pos, const motif_t *m,
                            const int rc, xrng_t *r) {
  if (rc) {
    /* RC PWM: output position i samples from forward column (size-1-i) with
       indices complemented (A↔T, C↔G). */
    for (uint64_t i = 0; i < m->size; i++) {
      const double *fwd = m->pwm_probs[m->size - 1 - i];
      const double rc_col[4] = { fwd[3], fwd[2], fwd[1], fwd[0] };
      seq[pos + i] = (unsigned char) index2dna[sample_base(rc_col, r)];
    }
  } else {
    for (uint64_t i = 0; i < m->size; i++) {
      seq[pos + i] = (unsigned char) index2dna[sample_base(m->pwm_probs[i], r)];
    }
  }
}

/* Poisson sample. Knuth's method below ~30; normal approximation above. */
static uint64_t poisson_sample(const double lambda, xrng_t *r) {
  if (lambda <= 0.0) return 0;
  if (lambda < 30.0) {
    double L = exp(-lambda);
    double p = 1.0;
    uint64_t k = 0;
    do {
      k++;
      p *= xrand_double(r);
    } while (p > L);
    return k - 1;
  } else {
    double u1 = xrand_double(r);
    double u2 = xrand_double(r);
    if (u1 < 1e-300) u1 = 1e-300;
    double z = sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
    double x = lambda + sqrt(lambda) * z;
    if (x < 0.0) return 0;
    return (uint64_t) (x + 0.5);
  }
}

typedef struct { uint64_t s, e; } ivl_t;

/* True iff [s, e) lies within `gap` of any existing interval. */
static int overlaps_any(const ivl_t *ivls, const uint64_t n_ivls,
                        const uint64_t s, const uint64_t e, const uint64_t gap) {
  for (uint64_t k = 0; k < n_ivls; k++) {
    /* clear iff existing ends ≥ gap before s, OR existing starts ≥ gap after e. */
    const int left_clear  = (ivls[k].e + gap <= s);
    const int right_clear = (e + gap <= ivls[k].s);
    if (!(left_clear || right_clear)) return 1;
  }
  return 0;
}

/* Per-sequence random mode. Tracks placed intervals; retries up to MAX_RETRIES
   per insertion on collision (collision = overlap or within args.min_spacing). */
static void do_random_mode(unsigned char *seq, const uint64_t L, const char *name, xrng_t *r) {
  const uint64_t N = poisson_sample(args.lambda * (double)L, r);
  if (N == 0) return;
  ivl_t *ivls = malloc(N * sizeof(*ivls));
  if (!ivls) {
    fprintf(stderr, "Error: Failed to allocate interval tracker for [%s].\n", name);
    badexit("");
  }
  uint64_t n_ivls = 0;
  for (uint64_t k = 0; k < N; k++) {
    int placed = 0;
    for (int attempt = 0; attempt < MAX_RETRIES && !placed; attempt++) {
      const uint64_t mi = xrand_r(r) % motif_info.n;
      motif_t *m = motifs[mi];
      if (m->size == 0 || m->size > L) break;
      const uint64_t pos = xrand_r(r) % (L - m->size + 1);
      const uint64_t end = pos + m->size;
      if (overlaps_any(ivls, n_ivls, pos, end, args.min_spacing)) continue;
      const int rc = (!args.no_rc) && (xrand_r(r) & 1);
      overwrite_motif(seq, pos, m, rc, r);
      if (files.O_open) {
        fprintf(files.O, "%s\t%" PRIu64 "\t%" PRIu64 "\t%s\t.\t%c\n",
          name, pos, end, m->name, rc ? '-' : '+');
      }
      ivls[n_ivls++] = (ivl_t){pos, end};
      total_insertions++;
      per_motif_count[mi]++;
      placed = 1;
    }
  }
  free(ivls);
}

/* ---- BED file (-x mode) ---- */

typedef struct bed_t {
  uint64_t  *starts;
  uint64_t  *ends;
  char      *strands;
  char     **seq_names;
  char     **range_names;   /* col 4 = motif name */
  uint64_t   n_regions;
  uint64_t   n_alloc;
} bed_t;

static bed_t bed = { .n_regions = 0, .n_alloc = 0 };

static khash_t(seq_str_h)   *bed_seq_hash   = NULL;
static khash_t(motif_str_h) *motif_name_hash = NULL;
static uint64_t             *bed_next       = NULL;  /* linked list of regions per seq */

static void free_bed(void) {
  if (bed.n_alloc) {
    if (bed.n_regions) {
      for (uint64_t i = 0; i < bed.n_regions; i++) {
        free(bed.seq_names[i]);
        free(bed.range_names[i]);
      }
    }
    free(bed.seq_names);
    free(bed.range_names);
    free(bed.starts);
    free(bed.ends);
    free(bed.strands);
    bed.n_alloc = 0;
    bed.n_regions = 0;
  }
  free(bed_next); bed_next = NULL;
  if (bed_seq_hash) { kh_destroy(seq_str_h, bed_seq_hash); bed_seq_hash = NULL; }
  if (motif_name_hash) { kh_destroy(motif_str_h, motif_name_hash); motif_name_hash = NULL; }
}

static uint64_t count_fields(const char *line) {
  uint64_t res = 1, i = 0;
  while (line[i] != '\0') { if (line[i] == '\t') res++; i++; }
  return res;
}

static uint64_t count_field_size(const char *line, const uint64_t k) {
  uint64_t res = 0, i = 0, n = 0;
  while (line[i] != '\0') {
    if (line[i] == '\t') { n++; }
    else if (n + 1 == k) res++;
    else if (n + 1 > k) break;
    i++;
  }
  return res;
}

static uint64_t field_start(const char *line, const uint64_t k) {
  uint64_t i = 0, n = 0;
  while (line[i] != '\0') {
    if (line[i] == '\t') n++;
    else if (n + 1 == k) break;
    i++;
  }
  return i;
}

static uint64_t field_end(const char *line, const uint64_t k) {
  uint64_t i = 0, n = 0;
  while (line[i] != '\0') {
    if (line[i] == '\t') n++;
    if (n == k) { i--; break; }
    i++;
  }
  return i;
}

static inline uint64_t parse_bed_field(const char *line, const uint64_t k, char *field, const int no_spaces) {
  uint64_t start_i = field_start(line, k);
  uint64_t end_i = field_end(line, k);
  uint64_t size_i = count_field_size(line, k);
  uint64_t field_it = 0;
  ERASE_ARRAY(field, BED_FIELD_MAX_CHAR);
  if (size_i > 0 && size_i < BED_FIELD_MAX_CHAR) {
    for (uint64_t i = start_i; i <= end_i; i++) {
      if (no_spaces && isspace((unsigned char) line[i])) size_i--;
      else { field[field_it++] = line[i]; }
    }
  }
  return size_i;
}

static void bed_grow_if_needed(void) {
  if (bed.n_regions + 1 <= bed.n_alloc) return;
  const uint64_t new_alloc = bed.n_alloc + BED_ALLOC_CHUNK_SIZE;
  char     **p1 = realloc(bed.seq_names,   sizeof(*bed.seq_names)   * new_alloc);
  uint64_t  *p2 = realloc(bed.starts,      sizeof(*bed.starts)      * new_alloc);
  uint64_t  *p3 = realloc(bed.ends,        sizeof(*bed.ends)        * new_alloc);
  char     **p4 = realloc(bed.range_names, sizeof(*bed.range_names) * new_alloc);
  char      *p5 = realloc(bed.strands,     sizeof(*bed.strands)     * new_alloc);
  if (!p1 || !p2 || !p3 || !p4 || !p5) {
    badexit("Error: Failed to grow BED storage.");
  }
  bed.seq_names = p1; bed.starts = p2; bed.ends = p3;
  bed.range_names = p4; bed.strands = p5;
  bed.n_alloc = new_alloc;
}

static void read_bed(void) {
  bed.seq_names   = malloc(sizeof(*bed.seq_names)   * BED_ALLOC_CHUNK_SIZE);
  bed.range_names = malloc(sizeof(*bed.range_names) * BED_ALLOC_CHUNK_SIZE);
  bed.starts      = malloc(sizeof(*bed.starts)      * BED_ALLOC_CHUNK_SIZE);
  bed.ends        = malloc(sizeof(*bed.ends)        * BED_ALLOC_CHUNK_SIZE);
  bed.strands     = malloc(sizeof(*bed.strands)     * BED_ALLOC_CHUNK_SIZE);
  if (!bed.seq_names || !bed.range_names || !bed.starts || !bed.ends || !bed.strands) {
    badexit("Error: Failed to allocate BED storage.");
  }
  bed.n_alloc = BED_ALLOC_CHUNK_SIZE;

  kstream_t *kbed = ks_init(files.x);
  kstring_t  line = { 0, 0, 0 };
  int        ret_val;
  uint64_t   line_num = 0;
  char       tmp_field[BED_FIELD_MAX_CHAR];

  while ((ret_val = ks_getuntil(kbed, '\n', &line, 0)) >= 0) {
    line_num++;
    if (count_nonempty_chars(line.s) == 0) continue;
    if (line.s[0] == '#') continue;
    if (line.l >= 7 && memcmp(line.s, "browser", 7) == 0) continue;
    if (line.l >= 5 && memcmp(line.s, "track", 5) == 0) continue;
    const uint64_t n_fields = count_fields(line.s);
    if (n_fields < 4) {
      ks_destroy(kbed);
      fprintf(stderr, "Error: BED line %" PRIu64 " has %" PRIu64 " fields; need ≥4 (col 4 = motif name).\n",
        line_num, n_fields);
      badexit("");
    }
    bed_grow_if_needed();

    /* Strand (col 6). Default to '.' if missing. */
    if (n_fields >= 6) {
      const uint64_t sz = parse_bed_field(line.s, 6, tmp_field, 1);
      if (sz != 1 || (tmp_field[0] != '+' && tmp_field[0] != '-' && tmp_field[0] != '.')) {
        ks_destroy(kbed);
        fprintf(stderr, "Error: BED line %" PRIu64 " strand field must be one of +/-/. (got '%s').\n",
          line_num, tmp_field);
        badexit("");
      }
      bed.strands[bed.n_regions] = tmp_field[0];
    } else {
      bed.strands[bed.n_regions] = '.';
    }

    /* Start (col 2) */
    if (parse_bed_field(line.s, 2, tmp_field, 1) == 0) {
      ks_destroy(kbed);
      fprintf(stderr, "Error: BED line %" PRIu64 " has empty start field.\n", line_num);
      badexit("");
    }
    uint64_t tmp_value;
    if (str_to_uint64_t(tmp_field, &tmp_value)) {
      ks_destroy(kbed);
      fprintf(stderr, "Error: Failed to parse BED start on line %" PRIu64 " ('%s').\n", line_num, tmp_field);
      badexit("");
    }
    bed.starts[bed.n_regions] = tmp_value;

    /* End (col 3) */
    if (parse_bed_field(line.s, 3, tmp_field, 1) == 0) {
      ks_destroy(kbed);
      fprintf(stderr, "Error: BED line %" PRIu64 " has empty end field.\n", line_num);
      badexit("");
    }
    if (str_to_uint64_t(tmp_field, &tmp_value)) {
      ks_destroy(kbed);
      fprintf(stderr, "Error: Failed to parse BED end on line %" PRIu64 " ('%s').\n", line_num, tmp_field);
      badexit("");
    }
    bed.ends[bed.n_regions] = tmp_value;
    if (bed.starts[bed.n_regions] >= bed.ends[bed.n_regions]) {
      ks_destroy(kbed);
      fprintf(stderr, "Error: BED line %" PRIu64 " has start >= end.\n", line_num);
      badexit("");
    }

    /* Motif name (col 4) */
    const uint64_t sz4 = parse_bed_field(line.s, 4, tmp_field, 0);
    if (sz4 == 0) {
      ks_destroy(kbed);
      fprintf(stderr, "Error: BED line %" PRIu64 " has empty motif-name field (col 4).\n", line_num);
      badexit("");
    }
    if (sz4 >= BED_FIELD_MAX_CHAR) {
      ks_destroy(kbed);
      fprintf(stderr, "Error: BED line %" PRIu64 " motif name too long (%" PRIu64 ">%" PRIu64 ").\n",
        line_num, sz4, BED_FIELD_MAX_CHAR - 1);
      badexit("");
    }
    bed.range_names[bed.n_regions] = malloc(sizeof(char) * (sz4 + 1));
    if (!bed.range_names[bed.n_regions]) { ks_destroy(kbed); badexit("Error: BED motif name alloc."); }
    memcpy(bed.range_names[bed.n_regions], tmp_field, sz4);
    bed.range_names[bed.n_regions][sz4] = '\0';
    /* Trim at first space (consistent with motif name trimming) */
    for (uint64_t i = 0; i < sz4; i++) {
      if (bed.range_names[bed.n_regions][i] == ' ') {
        bed.range_names[bed.n_regions][i] = '\0';
        break;
      }
    }

    /* Sequence name (col 1) */
    const uint64_t sz1 = parse_bed_field(line.s, 1, tmp_field, 0);
    if (sz1 == 0) {
      ks_destroy(kbed);
      free(bed.range_names[bed.n_regions]);
      fprintf(stderr, "Error: BED line %" PRIu64 " has empty sequence name.\n", line_num);
      badexit("");
    }
    if (sz1 >= BED_FIELD_MAX_CHAR) {
      ks_destroy(kbed);
      free(bed.range_names[bed.n_regions]);
      fprintf(stderr, "Error: BED line %" PRIu64 " sequence name too long.\n", line_num);
      badexit("");
    }
    bed.seq_names[bed.n_regions] = malloc(sizeof(char) * (sz1 + 1));
    if (!bed.seq_names[bed.n_regions]) {
      ks_destroy(kbed); free(bed.range_names[bed.n_regions]);
      badexit("Error: BED seq name alloc.");
    }
    memcpy(bed.seq_names[bed.n_regions], tmp_field, sz1);
    bed.seq_names[bed.n_regions][sz1] = '\0';
    for (uint64_t i = 0; i < sz1; i++) {
      if (bed.seq_names[bed.n_regions][i] == ' ') {
        bed.seq_names[bed.n_regions][i] = '\0';
        break;
      }
    }

    bed.n_regions++;
  }
  if (ret_val == -3) { ks_destroy(kbed); badexit("Error: Failed to read BED stream."); }
  if (!bed.n_regions) { ks_destroy(kbed); badexit("Error: BED has no usable records."); }
  ks_destroy(kbed);
  free(line.s);
  if (args.v) fprintf(stderr, "Read %" PRIu64 " BED region(s).\n", bed.n_regions);
}

static void build_bed_seq_hash(void) {
  bed_seq_hash = kh_init(seq_str_h);
  if (!bed_seq_hash) badexit("Error: Failed to init BED seq-name hash.");
  bed_next = malloc(sizeof(*bed_next) * bed.n_regions);
  if (!bed_next) badexit("Error: Failed to alloc BED next-region list.");
  for (uint64_t i = 0; i < bed.n_regions; i++) bed_next[i] = UINT64_MAX;
  int absent;
  khint_t k;
  for (uint64_t i = 0; i < bed.n_regions; i++) {
    k = kh_put(seq_str_h, bed_seq_hash, bed.seq_names[i], &absent);
    if (absent == -1) badexit("Error: Failed to hash BED sequence names.");
    if (absent == 0) {
      bed_next[i] = kh_val(bed_seq_hash, k);
      kh_val(bed_seq_hash, k) = i;
    } else {
      kh_val(bed_seq_hash, k) = i;
    }
  }
}

static void build_motif_name_hash(void) {
  motif_name_hash = kh_init(motif_str_h);
  if (!motif_name_hash) badexit("Error: Failed to init motif-name hash.");
  int absent;
  khint_t k;
  for (uint64_t i = 0; i < motif_info.n; i++) {
    k = kh_put(motif_str_h, motif_name_hash, motifs[i]->name, &absent);
    if (absent == -1) badexit("Error: Failed to hash motif names.");
    if (absent == 0 && args.v) {
      fprintf(stderr, "Warning: duplicate motif name '%s'; using last occurrence for BED lookup.\n",
        motifs[i]->name);
    }
    kh_val(motif_name_hash, k) = i;
  }
}

/* Collect all BED region indices belonging to a sequence into `out`; return count.
   `out` must be at least bed.n_regions long. */
static uint64_t collect_seq_regions(const char *seq_name, uint64_t *out) {
  if (!bed_seq_hash) return 0;
  khint_t k = kh_get(seq_str_h, bed_seq_hash, seq_name);
  if (k == kh_end(bed_seq_hash)) return 0;
  uint64_t n = 0;
  uint64_t idx = kh_val(bed_seq_hash, k);
  while (idx != UINT64_MAX) {
    out[n++] = idx;
    idx = bed_next[idx];
  }
  return n;
}

static int cmp_idx_by_start(const void *a, const void *b) {
  const uint64_t ia = *(const uint64_t *)a;
  const uint64_t ib = *(const uint64_t *)b;
  if (bed.starts[ia] < bed.starts[ib]) return -1;
  if (bed.starts[ia] > bed.starts[ib]) return  1;
  return 0;
}

static void do_bed_insertion(unsigned char *seq, const uint64_t L, const char *name,
                             const uint64_t region_i, xrng_t *r) {
  /* Look up motif by col-4 name */
  khint_t k = kh_get(motif_str_h, motif_name_hash, bed.range_names[region_i]);
  if (k == kh_end(motif_name_hash)) {
    fprintf(stderr, "Error: BED region '%s' (line refers to motif '%s') has no matching motif in -m.\n",
      name, bed.range_names[region_i]);
    badexit("");
  }
  const uint64_t mi = kh_val(motif_name_hash, k);
  motif_t *m = motifs[mi];
  const uint64_t w = m->size;
  if (w == 0) {
    fprintf(stderr, "Warning: motif '%s' is empty; skipping BED region on '%s'.\n",
      m->name, name);
    return;
  }
  if (L < w) {
    fprintf(stderr, "Warning: seq '%s' (len %" PRIu64 ") shorter than motif '%s' (width %" PRIu64 "); skipping.\n",
      name, L, m->name, w);
    return;
  }

  uint64_t bs = bed.starts[region_i];
  uint64_t be = bed.ends[region_i];
  uint64_t start;
  if ((be - bs) == w) {
    start = bs;
  } else {
    /* Center motif at BED midpoint. */
    fprintf(stderr,
      "Warning: BED region %s:%" PRIu64 "-%" PRIu64 " width %" PRIu64
      " != motif '%s' width %" PRIu64 "; centering at midpoint.\n",
      name, bs, be, be - bs, m->name, w);
    const uint64_t mid = (bs + be) / 2;
    const uint64_t half = w / 2;
    if (mid < half) {
      fprintf(stderr,
        "Warning: BED region too close to sequence start to center motif '%s'; skipping.\n", m->name);
      return;
    }
    start = mid - half;
    if (start + w > L) {
      fprintf(stderr,
        "Warning: BED region centered position exceeds seq '%s' length; skipping.\n", name);
      return;
    }
  }
  if (start + w > L) {
    fprintf(stderr,
      "Warning: BED start %" PRIu64 " + motif width %" PRIu64 " > seq '%s' length %" PRIu64 "; skipping.\n",
      start, w, name, L);
    return;
  }

  char strand = bed.strands[region_i];
  if (strand == '.') strand = '+';
  if (args.no_rc) strand = '+';
  const int rc = (strand == '-') ? 1 : 0;

  overwrite_motif(seq, start, m, rc, r);
  if (files.O_open) {
    fprintf(files.O, "%s\t%" PRIu64 "\t%" PRIu64 "\t%s\t.\t%c\n",
      name, start, start + w, m->name, strand);
  }
  total_insertions++;
  per_motif_count[mi]++;
}

/* ---- Usage ---- */

static void usage(void) {
  printf(
    "yamtk v%s  Copyright (C) %s  Benjamin Jean-Marie Tremblay\n"
    "Usage:  yamtk seed [options] -m motifs.txt -i seqs.fa[.gz]\n"
    "\n"
    " -m <str>   Motif file (MEME/JASPAR/HOMER/HOCOMOCO).\n"
    " -i <str>   Input FASTA/FASTQ ('-' = stdin). Sequence bases in seeded\n"
    "            regions are overwritten with samples from the motif PPM.\n"
    " -o <str>   Output FASTA (default: stdout).\n"
    " -O <str>   Write ground-truth BED of insertions (seq, start, end,\n"
    "            motif, '.', strand).\n"
    " -f <dbl>   Random mode: per-bp Poisson insertion rate. Excludes -x.\n"
    " -x <str>   BED mode: col-4 = motif name (must match a loaded motif),\n"
    "            col-6 = strand. If end-start != motif width, motif is\n"
    "            centered at the BED-range midpoint. Excludes -f.\n"
    " -M <int>   Minimum spacing (bp) between -f insertions (default: 0).\n"
    " -R         Disable reverse-strand sampling (always insert '+').\n"
    " -s <int>   RNG seed (default: %d).\n"
    " -r         Do not trim motif/sequence names to the first word.\n"
    " -g         Show progress bar.\n"
    " -v / -w / -h   Verbose / very-verbose / help.\n"
    , YAMTK_VERSION, YAMTK_YEAR, DEFAULT_SEED
  );
}

/* ---- Main ---- */

int main_seed(int argc, char **argv) {

  kseq_t *kseq;
  int opt;
  int use_stdout = 1;

  while ((opt = getopt(argc, argv, "m:i:o:O:f:x:M:Rrs:gvwh")) != -1) {
    switch (opt) {
      case 'm':
        if (files.m_open) badexit("Error: -m specified more than once.");
        files.m = fopen(optarg, "r");
        if (files.m == NULL) {
          fprintf(stderr, "Error: Failed to open motif file \"%s\" [%s]\n",
            optarg, strerror(errno));
          badexit("");
        }
        files.m_open = 1;
        break;
      case 'i':
        if (files.i_open) badexit("Error: -i specified more than once.");
        if (optarg[0] == '-' && optarg[1] == '\0') {
          files.i = gzdopen(fileno(stdin), "r");
        } else {
          files.i = gzopen(optarg, "r");
          if (files.i == NULL) {
            fprintf(stderr, "Error: Failed to open sequence file \"%s\" [%s]\n",
              optarg, strerror(errno));
            badexit("");
          }
        }
        files.i_open = 1;
        break;
      case 'o':
        if (files.o_open) badexit("Error: -o specified more than once.");
        use_stdout = 0;
        files.o = fopen(optarg, "w");
        if (files.o == NULL) {
          fprintf(stderr, "Error: Failed to create output file \"%s\" [%s]\n",
            optarg, strerror(errno));
          badexit("");
        }
        files.o_open = 1;
        break;
      case 'O':
        if (files.O_open) badexit("Error: -O specified more than once.");
        files.O = fopen(optarg, "w");
        if (files.O == NULL) {
          fprintf(stderr, "Error: Failed to create truth BED file \"%s\" [%s]\n",
            optarg, strerror(errno));
          badexit("");
        }
        files.O_open = 1;
        break;
      case 'f':
        if (str_to_double(optarg, &args.lambda)) {
          badexit("Error: Failed to parse -f value.");
        }
        if (args.lambda <= 0.0) {
          badexit("Error: -f must be a positive double.");
        }
        args.use_random = 1;
        break;
      case 'x':
        if (args.use_bed) badexit("Error: -x specified more than once.");
        files.x = gzopen(optarg, "r");
        if (files.x == NULL) {
          fprintf(stderr, "Error: Failed to open BED file \"%s\" [%s]\n",
            optarg, strerror(errno));
          badexit("");
        }
        files.x_open = 1;
        args.use_bed = 1;
        break;
      case 'M':
        if (str_to_uint64_t(optarg, &args.min_spacing)) {
          badexit("Error: Failed to parse -M value.");
        }
        break;
      case 'R':
        args.no_rc = 1;
        break;
      case 'r':
        args.trim_names = 0;
        break;
      case 's':
        if (str_to_int(optarg, &args.seed)) {
          badexit("Error: Failed to parse -s value.");
        }
        if (args.seed <= 0) {
          badexit("Error: -s must be a positive integer.");
        }
        break;
      case 'g':
        args.progress = 1;
        break;
      case 'w':
        args.w = 1;
        /* fallthrough */
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

  if (!files.m_open) badexit("Error: -m must be specified.");
  if (!files.i_open) badexit("Error: -i must be specified.");
  if (args.use_random && args.use_bed) {
    badexit("Error: -f and -x are mutually exclusive.");
  }
  if (!args.use_random && !args.use_bed) {
    badexit("Error: one of -f or -x must be specified.");
  }
  if (use_stdout) files.o = stdout;

  /* Parse motifs */
  load_motifs();

  /* Per-motif insertion counters (for -v / -w summary) */
  per_motif_count = calloc(motif_info.n, sizeof(*per_motif_count));
  if (!per_motif_count) badexit("Error: Failed to allocate per-motif counters.");

  /* BED setup */
  uint64_t *seq_regions_buf = NULL;
  if (args.use_bed) {
    read_bed();
    build_bed_seq_hash();
    build_motif_name_hash();
    seq_regions_buf = malloc(sizeof(*seq_regions_buf) * bed.n_regions);
    if (!seq_regions_buf) badexit("Error: Failed to alloc per-seq BED scratch.");
  }

  /* Seed the PRNG */
  sxrand_r(&xrng, (uint64_t) args.seed);

  /* Stream sequences */
  kseq = kseq_init(files.i);
  int ret_val;
  uint64_t n_seqs = 0;
  while ((ret_val = kseq_read(kseq)) >= 0) {
    n_seqs++;
    unsigned char *seq = (unsigned char *) kseq->seq.s;
    uint64_t L = kseq->seq.l;
    const char *name = kseq->name.s;
    const uint64_t before = total_insertions;
    if (args.use_random) {
      do_random_mode(seq, L, name, &xrng);
    } else if (args.use_bed) {
      const uint64_t nr = collect_seq_regions(name, seq_regions_buf);
      if (nr > 0) {
        qsort(seq_regions_buf, nr, sizeof(*seq_regions_buf), cmp_idx_by_start);
        for (uint64_t r_i = 0; r_i < nr; r_i++) {
          do_bed_insertion(seq, L, name, seq_regions_buf[r_i], &xrng);
        }
      }
    }
    if (args.w) {
      fprintf(stderr, "  %s (len %" PRIu64 "): %" PRIu64 " insertion(s)\n",
        name, L, total_insertions - before);
    }
    write_seq(seq, L, name);
    if (args.progress && (n_seqs % 100) == 0) {
      fprintf(stderr, "\rSequences processed: %" PRIu64, n_seqs);
      fflush(stderr);
    }
  }
  if (args.progress && n_seqs > 0) {
    fprintf(stderr, "\rSequences processed: %" PRIu64 "\n", n_seqs);
  }
  if (ret_val < -1) {
    fprintf(stderr, "Error: Failed to read input FASTA (kseq_read returned %d).\n", ret_val);
    kseq_destroy(kseq);
    free(seq_regions_buf);
    badexit("");
  }
  kseq_destroy(kseq);
  free(seq_regions_buf);

  if (args.v) {
    fprintf(stderr, "Processed %" PRIu64 " sequence(s); %" PRIu64 " insertion(s).\n",
      n_seqs, total_insertions);
  }
  if (args.w) {
    for (uint64_t i = 0; i < motif_info.n; i++) {
      fprintf(stderr, "  motif %s: %" PRIu64 " insertion(s)\n",
        motifs[i]->name, per_motif_count[i]);
    }
  }

  free(per_motif_count); per_motif_count = NULL;
  free_motifs();
  free_bed();
  close_files();
  return EXIT_SUCCESS;
}
