/*
 *   yamconv: Convert motifs between MEME/JASPAR/HOMER/HOCOMOCO formats
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
#include <stdint.h>
#include <inttypes.h>
#include <ctype.h>
#include "version.h"

/* ---- Constants ---- */

#define MAX_NAME_SIZE           ((uint64_t) 256)
#define MAX_MOTIF_WIDTH                       50
#define MOTIF_VALUE_MAX_CHAR    ((uint64_t) 256)
#define ALLOC_CHUNK_SIZE        ((uint64_t) 256)
#define DEFAULT_NSITES                      1000
#define DEFAULT_PSEUDOCOUNT                    0

#define ERASE_ARRAY(ARR, LEN) memset(ARR, 0, sizeof(ARR[0]) * (LEN))
#define VEC_DIV(VEC, X, VEC_LEN) \
  do { for (uint64_t Xi = 0; Xi < VEC_LEN; Xi++) VEC[Xi] /= X; } while (0)

/* ---- Motif struct (PPM-only; nsites preserved when known) ---- */

typedef struct motif_t {
  char     name[MAX_NAME_SIZE];
  uint64_t size;
  uint64_t file_line_num;
  uint64_t nsites;                          /* 0 = unknown */
  double   pwm_probs[MAX_MOTIF_WIDTH][4];   /* 0=A 1=C 2=G 3=T */
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
  int         to_fmt;        /* -t target format (FMT_*) */
  int         do_ic_trim;    /* -T sets to 1 */
  double      ic_min;        /* -T value (bits) */
  uint64_t    nsites_fb;     /* -N fallback nsites */
  int         pseudocount;   /* -p */
  double      bkg[4];        /* -b A,C,G,T */
  int         bkg_uniform;   /* 1 if bkg is the default uniform */
  int         trim_names;    /* 1 = trim (default); -r sets to 0 */
  int         v;
  int         w;
} args_t;

static args_t args = {
  .to_fmt      = 0,
  .do_ic_trim  = 0,
  .ic_min      = 0.5,
  .nsites_fb   = DEFAULT_NSITES,
  .pseudocount = DEFAULT_PSEUDOCOUNT,
  .bkg         = { 0.25, 0.25, 0.25, 0.25 },
  .bkg_uniform = 1,
  .trim_names  = 1,
  .v           = 0,
  .w           = 0
};

/* ---- Files ---- */

typedef struct files_t {
  int     m_open;
  int     o_open;
  FILE   *m;
  FILE   *o;
} files_t;

static files_t files = { .m_open = 0, .o_open = 0 };

/* ---- Cleanup ---- */

static void close_files(void) {
  if (files.m_open) fclose(files.m);
  if (files.o_open) fclose(files.o);
}

static void free_motifs(void) {
  if (!motifs) return;
  for (uint64_t i = 0; i < motif_info.n; i++) free(motifs[i]);
  free(motifs); motifs = NULL;
  motif_info.n = 0; motif_info.n_alloc = 0;
}

static void badexit(const char *msg) {
  if (msg && msg[0]) fprintf(stderr, "%s\n", msg);
  fprintf(stderr, "Run yamtk conv -h to see usage.\n");
  free_motifs();
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

/* Parse "nsites= 175" or "nsites=175" out of a MEME letter-probability matrix
   header line. Returns 0 (and sets *out) on success, 1 if no nsites token
   was found or it failed to parse. */
static int extract_meme_nsites(const char *line, uint64_t *out) {
  const char *p = strstr(line, "nsites");
  if (!p) return 1;
  p += 6;                                  /* past "nsites" */
  while (*p == ' ' || *p == '\t' || *p == '=') p++;
  if (!isdigit((unsigned char) *p)) return 1;
  char buf[32]; uint64_t k = 0;
  while (k < sizeof(buf) - 1 && isdigit((unsigned char) *p)) buf[k++] = *p++;
  buf[k] = '\0';
  return str_to_uint64_t(buf, out);
}

/* ---- Motif init/add ---- */

static void init_motif(motif_t *motif) {
  ERASE_ARRAY(motif->name, MAX_NAME_SIZE);
  motif->name[0]='m'; motif->name[1]='o'; motif->name[2]='t';
  motif->name[3]='i'; motif->name[4]='f'; motif->name[5]='\0';
  motif->size = 0;
  motif->file_line_num = 0;
  motif->nsites = 0;
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

/* ---- PPM/PCM column ingest ---- */

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

static int add_motif_pcm_column(motif_t *motif, const double counts[4], const uint64_t pos) {
  double col = counts[0] + counts[1] + counts[2] + counts[3];
  if (col <= 0.0) col = 1.0;
  if (pos >= MAX_MOTIF_WIDTH) {
    fprintf(stderr, "Error: Motif [%s] exceeds max width (%d).\n", motif->name, MAX_MOTIF_WIDTH);
    return 1;
  }
  /* Pseudocounted PPM (matches yamseed/yamscan); preserves raw nsites
     separately so PCM round-trips emit the original counts. */
  for (int i = 0; i < 4; i++) {
    double raw = counts[i] / col;
    double adj = (raw * DEFAULT_NSITES + (double) args.pseudocount / 4.0) /
                 (DEFAULT_NSITES + args.pseudocount);
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
            free(line); badexit("Error: yamconv cannot read HOCOMOCO PWMs.");
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
      uint64_t ns;
      if (extract_meme_nsites(line, &ns) == 0) motifs[motif_i]->nsites = ns;
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
  /* Capture nsites from first column's count sum (rounded). */
  if (motif->size > 0) {
    double s = counts[0][0] + counts[0][1] + counts[0][2] + counts[0][3];
    motif->nsites = (uint64_t) (s + 0.5);
  }
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
  if (pos == 0) motif->nsites = (uint64_t) (pcm_sum + 0.5);
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

/* ---- IC-based flank trimming (from yamref.c) ---- */

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

/* ---- HOMER consensus + threshold ---- */

/* Mask -> IUPAC: 4-bit (A=1,C=2,G=4,T=8). Single-base codes get the natural
   letter; ambiguous codes follow IUPAC convention. */
static const char mask2iupac[16] = {
  'N',   /* 0000 = empty (should never happen; caller falls back to argmax) */
  'A',   /* 0001 */
  'C',   /* 0010 */
  'M',   /* 0011 = A|C */
  'G',   /* 0100 */
  'R',   /* 0101 = A|G */
  'S',   /* 0110 = C|G */
  'V',   /* 0111 = A|C|G */
  'T',   /* 1000 */
  'W',   /* 1001 = A|T */
  'Y',   /* 1010 = C|T */
  'H',   /* 1011 = A|C|T */
  'K',   /* 1100 = G|T */
  'D',   /* 1101 = A|G|T */
  'B',   /* 1110 = C|G|T */
  'N'    /* 1111 = A|C|G|T */
};

static char column_consensus(const double probs[4]) {
  int mask = 0;
  for (int i = 0; i < 4; i++) {
    if (probs[i] >= 0.25) mask |= (1 << i);
  }
  if (mask == 0) {
    /* No base clears the 0.25 bar: fall back to argmax. */
    int amax = 0;
    for (int i = 1; i < 4; i++) if (probs[i] > probs[amax]) amax = i;
    mask = 1 << amax;
  }
  return mask2iupac[mask];
}

static void motif_consensus(const motif_t *m, char *out) {
  for (uint64_t i = 0; i < m->size; i++) out[i] = column_consensus(m->pwm_probs[i]);
  out[m->size] = '\0';
}

/* Max log-odds score = sum over columns of log2(max(p) / bkg). Pseudocount
   protection: the PWM has already been pseudocounted during PCM ingestion,
   but MEME/HOMER PPMs may legitimately contain zeros, so guard with a tiny
   epsilon to avoid -inf. */
static double motif_max_logodds(const motif_t *m) {
  double sum = 0.0;
  for (uint64_t i = 0; i < m->size; i++) {
    double best_lo = -(double)INFINITY;
    for (int b = 0; b < 4; b++) {
      double p = m->pwm_probs[i][b];
      double bg = args.bkg[b];
      if (p < 1e-12) p = 1e-12;
      if (bg < 1e-12) bg = 1e-12;
      double lo = log2(p / bg);
      if (lo > best_lo) best_lo = lo;
    }
    sum += best_lo;
  }
  return sum;
}

/* ---- Writers ---- */

static uint64_t emit_nsites(const motif_t *m) {
  return m->nsites > 0 ? m->nsites : args.nsites_fb;
}

static const char *fmt_name(int fmt) {
  switch (fmt) {
    case FMT_MEME:     return "MEME";
    case FMT_JASPAR:   return "JASPAR";
    case FMT_HOMER:    return "HOMER";
    case FMT_HOCOMOCO: return "HOCOMOCO";
  }
  return "?";
}

/* Print under -w: one line per motif summarizing what was parsed. */
static void report_motif_loaded(const motif_t *m) {
  char cons[MAX_MOTIF_WIDTH + 2];
  motif_consensus(m, cons);
  if (m->nsites > 0) {
    fprintf(stderr,
      "  motif '%s': width=%" PRIu64 " src_nsites=%" PRIu64 " IC=%.2f bits consensus=%s\n",
      m->name, m->size, m->nsites, motif_total_ic(m), cons);
  } else {
    fprintf(stderr,
      "  motif '%s': width=%" PRIu64 " src_nsites=unset IC=%.2f bits consensus=%s\n",
      m->name, m->size, motif_total_ic(m), cons);
  }
}

/* Print under -w: one line per motif when ic_trim_motif actually shrank it. */
static void report_motif_trimmed(const motif_t *m, uint64_t old_width) {
  fprintf(stderr,
    "  motif '%s': IC-trimmed %" PRIu64 " -> %" PRIu64 " columns (threshold %.2f bits)\n",
    m->name, old_width, m->size, args.ic_min);
}

/* Print under -w: one line per motif describing the emit nsites decision.
   Only meaningful for PCM-style targets. */
static void report_motif_emit(const motif_t *m) {
  uint64_t ns = emit_nsites(m);
  const char *src = (m->nsites > 0) ? "preserved from source" : "-N fallback";
  fprintf(stderr,
    "  emit '%s' -> %s with nsites=%" PRIu64 " (%s)\n",
    m->name, fmt_name(args.to_fmt), ns, src);
}

static void write_meme(void) {
  FILE *f = files.o;
  fprintf(f, "MEME version 4\n\n");
  fprintf(f, "ALPHABET= ACGT\n\n");
  fprintf(f, "strands: + -\n\n");
  fprintf(f, "Background letter frequencies\n");
  fprintf(f, "A %.6f C %.6f G %.6f T %.6f\n\n",
    args.bkg[0], args.bkg[1], args.bkg[2], args.bkg[3]);
  for (uint64_t mi = 0; mi < motif_info.n; mi++) {
    motif_t *m = motifs[mi];
    if (m->size == 0) continue;
    fprintf(f, "MOTIF %s\n", m->name);
    fprintf(f, "letter-probability matrix: alength= 4 w= %" PRIu64
               " nsites= %" PRIu64 " E= 0\n", m->size, emit_nsites(m));
    for (uint64_t i = 0; i < m->size; i++) {
      fprintf(f, " %.6f  %.6f  %.6f  %.6f\n",
        m->pwm_probs[i][0], m->pwm_probs[i][1],
        m->pwm_probs[i][2], m->pwm_probs[i][3]);
    }
    fprintf(f, "\n");
  }
}

static void write_jaspar(void) {
  FILE *f = files.o;
  const char *bases = "ACGT";
  for (uint64_t mi = 0; mi < motif_info.n; mi++) {
    motif_t *m = motifs[mi];
    if (m->size == 0) continue;
    uint64_t ns = emit_nsites(m);
    fprintf(f, ">%s\n", m->name);
    for (int b = 0; b < 4; b++) {
      fprintf(f, "%c [", bases[b]);
      for (uint64_t i = 0; i < m->size; i++) {
        long c = (long) (m->pwm_probs[i][b] * (double) ns + 0.5);
        fprintf(f, " %4ld", c);
      }
      fprintf(f, " ]\n");
    }
  }
}

static void write_homer(void) {
  FILE *f = files.o;
  char cons[MAX_MOTIF_WIDTH + 2];
  for (uint64_t mi = 0; mi < motif_info.n; mi++) {
    motif_t *m = motifs[mi];
    if (m->size == 0) continue;
    motif_consensus(m, cons);
    double thr = motif_max_logodds(m);
    fprintf(f, ">%s\t%s\t%.6f\n", cons, m->name, thr);
    for (uint64_t i = 0; i < m->size; i++) {
      fprintf(f, "%.3f\t%.3f\t%.3f\t%.3f\n",
        m->pwm_probs[i][0], m->pwm_probs[i][1],
        m->pwm_probs[i][2], m->pwm_probs[i][3]);
    }
  }
}

static void write_hocomoco(void) {
  FILE *f = files.o;
  for (uint64_t mi = 0; mi < motif_info.n; mi++) {
    motif_t *m = motifs[mi];
    if (m->size == 0) continue;
    uint64_t ns = emit_nsites(m);
    fprintf(f, ">%s\n", m->name);
    for (uint64_t i = 0; i < m->size; i++) {
      long c0 = (long) (m->pwm_probs[i][0] * (double) ns + 0.5);
      long c1 = (long) (m->pwm_probs[i][1] * (double) ns + 0.5);
      long c2 = (long) (m->pwm_probs[i][2] * (double) ns + 0.5);
      long c3 = (long) (m->pwm_probs[i][3] * (double) ns + 0.5);
      fprintf(f, "%ld %ld %ld %ld\n", c0, c1, c2, c3);
    }
  }
}

/* ---- Format-name parsing ---- */

static int parse_fmt_name(const char *s) {
  if (!s || !s[0]) return 0;
  char buf[16]; uint64_t i;
  for (i = 0; i < sizeof(buf) - 1 && s[i]; i++) buf[i] = (char) tolower((unsigned char) s[i]);
  buf[i] = '\0';
  if (strcmp(buf, "meme")     == 0) return FMT_MEME;
  if (strcmp(buf, "jaspar")   == 0) return FMT_JASPAR;
  if (strcmp(buf, "homer")    == 0) return FMT_HOMER;
  if (strcmp(buf, "hocomoco") == 0) return FMT_HOCOMOCO;
  return 0;
}

/* ---- Background parser (mirrors yamscan -b A,C,G,T) ---- */

static int parse_bkg(const char *s, double bkg[4]) {
  char buf[256];
  if (strlen(s) >= sizeof(buf)) return 1;
  strcpy(buf, s);
  char *tok = buf; int n = 0;
  while (*tok && n < 4) {
    char *comma = strchr(tok, ',');
    if (comma) *comma = '\0';
    if (str_to_double(tok, &bkg[n])) return 1;
    if (bkg[n] < 0.0 || bkg[n] > 1.0) return 1;
    n++;
    if (!comma) break;
    tok = comma + 1;
  }
  if (n != 4) return 1;
  double sum = bkg[0] + bkg[1] + bkg[2] + bkg[3];
  if (fabs(sum - 1.0) > 0.05) return 1;
  /* Normalize to exactly 1.0 to avoid drift. */
  for (int i = 0; i < 4; i++) bkg[i] /= sum;
  return 0;
}

/* ---- Usage ---- */

static void usage(void) {
  printf(
    "yamtk v%s  Copyright (C) %s  Benjamin Jean-Marie Tremblay\n"
    "Usage:  yamtk conv -m motifs.txt -t <fmt> [options]\n"
    "\n"
    " -m <str>   Input motif file (MEME/JASPAR/HOMER/HOCOMOCO; '-' = stdin).\n"
    " -t <str>   Target format: meme | jaspar | homer | hocomoco.\n"
    " -o <str>   Output file (default: stdout).\n"
    " -T <dbl>   IC threshold (bits) for flank trimming. When set, low-IC\n"
    "            flanking columns are trimmed before output. Suggested: 0.5.\n"
    " -N <int>   Fallback nsites for PCM output when the source has none\n"
    "            (HOMER input, or MEME without nsites=) (default: %d).\n"
    " -p <int>   Pseudocount applied during PCM->PPM ingestion of JASPAR/\n"
    "            HOCOMOCO sources (default: %d; pass >0 to soften zero\n"
    "            probabilities for downstream log-transformed scoring).\n"
    " -b A,C,G,T Background for the HOMER threshold (default: uniform).\n"
    " -r         Do not trim motif names to the first word.\n"
    " -v / -w / -h   Verbose / very-verbose / help.\n"
    , YAMTK_VERSION, YAMTK_YEAR, DEFAULT_NSITES, DEFAULT_PSEUDOCOUNT
  );
}

/* ---- Main ---- */

int main_conv(int argc, char **argv) {

  int opt;
  int use_stdout = 1;
  int use_stdin  = 0;

  while ((opt = getopt(argc, argv, "m:t:o:T:N:p:b:rvwh")) != -1) {
    switch (opt) {
      case 'm':
        if (files.m_open || use_stdin) badexit("Error: -m specified more than once.");
        if (optarg[0] == '-' && optarg[1] == '\0') {
          files.m = stdin;
          use_stdin = 1;
        } else {
          files.m = fopen(optarg, "r");
          if (files.m == NULL) {
            fprintf(stderr, "Error: Failed to open motif file \"%s\" [%s]\n",
              optarg, strerror(errno));
            badexit("");
          }
          files.m_open = 1;
        }
        break;
      case 't': {
        const int f = parse_fmt_name(optarg);
        if (!f) {
          fprintf(stderr, "Error: -t must be one of meme, jaspar, homer, hocomoco (got '%s').\n",
            optarg);
          badexit("");
        }
        args.to_fmt = f;
        break;
      }
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
      case 'T':
        if (str_to_double(optarg, &args.ic_min)) badexit("Error: Failed to parse -T value.");
        if (args.ic_min < 0.0 || args.ic_min > 2.0)
          badexit("Error: -T must be in [0, 2] bits.");
        args.do_ic_trim = 1;
        break;
      case 'N':
        if (str_to_uint64_t(optarg, &args.nsites_fb)) badexit("Error: Failed to parse -N value.");
        if (args.nsites_fb == 0) badexit("Error: -N must be > 0.");
        break;
      case 'p':
        if (str_to_int(optarg, &args.pseudocount)) badexit("Error: Failed to parse -p value.");
        if (args.pseudocount < 0) badexit("Error: -p must be >= 0.");
        break;
      case 'b':
        if (parse_bkg(optarg, args.bkg))
          badexit("Error: -b must be four comma-separated probabilities summing to 1.");
        args.bkg_uniform = 0;
        break;
      case 'r':
        args.trim_names = 0;
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

  if (!files.m_open && !use_stdin) badexit("Error: -m must be specified.");
  if (!args.to_fmt) badexit("Error: -t must be specified.");
  if (use_stdout) files.o = stdout;

  /* Load motifs (auto-detects input format) */
  load_motifs();

  if (args.w) {
    fprintf(stderr,
      "Converting from %s to %s; bkg=A:%.3f C:%.3f G:%.3f T:%.3f; pseudocount=%d; -N fallback=%" PRIu64 ".\n",
      fmt_name(motif_info.fmt), fmt_name(args.to_fmt),
      args.bkg[0], args.bkg[1], args.bkg[2], args.bkg[3],
      args.pseudocount, args.nsites_fb);
    for (uint64_t i = 0; i < motif_info.n; i++) report_motif_loaded(motifs[i]);
  }

  /* Optional IC-based flank trim */
  if (args.do_ic_trim) {
    for (uint64_t i = 0; i < motif_info.n; i++) {
      const uint64_t before = motifs[i]->size;
      ic_trim_motif(motifs[i], args.ic_min);
      if (args.w && motifs[i]->size != before) report_motif_trimmed(motifs[i], before);
    }
  }

  /* Per-motif emit info, useful for PCM targets to see nsites source. */
  if (args.w && (args.to_fmt == FMT_JASPAR || args.to_fmt == FMT_HOCOMOCO)) {
    for (uint64_t i = 0; i < motif_info.n; i++) report_motif_emit(motifs[i]);
  }

  /* Dispatch writer */
  switch (args.to_fmt) {
    case FMT_MEME:     write_meme();     break;
    case FMT_JASPAR:   write_jaspar();   break;
    case FMT_HOMER:    write_homer();    break;
    case FMT_HOCOMOCO: write_hocomoco(); break;
  }

  if (args.v) {
    fprintf(stderr, "Converted %" PRIu64 " motif(s).\n", motif_info.n);
  }

  free_motifs();
  close_files();
  return EXIT_SUCCESS;
}
