/*
 *   yamdedup: Deduplicate overlapping hits from yamscan.
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
#include <locale.h>
#include <string.h>
#include <errno.h>
#include <time.h>
#include <math.h>
#include <zlib.h>
#include "kseq.h"
#include "khash.h"

KSEQ_INIT(gzFile, gzread)
KHASH_MAP_INIT_STR(str_h, size_t)

#define YAMDEDUP_VERSION            "1.0"
#define YAMDEDUP_YEAR                2022

/* Maximum number of characters allowed for motif names */
#define MAX_NAME_SIZE             ((size_t) 256)
/* Maximum number of characters allowed for sequence names */
#define SEQ_NAME_MAX_CHAR         ((size_t) 512)
/* Maximum number of characters allowed for individual fields */
#define FIELD_MAX_CHAR            ((size_t) 512)
/* Maximum number of characters allowed for entire lines */
#define LINE_MAX_CHAR         ((size_t) 1048576)
/* Number of elements when calling malloc/realloc */
#define ALLOC_CHUNK_SIZE          ((size_t) 256)
/* Print a progress message (when -v is used) every N ranges */
#define PROGRESS_TRIGGER       ((size_t) 500000)
#define HASH_TABLE_SIZE          ((size_t) 1024)
#define HASH_TABLE_GROW          ((size_t) 1024)

#define ERASE_ARRAY(ARR, LEN) memset(ARR, 0, sizeof(ARR[0]) * (LEN))
#define MAX(x, y) (((x) > (y)) ? (x) : (y))

#define MALLOC_OR_RET1(OBJ, SIZE)                    \
  do {                                               \
    OBJ = malloc(SIZE);                              \
    if (OBJ == NULL)  return 1;                      \
  } while (0)

#define REALLOC_OR_RET1(OBJ, SIZE, TYPE)             \
  do {                                               \
    TYPE* TMP = realloc(OBJ, (SIZE) * sizeof(TYPE)); \
    if (TMP == NULL) return 1;                       \
    OBJ = TMP;                                       \
  } while (0)

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
long peak_mem(void) {
  return 0;
}
#else
#include <sys/resource.h>
long peak_mem(void) {
  struct rusage r_mem;
  getrusage(RUSAGE_SELF, &r_mem);
#ifdef __linux__
  return r_mem.ru_maxrss * 1024;
#else
  return r_mem.ru_maxrss;
#endif
}
#endif

#ifdef DEBUG
size_t malloc_count = 0;
size_t malloc_fun_counts = 0;
size_t realloc_count = 0;
size_t realloc_fun_counts = 0;
void print_mem_alloc_counts(int malloc, int realloc) {
  malloc_count += malloc;
  realloc_count += realloc;
  if (malloc) malloc_fun_counts += 1;
  if (realloc) realloc_fun_counts += 1;
  fprintf(stderr, "[DEBUG] malloc count: %'zu (%'zu)  realloc count: %'zu (%'zu)\n",
      malloc_fun_counts, malloc_count, realloc_fun_counts, realloc_count);
}
#endif

void print_peak_mb(void) {
  long bytes = peak_mem();
  if (bytes > (1 << 30)) {
    fprintf(stderr, "Approx. peak memory usage: %'.2f GB.\n",
      (((double) bytes / 1024.0) / 1024.0) / 1024.0);
  } else if (bytes > (1 << 20)) {
    fprintf(stderr, "Approx. peak memory usage: %'.2f MB.\n",
      ((double) bytes / 1024.0) / 1024.0);
  } else if (bytes) {
    fprintf(stderr, "Approx. peak memory usage: %'.2f KB.\n",
      (double) bytes / 1024.0);
  }
}

/* TODO: Give the option to specify what constitutes an overlap. */

void usage(void) {
  printf(
    "yamdedup v%s  Copyright (C) %d  Benjamin Jean-Marie Tremblay               \n"
    "                                                                              \n"
    "Usage:  yamdedup [options] -i [ results.txt | ranges.bed ]                    \n"
    "                                                                              \n"
    " -i <str>   Filename of yamscan results file or a tab-delimited BED file      \n"
    "            with at least six columns: (1) sequence name, (2) 0-based start,  \n"
    "            (3) end, (4) motif name, (5) motif score, and (6) the strand.     \n"
    "            If a yamscan results file is provided the P-value is used for     \n"
    "            deciding which overlapping hit(s) to discard, otherwise for BED   \n"
    "            the score column is used and lower scores are discarded. The input\n"
    "            is assumed to be partially sorted: different sequence/motif/strand\n"
    "            combinations can be interleaved, but the individual combinations  \n"
    "            themselves must be sorted by their start coordinates. The         \n"
    "            yamscan program outputs its results this way, so no additional    \n"
    "            sorting is needed. Can be gzipped. Use '-' for stdin.             \n"
    " -o <str>   Filename to output results. By default output goes to stdout.     \n"
    /*
    " -f <num>   Required minimum overlap of the smaller range. Use a value between\n"
    "            0 and 1 for proportions, an integer >=1 for absolute overlaps, and\n"
    "            -1 to require a complete overlap. If two ranges are the same size \n"
    "            then the value set with -f is used for both. The default is a     \n"
    "            minimum overlap of 1.                                             \n"
    " -F <num>   Required minimum overlap of the larger range. Use a value between \n"
    "            0 and 1 for proportions, an integer >=1 for absolute overlaps, and\n"
    "            -1 to require a complete overlap. The default is a minimum overlap\n"
    "            of 1.                                                             \n"
    */
    " -s         Ignore strand when finding overlapping ranges.                    \n"
    " -m         Ignore motif name when finding overlapping ranges.                \n"
    " -0         Ignore scores when removing overlapping ranges, causing yamdedup  \n"
    "            to simply remove overlapping ranges in the order they appear.     \n"
    " -S         Sort on range size instead of score/p-value (keeping larger ones).\n"
    " -r         Sort in opposite order (i.e., keep lower scores, higher p-values, \n"
    "            or smaller ranges).                                               \n"
    " -y         Force yamdedup to treat the input as a yamscan output file.       \n"
    " -b         Force yamdedup to treat the input as a BED file.                  \n"
    " -v         Verbose mode.                                                     \n"
    " -w         Very verbose mode.                                                \n"
    " -h         Print this help message.                                          \n"
    , YAMDEDUP_VERSION, YAMDEDUP_YEAR
  );
}

khash_t(str_h) *hash_tab;

typedef struct feat_tab_t {
  char  ***lines;
  char   **seq_name;
  char   **motif_name;
  char    *strand;
  size_t **starts;
  size_t **ends;
  double **scores;
  size_t  *n;
  size_t  *n_alloc;
  size_t   n_total;
  size_t   n_total_alloc;
} feat_tab_t;

feat_tab_t feat_tab = {
  .n_total = 0,
  .n_total_alloc = 0
};

void free_feat_tab(void) {
  if (feat_tab.n_total_alloc) {
    for (size_t i = 0; i < feat_tab.n_total; i++) {
      if (feat_tab.n_alloc[i]) {
        if (feat_tab.n[i]) {
          for (size_t j = 0; j < feat_tab.n[i]; j++) {
            free(feat_tab.lines[i][j]);
          }
        }
        free(feat_tab.lines[i]);
        free(feat_tab.seq_name[i]);
        free(feat_tab.motif_name[i]);
        free(feat_tab.starts[i]);
        free(feat_tab.ends[i]);
        free(feat_tab.scores[i]);
      }
    }
  }
  free(feat_tab.seq_name);
  free(feat_tab.motif_name);
  free(feat_tab.strand);
  free(feat_tab.starts);
  free(feat_tab.ends);
  free(feat_tab.scores);
  free(feat_tab.n);
  free(feat_tab.n_alloc);
  for (khint_t k = 0; k < kh_end(hash_tab); k++) {
    if (kh_exist(hash_tab, k)) free((char *) kh_key(hash_tab, k));
  }
  kh_destroy(str_h, hash_tab);
}

static inline int alloc_more_to_feat_one(const size_t i) {
  if (i < feat_tab.n_total) {
    REALLOC_OR_RET1(feat_tab.lines[i], ALLOC_CHUNK_SIZE + feat_tab.n_alloc[i], char *);
    REALLOC_OR_RET1(feat_tab.starts[i], ALLOC_CHUNK_SIZE + feat_tab.n_alloc[i], size_t);
    REALLOC_OR_RET1(feat_tab.ends[i], ALLOC_CHUNK_SIZE + feat_tab.n_alloc[i], size_t);
    REALLOC_OR_RET1(feat_tab.scores[i], ALLOC_CHUNK_SIZE + feat_tab.n_alloc[i], double);
  } else {
    feat_tab.lines[i] = malloc(sizeof(char *) * ALLOC_CHUNK_SIZE);
    if (feat_tab.lines[i] == NULL) {
      return 1;
    }
    feat_tab.seq_name[i] = malloc(sizeof(char) * SEQ_NAME_MAX_CHAR);
    if (feat_tab.seq_name[i] == NULL) {
      return 1;
    }
    feat_tab.motif_name[i] = malloc(sizeof(char) * MAX_NAME_SIZE);
    if (feat_tab.motif_name[i] == NULL) {
      free(feat_tab.seq_name[i]);
      return 1;
    }
    feat_tab.starts[i] = malloc(sizeof(size_t) * ALLOC_CHUNK_SIZE);
    if (feat_tab.starts[i] == NULL) {
      free(feat_tab.seq_name[i]);
      free(feat_tab.motif_name[i]);
      return 1;
    }
    feat_tab.ends[i] = malloc(sizeof(size_t) * ALLOC_CHUNK_SIZE);
    if (feat_tab.ends[i] == NULL) {
      free(feat_tab.seq_name[i]);
      free(feat_tab.motif_name[i]);
      free(feat_tab.starts[i]);
      return 1;
    }
    feat_tab.scores[i] = malloc(sizeof(double) * ALLOC_CHUNK_SIZE);
    if (feat_tab.scores[i] == NULL) {
      free(feat_tab.seq_name[i]);
      free(feat_tab.motif_name[i]);
      free(feat_tab.starts[i]);
      free(feat_tab.ends[i]);
      return 1;
    }
    feat_tab.n[i] = 0;
    feat_tab.n_alloc[i] = 0;
  }
  feat_tab.n_alloc[i] += ALLOC_CHUNK_SIZE;
  return 0;
}

static inline int alloc_more_to_feat_tab(void) {
  if (feat_tab.n_total_alloc) {
    REALLOC_OR_RET1(feat_tab.lines, ALLOC_CHUNK_SIZE + feat_tab.n_total_alloc, char **);
    REALLOC_OR_RET1(feat_tab.seq_name, ALLOC_CHUNK_SIZE + feat_tab.n_total_alloc, char *);
    REALLOC_OR_RET1(feat_tab.motif_name, ALLOC_CHUNK_SIZE + feat_tab.n_total_alloc, char *);
    REALLOC_OR_RET1(feat_tab.strand, ALLOC_CHUNK_SIZE + feat_tab.n_total_alloc, char);
    REALLOC_OR_RET1(feat_tab.starts, ALLOC_CHUNK_SIZE + feat_tab.n_total_alloc, size_t *);
    REALLOC_OR_RET1(feat_tab.ends, ALLOC_CHUNK_SIZE + feat_tab.n_total_alloc, size_t *);
    REALLOC_OR_RET1(feat_tab.scores, ALLOC_CHUNK_SIZE + feat_tab.n_total_alloc, double *);
    REALLOC_OR_RET1(feat_tab.n, ALLOC_CHUNK_SIZE + feat_tab.n_total_alloc, size_t);
    REALLOC_OR_RET1(feat_tab.n_alloc, ALLOC_CHUNK_SIZE + feat_tab.n_total_alloc, size_t);
  } else {
    MALLOC_OR_RET1(feat_tab.lines, sizeof(char **) * ALLOC_CHUNK_SIZE);
    MALLOC_OR_RET1(feat_tab.seq_name, sizeof(char *) * ALLOC_CHUNK_SIZE);
    MALLOC_OR_RET1(feat_tab.motif_name, sizeof(char *) * ALLOC_CHUNK_SIZE);
    MALLOC_OR_RET1(feat_tab.strand, sizeof(char) * ALLOC_CHUNK_SIZE);
    MALLOC_OR_RET1(feat_tab.starts, sizeof(size_t *) * ALLOC_CHUNK_SIZE);
    MALLOC_OR_RET1(feat_tab.ends, sizeof(size_t *) * ALLOC_CHUNK_SIZE);
    MALLOC_OR_RET1(feat_tab.scores, sizeof(double *) * ALLOC_CHUNK_SIZE);
    MALLOC_OR_RET1(feat_tab.n, sizeof(size_t) * ALLOC_CHUNK_SIZE);
    MALLOC_OR_RET1(feat_tab.n_alloc, sizeof(size_t) * ALLOC_CHUNK_SIZE);
  }
  feat_tab.n_total_alloc += ALLOC_CHUNK_SIZE;
  return 0;
}

typedef struct args_t {
  int ignore_strand : 1;
  int ignore_motif : 1;
  int ignore_score : 1;
  int reverse_sort : 1;
  int sort_sizes : 1;
  int override_is_yamscan : 1;
  int override_is_bed : 1;
  int v : 1;
  int w : 1;
} args_t;

args_t args = {
  .ignore_strand       = 0,
  .ignore_motif        = 0,
  .ignore_score        = 0,
  .reverse_sort        = 0,
  .sort_sizes          = 0,
  .override_is_yamscan = 0,
  .override_is_bed     = 0,
  .v                   = 0,
  .w                   = 0
};

typedef struct files_t {
  int     i_open : 1;
  int     o_open : 1;
  gzFile  i;
  FILE   *o;
} files_t;

files_t files = {
  .i_open = 0,
  .o_open = 0
};

void close_files(void) {
  if (files.i_open) gzclose(files.i);
  if (files.o_open) fclose(files.o);
}

typedef struct score2index_t {
  double score;
  size_t index;
  int active;
} score2index_t;

score2index_t *score2index;
size_t score2index_n = 0;
size_t score2index_n_alloc = 0;

void badexit(const char *msg) {
  fprintf(stderr, "%s\nRun yamdedup -h to see usage.\n", msg);
  close_files();
  free(score2index);
  free_feat_tab();
  exit(EXIT_FAILURE);
}

static inline void print_progress(const size_t n_ranges, const size_t n_discarded) {
  if (n_ranges % PROGRESS_TRIGGER == 0) {
    fprintf(stderr, "Processed %'zu ranges (%'zu discarded) ...\n",
        n_ranges, n_discarded);
  }
}

static inline size_t count_nonempty_chars(const char *line) {
  size_t total_chars = 0;
  for (size_t i = 0; i < LINE_MAX_CHAR; i++) {
    switch (line[i]) {
      case ' ':
      case '\t':
      case '\r':
      case '\v':
      case '\f':
      case '\n': break;
      case '\0': return total_chars;
      default: total_chars++;
    }
  }
  return total_chars;
}

static inline size_t count_fields(const char *line) {
  int res = 1;
  for (size_t i = 0; i < LINE_MAX_CHAR; i++) {
    if (line[i] == '\0') break;
    if (line[i] == '\t') res += 1;
  }
  return res;
}

int compare_scores(const void *a, const void *b) {
  const score2index_t *a_s = (score2index_t *) a;
  const score2index_t *b_s = (score2index_t *) b;
  if (a_s->score < b_s->score) {
    return 1;
  } else if (a_s->score > b_s->score) {
    return -1;
  } else {
    return 0;
  }
}

static inline int ranges_are_overlapping(const size_t start1, const size_t end1, const size_t start2, const size_t end2) {
  if ( 
      ((start1 >= start2) && (start1 <= end2)) ||
      ((end1   >= start2) && (end1   <= end2)) ||
      ((start1 <  start2) && (end1   >  end2))
     ) {
    return 1;
  } else {
    return 0;
  }
}

static inline size_t dedup_and_purge_feat(const size_t i) {
#ifdef DEBUG
    fprintf(stderr, "[DEBUG] Found non-overlapping range, triggering deduplication for:\n");
    fprintf(stderr, "[DEBUG]    %s(%c): %s\n", feat_tab.seq_name[i], feat_tab.strand[i], feat_tab.motif_name[i]);
#endif
  size_t n_discarded = 0;
  if (feat_tab.n[i] == 2) {
    if (feat_tab.scores[i][1] > feat_tab.scores[i][0]) {
      free(feat_tab.lines[i][0]);
      feat_tab.lines[i][0] = NULL;
    } else {
      free(feat_tab.lines[i][1]);
      feat_tab.lines[i][1] = NULL;
    }
    n_discarded = 1;
  } else if (feat_tab.n[i] > 2) {
    for (size_t j = 0; j < feat_tab.n[i]; j++) {
      score2index[j].score = feat_tab.scores[i][j];
      score2index[j].index = j;
      score2index[j].active = 1;
    }
    qsort(score2index, feat_tab.n[i], sizeof(score2index_t), compare_scores);
    for (size_t j = 1; j < feat_tab.n[i]; j++) {
      for (size_t k = 0; k < j; k++) {
        if (score2index[k].active &&
            ranges_are_overlapping(
              feat_tab.starts[i][score2index[j].index],
              feat_tab.ends[i][score2index[j].index],
              feat_tab.starts[i][score2index[k].index],
              feat_tab.ends[i][score2index[k].index])) {
          score2index[j].active = 0;
          free(feat_tab.lines[i][score2index[j].index]);
          feat_tab.lines[i][score2index[j].index] = NULL;
          n_discarded += 1;
          break;
        }
      }
    }
  }
  for (size_t j = 0; j < feat_tab.n[i]; j++) {
    if (feat_tab.lines[i][j] != NULL) {
      fprintf(files.o, "%s\n", feat_tab.lines[i][j]);
      free(feat_tab.lines[i][j]);
    }
  }
  feat_tab.n[i] = 0;
  return n_discarded;
}

static inline size_t find_matching_feat(const char *seq, const char *motif, const char strand) {
  size_t feat_index = 0, hash_key_i = 1, j, hash_key_n = 2 + MAX_NAME_SIZE + SEQ_NAME_MAX_CHAR;
  char hash_key[hash_key_n];
  hash_key[0] = strand;
  for (j = 0; hash_key_i < hash_key_n && j < MAX_NAME_SIZE; hash_key_i++, j++) {
    if (motif[j] == '\0') break;
    hash_key[hash_key_i] = motif[j];
  }
  for (j = 0; hash_key_i < hash_key_n && j < SEQ_NAME_MAX_CHAR; hash_key_i++, j++) {
    if (seq[j] == '\0') break;
    hash_key[hash_key_i] = seq[j];
  }
  hash_key[hash_key_i] = '\0';
#ifdef DEBUG
  fprintf(stderr, "HashKey:%s\n", hash_key);
#endif
  int absent;
  khint_t k = kh_put(str_h, hash_tab, hash_key, &absent);
  if (absent) {
    kh_key(hash_tab, k) = strdup(hash_key);
    kh_val(hash_tab, k) = feat_tab.n_total;
  } else {
    return kh_val(hash_tab, k);
  }
  if (args.w) {
    fprintf(stderr,
      "Found new seq/strand/motif combination (#%'zu):\n    %s(%c): %s\n",
      feat_tab.n_total + 1, seq, strand, motif);
  }
  if (feat_tab.n_total + 1 > feat_tab.n_total_alloc && alloc_more_to_feat_tab()) {
    return -1;
  }
  if (alloc_more_to_feat_one(feat_tab.n_total)) {
    return -1;
  }
  feat_tab.strand[feat_tab.n_total] = strand;
  for (size_t i = 0; i < SEQ_NAME_MAX_CHAR; i++) {
    feat_tab.seq_name[feat_tab.n_total][i] = seq[i];
    if (seq[i] == '\0') break;
  }
  for (size_t i = 0; i < MAX_NAME_SIZE; i++) {
    feat_tab.motif_name[feat_tab.n_total][i] = motif[i];
    if (motif[i] == '\0') break;
  }
  feat_index = feat_tab.n_total;
  feat_tab.n_total += 1;
  return feat_index;
}

static inline int push_feat_tab(const size_t i, char *line, const size_t start, const size_t end, const double score) {
  if (feat_tab.n_alloc[i] < feat_tab.n[i] + 1) {
    if (alloc_more_to_feat_one(i)) return 1;
  }
  feat_tab.scores[i][feat_tab.n[i]] = score;
  feat_tab.starts[i][feat_tab.n[i]] = start;
  feat_tab.ends[i][feat_tab.n[i]] = end;
  feat_tab.lines[i][feat_tab.n[i]] = line;
  feat_tab.n[i] += 1;
  return 0;
}

static inline size_t extract_field(const char *line, const size_t k, char *field) {
  size_t start = 0, end = 0, field_i = 1, size = 0;
  ERASE_ARRAY(field, FIELD_MAX_CHAR);
  for (size_t i = 0; i < LINE_MAX_CHAR; i++) {
    if (line[i] == '\0' || (line[i] == '\t' && field_i == k)) {
      end = i - 1;
      break;
    } else if (line[i] == '\t') {
      field_i += 1;
    } else if (field_i == k) {
      if (!size) start = i;
      size += 1;
    }
  }
  if (size > 0 && size < FIELD_MAX_CHAR) {
    for (size_t i = start, j = 0; i <= end; i++, j++) {
      field[j] = line[i];
    }
  }
  return size;
}

static inline int check_field_size(const size_t size, const size_t line_num) {
  if (!size) {
    fprintf(stderr, "Error: Found empty field on line #%'zu.", line_num);
    return 1;
  } else if (size > FIELD_MAX_CHAR) {
    fprintf(stderr, "Error: Field on line #%'zu exceeds max allowed size (%'zu>%'zu).",
      line_num, size, FIELD_MAX_CHAR);
    return 1;
  }
  return 0;
}

static inline int check_strand(const char *strand) {
  if ((strand[0] != '+' && strand[0] != '-' && strand[0] != '.') || strand[1] != '\0') {
    return 1;
  }
  return 0;
}

static inline size_t purge_remaining_feat_tab(void) {
  size_t n_discarded = 0;
  for (size_t i = 0; i < feat_tab.n_total; i++) {
    n_discarded += dedup_and_purge_feat(i);
  }
  return n_discarded;
}

static inline int safe_strtod(char *str, double *res) {
  char *tmp; errno = 0;
  *res = strtod(str, &tmp);
  if (str == tmp || errno != 0 || *tmp != '\0') {
    return 1;
  } else {
    return 0;
  }
}

static inline int safe_strtoull(char *str, size_t *res) {
  char *tmp; errno = 0;
  *res = (size_t) strtoull(str, &tmp, 10);
  if (str == tmp || errno != 0 || *tmp != '\0') {
    return 1;
  } else {
    return 0;
  }
}

void run_minidedup(void) {
  int ret_val;
  int is_yamscan = args.override_is_yamscan;
  int is_bed = args.override_is_bed;
  int is_yamscan_bed = 0;
  const int ignore_strand = args.ignore_strand;
  const int ignore_motif = args.ignore_motif;
  const int ignore_score = args.ignore_score;
  const int reverse_sort = args.reverse_sort;
  const int sort_sizes = args.sort_sizes;
  size_t n_ranges = 0, n_discarded = 0, n_comments = 0, n_empty = 0;
  size_t n_lines = 0, n_fields = 0, field_size = 0;
  size_t start, end;
  double score;
  char seq[FIELD_MAX_CHAR], motif[FIELD_MAX_CHAR], strand[FIELD_MAX_CHAR];
  char str_start[FIELD_MAX_CHAR], str_end[FIELD_MAX_CHAR], str_score[FIELD_MAX_CHAR];
  size_t start_loc, end_loc, seq_loc, motif_loc, strand_loc, score_loc;
  size_t feature_index;
  kstream_t *kinput = ks_init(files.i);
  kstring_t line = { 0, 0, 0 };
  while ((ret_val = ks_getuntil(kinput, '\n', &line, 0)) >= 0) {
    n_lines += 1;
    if (line.l > LINE_MAX_CHAR) {
      fprintf(stderr, "Error: Line #%'zu exceeded max allowed line length (%'zu>%'zu).\n",
        n_lines, line.l, LINE_MAX_CHAR);
      goto error_blank;
    }
    if (!is_bed && !feat_tab.n_total && line.l >= 9 &&
        line.s[0] == '#' && line.s[1] == '#' && line.s[2] == 'y' && line.s[3] == 'a' &&
        line.s[4] == 'm' && line.s[5] == 's' && line.s[6] == 'c' && line.s[7] == 'a' &&
        line.s[8] == 'n') {
      is_yamscan = 1;
      n_comments += 1;
      if (args.v && !args.override_is_yamscan) {
        fprintf(stderr, "Treating input as yamscan output.\n");
      }
      fprintf(files.o, "%s\n", line.s);
      continue;
    } else if (count_nonempty_chars(line.s) == 0) {
      n_empty += 1;
      continue;
    } else if (line.s[0] == '#') {
      n_comments += 1;
      fprintf(files.o, "%s\n", line.s);
      continue;
    } else if (!is_yamscan && line.l >= 7 && line.s[0] == 'b' && line.s[1] == 'r' &&
        line.s[2] == 'o' && line.s[3] == 'w' && line.s[4] == 's' && line.s[5] == 'e' &&
        line.s[6] == 'r') {
      n_comments += 1;
      fprintf(files.o, "%s\n", line.s);
      continue;
    } else if (!is_yamscan && line.l >= 5 && line.s[0] == 't' && line.s[1] == 'r' &&
        line.s[2] == 'a' && line.s[3] == 'c' && line.s[4] == 'k') {
      n_comments += 1;
      fprintf(files.o, "%s\n", line.s);
      continue;
    }
    n_fields = count_fields(line.s);
    if (is_yamscan && n_fields < 9) {
      fprintf(stderr, "Error: Found too few fields at line %'zu; expect 9-12 for yamscan output.",
        n_lines);
      goto error_blank;
    } else if (n_fields < 4) {
      fprintf(stderr, "Error: Found too few fields at line %'zu; expect at least 4 for BED or 9 for yamscan output.",
        n_lines);
      goto error_blank;
    }
    if (!is_yamscan && !is_bed) {
      is_bed = 1;
      if (args.v && !args.override_is_bed) {
        fprintf(stderr, "Treating input as BED.\n");
      }
    }
    if (is_yamscan) {
      if (n_fields == 9 || n_fields == 10) {
        if (is_yamscan_bed) {
          fprintf(stderr, "Error: Found 9-10 fields on line %'zu, but previous lines had 11-12.",
            n_lines);
          goto error_blank;
        }
        seq_loc = 1; start_loc = 2; end_loc = 3; strand_loc = 4; motif_loc = 5; score_loc = 6;
      } else if (n_fields == 11 || n_fields == 12) {
        if (!is_yamscan_bed && n_ranges) {
          fprintf(stderr, "Error: Found 11-12 fields on line %'zu, but previous lines had 9-10.",
            n_lines);
          goto error_blank;
        }
        seq_loc = 3; start_loc = 4; end_loc = 5; strand_loc = 6; motif_loc = 7; score_loc = 8;
      } else {
        fprintf(stderr, "Error: Found %'zu fields on line %'zu; expect 9-12 for yamscan output.",
          n_fields, n_lines);
        goto error_blank;
      }
    } else {
      seq_loc = 1; start_loc = 2; end_loc = 3; motif_loc = 4; score_loc = 5; strand_loc = 6;
    }
    if (check_field_size((field_size = extract_field(line.s, seq_loc, seq)), n_lines)) {
      goto error_blank;
    }
    if (check_field_size((field_size = extract_field(line.s, motif_loc, motif)), n_lines)) {
      goto error_blank;
    }
    if (check_field_size((field_size = extract_field(line.s, strand_loc, strand)), n_lines)) {
      goto error_blank;
    }
    if (check_field_size((field_size = extract_field(line.s, start_loc, str_start)), n_lines)) {
      goto error_blank;
    }
    if (check_field_size((field_size = extract_field(line.s, end_loc, str_end)), n_lines)) {
      goto error_blank;
    }
    if (check_field_size((field_size = extract_field(line.s, score_loc, str_score)), n_lines)) {
      goto error_blank;
    }
    if (safe_strtoull(str_start, &start)) {
      fprintf(stderr, "Error: Failed to parse number in start column on line %'zu; found: %s",
        n_lines, str_start);
      goto error_blank;
    }
    if (safe_strtoull(str_end, &end)) {
      fprintf(stderr, "Error: Failed to parse number in end column on line %'zu; found: %s",
        n_lines, str_end);
      goto error_blank;
    }
    if (safe_strtod(str_score, &score)) {
      fprintf(stderr, "Error: Failed to parse number in scores column on line %'zu; found: %s",
        n_lines, str_score);
      goto error_blank;
    }
    if (is_yamscan) {
      score = -log(score);
    } else {
      start += 1;
    }
    if (start > end) {
      fprintf(stderr, "Error: Incorrect start/end values on line %'zu (start: %'zu, end: %'zu)",
        n_lines, start, end);
      goto error_blank;
    }
    if (check_strand(strand)) {
      fprintf(stderr, "Error: Incorrect strand column on line %'zu; expect +/-/., found: %s",
        n_lines, strand);
      goto error_blank;
    }
    if (ignore_strand) strand[0] = '.';
    if (ignore_motif) {
      motif[0] = '.'; motif[1] = '\0';
    }
    if (sort_sizes) score = (double) (end - start + 1);
    if (ignore_score) score = 0;
    if (reverse_sort) score = -score;
    if ((feature_index = find_matching_feat(seq, motif, strand[0])) == -1) {
      goto error_mem;
    }
#ifdef DEBUG
    fprintf(stderr, "L:%zu\tI:%zu\tN:%zu\n", n_lines, feature_index, feat_tab.n[feature_index]);
#endif
    if (feat_tab.n[feature_index] > 0 &&
      start > feat_tab.ends[feature_index][feat_tab.n[feature_index] - 1]) {
      if (feat_tab.n[feature_index] > score2index_n_alloc) {
        score2index_t *tmp = realloc(score2index,
          sizeof(score2index) * (score2index_n_alloc + ALLOC_CHUNK_SIZE));
        if (score2index == NULL) goto error_mem;
        score2index = tmp;
        score2index_n_alloc += ALLOC_CHUNK_SIZE;
      }
      n_discarded += dedup_and_purge_feat(feature_index);
    } else if (feat_tab.n[feature_index] &&
        start < feat_tab.starts[feature_index][feat_tab.n[feature_index] - 1]) {
      fprintf(stderr, "Error: Input isn't properly sorted (line %'zu).\n", n_lines);
      fprintf(stderr, "Found the following order of ranges:\n%s\n%s",
        feat_tab.lines[feature_index][feat_tab.n[feature_index] - 1], line.s);
      goto error_blank;
    }
    if (push_feat_tab(feature_index, line.s, start, end, score)) {
      goto error_mem;
    }
    line.s = NULL; line.m = 0; line.l = 0;
    n_ranges += 1;
    if (args.v && n_ranges) print_progress(n_ranges, n_discarded);
  }
  for (size_t i = 0; i < feat_tab.n_total; i++) {
    if (feat_tab.n[i] > score2index_n_alloc) {
      score2index_t *tmp = realloc(score2index,
        sizeof(score2index) * (score2index_n_alloc + ALLOC_CHUNK_SIZE));
      if (score2index == NULL) goto error_mem;
      score2index = tmp;
      score2index_n_alloc += ALLOC_CHUNK_SIZE;
    }
  }
  n_discarded += purge_remaining_feat_tab();
  if (args.v) {
    fprintf(stderr,
      "Done. Total ranges: %'zu\nRemaining ranges: %'zu (%'.2f%%)\nDiscarded ranges: %'zu (%'.2f%%)\n",
      n_ranges, n_ranges - n_discarded,
      100.0 * ((double) (n_ranges - n_discarded) / (double) n_ranges),
      n_discarded,
      100.0 * ((double) n_discarded / (double) n_ranges));
  }
  free(line.s);
  ks_destroy(kinput);
  if (ret_val == -3) {
    badexit("Error: Failed to read file stream.");
  }
  return;

error_mem:
  fprintf(stderr, "Error: Failed memory allocation.");
error_blank:
  free(line.s);
  ks_destroy(kinput);
  badexit("");
}

void print_time(const size_t s, const char *what) {
  if (s > 7200) {
    fprintf(stderr, "Needed %'.2f hours to %s.\n", ((double) s / 60.0) / 60.0, what);
  } else if (s > 120) {
    fprintf(stderr, "Needed %'.2f minutes to %s.\n", (double) s / 60.0, what);
  } else if (s > 1) {
    fprintf(stderr, "Needed %'zu seconds to %s.\n", s, what);
  }
}

int main(int argc, char **argv) {

  int opt;
  int use_stdout = 1, has_input = 0;

  while ((opt = getopt(argc, argv, "i:o:smwvybr0Sh")) != -1) {
    switch (opt) {
      case 'i':
        has_input = 1;
        if (optarg[0] == '-' && optarg[1] == '\0') {
          files.i = gzdopen(fileno(stdin), "r");
        } else {
          files.i = gzopen(optarg, "r");
          if (files.i == NULL) {
            fprintf(stderr, "Error: Failed to open input file \"%s\" [%s]",
              optarg, strerror(errno));
            badexit("");
          }
        }
        files.i_open = 1;
        break;
      case 'o':
        use_stdout = 0;
        files.o = fopen(optarg, "w");
        if (files.o == NULL) {
          fprintf(stderr, "Error: Failed to create output file \"%s\" [%s]",
            optarg, strerror(errno));
          badexit("");
        }
        files.o_open = 1;
        break;
      case '0':
        args.ignore_score = 1;
        break;
      case 's':
        args.ignore_strand = 1;
        break;
      case 'm':
        args.ignore_motif = 1;
        break;
      case 'S':
        args.sort_sizes = 1;
        break;
      case 'r':
        args.reverse_sort = 1;
        break;
      case 'y':
        args.override_is_yamscan = 1;
        break;
      case 'b':
        args.override_is_bed = 1;
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

  if (setlocale(LC_NUMERIC, "en_US") == NULL && args.v) {
    fprintf(stderr, "Warning: setlocale(LC_NUMERIC, \"en_US\") failed.\n");
  }

  if (!has_input) {
    badexit("Error: Missing -i arg.");
  }

  if (args.override_is_bed && args.override_is_yamscan) {
    badexit("Error: Cannot use both -y and -b.");
  }

  if (args.reverse_sort && args.ignore_score) {
    badexit("Error: Cannot use both -r and -0.");
  }

  if (args.sort_sizes && args.ignore_score) {
    badexit("Error: Cannot use both -S and -0.");
  }

  if (use_stdout) {
    files.o = stdout;
    files.o_open = 1;
  }

  score2index = malloc(sizeof(score2index_t) * ALLOC_CHUNK_SIZE);
  if (score2index == NULL) {
    fprintf(stderr, "Error: Memory allocation failed");
    badexit("");
  }
  score2index_n_alloc += ALLOC_CHUNK_SIZE;

  hash_tab = kh_init(str_h);
  
  if (args.v) {
    if (args.ignore_strand) fprintf(stderr, "Ignoring strand column.\n");
    if (args.ignore_motif) fprintf(stderr, "Ignoring motif name column.\n");
    if (args.override_is_yamscan) fprintf(stderr, "Treating input as yamscan output.\n");
    if (args.override_is_bed) fprintf(stderr, "Treating input as BED.\n");
  }
  time_t time1 = time(NULL);
  run_minidedup();
  time_t time2 = time(NULL);
  if (args.v) {
    time_t time3 = difftime(time2, time1);
    print_time((size_t) time3, "deduplicate"); 
    print_peak_mb();
  }

  free(score2index);
  free_feat_tab();
  close_files();

  return EXIT_SUCCESS;

  if (0) {
    kseq_t *kseq = kseq_init(files.i);
    kseq_read(kseq);
    kseq_destroy(kseq);
  }

}

