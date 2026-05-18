/*
 *   yamshuf: A small super-fast higher-order DNA/RNA sequence shuffler
 *   Copyright (C) 2023  Benjamin Jean-Marie Tremblay
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

/* A complete k-mer table is constructed, meaning an exponential increase in
 * memory usage with increasing k. Since counts are stored using uint64_t,
 * this means for k=10:
 *
 * 8 [bytes] * 5 [# of letters, ACGTN] ^ 10 [k] = 74.51 MB
 *
 * In addition, this is accompanied by a k-1 edge graph:
 *
 * 8 [bytes] * 5 [# of letters, ACGTN] ^ 9 [k - 1] = 14.90 MB
 *
 * Final memory usage for a k=10 scenario: 74.51 + 14.90 = 89.41 MB
 *
 * Other approx. memory footprints:
 *
 * k<5  --> <0.01 MB
 * k=5  -->  0.03 MB
 * k=6  -->  0.14 MB
 * k=7  -->  0.72 MB
 * k=8  -->  3.58 MB
 * k=9  --> 17.88 MB
 * k=10 --> 89.41 MB
 * ...
 * k=12 --> 2.18 GB
 * ...
 * k=16 --> 1.33 TB
 *
 * A more effcient way of storing k-mer counts would be to use hash tables,
 * thus only ever storing existing k-mers in memory. However since sequence
 * shuffling tasks only ever need k=1-3 99.9% of the time, yamshuf just uses
 * dumb complete k-mer tables.
 */

/* Keeping gaps with -g
 * - what about making an array to keep track of location + length of gaps,
 *   then deleting them? afterwhich the modified sequence can be shuffled as
 *   normal, and gaps inserted again when writing. this seems simpler than
 *   embedding if checks everywhere in the code for gap chars.
 * - this would lead to an incorrect number of each k-mer afterwards though..
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
#include "khash.h"
#include "version.h"

KSEQ_INIT(gzFile, gzread)
KHASH_MAP_INIT_STR(seq_str_h, uint64_t)

#define YAMSHUF_VERSION             "1.3"
#define YAMSHUF_YEAR                 2023

/* ChangeLog
 *
 * v1.3 (20 Nov 2023)
 * - Reduce branching for Euler/Markov shuffling, slight speed ups
 *
 * v1.2 (18 Aug 2023)
 * - Use xoroshiro128++ intead of xoroshiro128+
 *
 * v1.1 (15 Aug 2023)
 * - Add -p flag to print k-mer counts
 *
 */

// Maybe test which one is faster?
/* Using a complete k-mer table: */
/* #define MAX_K_TABLE                            8 */
/* Using hashes to only represent available k-mers: */
/* #define MAX_K_HASH                            27 */

/* HARD LIMIT FOR 64-BIT: 5^27 */
#define MAX_K                                  9

#define FASTA_LINE_LEN                        60

/* These must be positive integers: */
#define DEFAULT_K                              3
#define DEFAULT_SEED                           4

/* Which characters should be seen as gaps? */
/* #define GAP_CHARS                           ".-" */
/* When writing the gaps, use which character? */
/* #define FINAL_GAP                            '-' */

#define ERASE_ARRAY(ARR, LEN) memset(ARR, 0, sizeof(ARR[0]) * (LEN))

#define LIKELY(COND) __builtin_expect(COND, 1)
#define UNLIKELY(COND) __builtin_expect(COND, 0)

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
static long peak_mem(void) {
  return 0;
}
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

/* ---- RNG ---- */

// Modified krng.h code from Heng Li to use xoroshiro128++ 1.0 instead

typedef struct {
  uint64_t s[2];
} xrng_t;

static inline uint64_t splitmix64(uint64_t x) {
  uint64_t z = (x += 0x9E3779B97F4A7C15ULL);
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
  z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
  return z ^ (z >> 31);
}

#define ROTL(X, K) (((X) << (K)) | ((X) >> (64 - (K))))

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

// rng code end

static void print_peak_mb(void) {
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

static void print_time(const uint64_t s, const char *what) {
  if (s > 7200) {
    fprintf(stderr, "Needed %'.2f hours to %s.\n", ((double) s / 60.0) / 60.0, what);
  } else if (s > 120) {
    fprintf(stderr, "Needed %'.2f minutes to %s.\n", (double) s / 60.0, what);
  } else if (s > 1) {
    fprintf(stderr, "Needed %'" PRIu64 " seconds to %s.\n", s, what);
  }
}

/* ---- Usage ---- */

static void usage(void) {
  printf(
    "yamtk v%s  Copyright (C) %s  Benjamin Jean-Marie Tremblay\n"
    "Usage:  yamtk shuf [options] -i sequences.fa[.gz]\n"
    "\n"
    " -i <str>   Input FASTA/FASTQ ('-' = stdin). Non-ACGTU chars become N\n"
    "            (except with -l or -k 1).\n"
    " -k <int>   k-mer size for shuffling (default: %d). k=1 uses Fisher-Yates;\n"
    "            max k for Euler/Markov: %d.\n"
    " -o <str>   Output file (default: stdout).\n"
    " -s <int>   RNG seed (default: time-seeded).\n"
    " -m         Markov shuffling: generates sequences with similar k-mer\n"
    "            frequencies. Best for large sequences.\n"
    " -l         Linear k-mer shuffle (fast Fisher-Yates over k-mer blocks).\n"
    " -r <int>   Repeat shuffle N times per sequence; index appended to name.\n"
    " -R         Reset RNG to seed before each sequence instead of just once.\n"
    " -x <str>   BED file of ranges to restrict shuffling to. Bases inside\n"
    "            ranges are shuffled in place; bases outside pass through\n"
    "            unchanged. Incompatible with -p.\n"
    " -n         Output RNA instead of DNA. Only applies when k > 1 and -l\n"
    "            is not used.\n"
    " -p         Print k-mer counts instead of shuffling (-i, -k, -o only).\n"
    " -v / -w / -h   Verbose / very-verbose / help.\n"
    , YAMTK_VERSION, YAMTK_YEAR, DEFAULT_K, MAX_K
  );
}

/* ---- Args ---- */

typedef struct args_t {
  int      k;
  int      seed;
  int      reset_seed;
  int      use_seed;
  int      use_markov;
  int      use_linear;
  /* int      leave_gaps; */
  int      rna_out;
  int      v;
  int      w;
  int      print_kmers;
  int      shuf_repeats;
  int      use_bed;
  uint64_t window_step;
  uint64_t window_overlap;
} args_t;

static args_t args = {
  .k              = DEFAULT_K,
  .seed           = DEFAULT_SEED,
  .reset_seed     = 0,
  .use_seed       = 0,
  .use_markov     = 0,
  .use_linear     = 0,
  /* .leave_gaps     = 0, */
  .rna_out        = 0,
  .v              = 0,
  .w              = 0,
  .print_kmers    = 0,
  .shuf_repeats   = 0,
  .use_bed        = 0,
  .window_step    = 0,
  .window_overlap = 0
};

/* ---- Files ---- */

typedef struct files_t {
  int    s_open;
  int    o_open;
  int    b_open;
  gzFile s;
  FILE  *o;
  gzFile b;
} files_t;

static files_t files = {
  .s_open = 0,
  .o_open = 0,
  .b_open = 0
};

void close_files(void) {
  if (files.s_open) gzclose(files.s);
  if (files.o_open) fclose(files.o);
  if (files.b_open) gzclose(files.b);
}

/* ---- Lookup tables ---- */

/* aA = 0, cC = 1, gG = 2, tTuU = 3 */
static const unsigned char char2index[256] = { /* 16 x 16 */
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

static const char index2dna[6] = "ACGTN";
static const char index2rna[6] = "ACGUN";

/* unsigned char gap2bool[256] = { 0 }; */
/* const char gapchars[] = GAP_CHARS; */

static uint64_t char_counts[256];

static const uint64_t pow5[16] = {
           1 ,           5 ,         25 ,         125
,        625 ,        3125 ,      15625 ,       78125
,     390625 ,     1953125 ,    9765625 ,    48828125
,  244140625 ,  1220703125 , 6103515625 , 30517578125
};

static xrng_t xrng;

/* ---- Helpers ---- */

static khash_t(seq_str_h) *bed_name_hash;
static uint64_t            *bed_next;
static void free_bed(void);

static void badexit(const char *msg) {
  fprintf(stderr, "%s\nRun yamshuf -h to see usage.\n", msg);
  if (bed_name_hash) { kh_destroy(seq_str_h, bed_name_hash); bed_name_hash = NULL; }
  free(bed_next); bed_next = NULL;
  free_bed();
  close_files();
  exit(EXIT_FAILURE);
}

static inline int str_to_int(char *str, int *res) {
  char *tmp; errno = 0;
  long int res_long = strtol(str, &tmp, 10); 
  if (res_long > INT_MAX) return 1;
  *res = (int) res_long;
  if (str == tmp || errno != 0 || *tmp != '\0') {
    return 1;
  } else {
    return 0;
  }
}

/* static inline int str_to_uint64_t(char *str, uint64_t *res) { */
/*   char *tmp; errno = 0; */
/*   *res = (uint64_t) strtoull(str, &tmp, 10); */
/*   if (str == tmp || errno != 0 || *tmp != '\0') { */
/*     return 1; */
/*   } else { */
/*     return 0; */
/*   } */
/* } */

static inline void count_bases(const unsigned char *seq, const uint64_t size) {
  for (uint64_t i = 0; i < size; i++) {
    char_counts[seq[i]]++;
  }
}

static inline void swap(unsigned char *seq, const uint64_t i, const uint64_t j) {
  const unsigned char tmp = seq[i]; seq[i] = seq[j]; seq[j] = tmp;
}

/* ---- Shuffle algorithms ---- */

static int shuffle_fisher_yates(unsigned char *seq, const uint64_t len) {
  for (uint64_t i = 0, l = len - 1; i < l; i++) {
    swap(seq, i, i + xrand_r(&xrng) % (l - i));
  }
  return 0;
}

static inline void swap_k(unsigned char *seq, const uint64_t i, const uint64_t j, const uint64_t k) {
  for (uint64_t a = 0; a < k; a++) {
    swap(seq, i + a, j + a);
  }
}

static int shuffle_linear(unsigned char *seq, const uint64_t size, const uint64_t k) {
  const uint64_t n_blocks = size / k;
  for (uint64_t r, i = 0; i + 1 < n_blocks; i++) {
    r = i + 1 + xrand_r(&xrng) % (n_blocks - i - 1);
    swap_k(seq, i * k, r * k, k);
  }
  return 0;
}

static inline uint64_t chars2kmer(const unsigned char *seq, const uint64_t k, const uint64_t offset) {
  uint64_t kmer = 0;
  for (uint64_t j = 0, i = k - 1; i < -1; j++, i--) {
    kmer += pow5[i] * char2index[seq[offset + j]];
  }
  return kmer;
}

static inline void count_kmers(const unsigned char *seq, const uint64_t size, uint64_t *kmer_tab, const uint64_t k) {
  for (uint64_t i = 0; i < size - k + 1; i++) {
    kmer_tab[chars2kmer(seq, k, i)]++;
  }
}

/* ---- BED ---- */

#define BED_FIELD_MAX_CHAR    ((uint64_t) 256)
#define BED_NAME_MAX_CHAR     ((uint64_t) 512)
#define BED_ALLOC_CHUNK_SIZE  ((uint64_t) 256)

typedef struct bed_t {
  uint64_t  *seq_indices;
  uint64_t  *starts;
  uint64_t  *ends;
  char      *strands;
  char     **seq_names;
  char     **range_names;
  uint64_t   n_regions;
  uint64_t   n_comments;
  uint64_t   n_lines;
  uint64_t   n_empty;
  uint64_t   n_seqs;
  uint64_t   n_alloc;
  uint64_t   indices_are_filled;
} bed_t;

static bed_t bed = {
  .n_regions          = 0,
  .n_comments         = 0,
  .n_lines            = 0,
  .n_empty            = 0,
  .n_seqs             = 0,
  .n_alloc            = 0,
  .indices_are_filled = 0
};

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
    if (bed.indices_are_filled) {
      free(bed.seq_indices);
    }
    bed.n_alloc = 0;
  }
}

static uint64_t count_nonempty_chars(const char *line) {
  uint64_t total_chars = 0, i = 0;
  for (;;) {
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
    i++;
  }
  return total_chars;
}

static inline int str_to_uint64_t(char *str, uint64_t *res) {
  char *tmp; errno = 0;
  *res = (uint64_t) strtoull(str, &tmp, 10);
  if (str == tmp || errno != 0 || *tmp != '\0') {
    return 1;
  } else {
    return 0;
  }
}

static uint64_t count_fields(const char *line) {
  int res = 1, i = 0;
  for (;;) {
    if (line[i] == '\0') break;
    if (line[i] == '\t') res += 1;
    i++;
  }
  return res;
}

static uint64_t count_field_size(const char *line, const uint64_t k) {
  int res = 0, i = 0, n = 0;
  for (;;) {
    if (line[i] == '\0') break;
    if (line[i] == '\t') {
      n++;
    } else if (n + 1 == (int) k) {
      res++;
    } else if (n + 1 > (int) k) {
      break;
    }
    i++;
  }
  return res;
}

static uint64_t field_start(const char *line, const uint64_t k) {
  int i = 0, n = 0;
  for (;;) {
    if (line[i] == '\0') break;
    if (line[i] == '\t') {
      n++;
    } else if (n + 1 == (int) k) {
      break;
    }
    i++;
  }
  return i;
}

static uint64_t field_end(const char *line, const uint64_t k) {
  int i = 0, n = 0;
  for (;;) {
    if (line[i] == '\0') break;
    if (line[i] == '\t') {
      n++;
    }
    if (n == (int) k) {
      i--;
      break;
    }
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
      if (no_spaces && isspace(line[i])) {
        size_i--;
      } else {
        field[field_it] = line[i];
        field_it++;
      }
    }
  }
  return size_i;
}

static void read_bed(void) {
  bed.seq_names = malloc(sizeof(*bed.seq_names) * BED_ALLOC_CHUNK_SIZE);
  if (bed.seq_names == NULL) {
    badexit("Error: Failed to allocate memory for bed sequence names.");
  }
  bed.range_names = malloc(sizeof(*bed.range_names) * BED_ALLOC_CHUNK_SIZE);
  if (bed.range_names == NULL) {
    free(bed.seq_names);
    badexit("Error: Failed to allocate memory for bed range names.");
  }
  bed.strands = malloc(sizeof(*bed.strands) * BED_ALLOC_CHUNK_SIZE);
  if (bed.strands == NULL) {
    free(bed.range_names);
    free(bed.seq_names);
    badexit("Error: Failed to allocate memory for bed strand.");
  }
  bed.starts = malloc(sizeof(*bed.starts) * BED_ALLOC_CHUNK_SIZE);
  if (bed.starts == NULL) {
    free(bed.range_names);
    free(bed.strands);
    free(bed.seq_names);
    badexit("Error: Failed to allocate memory for bed starts.");
  }
  bed.ends = malloc(sizeof(*bed.ends) * BED_ALLOC_CHUNK_SIZE);
  if (bed.ends == NULL) {
    free(bed.range_names);
    free(bed.strands);
    free(bed.seq_names);
    free(bed.starts);
    badexit("Error: Failed to allocate memory for bed ends.");
  }
  bed.n_alloc = BED_ALLOC_CHUNK_SIZE;
  int ret_val;
  uint64_t line_num = 0, empty_lines = 0, comment_lines = 0;
  uint64_t n_fields = 0, field_size = 0, tmp_value = 0;
  char tmp_field[BED_FIELD_MAX_CHAR];
  kstream_t *kbed = ks_init(files.b);
  kstring_t line = { 0, 0, 0 };
  while ((ret_val = ks_getuntil(kbed, '\n', &line, 0)) >= 0) {
    line_num += 1;
    if (bed.n_regions + 1 > bed.n_alloc) {
      char **tmp_ptr1 = realloc(bed.seq_names,
        sizeof(*bed.seq_names) * bed.n_alloc + sizeof(*bed.seq_names) * BED_ALLOC_CHUNK_SIZE);
      if (tmp_ptr1 == NULL) {
        ks_destroy(kbed);
        badexit("Error: Failed to allocate more memory for bed sequence names.");
      } else {
        bed.seq_names = tmp_ptr1;
      }
      uint64_t *tmp_ptr2 = realloc(bed.starts,
        sizeof(*bed.starts) * bed.n_alloc + sizeof(*bed.starts) * BED_ALLOC_CHUNK_SIZE);
      if (tmp_ptr2 == NULL) {
        ks_destroy(kbed);
        badexit("Error: Failed to allocate more memory for bed starts.");
      } else {
        bed.starts = tmp_ptr2;
      }
      uint64_t *tmp_ptr3 = realloc(bed.ends,
        sizeof(*bed.ends) * bed.n_alloc + sizeof(*bed.ends) * BED_ALLOC_CHUNK_SIZE);
      if (tmp_ptr3 == NULL) {
        ks_destroy(kbed);
        badexit("Error: Failed to allocate more memory for bed ends.");
      } else {
        bed.ends = tmp_ptr3;
      }
      char **tmp_ptr4 = realloc(bed.range_names,
        sizeof(*bed.range_names) * bed.n_alloc + sizeof(*bed.range_names) * BED_ALLOC_CHUNK_SIZE);
      if (tmp_ptr4 == NULL) {
        ks_destroy(kbed);
        badexit("Error: Failed to allocate more memory for bed range names.");
      } else {
        bed.range_names = tmp_ptr4;
      }
      char *tmp_ptr5 = realloc(bed.strands,
        sizeof(*bed.strands) * bed.n_alloc + sizeof(*bed.strands) * BED_ALLOC_CHUNK_SIZE);
      if (tmp_ptr5 == NULL) {
        ks_destroy(kbed);
        badexit("Error: Failed to allocate more memory for bed strands.");
      } else {
        bed.strands = tmp_ptr5;
      }
      bed.n_alloc += BED_ALLOC_CHUNK_SIZE;
    }
    n_fields = count_fields(line.s);
    if (count_nonempty_chars(line.s) == 0) {
      empty_lines += 1;
      continue;
    } else if (line.s[0] == '#') {
      comment_lines += 1;
      continue;
    } else if (line.l >= 7 && line.s[0] == 'b' && line.s[1] == 'r' && line.s[2] == 'o' &&
        line.s[3] == 'w' && line.s[4] == 's' && line.s[5] == 'e' && line.s[6] == 'r') {
      comment_lines += 1;
      continue;
    } else if (line.l >= 5 && line.s[0] == 't' && line.s[1] == 'r' && line.s[2] == 'a' &&
        line.s[3] == 'c' && line.s[4] == 'k') {
      comment_lines += 1;
      continue;
    } else if (n_fields < 3) {
      ks_destroy(kbed);
      fprintf(stderr, "Line %'" PRIu64 " has %'" PRIu64 " fields and %'" PRIu64 " non-whitespace characters.\n",
          line_num, count_fields(line.s), count_nonempty_chars(line.s));
      badexit("Error: Encountered line in bed with fewer than 3 tab-separated fields.");
    }
    if (n_fields >= 6) {
      if ((field_size = parse_bed_field(line.s, 6, tmp_field, 1)) != 1) {
        ks_destroy(kbed);
        fprintf(stderr, "Error: Line %'" PRIu64 " in bed does not have a single character in the strand field (found %" PRIu64 ": '%s').",
          line_num, field_size, tmp_field);
        badexit("");
      } else if (tmp_field[0] != '+' && tmp_field[0] != '-' && tmp_field[0] != '.') {
        ks_destroy(kbed);
        fprintf(stderr, "Error: Line %'" PRIu64 " in bed has an incorrect strand character (found '%s', need +/-/.).",
          line_num, tmp_field);
        badexit("");
      }
      bed.strands[bed.n_regions] = tmp_field[0];
    } else {
      bed.strands[bed.n_regions] = '.';
    }
    if (parse_bed_field(line.s, 2, tmp_field, 1) == 0) {
      ks_destroy(kbed);
      fprintf(stderr, "Error: Line %'" PRIu64 " in bed has an empty start field.", line_num);
      badexit("");
    }
    if (str_to_uint64_t(tmp_field, &tmp_value)) {
      ks_destroy(kbed);
      fprintf(stderr, "Error: Failed to parse bed start value on line %'" PRIu64 ".\n", line_num);
      fprintf(stderr, "  Line: %s\n  Bad value: '%s'", line.s, tmp_field);
      badexit("");
    }
    bed.starts[bed.n_regions] = tmp_value;
    if (parse_bed_field(line.s, 3, tmp_field, 1) == 0) {
      ks_destroy(kbed);
      fprintf(stderr, "Error: Line %'" PRIu64 " in bed has an empty end field.", line_num);
      badexit("");
    }
    if (str_to_uint64_t(tmp_field, &tmp_value)) {
      ks_destroy(kbed);
      fprintf(stderr, "Error: Failed to parse bed end value on line %'" PRIu64 ".\n", line_num);
      fprintf(stderr, "  Line: %s\n  Bad value: '%s'", line.s, tmp_field);
      badexit("");
    }
    bed.ends[bed.n_regions] = tmp_value;
    if (bed.starts[bed.n_regions] >= bed.ends[bed.n_regions]) {
      ks_destroy(kbed);
      fprintf(stderr, "Error: Line %'" PRIu64 " in bed has a start >= end value.", line_num);
      badexit("");
    }
    if (n_fields >= 4) {
      if ((field_size = parse_bed_field(line.s, 4, tmp_field, 0)) == 0) {
        ks_destroy(kbed);
        fprintf(stderr, "Error: Line %'" PRIu64 " in bed has an empty range name.", line_num);
        badexit("");
      }
      if (field_size >= BED_FIELD_MAX_CHAR) {
        ks_destroy(kbed);
        fprintf(stderr, "Error: Range name in bed on line  %'" PRIu64 " is too large (%" PRIu64 ">%" PRIu64 ").",
          line_num, field_size, BED_FIELD_MAX_CHAR - 1);
        badexit("");
      }
      bed.range_names[bed.n_regions] = malloc(sizeof(char) * (field_size + 1));
      if (bed.range_names[bed.n_regions] == NULL) {
        ks_destroy(kbed);
        badexit("Error: Failed to allocate memory for bed range name.");
      }
      for (uint64_t i = 0; i < field_size; i++) {
        bed.range_names[bed.n_regions][i] = tmp_field[i];
      }
      bed.range_names[bed.n_regions][field_size] = '\0';
      for (uint64_t i = 0; i < BED_NAME_MAX_CHAR; i++) {
        if (bed.range_names[bed.n_regions][i] == ' ') {
          bed.range_names[bed.n_regions][i] = '\0';
          break;
        } else if (bed.range_names[bed.n_regions][i] == '\0') {
          break;
        }
      }
    } else {
      bed.range_names[bed.n_regions] = malloc(sizeof(char) * 2);
      if (bed.range_names[bed.n_regions] == NULL) {
        ks_destroy(kbed);
        badexit("Error: Failed to allocate memory for bed range name.");
      }
      bed.range_names[bed.n_regions][0] = '.';
      bed.range_names[bed.n_regions][1] = '\0';
    }
    if ((field_size = parse_bed_field(line.s, 1, tmp_field, 0)) == 0) {
      ks_destroy(kbed);
      free(bed.range_names[bed.n_regions]);
      fprintf(stderr, "Error: Line %'" PRIu64 " in bed has an empty sequence name.", line_num);
      badexit("");
    }
    if (field_size >= BED_FIELD_MAX_CHAR) {
      ks_destroy(kbed);
      free(bed.range_names[bed.n_regions]);
      fprintf(stderr, "Error: Sequence name in bed on line  %'" PRIu64 " is too large (%" PRIu64 ">%" PRIu64 ").",
        line_num, field_size, BED_FIELD_MAX_CHAR - 1);
      badexit("");
    }
    bed.seq_names[bed.n_regions] = malloc(sizeof(char) * (field_size + 1));
    if (bed.seq_names[bed.n_regions] == NULL) {
      ks_destroy(kbed);
      free(bed.range_names[bed.n_regions]);
      badexit("Error: Failed to allocate memory for bed sequence name.");
    }
    for (uint64_t i = 0; i < field_size; i++) {
      bed.seq_names[bed.n_regions][i] = tmp_field[i];
    }
    bed.seq_names[bed.n_regions][field_size] = '\0';
    for (uint64_t i = 0; i < BED_NAME_MAX_CHAR; i++) {
      if (bed.seq_names[bed.n_regions][i] == ' ') {
        bed.seq_names[bed.n_regions][i] = '\0';
        break;
      } else if (bed.seq_names[bed.n_regions][i] == '\0') {
        break;
      }
    }
    bed.n_regions += 1;
  }
  if (ret_val == -3) {
    ks_destroy(kbed);
    badexit("Error: Failed to read file stream.");
  } else if (!bed.n_regions) {
    ks_destroy(kbed);
    badexit("Error: Failed to read any records in bed file.");
  }
  ks_destroy(kbed);
  free(line.s);
  bed.n_lines = line_num;
  bed.n_comments = comment_lines;
  bed.n_empty = empty_lines;
}

static void build_bed_name_hash(void) {
  bed_next = malloc(sizeof(uint64_t) * bed.n_regions);
  if (!bed_next) badexit("Error: Failed to allocate BED next-range array.");
  for (uint64_t i = 0; i < bed.n_regions; i++) bed_next[i] = UINT64_MAX;
  int absent;
  khint64_t k;
  for (uint64_t i = 0; i < bed.n_regions; i++) {
    k = kh_put(seq_str_h, bed_name_hash, bed.seq_names[i], &absent);
    if (absent == -1) badexit("Error: Failed to hash BED sequence names.");
    if (absent == 0) {
      bed_next[i] = kh_val(bed_name_hash, k);
      kh_val(bed_name_hash, k) = i;
    } else {
      kh_val(bed_name_hash, k) = i;
    }
  }
}

static inline uint64_t cumsum_and_pick_next_letter(const uint64_t *kmers) {
  const uint64_t k0 = kmers[0];
  const uint64_t k1 = kmers[0] + kmers[1];
  const uint64_t k2 = kmers[0] + kmers[1] + kmers[2];
  const uint64_t k3 = kmers[0] + kmers[1] + kmers[2] + kmers[3];
  const uint64_t k4 = kmers[0] + kmers[1] + kmers[2] + kmers[3] + kmers[4];
  const uint64_t r = xrand_r(&xrng) % k4; 
  return 4 - ((r < k0) + (r < k1) + (r < k2) + (r < k3));
}

static inline uint64_t pick_next_letter(const uint64_t *kmers) {
  if (UNLIKELY(!kmers[4])) { // Can't think of how to remove this branch
    return (uint64_t) xrand_r(&xrng) % 4;
  } else {
    const uint64_t r = xrand_r(&xrng) % kmers[4]; 
    return 4 - ((r < kmers[0]) + (r < kmers[1]) + (r < kmers[2]) + (r < kmers[3]));
  }
}

static int shuffle_markov(unsigned char *seq, const uint64_t size, const uint64_t k, uint64_t *kmer_tab, const int is_dna) {
  const char *index2xna = is_dna ? index2dna : index2rna;
  for (uint64_t i = 0; i < pow5[k]; i += 5) {
    kmer_tab[i + 1] += kmer_tab[i];
    kmer_tab[i + 2] += kmer_tab[i + 1];
    kmer_tab[i + 3] += kmer_tab[i + 2];
    kmer_tab[i + 4] += kmer_tab[i + 3];
  }
  for (uint64_t i = 0; i < k - 1; i++) {
    seq[i] = index2xna[char2index[seq[i]]];
  }
  for (uint64_t previous = 0, i = k - 1; i < size; previous = 0, i++) {
    for (uint64_t j = k - 1; j > 0; j--) {
      previous += char2index[seq[i - j]] * pow5[j];
    }
    seq[i] = index2xna[pick_next_letter(kmer_tab + previous)];
  }
  return 0;
}

#define COUNT_EDGES(OFFSET, TABLE) (TABLE[OFFSET]+TABLE[OFFSET+1]+TABLE[OFFSET+2]+TABLE[OFFSET+3]+TABLE[OFFSET+4])

static int shuffle_euler(unsigned char *seq, const uint64_t size, const uint64_t k, uint64_t *kmer_tab, const int is_dna, unsigned char *invalid_vertex, uint64_t *euler_path, uint64_t *next_index) {

  const char *index2xna = is_dna ? index2dna : index2rna;

  // Initialize new sequence with starting vertex and final edge.

  for (uint64_t i = 0; i < k - 1; i++) {
    seq[i] = index2xna[char2index[seq[i]]];
  }
  seq[size - 1] = index2xna[char2index[seq[size - 1]]];

  // Remove the final edge from the edge pool.

  uint64_t last_edge = chars2kmer(seq, k, size - k);
  kmer_tab[last_edge]--;

  // Initialize unavailable vertices.

  for (uint64_t i = 0, j = 0; i < pow5[k - 1]; i++, j+= 5) {
    if (!(COUNT_EDGES(j, kmer_tab))) invalid_vertex[i] = 1;
  }

  // Reserve the final vertex.

  const uint64_t final_vertex = chars2kmer(seq, k - 1, size - k + 1);
  invalid_vertex[final_vertex] = 1;

  // For k > 2, prefill an array containing the indices to the next vertex from an edge.
  // (For k = 2, it is just 0 as the edge is already the correct index.)

  if (k > 2) {
    for (uint64_t i = 0, j = 0, j_max = pow5[k - 2]; i < pow5[k - 1]; i++, j++) {
      if (j == j_max) j = 0;
      next_index[i] = j * 5;
    }
  }

  // Find a random Eulerian path through all available vertices.

  for (uint64_t u, i = 0; i < pow5[k - 1]; i++) {
    u = i;
    while (!invalid_vertex[u]) {
      euler_path[u] = cumsum_and_pick_next_letter(kmer_tab + u * 5);
      u = euler_path[u] + next_index[u];
    }
    u = i;
    while (!invalid_vertex[u]) {
      invalid_vertex[u] = 1;
      u = euler_path[u] + next_index[u];
    }
  }

  // Remove reserved edges from pool.

  for (uint64_t i = 0, edge; i < pow5[k - 1]; i++) {
    if (i == final_vertex) continue;
    edge = i * 5 + euler_path[i];
    if (edge != last_edge && kmer_tab[edge]) kmer_tab[edge]--;
  }

  // Walk through Eulerian path, using up all available edges.

  for (uint64_t current_vertex, next_edge, kmer_index, i = k - 2; i < size - 2; i++) {
    current_vertex = chars2kmer(seq, k - 1, (i + 2) - k);
    kmer_index = current_vertex * 5;
    if (LIKELY(COUNT_EDGES(kmer_index, kmer_tab))) {
      // Use up all available edges per vertex.
      next_edge = cumsum_and_pick_next_letter(kmer_tab + kmer_index);
      kmer_tab[next_edge + kmer_index]--;
    } else {
      // Once exhausted, walk the Eulerian path to the next vertex.
      next_edge = euler_path[current_vertex];
    }
    seq[i + 1] = index2xna[next_edge];
  }

  return 0;

}

/* ---- Output ---- */

static void write_seq(const unsigned char *seq, const uint64_t size, const char *name, const char *comment, const uint64_t comment_l, const uint64_t n) {
  if (comment_l && n) {
    fprintf(files.o, ">%s %s-%" PRIu64 "\n", name, comment, n);
  } else if (comment_l) {
    fprintf(files.o, ">%s %s\n", name, comment);
  } else if (n) {
    fprintf(files.o, ">%s-%" PRIu64 "\n", name, n);
  } else {
    fprintf(files.o, ">%s\n", name);
  }
  for (uint64_t i = 0; i < size; i += FASTA_LINE_LEN) {
    fprintf(files.o, "%.*s\n", FASTA_LINE_LEN, seq + i);
  }
}

static void print_kmer_table_header(const uint64_t k, const int is_dna) {
  const char *index2xna = is_dna ? index2dna : index2rna;
  fputs("seq", files.o);
  if (k == 1) {
    fprintf(files.o, "\tA\tC\tG\t%c\tN", index2xna[3]);
  } else {
    // This is way more complicated than it needs to be, but I wish to
    // keep this code around in case I ever encounter a situation where
    // it would be needed.
    char kmer[MAX_K + 1];
    uint64_t decomposed[MAX_K];
    for (uint64_t j = 0; j < pow5[k]; j++) {
      fputc('\t', files.o);
      decomposed[k - 1] = j;
      for (uint64_t i = k - 2; i < -1; i--) {
        decomposed[i] = decomposed[i + 1] / 5;
      }
      kmer[0] = index2xna[decomposed[0]];
      kmer[k] = '\0';
      for (uint64_t i = 1; i < k; i++) {
        kmer[i] = index2xna[5 + (decomposed[i] - (decomposed[i - 1] + 1) * 5)];
      }
      fputs(kmer, files.o);
    }
  }
  fputc('\n', files.o);
}

/* ---- Main ---- */

int main_shuf(int argc, char **argv) {

  struct timespec ts_program;
  clock_gettime(CLOCK_MONOTONIC, &ts_program);

  kseq_t *kseq;
  int opt;
  int use_stdout = 1;

  while ((opt = getopt(argc, argv, "i:k:o:s:mlr:Rnpvwhx:")) != -1) {
    switch (opt) {
      case 'i':
        if (files.s_open) badexit("Error: -i specified more than once.");
        if (optarg[0] == '-' && optarg[1] == '\0') {
          files.s = gzdopen(fileno(stdin), "r");
        } else {
          files.s = gzopen(optarg, "r");
          if (files.s == NULL) {
            fprintf(stderr, "Error: Failed to open sequence file \"%s\" [%s]",
              optarg, strerror(errno));
            badexit("");
          }
        }
        files.s_open = 1;
        break;
      case 'o':
        if (files.o_open) badexit("Error: -o specified more than once.");
        use_stdout = 0;
        files.o = fopen(optarg, "w");
        if (files.o == NULL) {
          fprintf(stderr, "Error: Failed to create output file \"%s\" [%s]",
            optarg, strerror(errno));
          badexit("");
        }
        files.o_open = 1;
        break;
      case 'k':
        if (str_to_int(optarg, &args.k)) {
          badexit("Error: Failed to parse -k value.");
        }
        if (args.k <= 0) {
          badexit("Error: -k must be a positive integer.");
        }
        break;
      case 'm':
        args.use_markov = 1;
        break;
      case 'l':
        args.use_linear = 1;
        break;
      case 's':
        if (str_to_int(optarg, &args.seed)) {
          badexit("Error: Failed to parse -s value.");
        }
        if (args.seed <= 0) {
          badexit("Error: -s must be a positive integer.");
        }
        args.use_seed = 1;
        break;
      case 'r':
        if (str_to_int(optarg, &args.shuf_repeats)) {
          badexit("Error: Failed to parse -r value.");
        }
        if (args.shuf_repeats <= 0 || args.shuf_repeats == INT_MAX) {
          badexit("Error: -r must be a positive integer less than INT_MAX.");
        }
        break;
      /* case 'W': */
      /*   if (str_to_uint64_t(optarg, &args.window_step)) { */
      /*     badexit("Error: Failed to parse -w value."); */
      /*   } */
      /*   if (!args.window_step) { */
      /*     badexit("Error: -w must be a positive integer."); */
      /*   } */
      /*   break; */
      /* case 'S': */
      /*   if (str_to_uint64_t(optarg, &args.window_overlap)) { */
      /*     badexit("Error: Failed to parse -S value."); */
      /*   } */
      /*   if (!args.window_overlap) { */
      /*     badexit("Error: -S must be a positive integer."); */
      /*   } */
      /*   break; */
      case 'R':
        args.reset_seed = 1;
        break;
      /* case 'g': */
      /*   args.leave_gaps = 1; */
      /*   break; */
      case 'n':
        args.rna_out = 1;
        break;
      case 'p':
        args.print_kmers = 1;
        break;
      case 'x':
        if (args.use_bed) badexit("Error: -x specified more than once.");
        files.b = gzopen(optarg, "r");
        if (files.b == NULL) {
          fprintf(stderr, "Error: Failed to open BED file \"%s\" [%s]", optarg, strerror(errno));
          badexit("");
        }
        files.b_open = 1;
        args.use_bed = 1;
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

  if (!files.s_open) {
    badexit("Error: -i must be specified.");
  }
  if (args.use_bed && args.print_kmers) {
    badexit("Error: -p and -x cannot be combined.");
  }
  if (args.use_linear && args.use_markov && !args.print_kmers) {
    badexit("Error: Cannot use both -m and -l.");
  }
  if (!args.use_linear && args.k > MAX_K && !args.print_kmers) {
    fprintf(stderr, "Error: -k%d exceeds allowed max for Euler/Markov [MAX_K=%d]", args.k, MAX_K);
    badexit("");
  } else if (args.k > MAX_K && args.print_kmers) {
    fprintf(stderr, "Error: -k%d exceeds allowed max for k-mer printing [MAX_K=%d]", args.k, MAX_K);
    badexit("");
  }

  if (args.print_kmers) {
    args.v = 0;
    args.w = 0;
    args.use_markov = 0;
    args.use_linear = 0;
  }

  if (setlocale(LC_NUMERIC, "en_US") == NULL && args.v) {
    fprintf(stderr, "Warning: setlocale(LC_NUMERIC, \"en_US\") failed.\n");
  }

  if (args.v && args.rna_out && (args.k == 1 || args.use_linear) && !args.print_kmers) {
    fprintf(stderr, "Warning: The -n flag is ignored when -k is 1 or -l is used.\n");
  }

  if (use_stdout) {
    files.o = stdout;
    files.o_open = 1;
  }

  if (args.window_step && !args.window_overlap && !args.print_kmers) {
    args.window_overlap = args.window_step;
  }

  /*
  if (args.leave_gaps && !args.print_kmers) {
    const unsigned char s_len = (unsigned char) strlen(gapchars);
    for (unsigned char i = 0; i < s_len; i++) {
      gap2bool[(unsigned char) gapchars[i]] = 1;
    }
  }
  */

  time_t time1 = time(NULL);
  unsigned char *seq;
  /* uint64_t gaps; */
  uint64_t size, unknowns, seq_comment_l;
  uint64_t n_seqs = 0;
  uint64_t *kmer_tab, *euler_path, *next_index;
  unsigned char *invalid_vertex;
  char *seq_name, *seq_comment;
  double gc_pct;
  const int total_reps = args.shuf_repeats + 1, is_dna = 1 - args.rna_out;
  const uint64_t k = (uint64_t) args.k;
  int ret_val, rep_counter;
  unsigned char *scratch = NULL;
  uint64_t scratch_size = 0;
  uint8_t *bed_range_seen = NULL;
  kseq = kseq_init(files.s);
  const uint64_t actual_seed = args.use_seed
    ? (uint64_t) args.seed : (uint64_t) time(NULL);
  sxrand_r(&xrng, actual_seed);
  if (args.v)
    fprintf(stderr, "RNG seed: %" PRIu64 "\n", actual_seed);

  if (args.use_bed) {
    read_bed();
    bed_name_hash = kh_init(seq_str_h);
    if (!bed_name_hash) badexit("Error: Failed to initialize BED name hash.");
    build_bed_name_hash();
    bed_range_seen = calloc(bed.n_regions, sizeof(uint8_t));
    if (!bed_range_seen) badexit("Error: Failed to allocate BED range seen array.");
    if (args.v) {
      fprintf(stderr, "Found %'" PRIu64 " range(s) in BED file.\n", bed.n_regions);
    }
  }

  if (k > 1 && !args.use_linear) {
    kmer_tab = malloc(sizeof(uint64_t) * pow5[k]);
    if (kmer_tab == NULL) {
      badexit("Failed to allocate memory for k-mer table.");
    }
  }

  if (k > 1 && !args.use_linear && !args.use_markov && !args.print_kmers) {
    invalid_vertex = malloc(sizeof(unsigned char) * pow5[k - 1]);
    if (invalid_vertex == NULL) {
      badexit("Failed to allocate memory for vertex table.");
    }
    euler_path = malloc(sizeof(uint64_t) * pow5[k - 1]);
    if (euler_path == NULL) {
      badexit("Failed to allocate memory for Eulerian path vector.");
    }
    next_index = malloc(sizeof(uint64_t) * pow5[k - 1]);
    if (next_index == NULL) {
      badexit("Failed to allocate memory for Eulerian index vector.");
    }
  }

  if (args.print_kmers) {
    print_kmer_table_header(k, is_dna);
  }

  int markov_warning_has_been_emitted = 0;

  while ((ret_val = kseq_read(kseq)) >= 0) {

    n_seqs++;
    rep_counter = total_reps;
    seq = (unsigned char *) kseq->seq.s;
    size = kseq->seq.l;
    seq_name = kseq->name.s;
    seq_comment = kseq->comment.s;
    seq_comment_l = kseq->comment.l;

    if (args.v) {
      if (seq_comment_l) {
        fprintf(stderr, "Shuffling sequence #%'" PRIu64 ": %s %s\n", n_seqs, seq_name, seq_comment);
      } else {
        fprintf(stderr, "Shuffling sequence #%'" PRIu64 ": %s\n", n_seqs, seq_name);
      }
    }
    if (args.w || args.print_kmers) {
      ERASE_ARRAY(char_counts, 256);
      count_bases(seq, size);
      /* gaps = 0; */
      /* for (size_t i = 0; i < strlen(gapchars); i++) { */
      /*   gaps += char_counts[(unsigned char) gapchars[i]]; */
      /* } */
      unknowns = size - //gaps -
        char_counts['A'] - char_counts['a'] -
        char_counts['C'] - char_counts['c'] -
        char_counts['G'] - char_counts['g'] -
        char_counts['T'] - char_counts['t'] - char_counts['U'] - char_counts['u'];
      gc_pct = (double) (char_counts['G'] + char_counts['C'] + char_counts['g'] + char_counts['c']);
      gc_pct /= (double) (size - unknowns);// - gaps);
      gc_pct *= 100.0;
    }
    if (args.w) {
      fprintf(stderr, "  Sequence size: %'" PRIu64 " (%.2f%% non-standard)\n", size,
        100.0 * (double) unknowns / (double) size);
      fprintf(stderr, "  GC content: %.2f%%\n", gc_pct); 
    }

    if (args.reset_seed) {
      sxrand_r(&xrng, actual_seed);
    }

    if (args.print_kmers) {
      fputs(seq_name, files.o);
      if (seq_comment_l) {
        for (uint64_t i = 0; i < seq_comment_l; i++) {
          // Since the output is a TSV, it makes no sense to allow tabs in the column
          if (seq_comment[i] == '\t') seq_comment[i] = ' ';
        }
        fprintf(files.o, " %s", seq_comment);
      }
      if (k == 1) {
        fprintf(files.o, "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "",
          char_counts['a'] + char_counts['A'],
          char_counts['c'] + char_counts['C'],
          char_counts['g'] + char_counts['G'],
          char_counts['t'] + char_counts['T'] + char_counts['u'] + char_counts['U'],
          unknowns);// + gaps);
      } else if (size >= k) {
        ERASE_ARRAY(kmer_tab, pow5[k]);
        count_kmers(seq, size, kmer_tab, k);
        for (uint64_t i = 0; i < pow5[k]; i++) {
          fprintf(files.o, "\t%" PRIu64 "", kmer_tab[i]);
        }
      } else {
        for (uint64_t i = 0; i < pow5[k]; i++) {
          fputs("\t0", files.o);
        }
      }
      fputc('\n', files.o);
    } else if (args.use_bed) {
      khint64_t bed_hit = kh_get(seq_str_h, bed_name_hash, seq_name);
      if (bed_hit == kh_end(bed_name_hash)) {
        while (rep_counter) {
          write_seq(seq, size, seq_name, seq_comment, seq_comment_l, total_reps - rep_counter);
          rep_counter--;
        }
      } else {
        if (size > scratch_size) {
          unsigned char *tmp = realloc(scratch, size);
          if (!tmp) badexit("Error: Failed to allocate scratch buffer.");
          scratch = tmp;
          scratch_size = size;
        }
        while (rep_counter) {
          memcpy(scratch, seq, size);
          for (uint64_t r = kh_val(bed_name_hash, bed_hit); r != UINT64_MAX; r = bed_next[r]) {
            uint64_t rstart = bed.starts[r];
            uint64_t rend   = bed.ends[r];
            if (rstart >= size) {
              if (args.v) {
                fprintf(stderr, "! Warning: Range #%'" PRIu64 " (%s:%'" PRIu64 "-%'" PRIu64 ") starts beyond end of %s (size=%'" PRIu64 "), skipping.\n",
                  r + 1, seq_name, rstart + 1, rend, seq_name, size);
              }
              bed_range_seen[r] = 1;
              continue;
            }
            if (rend > size) {
              if (args.v) {
                fprintf(stderr, "! Warning: Trimming range #%'" PRIu64 " on %s to %" PRIu64 ".\n",
                  r + 1, seq_name, size);
              }
              rend = size;
            }
            uint64_t rlen = rend - rstart;
            bed_range_seen[r] = 1;
            if (args.k == 1) {
              if (rlen < 2) {
                if (args.v) {
                  fprintf(stderr, "! Warning: Range #%'" PRIu64 " on %s too short (len=%" PRIu64 "), skipping.\n",
                    r + 1, seq_name, rlen);
                }
                continue;
              }
              if (shuffle_fisher_yates(scratch + rstart, rlen)) goto error_shuffle;
            } else if (args.use_linear) {
              if (rlen < k * 2) {
                if (args.v) {
                  fprintf(stderr, "! Warning: Range #%'" PRIu64 " on %s too short for -l (len=%" PRIu64 " < 2k=%" PRIu64 "), skipping.\n",
                    r + 1, seq_name, rlen, k * 2);
                }
                continue;
              }
              if (shuffle_linear(scratch + rstart, rlen, k)) goto error_shuffle;
            } else if (args.use_markov) {
              if (rlen < k * 2) {
                if (args.v) {
                  fprintf(stderr, "! Warning: Range #%'" PRIu64 " on %s too short for Markov (len=%" PRIu64 " < 2k=%" PRIu64 "), skipping.\n",
                    r + 1, seq_name, rlen, k * 2);
                }
                continue;
              }
              ERASE_ARRAY(kmer_tab, pow5[k]);
              count_kmers(scratch + rstart, rlen, kmer_tab, k);
              if (shuffle_markov(scratch + rstart, rlen, k, kmer_tab, is_dna)) goto error_shuffle;
            } else {
              if (rlen < k * 2) {
                if (args.v) {
                  fprintf(stderr, "! Warning: Range #%'" PRIu64 " on %s too short for Euler (len=%" PRIu64 " < 2k=%" PRIu64 "), skipping.\n",
                    r + 1, seq_name, rlen, k * 2);
                }
                continue;
              }
              ERASE_ARRAY(invalid_vertex, pow5[k - 1]);
              ERASE_ARRAY(euler_path, pow5[k - 1]);
              ERASE_ARRAY(next_index, pow5[k - 1]);
              ERASE_ARRAY(kmer_tab, pow5[k]);
              count_kmers(scratch + rstart, rlen, kmer_tab, k);
              if (shuffle_euler(scratch + rstart, rlen, k, kmer_tab, is_dna, invalid_vertex, euler_path, next_index)) goto error_shuffle;
            }
          }
          write_seq(scratch, size, seq_name, seq_comment, seq_comment_l, total_reps - rep_counter);
          rep_counter--;
        }
      }
    } else if (size < k * 2) {
      if (args.v) {
        fprintf(stderr, "! Warning: Sequence too short to shuffle (size = %'" PRIu64 ", k = %" PRIu64 ")\n",
          size, k);
      }
    } else if (args.k == 1) {
      while (rep_counter) {
        if (shuffle_fisher_yates(seq, size)) goto error_shuffle;
        write_seq(seq, size, seq_name, seq_comment, seq_comment_l, total_reps - rep_counter);
        rep_counter--;
      }
    } else if (!args.use_linear && !args.use_markov) {
      while (rep_counter) {
        ERASE_ARRAY(invalid_vertex, pow5[k - 1]);
        ERASE_ARRAY(euler_path, pow5[k - 1]);
        ERASE_ARRAY(next_index, pow5[k - 1]);
        ERASE_ARRAY(kmer_tab, pow5[k]);
        count_kmers(seq, size, kmer_tab, k);
        if (shuffle_euler(seq, size, k, kmer_tab, is_dna, invalid_vertex, euler_path, next_index)) goto error_shuffle;
        write_seq(seq, size, seq_name, seq_comment, seq_comment_l, total_reps - rep_counter);
        rep_counter--;
      }
    } else if (args.use_linear) {
      while (rep_counter) {
        if (shuffle_linear(seq, size, k)) goto error_shuffle;
        write_seq(seq, size, seq_name, seq_comment, seq_comment_l, total_reps - rep_counter);
        rep_counter--;
      }
    } else if (args.use_markov) {
      if (size < 100 && args.v && !markov_warning_has_been_emitted) {
        fprintf(stderr, "! Warning: Markov shuffling of small sequences may generate homopolymer repeats\n");
        markov_warning_has_been_emitted = 1;
      }
      while (rep_counter) {
        ERASE_ARRAY(kmer_tab, pow5[k]);
        count_kmers(seq, size, kmer_tab, k);
        if (shuffle_markov(seq, size, k, kmer_tab, is_dna)) goto error_shuffle;
        write_seq(seq, size, seq_name, seq_comment, seq_comment_l, total_reps - rep_counter);
        rep_counter--;
      }
    }

    // TODO: send to gap-aware handler or not
    // - for markov: just skip over gaps
    // - for linear: don't swap chars when a gap is present
    // - for euler: shuffle between gaps (only way to preserve kmers)

    // TODO: send to window-aware handler or not

  }

  if (k > 1 && !args.use_linear) free(kmer_tab);
  if (k > 1 && !args.use_linear && !args.use_markov && !args.print_kmers) {
    free(invalid_vertex);
    free(euler_path);
    free(next_index);
  }
  kseq_destroy(kseq);
  if (ret_val == -2) {
    badexit("Error: Failed to parse FASTQ qualities.");
  } else if (ret_val < -2) {
    badexit("Error: Failed to read input.");
  } else if (!n_seqs) {
    badexit("Error: Failed to read any sequences from input.");
  }

  if (args.use_bed) {
    for (uint64_t i = 0; i < bed.n_regions; i++) {
      if (!bed_range_seen[i]) {
        fprintf(stderr, "Error: BED range #%'" PRIu64 " references sequence not in input (%s).\n",
          i + 1, bed.seq_names[i]);
        free(bed_range_seen);
        free(scratch);
        badexit("");
      }
    }
    free(bed_range_seen);
    free(scratch);
    if (bed_name_hash) { kh_destroy(seq_str_h, bed_name_hash); bed_name_hash = NULL; }
    free(bed_next); bed_next = NULL;
    free_bed();
  }

  close_files();

  time_t time2 = time(NULL);
  if (args.v) {
    struct timespec ts_end;
    clock_gettime(CLOCK_MONOTONIC, &ts_end);
    double elapsed = (double)(ts_end.tv_sec - ts_program.tv_sec)
                   + (double)(ts_end.tv_nsec - ts_program.tv_nsec) / 1e9;
    fprintf(stderr, "Done.\n");
    time_t time3 = difftime(time2, time1);
    print_time((uint64_t) time3, "shuffle");
    fprintf(stderr, "Total runtime: %.3fs\n", elapsed);
    print_peak_mb();
  }

  return EXIT_SUCCESS;

error_shuffle:
  if (k > 1 && !args.use_linear) free(kmer_tab);
  if (k > 1 && !args.use_linear && !args.use_markov) {
    free(invalid_vertex);
    free(euler_path);
    free(next_index);
  }
  free(bed_range_seen);
  free(scratch);
  kseq_destroy(kseq);
  badexit("");
  return EXIT_FAILURE;

}

