/*
 *   yamseq: Sequence manipulation (stats, rc, rna/dna, dup, subset, mask, subsample)
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
#include <limits.h>
#include <stdint.h>
#include <inttypes.h>
#include <ctype.h>
#include <math.h>
#include <time.h>
#include <locale.h>
#include <zlib.h>
#include "kseq.h"
#include "khash.h"
#include "version.h"

KSEQ_INIT(gzFile, gzread)
KHASH_MAP_INIT_STR(seq_str_h, uint64_t)

#define FASTA_LINE_LEN          60
#define BED_FIELD_MAX_CHAR      ((uint64_t) 256)
#define BED_NAME_MAX_CHAR       ((uint64_t) 512)
#define BED_ALLOC_CHUNK_SIZE    ((uint64_t) 256)
#define DEFAULT_GC_HIST_STEP    0.05

#define ERASE_ARRAY(ARR, LEN) memset(ARR, 0, sizeof(ARR[0]) * (LEN))

/* ---- Character classification ----
   0=A 1=C 2=G 3=T/U 4=N or other. Recognises both cases. */
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

/* ---- Action enum ---- */

typedef enum {
  ACTION_NONE = 0,
  ACTION_STATS,
  ACTION_RC,
  ACTION_RNA,
  ACTION_DNA,
  ACTION_DUP,
  ACTION_SUBSET,
  ACTION_MASK,
  ACTION_SUBSAMPLE,
  ACTION_HIST,
  ACTION_FORMAT
} action_t;

static action_t parse_action(const char *s) {
  if (!s) return ACTION_NONE;
  if (strcmp(s, "stats")     == 0) return ACTION_STATS;
  if (strcmp(s, "rc")        == 0) return ACTION_RC;
  if (strcmp(s, "rna")       == 0) return ACTION_RNA;
  if (strcmp(s, "dna")       == 0) return ACTION_DNA;
  if (strcmp(s, "dup")       == 0) return ACTION_DUP;
  if (strcmp(s, "subset")    == 0) return ACTION_SUBSET;
  if (strcmp(s, "mask")      == 0) return ACTION_MASK;
  if (strcmp(s, "subsample") == 0) return ACTION_SUBSAMPLE;
  if (strcmp(s, "hist")      == 0) return ACTION_HIST;
  if (strcmp(s, "format")    == 0) return ACTION_FORMAT;
  return ACTION_NONE;
}

static const char *action_name(action_t a) {
  switch (a) {
    case ACTION_STATS:     return "stats";
    case ACTION_RC:        return "rc";
    case ACTION_RNA:       return "rna";
    case ACTION_DNA:       return "dna";
    case ACTION_DUP:       return "dup";
    case ACTION_SUBSET:    return "subset";
    case ACTION_MASK:      return "mask";
    case ACTION_SUBSAMPLE: return "subsample";
    case ACTION_HIST:      return "hist";
    case ACTION_FORMAT:    return "format";
    default:               return "<none>";
  }
}

/* ---- Args ---- */

typedef struct args_t {
  action_t action;
  uint64_t n_arg;        /* -n N (used by dup count, subsample reservoir size) */
  double   sub_f;        /* -f p Bernoulli probability for subsample; <0 = unset */
  double   bin_step;     /* -b step for hist (fraction in (0, 1]) */
  int      use_bin_step;
  uint64_t line_len;     /* -l output FASTA line width for format (0 = unwrapped) */
  int      use_line_len;
  uint64_t seed;
  int      use_seed;
  int      mask_hard;
  int      use_bed;
  int      trim_names;
  int      progress;
  int      v;
  int      w;
} args_t;

static args_t args = {
  .action       = ACTION_NONE,
  .n_arg        = 0,
  .sub_f        = -1.0,
  .bin_step     = DEFAULT_GC_HIST_STEP,
  .use_bin_step = 0,
  .line_len     = FASTA_LINE_LEN,
  .use_line_len = 0,
  .seed         = 0,
  .use_seed     = 0,
  .mask_hard    = 0,
  .use_bed      = 0,
  .trim_names   = 1,
  .progress     = 0,
  .v            = 0,
  .w            = 0
};

/* ---- Files ---- */

typedef struct files_t {
  int     i_open, o_open, x_open;
  gzFile  i;
  FILE   *o;
  gzFile  x;
} files_t;

static files_t files = { .i_open = 0, .o_open = 0, .x_open = 0 };

static void close_files(void) {
  if (files.i_open) gzclose(files.i);
  if (files.o_open) fclose(files.o);
  if (files.x_open) gzclose(files.x);
}

static void free_bed(void);

static void badexit(const char *msg) {
  if (msg && msg[0]) fprintf(stderr, "%s\n", msg);
  fprintf(stderr, "Run yamtk seq -h to see usage.\n");
  free_bed();
  close_files();
  exit(EXIT_FAILURE);
}

/* ---- String helpers ---- */

static inline int str_to_uint64_t(char *str, uint64_t *res) {
  char *tmp; errno = 0;
  *res = (uint64_t) strtoull(str, &tmp, 10);
  if (str == tmp || errno != 0 || *tmp != '\0') return 1;
  return 0;
}

/* ---- Usage ---- */

static void usage(void) {
  printf(
    "yamtk v%s  Copyright (C) %s  Benjamin Jean-Marie Tremblay\n"
    "Usage:  yamtk seq -a <action> [options] -i seqs.fa[.gz]\n"
    "\n"
    "Actions (-a <name>):\n"
    "  stats    Per-sequence TSV: name, size, gc_pct, n_count.\n"
    "  rc       Reverse-complement each sequence (IUPAC-aware, case preserved).\n"
    "  rna      Convert T->U (case preserved).\n"
    "  dna      Convert U->T (case preserved).\n"
    "  dup      Emit each input sequence N times (requires -n N).\n"
    "  subset   Extract substrings defined by BED ranges (requires -x).\n"
    "           BED col-4 = output name (else seq:start-end), col-6 = strand.\n"
    "  mask     Soft-mask (lowercase) BED regions; -N for hard mask (N).\n"
    "           Requires -x; strand ignored.\n"
    "  subsample  Random subset of input sequences. Requires either -n N\n"
    "           (reservoir, exact count) or -f p (Bernoulli, 0<p<1).\n"
    "           Input order is preserved in the output.\n"
    "  hist     Per-bin TSV histogram of per-sequence GC fractions, with\n"
    "           seq_count and bp_count columns. Bin step set by -b.\n"
    "  format   Re-emit input FASTA: rewrap lines, trim names to the first\n"
    "           word. Width set by -l (0 = unwrapped). -v/-w report the\n"
    "           character composition.\n"
    "\n"
    " -i <str>   Input FASTA/FASTQ ('-' = stdin).\n"
    " -o <str>   Output (default: stdout).\n"
    " -x <str>   BED file (subset/mask only). Can be gzipped.\n"
    " -n <int>   Count: copies per input (dup) or reservoir size (subsample).\n"
    " -f <dbl>   Per-sequence keep probability for subsample (0,1).\n"
    " -b <dbl>   GC bin step for hist, in (0, 1] (default: %.2f).\n"
    " -l <int>   Output FASTA line width for format (0 = unwrapped; default: %d).\n"
    " -s <uint>  RNG seed for subsample (default: time-seeded).\n"
    " -N         Hard-mask (replace with N) instead of soft-mask. mask only.\n"
    " -r         Do not trim sequence names to the first word.\n"
    " -g         Show progress bar.\n"
    " -v / -w / -h   Verbose / very-verbose / help.\n"
    , YAMTK_VERSION, YAMTK_YEAR, DEFAULT_GC_HIST_STEP, FASTA_LINE_LEN
  );
}

/* ---- IUPAC complement ---- */

static inline unsigned char iupac_comp(unsigned char c) {
  switch (c) {
    case 'A':           return 'T'; case 'a':           return 't';
    case 'T': case 'U': return 'A'; case 't': case 'u': return 'a';
    case 'C':           return 'G'; case 'c':           return 'g';
    case 'G':           return 'C'; case 'g':           return 'c';
    case 'R':           return 'Y'; case 'r':           return 'y';
    case 'Y':           return 'R'; case 'y':           return 'r';
    case 'K':           return 'M'; case 'k':           return 'm';
    case 'M':           return 'K'; case 'm':           return 'k';
    case 'B':           return 'V'; case 'b':           return 'v';
    case 'V':           return 'B'; case 'v':           return 'b';
    case 'D':           return 'H'; case 'd':           return 'h';
    case 'H':           return 'D'; case 'h':           return 'd';
    /* W, S, N are self-complement; anything else passes through (gaps, etc). */
    default: return c;
  }
}

static void reverse_complement(unsigned char *seq, const uint64_t L) {
  if (L == 0) return;
  for (uint64_t i = 0, j = L - 1; i < j; i++, j--) {
    const unsigned char a = iupac_comp(seq[i]);
    const unsigned char b = iupac_comp(seq[j]);
    seq[i] = b; seq[j] = a;
  }
  if (L % 2 == 1) seq[L / 2] = iupac_comp(seq[L / 2]);
}

/* ---- T<->U toggle ---- */

static void toggle_t_to_u(unsigned char *seq, const uint64_t L) {
  for (uint64_t i = 0; i < L; i++) {
    if      (seq[i] == 'T') seq[i] = 'U';
    else if (seq[i] == 't') seq[i] = 'u';
  }
}

static void toggle_u_to_t(unsigned char *seq, const uint64_t L) {
  for (uint64_t i = 0; i < L; i++) {
    if      (seq[i] == 'U') seq[i] = 'T';
    else if (seq[i] == 'u') seq[i] = 't';
  }
}

/* ---- PRNG (xoroshiro128++, seeded via splitmix64). Shared with subsample. ---- */

#define ROTL(x, k) (((x) << (k)) | ((x) >> (64 - (k))))

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

/* Uniform integer in [0, k). Unbiased rejection sampling. */
static inline uint64_t xrand_below(xrng_t *r, uint64_t k) {
  if (k == 0) return 0;
  const uint64_t lim = UINT64_MAX - (UINT64_MAX % k);
  uint64_t x;
  do { x = xrand_r(r); } while (x >= lim);
  return x % k;
}

static xrng_t xrng;

/* ---- Peak memory / elapsed-time reporting (shared with other subcommands) ---- */

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

/* ---- FASTA output ---- */

static void write_seq(const unsigned char *seq, const uint64_t L, const char *name,
                      const char *comment, const uint64_t comment_l) {
  if (!args.trim_names && comment_l) {
    fprintf(files.o, ">%s %s\n", name, comment);
  } else {
    fprintf(files.o, ">%s\n", name);
  }
  if (L == 0) return;
  if (args.line_len == 0) {
    /* Unwrapped: whole sequence on one line. */
    fwrite(seq, 1, L, files.o);
    fputc('\n', files.o);
    return;
  }
  for (uint64_t i = 0; i < L; i += args.line_len) {
    const uint64_t chunk = (L - i < args.line_len) ? L - i : args.line_len;
    fwrite(seq + i, 1, chunk, files.o);
    fputc('\n', files.o);
  }
}

/* ---- BED ---- */

typedef struct bed_t {
  uint64_t  *starts;
  uint64_t  *ends;
  char      *strands;
  char     **seq_names;
  char     **range_names;
  uint64_t   n_regions;
  uint64_t   n_alloc;
} bed_t;

static bed_t bed = { .n_regions = 0, .n_alloc = 0 };

static khash_t(seq_str_h) *bed_seq_hash = NULL;
static uint64_t           *bed_next     = NULL;

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
    bed.n_alloc = 0; bed.n_regions = 0;
  }
  free(bed_next); bed_next = NULL;
  if (bed_seq_hash) { kh_destroy(seq_str_h, bed_seq_hash); bed_seq_hash = NULL; }
}

static uint64_t count_nonempty_chars(const char *line) {
  uint64_t n = 0;
  for (uint64_t i = 0; line[i] != '\0'; i++) {
    if (line[i] != ' ' && line[i] != '\t' && line[i] != '\r' && line[i] != '\n') n++;
  }
  return n;
}

static uint64_t count_fields(const char *line) {
  uint64_t res = 1, i = 0;
  while (line[i] != '\0') { if (line[i] == '\t') res++; i++; }
  return res;
}

static uint64_t count_field_size(const char *line, const uint64_t k) {
  uint64_t res = 0, i = 0, n = 0;
  while (line[i] != '\0') {
    if (line[i] == '\t') n++;
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
  if (!p1 || !p2 || !p3 || !p4 || !p5) badexit("Error: Failed to grow BED storage.");
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
    if (n_fields < 3) {
      ks_destroy(kbed);
      fprintf(stderr, "Error: BED line %" PRIu64 " has %" PRIu64 " fields; need >=3.\n",
        line_num, n_fields);
      badexit("");
    }
    bed_grow_if_needed();

    if (n_fields >= 6) {
      const uint64_t sz = parse_bed_field(line.s, 6, tmp_field, 1);
      if (sz != 1 || (tmp_field[0] != '+' && tmp_field[0] != '-' && tmp_field[0] != '.')) {
        ks_destroy(kbed);
        fprintf(stderr, "Error: BED line %" PRIu64 " strand must be +/-/. (got '%s').\n",
          line_num, tmp_field);
        badexit("");
      }
      bed.strands[bed.n_regions] = tmp_field[0];
    } else {
      bed.strands[bed.n_regions] = '.';
    }

    if (parse_bed_field(line.s, 2, tmp_field, 1) == 0) {
      ks_destroy(kbed);
      fprintf(stderr, "Error: BED line %" PRIu64 " empty start.\n", line_num);
      badexit("");
    }
    uint64_t tmp_value;
    if (str_to_uint64_t(tmp_field, &tmp_value)) {
      ks_destroy(kbed);
      fprintf(stderr, "Error: Failed to parse BED start on line %" PRIu64 " ('%s').\n", line_num, tmp_field);
      badexit("");
    }
    bed.starts[bed.n_regions] = tmp_value;

    if (parse_bed_field(line.s, 3, tmp_field, 1) == 0) {
      ks_destroy(kbed);
      fprintf(stderr, "Error: BED line %" PRIu64 " empty end.\n", line_num);
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

    /* Col 4 (optional): output name for subset. Default to '.' */
    if (n_fields >= 4) {
      const uint64_t sz4 = parse_bed_field(line.s, 4, tmp_field, 0);
      if (sz4 == 0 || sz4 >= BED_FIELD_MAX_CHAR) {
        bed.range_names[bed.n_regions] = malloc(2);
        if (!bed.range_names[bed.n_regions]) { ks_destroy(kbed); badexit("Error: alloc."); }
        bed.range_names[bed.n_regions][0] = '.';
        bed.range_names[bed.n_regions][1] = '\0';
      } else {
        bed.range_names[bed.n_regions] = malloc(sz4 + 1);
        if (!bed.range_names[bed.n_regions]) { ks_destroy(kbed); badexit("Error: alloc."); }
        memcpy(bed.range_names[bed.n_regions], tmp_field, sz4);
        bed.range_names[bed.n_regions][sz4] = '\0';
        for (uint64_t i = 0; i < sz4; i++) {
          if (bed.range_names[bed.n_regions][i] == ' ') {
            bed.range_names[bed.n_regions][i] = '\0';
            break;
          }
        }
      }
    } else {
      bed.range_names[bed.n_regions] = malloc(2);
      if (!bed.range_names[bed.n_regions]) { ks_destroy(kbed); badexit("Error: alloc."); }
      bed.range_names[bed.n_regions][0] = '.';
      bed.range_names[bed.n_regions][1] = '\0';
    }

    const uint64_t sz1 = parse_bed_field(line.s, 1, tmp_field, 0);
    if (sz1 == 0 || sz1 >= BED_FIELD_MAX_CHAR) {
      ks_destroy(kbed);
      free(bed.range_names[bed.n_regions]);
      fprintf(stderr, "Error: BED line %" PRIu64 " has empty or too-long seq name.\n", line_num);
      badexit("");
    }
    bed.seq_names[bed.n_regions] = malloc(sz1 + 1);
    if (!bed.seq_names[bed.n_regions]) {
      ks_destroy(kbed); free(bed.range_names[bed.n_regions]);
      badexit("Error: alloc.");
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

/* ---- Action: subset ---- */

static void do_subset(const unsigned char *seq, const uint64_t L, const char *seq_name,
                      const uint64_t region_i) {
  const uint64_t s = bed.starts[region_i];
  const uint64_t e = bed.ends[region_i];
  if (s >= L) {
    fprintf(stderr, "Warning: BED start %" PRIu64 " >= seq '%s' length %" PRIu64 "; skipping.\n",
      s, seq_name, L);
    return;
  }
  const uint64_t end = (e > L) ? L : e;
  const uint64_t n = end - s;
  unsigned char *buf = malloc(n);
  if (!buf) badexit("Error: Failed to alloc subset buffer.");
  memcpy(buf, seq + s, n);
  const char strand = bed.strands[region_i];
  if (strand == '-') reverse_complement(buf, n);

  /* Output name: BED col-4 if not '.', else fallback "seq:start-end" */
  const char *range_name = bed.range_names[region_i];
  char fallback[BED_NAME_MAX_CHAR];
  if (range_name[0] == '.' && range_name[1] == '\0') {
    snprintf(fallback, sizeof(fallback), "%s:%" PRIu64 "-%" PRIu64, seq_name, s, end);
    range_name = fallback;
  }
  write_seq(buf, n, range_name, NULL, 0);
  free(buf);
}

/* ---- Action: mask ---- */

static void do_mask(unsigned char *seq, const uint64_t L, const char *seq_name,
                    const uint64_t region_i) {
  const uint64_t s = bed.starts[region_i];
  if (s >= L) {
    fprintf(stderr, "Warning: BED start %" PRIu64 " >= seq '%s' length %" PRIu64 "; skipping.\n",
      s, seq_name, L);
    return;
  }
  uint64_t e = bed.ends[region_i];
  if (e > L) e = L;
  if (args.mask_hard) {
    for (uint64_t i = s; i < e; i++) seq[i] = 'N';
  } else {
    for (uint64_t i = s; i < e; i++) seq[i] = (unsigned char) tolower(seq[i]);
  }
}

/* ---- Action: stats ---- */

static void emit_stats_header(void) {
  fprintf(files.o, "##seq_num\tseq_name\tsize\tgc_pct\tn_count\n");
}

static void emit_stats_row(const unsigned char *seq, const uint64_t L, const char *name,
                           const uint64_t idx) {
  uint64_t n_a = 0, n_c = 0, n_g = 0, n_t = 0, n_n = 0;
  for (uint64_t i = 0; i < L; i++) {
    switch (char2index[seq[i]]) {
      case 0: n_a++; break;
      case 1: n_c++; break;
      case 2: n_g++; break;
      case 3: n_t++; break;
      default: n_n++; break;
    }
  }
  const uint64_t known = n_a + n_c + n_g + n_t;
  const double gc_pct = known ? ((double)(n_c + n_g) / (double)known) * 100.0 : 0.0;
  fprintf(files.o, "%" PRIu64 "\t%s\t%" PRIu64 "\t%.2f\t%" PRIu64 "\n",
    idx, name, L, gc_pct, n_n);
}

/* ---- Action: format ----
   Re-emit the input FASTA, rewrapping lines (write_seq honours args.line_len)
   and trimming names to the first word unless -r. With -v/-w, accumulate a
   per-character frequency table and print a composition report to stderr. */

static uint64_t comp_freq[256] = {0};   /* file-wide accumulator */

static void tally_into(uint64_t freq[256], const unsigned char *seq, const uint64_t L) {
  for (uint64_t i = 0; i < L; i++) freq[seq[i]]++;
}

static void report_composition(FILE *out, const uint64_t freq[256],
                               const uint64_t nseq, const char *label) {
  uint64_t total = 0;
  for (int c = 0; c < 256; c++) total += freq[c];

  const uint64_t n_a = freq['A'] + freq['a'];
  const uint64_t n_c = freq['C'] + freq['c'];
  const uint64_t n_g = freq['G'] + freq['g'];
  const uint64_t n_t = freq['T'] + freq['t'];
  const uint64_t n_u = freq['U'] + freq['u'];
  const uint64_t acgtu = n_a + n_c + n_g + n_t + n_u;
  const uint64_t other = total - acgtu;

  /* Soft mask = lowercase letters; hard mask = N/n/-/. */
  uint64_t soft = 0;
  for (int c = 'a'; c <= 'z'; c++) soft += freq[c];
  const uint64_t hard = freq['N'] + freq['n'] + freq['-'] + freq['.'];

  /* Most frequent character that is not A/C/G/T/U (either case). */
  int top_c = -1;
  uint64_t top_n = 0;
  for (int c = 0; c < 256; c++) {
    if (c=='A'||c=='a'||c=='C'||c=='c'||c=='G'||c=='g'||
        c=='T'||c=='t'||c=='U'||c=='u') continue;
    if (freq[c] > top_n) { top_n = freq[c]; top_c = c; }
  }

  const double denom = total ? (double) total : 1.0;

  /* Header line, then an indented block of one item per line. */
  if (label)
    fprintf(out, "[%s] %" PRIu64 " bp:\n", label, total);
  else
    fprintf(out, "Composition (%" PRIu64 " sequence(s), %" PRIu64 " bp):\n",
      nseq, total);

  fprintf(out, "  bases:        A=%" PRIu64 " C=%" PRIu64 " G=%" PRIu64
    " T=%" PRIu64 " U=%" PRIu64 "\n", n_a, n_c, n_g, n_t, n_u);
  fprintf(out, "  other:        %" PRIu64 " (%.1f%%)\n",
    other, 100.0 * (double) other / denom);
  fprintf(out, "  soft-masked:  %" PRIu64 " (%.1f%%)\n",
    soft,  100.0 * (double) soft  / denom);
  fprintf(out, "  hard-masked:  %" PRIu64 " (%.1f%%)\n",
    hard,  100.0 * (double) hard  / denom);

  /* Call out a dominant unknown character (more than 1% of all bases). */
  if (top_c >= 0 && total && (double) top_n / denom > 0.01) {
    if (top_c >= 32 && top_c < 127)
      fprintf(out, "  top unknown:  '%c' = %" PRIu64 " (%.1f%%)\n",
        (char) top_c, top_n, 100.0 * (double) top_n / denom);
    else
      fprintf(out, "  top unknown:  '\\x%02X' = %" PRIu64 " (%.1f%%)\n",
        (unsigned) top_c, top_n, 100.0 * (double) top_n / denom);
  }
}

/* ---- Action: hist ----
   Per-sequence GC fraction binned into uniform-width bins of width
   args.bin_step. All-N / empty sequences contribute to n_skipped
   (not bin 0). */

static int gc_hist_n_bins(double step) {
  int n = (int)ceil(1.0 / step);
  if (n < 1) n = 1;
  return n;
}

static int gc_hist_bin(double gc_frac, double step, int nb) {
  int b = (int)floor(gc_frac / step);
  if (b < 0) b = 0;
  if (b >= nb) b = nb - 1;
  return b;
}

/* ---- Action: subsample ----

   Two modes (mutually exclusive):
     -n N : reservoir sampling, exact output count (or all if input < N).
     -f p : independent Bernoulli per sequence, streaming, approximate count.

   In both modes the output preserves input order. For -f that's free (we
   decide on the fly); for -n we sort the reservoir by original index
   before emitting. */

typedef struct {
  uint64_t       idx;       /* 1-based original input index */
  uint64_t       L;
  uint64_t       comment_l;
  char          *name;
  unsigned char *seq;
  char          *comment;   /* NULL if no comment or if !args.trim_names==0 */
} reservoir_slot_t;

static reservoir_slot_t *reservoir = NULL;
static uint64_t          reservoir_filled = 0;

static int slot_copy(reservoir_slot_t *dst, const kseq_t *kseq, const uint64_t idx) {
  dst->idx = idx;
  dst->L = kseq->seq.l;
  dst->name = strdup(kseq->name.s ? kseq->name.s : "");
  if (!dst->name) return 1;
  if (dst->L) {
    dst->seq = malloc(dst->L + 1);
    if (!dst->seq) { free(dst->name); dst->name = NULL; return 1; }
    memcpy(dst->seq, kseq->seq.s, dst->L);
    dst->seq[dst->L] = '\0';
  } else {
    dst->seq = NULL;
  }
  dst->comment_l = kseq->comment.l;
  if (kseq->comment.l && kseq->comment.s) {
    dst->comment = strdup(kseq->comment.s);
    if (!dst->comment) { free(dst->name); free(dst->seq); dst->name=NULL; dst->seq=NULL; return 1; }
  } else {
    dst->comment = NULL;
  }
  return 0;
}

static void slot_free(reservoir_slot_t *s) {
  free(s->name); free(s->seq); free(s->comment);
  s->name = NULL; s->seq = NULL; s->comment = NULL;
}

static int cmp_slot_by_idx(const void *a, const void *b) {
  const reservoir_slot_t *sa = a, *sb = b;
  if (sa->idx < sb->idx) return -1;
  if (sa->idx > sb->idx) return  1;
  return 0;
}

/* ---- Validation of action + flag combos ---- */

static void validate_args(void) {
  if (args.action == ACTION_NONE) {
    badexit("Error: -a <action> must be specified.");
  }
  /* -x only with subset/mask, required there */
  if (args.use_bed && args.action != ACTION_SUBSET && args.action != ACTION_MASK) {
    fprintf(stderr, "Error: -x is only valid with -a subset or -a mask (got -a %s).\n",
      action_name(args.action));
    badexit("");
  }
  if ((args.action == ACTION_SUBSET || args.action == ACTION_MASK) && !args.use_bed) {
    fprintf(stderr, "Error: -a %s requires -x <bed>.\n", action_name(args.action));
    badexit("");
  }
  /* -n only with dup/subsample, required for dup */
  if (args.n_arg > 0 && args.action != ACTION_DUP && args.action != ACTION_SUBSAMPLE) {
    fprintf(stderr, "Error: -n is only valid with -a dup or -a subsample (got -a %s).\n",
      action_name(args.action));
    badexit("");
  }
  if (args.action == ACTION_DUP && args.n_arg == 0) {
    badexit("Error: -a dup requires -n <int>=1.");
  }
  /* -f only with subsample */
  if (args.sub_f >= 0.0 && args.action != ACTION_SUBSAMPLE) {
    fprintf(stderr, "Error: -f is only valid with -a subsample (got -a %s).\n",
      action_name(args.action));
    badexit("");
  }
  if (args.action == ACTION_SUBSAMPLE) {
    int has_n = (args.n_arg > 0);
    int has_f = (args.sub_f >= 0.0);
    if (!has_n && !has_f) badexit("Error: -a subsample requires either -n <int> or -f <dbl>.");
    if ( has_n &&  has_f) badexit("Error: -a subsample: -n and -f are mutually exclusive.");
    if (has_f && (args.sub_f <= 0.0 || args.sub_f >= 1.0))
      badexit("Error: -f must be in (0,1).");
  }
  /* -s only with subsample */
  if (args.use_seed && args.action != ACTION_SUBSAMPLE) {
    fprintf(stderr, "Error: -s is only valid with -a subsample (got -a %s).\n",
      action_name(args.action));
    badexit("");
  }
  /* -b only with hist */
  if (args.use_bin_step && args.action != ACTION_HIST) {
    fprintf(stderr, "Error: -b is only valid with -a hist (got -a %s).\n",
      action_name(args.action));
    badexit("");
  }
  /* -l only with format */
  if (args.use_line_len && args.action != ACTION_FORMAT) {
    fprintf(stderr, "Error: -l is only valid with -a format (got -a %s).\n",
      action_name(args.action));
    badexit("");
  }
  /* -N only with mask */
  if (args.mask_hard && args.action != ACTION_MASK) {
    fprintf(stderr, "Error: -N is only valid with -a mask (got -a %s).\n",
      action_name(args.action));
    badexit("");
  }
}

/* ---- Main ---- */

int main_seq(int argc, char **argv) {

  struct timespec ts_program;
  clock_gettime(CLOCK_MONOTONIC, &ts_program);

  int opt;
  int use_stdout = 1;

  while ((opt = getopt(argc, argv, "a:i:o:x:n:f:b:l:s:Nrgvwh")) != -1) {
    switch (opt) {
      case 'a':
        if (args.action != ACTION_NONE) badexit("Error: -a specified more than once.");
        args.action = parse_action(optarg);
        if (args.action == ACTION_NONE) {
          fprintf(stderr, "Error: Unknown action '%s'.\n", optarg);
          badexit("");
        }
        break;
      case 'i':
        if (files.i_open) badexit("Error: -i specified more than once.");
        if (optarg[0] == '-' && optarg[1] == '\0') {
          files.i = gzdopen(fileno(stdin), "r");
        } else {
          files.i = gzopen(optarg, "r");
          if (files.i == NULL) {
            fprintf(stderr, "Error: Failed to open input file \"%s\" [%s]\n",
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
      case 'n':
        if (str_to_uint64_t(optarg, &args.n_arg)) {
          badexit("Error: Failed to parse -n value.");
        }
        if (args.n_arg == 0) {
          badexit("Error: -n must be a positive integer.");
        }
        break;
      case 'f': {
        char *tmp; errno = 0;
        args.sub_f = strtod(optarg, &tmp);
        if (optarg == tmp || errno != 0 || *tmp != '\0') {
          badexit("Error: Failed to parse -f value.");
        }
        break;
      }
      case 'b': {
        char *tmp; errno = 0;
        args.bin_step = strtod(optarg, &tmp);
        if (optarg == tmp || errno != 0 || *tmp != '\0') {
          badexit("Error: Failed to parse -b value.");
        }
        if (args.bin_step <= 0.0 || args.bin_step > 1.0) {
          badexit("Error: -b must be in (0, 1].");
        }
        args.use_bin_step = 1;
        break;
      }
      case 'l':
        /* Permit 0 (unwrapped); str_to_uint64_t accepts it. */
        if (str_to_uint64_t(optarg, &args.line_len)) {
          badexit("Error: Failed to parse -l value.");
        }
        args.use_line_len = 1;
        break;
      case 's':
        if (str_to_uint64_t(optarg, &args.seed)) {
          badexit("Error: Failed to parse -s value.");
        }
        args.use_seed = 1;
        break;
      case 'N':
        args.mask_hard = 1;
        break;
      case 'r':
        args.trim_names = 0;
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

  if (!files.i_open) badexit("Error: -i must be specified.");
  validate_args();
  if (use_stdout) files.o = stdout;

  if (setlocale(LC_NUMERIC, "en_US") == NULL && args.v) {
    fprintf(stderr, "Warning: setlocale(LC_NUMERIC, \"en_US\") failed.\n");
  }

  /* BED setup for subset / mask */
  uint64_t *seq_regions_buf = NULL;
  if (args.use_bed) {
    read_bed();
    build_bed_seq_hash();
    seq_regions_buf = malloc(sizeof(*seq_regions_buf) * bed.n_regions);
    if (!seq_regions_buf) badexit("Error: Failed to alloc per-seq BED scratch.");
  }

  /* Per-action one-time header */
  if (args.action == ACTION_STATS) emit_stats_header();

  /* Subsample setup: seed RNG, allocate reservoir if -n. */
  if (args.action == ACTION_SUBSAMPLE) {
    const uint64_t actual_seed = args.use_seed
      ? args.seed : (uint64_t)time(NULL);
    sxrand_r(&xrng, actual_seed);
    if (args.v) fprintf(stderr, "RNG seed: %" PRIu64 "\n", actual_seed);
    if (args.n_arg > 0) {
      reservoir = calloc(args.n_arg, sizeof(*reservoir));
      if (!reservoir) badexit("Error: Failed to alloc reservoir.");
    }
  }

  /* Hist setup: allocate per-bin counters. */
  int       hist_nb       = 0;
  uint64_t *hist_seq_cnt  = NULL;
  uint64_t *hist_bp_cnt   = NULL;
  uint64_t  hist_skipped  = 0;
  if (args.action == ACTION_HIST) {
    hist_nb = gc_hist_n_bins(args.bin_step);
    hist_seq_cnt = calloc((size_t)hist_nb, sizeof(*hist_seq_cnt));
    hist_bp_cnt  = calloc((size_t)hist_nb, sizeof(*hist_bp_cnt));
    if (!hist_seq_cnt || !hist_bp_cnt) {
      free(hist_seq_cnt); free(hist_bp_cnt);
      badexit("Error: Failed to alloc hist bins.");
    }
  }

  /* Stream sequences */
  time_t time1 = time(NULL);
  kseq_t *kseq = kseq_init(files.i);
  int ret_val;
  uint64_t idx = 0;
  uint64_t kept = 0;  /* for subsample summary */
  while ((ret_val = kseq_read(kseq)) >= 0) {
    idx++;
    unsigned char *seq = (unsigned char *) kseq->seq.s;
    const uint64_t L = kseq->seq.l;
    const char *name = kseq->name.s;

    switch (args.action) {
      case ACTION_STATS:
        emit_stats_row(seq, L, name, idx);
        break;
      case ACTION_RC:
        reverse_complement(seq, L);
        write_seq(seq, L, name, kseq->comment.s, kseq->comment.l);
        break;
      case ACTION_RNA:
        toggle_t_to_u(seq, L);
        write_seq(seq, L, name, kseq->comment.s, kseq->comment.l);
        break;
      case ACTION_DNA:
        toggle_u_to_t(seq, L);
        write_seq(seq, L, name, kseq->comment.s, kseq->comment.l);
        break;
      case ACTION_DUP: {
        char dup_name[1024];
        for (uint64_t i = 1; i <= args.n_arg; i++) {
          snprintf(dup_name, sizeof(dup_name), "%s-%" PRIu64, name, i);
          write_seq(seq, L, dup_name, NULL, 0);
        }
        break;
      }
      case ACTION_SUBSET: {
        const uint64_t nr = collect_seq_regions(name, seq_regions_buf);
        if (nr > 0) {
          qsort(seq_regions_buf, nr, sizeof(*seq_regions_buf), cmp_idx_by_start);
          for (uint64_t r_i = 0; r_i < nr; r_i++) {
            do_subset(seq, L, name, seq_regions_buf[r_i]);
          }
        }
        break;
      }
      case ACTION_MASK: {
        const uint64_t nr = collect_seq_regions(name, seq_regions_buf);
        for (uint64_t r_i = 0; r_i < nr; r_i++) {
          do_mask(seq, L, name, seq_regions_buf[r_i]);
        }
        write_seq(seq, L, name, kseq->comment.s, kseq->comment.l);
        break;
      }
      case ACTION_HIST: {
        uint64_t n_cg = 0, n_known = 0;
        for (uint64_t i = 0; i < L; i++) {
          const unsigned char idx_b = char2index[seq[i]];
          if (idx_b == 4) continue;
          n_known++;
          if (idx_b == 1 || idx_b == 2) n_cg++;
        }
        if (n_known == 0) {
          hist_skipped++;
        } else {
          const double gc_frac = (double)n_cg / (double)n_known;
          const int b = gc_hist_bin(gc_frac, args.bin_step, hist_nb);
          hist_seq_cnt[b]++;
          hist_bp_cnt[b]  += L;
        }
        break;
      }
      case ACTION_SUBSAMPLE:
        if (args.n_arg > 0) {
          /* Reservoir sampling (Vitter algorithm R). */
          if (idx <= args.n_arg) {
            if (slot_copy(&reservoir[idx - 1], kseq, idx))
              badexit("Error: Failed to copy sequence into reservoir.");
            reservoir_filled = idx;
          } else {
            const uint64_t j = xrand_below(&xrng, idx);
            if (j < args.n_arg) {
              slot_free(&reservoir[j]);
              if (slot_copy(&reservoir[j], kseq, idx))
                badexit("Error: Failed to copy sequence into reservoir.");
            }
          }
        } else {
          /* Bernoulli: emit immediately if u < p. Streams in input order. */
          if (xrand_double(&xrng) < args.sub_f) {
            write_seq(seq, L, name, kseq->comment.s, kseq->comment.l);
            kept++;
          }
        }
        break;
      case ACTION_FORMAT:
        if (args.v) {
          if (args.w) {
            uint64_t f[256] = {0};
            tally_into(f, seq, L);
            report_composition(stderr, f, 1, name);   /* per-seq line */
            for (int c = 0; c < 256; c++) comp_freq[c] += f[c];
          } else {
            tally_into(comp_freq, seq, L);
          }
        }
        write_seq(seq, L, name, kseq->comment.s, kseq->comment.l);
        break;
      default:
        /* Other actions wired in subsequent commits */
        break;
    }
  }
  if (ret_val < -1) {
    fprintf(stderr, "Error: Failed to read input FASTA (kseq_read returned %d).\n", ret_val);
    kseq_destroy(kseq);
    badexit("");
  }
  kseq_destroy(kseq);
  free(seq_regions_buf);

  /* Emit hist TSV. */
  if (args.action == ACTION_HIST) {
    fprintf(files.o, "##yamseq hist gc_step=%g\n", args.bin_step);
    fprintf(files.o, "#bin_lo\tbin_hi\tseq_count\tbp_count\n");
    for (int b = 0; b < hist_nb; b++) {
      double lo = (double)b * args.bin_step;
      double hi = lo + args.bin_step;
      if (hi > 1.0) hi = 1.0;
      fprintf(files.o, "%g\t%g\t%" PRIu64 "\t%" PRIu64 "\n",
        lo, hi, hist_seq_cnt[b], hist_bp_cnt[b]);
    }
    if (args.v && hist_skipped) {
      fprintf(stderr, "Skipped %" PRIu64 " sequence(s) with no A/C/G/T bases.\n",
        hist_skipped);
    }
    free(hist_seq_cnt); free(hist_bp_cnt);
    hist_seq_cnt = NULL; hist_bp_cnt = NULL;
  }

  /* Flush reservoir (subsample -n N): sort by original index, write, free. */
  if (args.action == ACTION_SUBSAMPLE && args.n_arg > 0) {
    if (reservoir_filled > 1)
      qsort(reservoir, reservoir_filled, sizeof(*reservoir), cmp_slot_by_idx);
    for (uint64_t i = 0; i < reservoir_filled; i++) {
      reservoir_slot_t *s = &reservoir[i];
      write_seq(s->seq, s->L, s->name, s->comment, s->comment_l);
      slot_free(s);
    }
    free(reservoir); reservoir = NULL;
    kept = reservoir_filled;
  }

  /* File-wide composition report for format. */
  if (args.action == ACTION_FORMAT && args.v) {
    report_composition(stderr, comp_freq, idx, NULL);
  }

  time_t time2 = time(NULL);
  if (args.v) {
    if (args.action == ACTION_SUBSAMPLE)
      fprintf(stderr, "Processed %" PRIu64 " sequence(s); kept %" PRIu64 ".\n", idx, kept);
    else
      fprintf(stderr, "Processed %" PRIu64 " sequence(s).\n", idx);
    struct timespec ts_end;
    clock_gettime(CLOCK_MONOTONIC, &ts_end);
    double elapsed = (double)(ts_end.tv_sec - ts_program.tv_sec)
                   + (double)(ts_end.tv_nsec - ts_program.tv_nsec) / 1e9;
    print_time((uint64_t) difftime(time2, time1), "process sequences");
    fprintf(stderr, "Total runtime: %.3fs\n", elapsed);
    print_peak_mb();
  }

  free_bed();
  close_files();
  return EXIT_SUCCESS;
}
