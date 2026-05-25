/*
 *   yambkg: GC/length-matched background sequence sampling
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

#define ALLOC_CHUNK_SIZE        ((uint64_t) 256)
#define FASTA_LINE_LEN                        60
#define SEQ_NAME_MAX_CHAR       ((uint64_t) 512)
#define DEFAULT_GC_STEP                     0.05
#define DEFAULT_LEN_TOL                       50
#define DEFAULT_N_MAX                       0.10
#define DEFAULT_MAX_ATTEMPTS                1000
#define DEFAULT_N_PER_TARGET                   1

#define ERASE_ARRAY(ARR, LEN) memset(ARR, 0, sizeof(ARR[0]) * (LEN))

/* ---- Peak memory / elapsed-time reporting (copied from yamshuf.c) ---- */

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

static void print_time_h(const uint64_t s, const char *what) {
  if (s > 7200)      fprintf(stderr, "Needed %'.2f hours to %s.\n",   ((double)s / 60.0) / 60.0, what);
  else if (s > 120)  fprintf(stderr, "Needed %'.2f minutes to %s.\n", (double)s / 60.0, what);
  else if (s > 1)    fprintf(stderr, "Needed %'" PRIu64 " seconds to %s.\n", s, what);
}

/* ---- Character classification (independent of yamseq/yamenr tables) ---- */

/* 1 if A/C/G/T/U (either case), else 0. */
static const unsigned char is_std[256] = {
  ['A']=1, ['a']=1, ['C']=1, ['c']=1, ['G']=1, ['g']=1,
  ['T']=1, ['t']=1, ['U']=1, ['u']=1
};
/* 1 if G/C (either case), else 0. */
static const unsigned char is_gc[256] = {
  ['C']=1, ['c']=1, ['G']=1, ['g']=1
};

/* Reverse-complement lookup. Non-ACGTU map to 0 here; the caller checks
   and passes the original character through (so N stays N etc). */
static const unsigned char rc_tab[256] = {
  ['A']='T', ['a']='t', ['C']='G', ['c']='g',
  ['G']='C', ['g']='c', ['T']='A', ['t']='a',
  ['U']='A', ['u']='a'
};

/* ---- Args ---- */

typedef struct args_t {
  double   gc_step;
  int      len_tol;
  int      n_per_target;
  double   n_max;
  int      max_attempts;
  int      randomize_strand;
  int      no_replace;
  uint64_t seed;
  int      use_seed;
  int      v;
  int      w;
} args_t;

static args_t args = {
  .gc_step          = DEFAULT_GC_STEP,
  .len_tol          = DEFAULT_LEN_TOL,
  .n_per_target     = DEFAULT_N_PER_TARGET,
  .n_max            = DEFAULT_N_MAX,
  .max_attempts     = DEFAULT_MAX_ATTEMPTS,
  .randomize_strand = 0,
  .no_replace       = 0,
  .seed             = 0,
  .use_seed         = 0,
  .v                = 0,
  .w                = 0
};

/* ---- Files ---- */

typedef struct files_t {
  int     i_open, p_open, g_open, o_open, B_open;
  gzFile  i;          /* targets */
  gzFile  p;          /* pool */
  gzFile  g;          /* genome */
  FILE   *o;          /* output FASTA */
  FILE   *B;          /* output BED (genome mode) */
} files_t;

static files_t files = {
  .i_open=0, .p_open=0, .g_open=0, .o_open=0, .B_open=0
};

static void close_files(void) {
  if (files.i_open) gzclose(files.i);
  if (files.p_open) gzclose(files.p);
  if (files.g_open) gzclose(files.g);
  if (files.o_open && files.o != stdout) fclose(files.o);
  if (files.B_open) fclose(files.B);
}

/* ---- PRNG (xoroshiro128++, seeded via splitmix64) ---- */

#define ROTL(x, k) (((x) << (k)) | ((x) >> (64 - (k))))

typedef struct { uint64_t s[2]; } xrng_t;

static inline uint64_t splitmix64(uint64_t x) {
  uint64_t z = (x += 0x9e3779b97f4a7c15ULL);
  z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
  z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
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

static inline uint64_t xrand_below(xrng_t *r, uint64_t k) {
  if (k == 0) return 0;
  const uint64_t lim = UINT64_MAX - (UINT64_MAX % k);
  uint64_t x;
  do { x = xrand_r(r); } while (x >= lim);
  return x % k;
}

static xrng_t xrng;

/* ---- Sequence sets ---- */

typedef struct {
  char           **names;
  unsigned char  **seqs;
  uint64_t        *sizes;
  uint64_t         n;
  uint64_t         n_alloc;
} seq_set_t;

static seq_set_t targets = {0};
static seq_set_t pool    = {0};
static seq_set_t genome  = {0};

/* Parallel arrays sized to targets.n / pool.n holding the GC bin of each
   sequence (computed once at load time so we don't repeat the O(L) scan). */
static int *target_bin = NULL;
static int *pool_bin   = NULL;

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

static void badexit(const char *msg) {
  if (msg && *msg) fprintf(stderr, "%s\n", msg);
  fprintf(stderr, "Run yamtk bkg -h to see usage.\n");
  close_files();
  free_seq_set(&targets);
  free_seq_set(&pool);
  free_seq_set(&genome);
  free(target_bin);
  free(pool_bin);
  exit(EXIT_FAILURE);
}

/* Slurp a kseq stream into a seq_set_t. Names take the FASTA id only
   (comment dropped) to keep output names stable and uniquely usable
   downstream. */
static void load_seq_set(kseq_t *kseq, seq_set_t *set, const char *label) {
  int ret;
  while ((ret = kseq_read(kseq)) >= 0) {
    if (set->n + 1 > set->n_alloc) {
      uint64_t na = set->n_alloc + ALLOC_CHUNK_SIZE;
      char **t1 = realloc(set->names, sizeof(*set->names) * na);
      if (!t1) { kseq_destroy(kseq); badexit("Error: Failed to alloc seq names."); }
      set->names = t1;
      unsigned char **t2 = realloc(set->seqs, sizeof(*set->seqs) * na);
      if (!t2) { kseq_destroy(kseq); badexit("Error: Failed to alloc seqs."); }
      set->seqs = t2;
      uint64_t *t3 = realloc(set->sizes, sizeof(*set->sizes) * na);
      if (!t3) { kseq_destroy(kseq); badexit("Error: Failed to alloc seq sizes."); }
      set->sizes = t3;
      set->n_alloc = na;
    }
    if (kseq->name.l > SEQ_NAME_MAX_CHAR) {
      kseq_destroy(kseq); badexit("Error: Seq name too large.");
    }
    set->seqs[set->n]  = (unsigned char *)kseq->seq.s;
    kseq->seq.s = NULL;
    set->sizes[set->n] = kseq->seq.l;
    set->names[set->n] = malloc(kseq->name.l + 1);
    if (!set->names[set->n]) { kseq_destroy(kseq); badexit("Error: Failed to alloc name."); }
    memcpy(set->names[set->n], kseq->name.s, kseq->name.l);
    set->names[set->n][kseq->name.l] = '\0';
    set->n++;
  }
  if (ret == -2) { kseq_destroy(kseq); badexit("Error: Failed to parse FASTQ qualities."); }
  else if (ret < -2) { kseq_destroy(kseq); badexit("Error: Failed to read input."); }
  else if (!set->n) { kseq_destroy(kseq); badexit("Error: No sequences read from input."); }
  kseq_destroy(kseq);

  uint64_t total = 0;
  for (uint64_t i = 0; i < set->n; i++) total += set->sizes[i];
  if (args.v)
    fprintf(stderr, "%s: %'" PRIu64 " seq(s), %'" PRIu64 " bp\n",
      label, set->n, total);
}

/* ---- GC / window statistics ----
   Iterate over a window once, returning the GC fraction (over A/C/G/T
   bases — N excluded from the denominator) and the N fraction (1 - std/L).
   Returns 0 if the window has zero standard bases (caller should treat
   as "skip"). */
static int window_stats(const unsigned char *seq, uint64_t s, uint64_t e,
                        double *gc_out, double *n_frac_out) {
  uint64_t std = 0, gc = 0;
  for (uint64_t i = s; i < e; i++) {
    const unsigned char c = seq[i];
    std += is_std[c];
    gc  += is_gc[c];
  }
  const uint64_t L = e - s;
  if (std == 0) { *gc_out = 0.0; *n_frac_out = 1.0; return 0; }
  *gc_out     = (double)gc / (double)std;
  *n_frac_out = 1.0 - (double)std / (double)L;
  return 1;
}

/* Number of GC bins given step. Step is validated > 0 and <= 1 at parse. */
static int n_gc_bins(void) {
  int n = (int)ceil(1.0 / args.gc_step);
  if (n < 1) n = 1;
  return n;
}

static int gc_to_bin(double gc_frac) {
  int b = (int)floor(gc_frac / args.gc_step);
  const int nb = n_gc_bins();
  if (b < 0) b = 0;
  if (b >= nb) b = nb - 1;
  return b;
}

/* ---- Pool indexing ----
   Bin pool seqs by GC bin; within each bin keep an array of indices sorted
   by length so we can binary-search the length-tolerance slice. */

typedef struct {
  uint64_t *idx;   /* sorted by pool.sizes[idx[k]] ascending */
  uint64_t  n;
} pool_bin_idx_t;

static pool_bin_idx_t *pool_bins = NULL;
static int n_bins_global = 0;

static int cmp_pool_idx_by_len(const void *a, const void *b) {
  const uint64_t ia = *(const uint64_t *)a;
  const uint64_t ib = *(const uint64_t *)b;
  const uint64_t la = pool.sizes[ia];
  const uint64_t lb = pool.sizes[ib];
  if (la < lb) return -1;
  if (la > lb) return 1;
  return 0;
}

static void build_pool_bins(void) {
  n_bins_global = n_gc_bins();
  pool_bins = calloc(n_bins_global, sizeof(*pool_bins));
  if (!pool_bins) badexit("Error: Failed to alloc pool bins.");
  /* count per bin */
  for (uint64_t i = 0; i < pool.n; i++) pool_bins[pool_bin[i]].n++;
  for (int b = 0; b < n_bins_global; b++) {
    if (pool_bins[b].n) {
      pool_bins[b].idx = malloc(sizeof(uint64_t) * pool_bins[b].n);
      if (!pool_bins[b].idx) badexit("Error: Failed to alloc pool bin slot.");
      pool_bins[b].n = 0; /* reuse as fill pointer */
    }
  }
  for (uint64_t i = 0; i < pool.n; i++) {
    const int b = pool_bin[i];
    pool_bins[b].idx[pool_bins[b].n++] = i;
  }
  for (int b = 0; b < n_bins_global; b++) {
    if (pool_bins[b].n > 1) {
      qsort(pool_bins[b].idx, pool_bins[b].n, sizeof(uint64_t), cmp_pool_idx_by_len);
    }
  }
}

static void free_pool_bins(void) {
  if (!pool_bins) return;
  for (int b = 0; b < n_bins_global; b++) free(pool_bins[b].idx);
  free(pool_bins);
  pool_bins = NULL;
}

/* Lower bound: first k in bin with pool.sizes[idx[k]] >= want. */
static uint64_t bin_lower_bound(const pool_bin_idx_t *pb, uint64_t want) {
  uint64_t lo = 0, hi = pb->n;
  while (lo < hi) {
    uint64_t m = lo + (hi - lo) / 2;
    if (pool.sizes[pb->idx[m]] < want) lo = m + 1;
    else hi = m;
  }
  return lo;
}

/* Upper bound: first k in bin with pool.sizes[idx[k]] > want. */
static uint64_t bin_upper_bound(const pool_bin_idx_t *pb, uint64_t want) {
  uint64_t lo = 0, hi = pb->n;
  while (lo < hi) {
    uint64_t m = lo + (hi - lo) / 2;
    if (pool.sizes[pb->idx[m]] <= want) lo = m + 1;
    else hi = m;
  }
  return lo;
}

/* Pick a pool seq from a specific bin matching the length tolerance.
   With replacement: pick uniformly from the slice. Without replacement:
   collect unused indices in the slice into a scratch buffer and pick from
   that. Returns 1 on success (stores pool index in *out), 0 if no match. */
static int pick_pool_in_bin(int bin, uint64_t target_len,
                            uint8_t *used, uint64_t *scratch,
                            uint64_t *out) {
  if (bin < 0 || bin >= n_bins_global) return 0;
  const pool_bin_idx_t *pb = &pool_bins[bin];
  if (pb->n == 0) return 0;
  const uint64_t want_lo = (target_len > (uint64_t)args.len_tol)
                             ? target_len - (uint64_t)args.len_tol : 0;
  const uint64_t want_hi = target_len + (uint64_t)args.len_tol;
  const uint64_t lo = bin_lower_bound(pb, want_lo);
  const uint64_t hi = bin_upper_bound(pb, want_hi);
  if (lo >= hi) return 0;

  if (!args.no_replace) {
    const uint64_t k = lo + xrand_below(&xrng, hi - lo);
    *out = pb->idx[k];
    return 1;
  }
  /* without replacement: filter to unused */
  uint64_t n_avail = 0;
  for (uint64_t k = lo; k < hi; k++) {
    if (!used[pb->idx[k]]) scratch[n_avail++] = pb->idx[k];
  }
  if (n_avail == 0) return 0;
  const uint64_t k = xrand_below(&xrng, n_avail);
  *out = scratch[k];
  return 1;
}

/* Search outward from target_bin for the nearest non-empty bin that can
   serve a target of the given length. Deterministic tiebreak: prefer the
   higher-GC bin (target_bin + d) before the lower-GC one (target_bin - d).
   Returns 1 + bin (so 0 = none). */
static int find_fallback_pool_bin(int target_bin, uint64_t target_len,
                                  uint8_t *used, uint64_t *scratch,
                                  uint64_t *out) {
  for (int d = 1; d < n_bins_global; d++) {
    if (pick_pool_in_bin(target_bin + d, target_len, used, scratch, out))
      return target_bin + d + 1;
    if (pick_pool_in_bin(target_bin - d, target_len, used, scratch, out))
      return target_bin - d + 1;
  }
  return 0;
}

/* Render the GC-percentage range covered by `bin` into `out` (e.g. "50-60%").
   Bin ranges are half-open [lo, hi); the last bin is clamped to 100%. */
static void format_bin_range(int bin, char *out, size_t outsz) {
  double lo = (double)bin * args.gc_step * 100.0;
  double hi = lo + args.gc_step * 100.0;
  if (lo < 0.0)   lo = 0.0;
  if (hi > 100.0) hi = 100.0;
  snprintf(out, outsz, "%g-%g%%", lo, hi);
}

/* ---- Bin distribution dump (under -w) ---- */

static void print_bin_histogram(const int *bins, uint64_t n, const char *label) {
  if (n == 0 || n_bins_global <= 0) return;
  uint64_t *hist = calloc((size_t)n_bins_global, sizeof(uint64_t));
  if (!hist) return; /* best-effort; histogram is diagnostics only */
  for (uint64_t i = 0; i < n; i++) {
    const int b = bins[i];
    if (b >= 0 && b < n_bins_global) hist[b]++;
  }
  fprintf(stderr, "%s GC bins (step=%g%%):", label, args.gc_step * 100.0);
  int any = 0;
  char range[32];
  for (int b = 0; b < n_bins_global; b++) {
    if (hist[b] == 0) continue;
    format_bin_range(b, range, sizeof(range));
    fprintf(stderr, " %s=%" PRIu64, range, hist[b]);
    any = 1;
  }
  if (!any) fprintf(stderr, " (empty)");
  fputc('\n', stderr);
  free(hist);
}

/* ---- FASTA output ---- */

static void write_fasta(FILE *o, const char *name, const unsigned char *seq, uint64_t L) {
  fprintf(o, ">%s\n", name);
  for (uint64_t i = 0; i < L; i += FASTA_LINE_LEN) {
    const uint64_t take = (L - i < FASTA_LINE_LEN) ? L - i : FASTA_LINE_LEN;
    fwrite(seq + i, 1, take, o);
    fputc('\n', o);
  }
}

/* ---- Reverse complement (in place on a writable copy) ---- */
static void rc_inplace(unsigned char *seq, uint64_t L) {
  for (uint64_t i = 0, j = L - 1; i < j; i++, j--) {
    const unsigned char a = seq[i], b = seq[j];
    seq[i] = rc_tab[b] ? rc_tab[b] : b;
    seq[j] = rc_tab[a] ? rc_tab[a] : a;
  }
  if (L & 1) {
    const uint64_t m = L / 2;
    if (rc_tab[seq[m]]) seq[m] = rc_tab[seq[m]];
  }
}

/* ---- Pool mode ---- */

static int run_pool_mode(void) {
  if (args.v) fprintf(stderr, "Mode: pool subsample.\n");

  /* Compute GC bin for each pool sequence. */
  pool_bin = malloc(sizeof(*pool_bin) * pool.n);
  if (!pool_bin) badexit("Error: Failed to alloc pool_bin.");
  for (uint64_t i = 0; i < pool.n; i++) {
    double gc, nf;
    if (!window_stats(pool.seqs[i], 0, pool.sizes[i], &gc, &nf)) {
      pool_bin[i] = 0;
    } else {
      pool_bin[i] = gc_to_bin(gc);
    }
  }
  build_pool_bins();

  if (args.w) {
    print_bin_histogram(target_bin, targets.n, "Targets");
    print_bin_histogram(pool_bin,    pool.n,    "Pool   ");
  }

  /* Scratch arrays for the no-replace path. */
  uint8_t  *used    = args.no_replace ? calloc(pool.n, sizeof(uint8_t)) : NULL;
  uint64_t *scratch = args.no_replace ? malloc(sizeof(uint64_t) * pool.n) : NULL;
  if (args.no_replace && (!used || !scratch))
    badexit("Error: Failed to alloc pool bookkeeping.");

  uint64_t n_emitted = 0, n_fallback = 0, n_skipped = 0;

  for (uint64_t i = 0; i < targets.n; i++) {
    const uint64_t L_i = targets.sizes[i];
    const int      b_i = target_bin[i];
    for (int k = 0; k < args.n_per_target; k++) {
      uint64_t pick = 0;
      int is_fallback = 0;
      int got = pick_pool_in_bin(b_i, L_i, used, scratch, &pick);
      if (!got) {
        if (find_fallback_pool_bin(b_i, L_i, used, scratch, &pick)) {
          got = 1;
          is_fallback = 1;
          n_fallback++;
        }
      }
      if (!got) {
        char range[32]; format_bin_range(b_i, range, sizeof(range));
        n_skipped++;
        fprintf(stderr,
          "Warning: no pool match for target [%s] (len=%" PRIu64 ", gc=%s); skipping.\n",
          targets.names[i], L_i, range);
        continue;
      }
      if (args.no_replace) used[pick] = 1;
      if (args.w) {
        const int b_pick = pool_bin[pick];
        char r_tgt[32], r_pick[32];
        format_bin_range(b_i,    r_tgt,  sizeof(r_tgt));
        format_bin_range(b_pick, r_pick, sizeof(r_pick));
        fprintf(stderr,
          "  [%" PRIu64 "/%" PRIu64 "] tgt=[%s] L=%" PRIu64 " gc=%s -> "
          "pool=[%s] L=%" PRIu64 " gc=%s%s\n",
          i + 1, targets.n, targets.names[i], L_i, r_tgt,
          pool.names[pick], pool.sizes[pick], r_pick,
          is_fallback ? " [fallback]" : "");
      }
      write_fasta(files.o, pool.names[pick], pool.seqs[pick], pool.sizes[pick]);
      n_emitted++;
    }
  }

  free(used); free(scratch);
  free_pool_bins();
  free(pool_bin); pool_bin = NULL;

  if (args.v) {
    fprintf(stderr,
      "Emitted %'" PRIu64 " bkg seq(s); %'" PRIu64 " fell back to nearest bin; %'" PRIu64 " skipped.\n",
      n_emitted, n_fallback, n_skipped);
  }
  if (n_emitted == 0) {
    fprintf(stderr, "Error: No background sequences could be sampled.\n");
    return 1;
  }
  return 0;
}

/* ---- Genome mode ---- */

/* One attempt: pick a uniformly random valid (chrom, start) pair such that
   chrom_len[c] >= L_i. Returns 1 and stores result on success; 0 if no
   chrom is long enough (i.e. target longer than every chrom). */
static int sample_random_window(uint64_t L_i, uint64_t *chrom_out, uint64_t *start_out) {
  uint64_t total = 0;
  for (uint64_t c = 0; c < genome.n; c++) {
    if (genome.sizes[c] >= L_i) total += genome.sizes[c] - L_i + 1;
  }
  if (total == 0) return 0;
  uint64_t r = xrand_below(&xrng, total);
  for (uint64_t c = 0; c < genome.n; c++) {
    if (genome.sizes[c] < L_i) continue;
    const uint64_t span = genome.sizes[c] - L_i + 1;
    if (r < span) {
      *chrom_out = c;
      *start_out = r;
      return 1;
    }
    r -= span;
  }
  /* unreachable if total computed consistently */
  return 0;
}

/* Try `attempts` random windows looking for one whose GC bin matches
   `want_bin` (or any bin if want_bin == -1) and whose N fraction is
   within args.n_max. Returns 1 on success. `attempts_used` is set to the
   actual number of rejection-sampling iterations performed (1..attempts on
   success, == attempts on failure). */
static int try_window(uint64_t L_i, int want_bin, int attempts,
                      uint64_t *chrom_out, uint64_t *start_out,
                      int *attempts_used) {
  for (int a = 0; a < attempts; a++) {
    uint64_t c, s;
    if (!sample_random_window(L_i, &c, &s)) { *attempts_used = a + 1; return 0; }
    double gc, nf;
    window_stats(genome.seqs[c], s, s + L_i, &gc, &nf);
    if (nf > args.n_max) continue;
    if (want_bin >= 0) {
      if (gc_to_bin(gc) != want_bin) continue;
    }
    *chrom_out = c;
    *start_out = s;
    *attempts_used = a + 1;
    return 1;
  }
  *attempts_used = attempts;
  return 0;
}

static int run_genome_mode(void) {
  if (args.v) fprintf(stderr, "Mode: genome window sampling.\n");

  /* Sanity: at least one chrom must be at least as long as the shortest
     target. We don't preflight every target because the per-attempt loop
     handles the "no eligible chrom" case naturally. */
  uint64_t max_chrom = 0;
  for (uint64_t c = 0; c < genome.n; c++)
    if (genome.sizes[c] > max_chrom) max_chrom = genome.sizes[c];
  uint64_t min_target = UINT64_MAX;
  for (uint64_t i = 0; i < targets.n; i++)
    if (targets.sizes[i] < min_target) min_target = targets.sizes[i];
  if (max_chrom < min_target) {
    fprintf(stderr,
      "Error: All target sequences are longer than the longest chromosome "
      "(max chrom=%" PRIu64 " bp, min target=%" PRIu64 " bp).\n",
      max_chrom, min_target);
    return 1;
  }

  if (args.w) print_bin_histogram(target_bin, targets.n, "Targets");

  uint64_t n_emitted = 0, n_fallback = 0, n_skipped = 0;
  uint64_t total_attempts = 0, max_attempts_seen = 0;
  char namebuf[SEQ_NAME_MAX_CHAR + 64];

  /* Used for in-place reverse-complement; one re-usable buffer per draw. */
  unsigned char *rc_buf = NULL;
  uint64_t rc_buf_cap = 0;

  for (uint64_t i = 0; i < targets.n; i++) {
    const uint64_t L_i = targets.sizes[i];
    const int      b_i = target_bin[i];
    for (int k = 0; k < args.n_per_target; k++) {
      uint64_t c = 0, s = 0;
      int placed_bin = -2;   /* -2 = not placed; -1 = any-bin; else exact b */
      uint64_t pick_attempts = 0;
      int au = 0;

      int placed = try_window(L_i, b_i, args.max_attempts, &c, &s, &au);
      pick_attempts += (uint64_t)au;
      if (placed) placed_bin = b_i;

      if (!placed) {
        /* Nearest-bin fallback. Quarter-budget per bin so we don't blow up
           on heavily skewed genomes. */
        const int fb_attempts = (args.max_attempts / 4 > 1) ? args.max_attempts / 4 : 1;
        for (int d = 1; d < n_bins_global && !placed; d++) {
          if (b_i + d < n_bins_global) {
            if (try_window(L_i, b_i + d, fb_attempts, &c, &s, &au)) {
              placed = 1; placed_bin = b_i + d;
            }
            pick_attempts += (uint64_t)au;
            if (placed) break;
          }
          if (b_i - d >= 0) {
            if (try_window(L_i, b_i - d, fb_attempts, &c, &s, &au)) {
              placed = 1; placed_bin = b_i - d;
            }
            pick_attempts += (uint64_t)au;
            if (placed) break;
          }
        }
        if (placed) n_fallback++;
      }
      if (!placed) {
        /* Last-ditch: any bin at all (respecting only the N filter). */
        if (try_window(L_i, -1, args.max_attempts, &c, &s, &au)) {
          placed = 1; placed_bin = -1; n_fallback++;
        }
        pick_attempts += (uint64_t)au;
      }
      total_attempts += pick_attempts;
      if (pick_attempts > max_attempts_seen) max_attempts_seen = pick_attempts;

      if (!placed) {
        char range[32]; format_bin_range(b_i, range, sizeof(range));
        n_skipped++;
        fprintf(stderr,
          "Warning: failed to place target [%s] (len=%" PRIu64 ", gc=%s, attempts=%" PRIu64 "); skipping.\n",
          targets.names[i], L_i, range, pick_attempts);
        continue;
      }
      const int rc = args.randomize_strand ? (int)(xrand_r(&xrng) & 1ULL) : 0;
      const unsigned char *src = genome.seqs[c] + s;
      const unsigned char *emit_src = src;
      if (rc) {
        if (L_i > rc_buf_cap) {
          unsigned char *t = realloc(rc_buf, L_i);
          if (!t) { free(rc_buf); badexit("Error: Failed to alloc RC buffer."); }
          rc_buf = t; rc_buf_cap = L_i;
        }
        memcpy(rc_buf, src, L_i);
        rc_inplace(rc_buf, L_i);
        emit_src = rc_buf;
      }
      snprintf(namebuf, sizeof(namebuf),
               "%s:%" PRIu64 "-%" PRIu64 "(%c)",
               genome.names[c], s, s + L_i, rc ? '-' : '+');
      if (args.w) {
        const char *tag = (placed_bin == b_i) ? ""
                        : (placed_bin == -1)  ? " [fallback any-bin]"
                                              : " [fallback]";
        char r_tgt[32], r_pick[32];
        format_bin_range(b_i, r_tgt, sizeof(r_tgt));
        if (placed_bin >= 0) format_bin_range(placed_bin, r_pick, sizeof(r_pick));
        else                 snprintf(r_pick, sizeof(r_pick), "any");
        fprintf(stderr,
          "  [%" PRIu64 "/%" PRIu64 "] tgt=[%s] L=%" PRIu64 " gc=%s -> "
          "%s:%" PRIu64 "-%" PRIu64 "(%c) gc=%s attempts=%" PRIu64 "%s\n",
          i + 1, targets.n, targets.names[i], L_i, r_tgt,
          genome.names[c], s, s + L_i, rc ? '-' : '+',
          r_pick, pick_attempts, tag);
      }
      write_fasta(files.o, namebuf, emit_src, L_i);
      if (files.B_open) {
        fprintf(files.B,
                "%s\t%" PRIu64 "\t%" PRIu64 "\t%s\t.\t%c\n",
                genome.names[c], s, s + L_i, targets.names[i],
                rc ? '-' : '+');
      }
      n_emitted++;
    }
  }
  free(rc_buf);

  if (args.v) {
    fprintf(stderr,
      "Emitted %'" PRIu64 " bkg seq(s); %'" PRIu64 " fell back to nearest bin; %'" PRIu64 " skipped.\n",
      n_emitted, n_fallback, n_skipped);
  }
  if (args.w) {
    const uint64_t n_picks = n_emitted + n_skipped;
    const double mean_att = n_picks ? (double)total_attempts / (double)n_picks : 0.0;
    fprintf(stderr,
      "Attempts: total=%'" PRIu64 ", mean=%.1f/pick, max=%'" PRIu64 "/pick (limit -A %d).\n",
      total_attempts, mean_att, max_attempts_seen, args.max_attempts);
  }
  if (n_emitted == 0) {
    fprintf(stderr, "Error: No background sequences could be sampled.\n");
    return 1;
  }
  return 0;
}

/* ---- Usage / parsing ---- */

static void usage(void) {
  printf(
    "yamtk v%s  Copyright (C) %s  Benjamin Jean-Marie Tremblay\n"
    "Usage:  yamtk bkg [options] -i targets.fa[.gz] { -p pool.fa | -g genome.fa }\n"
    "\n"
    " -i <str>   Target FASTA/FASTQ ('-' = stdin).\n"
    " -p <str>   Candidate pool FASTA. Sample matching sequences from this set.\n"
    " -g <str>   Genome FASTA. Randomly sample matching windows from these seqs.\n"
    "            Exactly one of -p or -g is required.\n"
    " -o <str>   Output FASTA file (default: stdout).\n"
    " -B <str>   Output BED of sampled coordinates (genome mode only).\n"
    " -G <dbl>   GC bin step, 0 < step <= 1 (default: %.2f).\n"
    " -T <int>   Length tolerance in bp (pool mode only; default: %d).\n"
    " -n <int>   Background sequences per target (default: %d).\n"
    " -N <dbl>   Max N-fraction allowed in a window, in [0, 1] (default: %.2f).\n"
    " -A <int>   Max sampling attempts per target before fallback (default: %d).\n"
    " -r         Randomize strand (genome mode; reverse-complement on emit).\n"
    " -u         Sample pool without replacement (default: with replacement).\n"
    " -s <uint>  RNG seed (default: time-seeded).\n"
    " -v / -w / -h   Verbose / very-verbose / help.\n"
    , YAMTK_VERSION, YAMTK_YEAR,
    DEFAULT_GC_STEP, DEFAULT_LEN_TOL, DEFAULT_N_PER_TARGET,
    DEFAULT_N_MAX, DEFAULT_MAX_ATTEMPTS
  );
}

static inline int str_to_double(const char *s, double *res) {
  char *tmp; errno = 0;
  *res = strtod(s, &tmp);
  return (s == tmp || errno != 0 || *tmp != '\0');
}
static inline int str_to_int(const char *s, int *res) {
  char *tmp; errno = 0;
  long r = strtol(s, &tmp, 10);
  if (r > INT_MAX || r < INT_MIN) return 1;
  *res = (int)r;
  return (s == tmp || errno != 0 || *tmp != '\0');
}
static inline int str_to_uint64(const char *s, uint64_t *res) {
  char *tmp; errno = 0;
  *res = (uint64_t)strtoull(s, &tmp, 10);
  return (s == tmp || errno != 0 || *tmp != '\0');
}

int main_bkg(int argc, char **argv) {
  int has_target = 0, has_pool = 0, has_genome = 0, use_stdout = 1;
  int opt;

  struct timespec ts_program;
  clock_gettime(CLOCK_MONOTONIC, &ts_program);
  time_t t_start = time(NULL);

  while ((opt = getopt(argc, argv, "i:p:g:o:B:G:T:n:N:A:rus:vwh")) != -1) {
    switch (opt) {
      case 'i':
        if (files.i_open) badexit("Error: -i specified more than once.");
        has_target = 1;
        if (optarg[0] == '-' && optarg[1] == '\0') {
          files.i = gzdopen(fileno(stdin), "r");
        } else {
          files.i = gzopen(optarg, "r");
        }
        if (!files.i) {
          fprintf(stderr, "Error: Failed to open targets \"%s\" [%s]\n",
            optarg, strerror(errno)); badexit("");
        }
        files.i_open = 1;
        break;
      case 'p':
        if (files.p_open) badexit("Error: -p specified more than once.");
        has_pool = 1;
        if (optarg[0] == '-' && optarg[1] == '\0') {
          files.p = gzdopen(fileno(stdin), "r");
        } else {
          files.p = gzopen(optarg, "r");
        }
        if (!files.p) {
          fprintf(stderr, "Error: Failed to open pool \"%s\" [%s]\n",
            optarg, strerror(errno)); badexit("");
        }
        files.p_open = 1;
        break;
      case 'g':
        if (files.g_open) badexit("Error: -g specified more than once.");
        has_genome = 1;
        if (optarg[0] == '-' && optarg[1] == '\0') {
          files.g = gzdopen(fileno(stdin), "r");
        } else {
          files.g = gzopen(optarg, "r");
        }
        if (!files.g) {
          fprintf(stderr, "Error: Failed to open genome \"%s\" [%s]\n",
            optarg, strerror(errno)); badexit("");
        }
        files.g_open = 1;
        break;
      case 'o':
        if (files.o_open) badexit("Error: -o specified more than once.");
        use_stdout = 0;
        files.o = fopen(optarg, "w");
        if (!files.o) {
          fprintf(stderr, "Error: Failed to create output \"%s\" [%s]\n",
            optarg, strerror(errno)); badexit("");
        }
        files.o_open = 1;
        break;
      case 'B':
        if (files.B_open) badexit("Error: -B specified more than once.");
        files.B = fopen(optarg, "w");
        if (!files.B) {
          fprintf(stderr, "Error: Failed to create BED \"%s\" [%s]\n",
            optarg, strerror(errno)); badexit("");
        }
        files.B_open = 1;
        break;
      case 'G':
        if (str_to_double(optarg, &args.gc_step) || args.gc_step <= 0.0 || args.gc_step > 1.0)
          badexit("Error: -G must be a number in (0, 1].");
        break;
      case 'T':
        if (str_to_int(optarg, &args.len_tol) || args.len_tol < 0)
          badexit("Error: -T must be a non-negative integer.");
        break;
      case 'n':
        if (str_to_int(optarg, &args.n_per_target) || args.n_per_target < 1)
          badexit("Error: -n must be a positive integer.");
        break;
      case 'N':
        if (str_to_double(optarg, &args.n_max) || args.n_max < 0.0 || args.n_max > 1.0)
          badexit("Error: -N must be in [0, 1].");
        break;
      case 'A':
        if (str_to_int(optarg, &args.max_attempts) || args.max_attempts < 1)
          badexit("Error: -A must be a positive integer.");
        break;
      case 'r':
        args.randomize_strand = 1;
        break;
      case 'u':
        args.no_replace = 1;
        break;
      case 's':
        if (str_to_uint64(optarg, &args.seed))
          badexit("Error: -s must be a non-negative integer.");
        args.use_seed = 1;
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

  if (setlocale(LC_NUMERIC, "en_US") == NULL && args.v) {
    fprintf(stderr, "Warning: setlocale(LC_NUMERIC, \"en_US\") failed.\n");
  }

  if (!has_target)             badexit("Error: Missing -i (targets).");
  if (!has_pool && !has_genome) badexit("Error: One of -p or -g is required.");
  if (has_pool && has_genome)   badexit("Error: -p and -g are mutually exclusive.");
  if (files.B_open && has_pool) badexit("Error: -B is only valid with -g (genome mode).");
  if (args.randomize_strand && has_pool && args.v)
    fprintf(stderr, "Warning: -r has no effect in pool mode; ignoring.\n");

  if (use_stdout) {
    files.o = stdout;
    files.o_open = 1;
  }

  /* Seed PRNG (uses CLI seed if provided, else clock). */
  const uint64_t actual_seed = args.use_seed ? args.seed : (uint64_t)time(NULL);
  sxrand_r(&xrng, actual_seed);
  if (args.v) fprintf(stderr, "Using seed: %" PRIu64 "\n", actual_seed);

  /* Load targets. */
  {
    kseq_t *k = kseq_init(files.i);
    load_seq_set(k, &targets, "Targets");
  }

  /* Compute per-target GC bin. */
  n_bins_global = n_gc_bins();
  target_bin = malloc(sizeof(*target_bin) * targets.n);
  if (!target_bin) badexit("Error: Failed to alloc target_bin.");
  for (uint64_t i = 0; i < targets.n; i++) {
    double gc, nf;
    if (!window_stats(targets.seqs[i], 0, targets.sizes[i], &gc, &nf)) {
      target_bin[i] = 0;
    } else {
      target_bin[i] = gc_to_bin(gc);
    }
  }

  int rc = 0;
  if (has_pool) {
    kseq_t *k = kseq_init(files.p);
    load_seq_set(k, &pool, "Pool");
    rc = run_pool_mode();
  } else {
    kseq_t *k = kseq_init(files.g);
    load_seq_set(k, &genome, "Genome");
    rc = run_genome_mode();
  }

  if (files.o_open) fflush(files.o);
  if (files.B_open) fflush(files.B);

  close_files();
  free_seq_set(&targets);
  free_seq_set(&pool);
  free_seq_set(&genome);
  free(target_bin); target_bin = NULL;

  if (args.v) {
    struct timespec ts_end;
    clock_gettime(CLOCK_MONOTONIC, &ts_end);
    double elapsed = (double)(ts_end.tv_sec - ts_program.tv_sec)
                   + (double)(ts_end.tv_nsec - ts_program.tv_nsec) / 1e9;
    print_time_h((uint64_t)difftime(time(NULL), t_start), "sample background");
    fprintf(stderr, "Total runtime: %.3fs\n", elapsed);
    print_peak_mb();
  }

  return rc ? EXIT_FAILURE : EXIT_SUCCESS;
}
