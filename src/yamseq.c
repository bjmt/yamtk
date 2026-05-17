/*
 *   yamseq: Sequence manipulation (stats, rc, rna/dna, dup, subset, mask)
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
#include <zlib.h>
#include "kseq.h"
#include "version.h"

KSEQ_INIT(gzFile, gzread)

#define FASTA_LINE_LEN  60

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
  ACTION_MASK
} action_t;

static action_t parse_action(const char *s) {
  if (!s) return ACTION_NONE;
  if (strcmp(s, "stats")  == 0) return ACTION_STATS;
  if (strcmp(s, "rc")     == 0) return ACTION_RC;
  if (strcmp(s, "rna")    == 0) return ACTION_RNA;
  if (strcmp(s, "dna")    == 0) return ACTION_DNA;
  if (strcmp(s, "dup")    == 0) return ACTION_DUP;
  if (strcmp(s, "subset") == 0) return ACTION_SUBSET;
  if (strcmp(s, "mask")   == 0) return ACTION_MASK;
  return ACTION_NONE;
}

static const char *action_name(action_t a) {
  switch (a) {
    case ACTION_STATS:  return "stats";
    case ACTION_RC:     return "rc";
    case ACTION_RNA:    return "rna";
    case ACTION_DNA:    return "dna";
    case ACTION_DUP:    return "dup";
    case ACTION_SUBSET: return "subset";
    case ACTION_MASK:   return "mask";
    default:            return "<none>";
  }
}

/* ---- Args ---- */

typedef struct args_t {
  action_t action;
  uint64_t dup_n;
  int      mask_hard;
  int      use_bed;
  int      trim_names;
  int      progress;
  int      v;
  int      w;
} args_t;

static args_t args = {
  .action     = ACTION_NONE,
  .dup_n      = 0,
  .mask_hard  = 0,
  .use_bed    = 0,
  .trim_names = 1,
  .progress   = 0,
  .v          = 0,
  .w          = 0
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

static void badexit(const char *msg) {
  if (msg && msg[0]) fprintf(stderr, "%s\n", msg);
  fprintf(stderr, "Run yamtk seq -h to see usage.\n");
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
    "\n"
    " -i <str>   Input FASTA/FASTQ ('-' = stdin).\n"
    " -o <str>   Output (default: stdout).\n"
    " -x <str>   BED file (subset/mask only). Can be gzipped.\n"
    " -n <int>   Repeat count for dup.\n"
    " -N         Hard-mask (replace with N) instead of soft-mask. mask only.\n"
    " -r         Do not trim sequence names to the first word.\n"
    " -g         Show progress bar.\n"
    " -v / -w / -h   Verbose / very-verbose / help.\n"
    , YAMTK_VERSION, YAMTK_YEAR
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

/* ---- FASTA output ---- */

static void write_seq(const unsigned char *seq, const uint64_t L, const char *name,
                      const char *comment, const uint64_t comment_l) {
  if (!args.trim_names && comment_l) {
    fprintf(files.o, ">%s %s\n", name, comment);
  } else {
    fprintf(files.o, ">%s\n", name);
  }
  for (uint64_t i = 0; i < L; i += FASTA_LINE_LEN) {
    fprintf(files.o, "%.*s\n", FASTA_LINE_LEN, seq + i);
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
  /* -n only with dup, required there */
  if (args.dup_n > 0 && args.action != ACTION_DUP) {
    fprintf(stderr, "Error: -n is only valid with -a dup (got -a %s).\n",
      action_name(args.action));
    badexit("");
  }
  if (args.action == ACTION_DUP && args.dup_n == 0) {
    badexit("Error: -a dup requires -n <int>=1.");
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

  int opt;
  int use_stdout = 1;

  while ((opt = getopt(argc, argv, "a:i:o:x:n:Nrgvwh")) != -1) {
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
        if (str_to_uint64_t(optarg, &args.dup_n)) {
          badexit("Error: Failed to parse -n value.");
        }
        if (args.dup_n == 0) {
          badexit("Error: -n must be a positive integer.");
        }
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

  /* Per-action one-time header */
  if (args.action == ACTION_STATS) emit_stats_header();

  /* Stream sequences */
  kseq_t *kseq = kseq_init(files.i);
  int ret_val;
  uint64_t idx = 0;
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

  if (args.v) fprintf(stderr, "Processed %" PRIu64 " sequence(s).\n", idx);

  close_files();
  return EXIT_SUCCESS;
}
