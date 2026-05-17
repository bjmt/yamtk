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
#include <limits.h>
#include <time.h>
#include <stdint.h>
#include <inttypes.h>
#include <zlib.h>
#include "version.h"

#define DEFAULT_SEED                   4

/* ---- Args ---- */

typedef struct args_t {
  double   lambda;        /* -f, per-bp Poisson rate (random mode) */
  uint64_t min_spacing;   /* -M */
  int      seed;          /* -s */
  int      trim_names;    /* -r toggles preservation; default trims */
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

static void close_files(void) {
  if (files.m_open) fclose(files.m);
  if (files.i_open) gzclose(files.i);
  if (files.o_open) fclose(files.o);
  if (files.x_open) gzclose(files.x);
  if (files.O_open) fclose(files.O);
}

static void badexit(const char *msg) {
  if (msg && msg[0]) fprintf(stderr, "%s\n", msg);
  fprintf(stderr, "Run yamtk seed -h to see usage.\n");
  close_files();
  exit(EXIT_FAILURE);
}

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

  fprintf(stderr, "yamtk seed: not implemented yet (skeleton).\n");

  close_files();
  return EXIT_SUCCESS;
}
