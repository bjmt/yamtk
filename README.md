# minimotif: A small super-fast DNA/RNA motif scanner

## Installation

```sh
git clone https://github.com/bjmt/minimotif
cd minimotif
make
```

This will create the final binary as `bin/minimotif` within the project folder.

## Motivation

I occasionally find myself needing to scan motifs against the Arabidopsis
genome from the command line. Usually having the MEME suite installed locally,
I resort to using [fimo](https://meme-suite.org/meme/tools/fimo). Inevitably
the wait time exceeds my limited patience. Eventually I decided to perform my
own re-write of fimo, only this time I would only include features I need in
addition to making sure it could run fast enough to satisfy me. This effort led
to this small utility, minimotif.

## Usage

```
minimotif v1.2  Copyright (C) 2022  Benjamin Jean-Marie Tremblay

Usage:  minimotif [options] [ -m motifs.txt | -1 CONSENSUS ] -s sequences.fa

 -m <str>   Filename of text file containing motifs. Acceptable formats: MEME,
            JASPAR, HOMER, HOCOMOCO (PCM). Must be 1-50 bases wide.
 -1 <str>   Instead of -m, scan a single consensus sequence. Ambiguity letters
            are allowed. Must be 1-50 bases wide. The -b, -t, -0, -p, and -n
            flags are unused.
 -s <str>   Filename of fast(a|q)-formatted file containing DNA/RNA sequences
            to scan. Can be gzipped. Use '-' for stdin. Omitting -s will cause
            minimotif to print the parsed motifs instead of scanning.
            Alternatively, solely providing -s and not -m/-1 will cause
            minimotif to return sequence stats. Non-standard characters (i.e.
            other than ACGTU) will be read but are treated as gaps during
            scanning.
 -x <str>   Filename of a BED-formatted file containing ranges within
            sequences which scanning will be restricted to. Must have at least
            three tab-separated columns. If a fourth column is present it will
            be used as the range name. If a sixth strand column is present
            scanning will be restricted to the indicated strand. Note that -f
            is disabled when -x is used. It is recommended the BED be sorted
            for speed. Overlapping ranges are allowed, but be warned that they
            will individually scanned thus potentially introducing duplicate
            hits. The file can be gzipped.
 -o <str>   Filename to output results. By default output goes to stdout.
 -b <dbl,   Comma-separated background probabilities for A,C,G,T|U. By default
     dbl,   the background probability values from the motif file (MEME only)
     dbl,   are used, or a uniform background is assumed. Used in PWM
     dbl>   generation.
 -f         Only scan the forward strand.
 -t <dbl>   Threshold P-value. Default: 0.0001.
 -0         Instead of using a threshold, simply report all hits with a score
            of zero or greater. Useful for manual filtering.
 -p <int>   Pseudocount for PWM generation. Default: 1. Must be a positive
            integer.
 -n <int>   Number of motif sites used in PWM generation. Default: 1000.
 -d         Deduplicate motif/sequence names. Default: abort. Duplicates will
            have the motif/sequence numbers appended.
 -r         Don't trim motif (HOCOMOCO/JASPAR only) and sequence names to the
            first word.
 -l         Deactivate low memory mode. Normally only a single sequence is
            stored in memory at a time. Setting this flag allows the program
            to instead store the entire input into memory, which can help with
            performance in cases of slow disk access or gzipped files. Note
            that this flag is automatically set when reading sequences from
            stdin, and when multithreading is enabled.
 -j <int>   Number of threads minimotif can use to scan. Default: 1. Note that
            increasing this number will also increase memory usage slightly.
            The number of threads is limited by the number of motifs being
            scanned.
 -g         Print a progress bar during scanning. This turns off some of the
            messages printed by -w. Note that it's only useful if there is
            more than one input motif.
 -v         Verbose mode.
 -w         Very verbose mode.
 -h         Print this help message.
```

## Output

minimotif reports basic information about matches, including coordinates,
scores, P-values, percent of the scores from the max, and the actual match.
Additional information about the motifs and sequences is included in the
header, which can be used for calculating Q-values after the fact. The
coordinates are 1-based.

Example output:

```
##minimotif v1.2 [ -t 0.04 -m test/motif.jaspar -s test/dna.fa ]
##MotifCount=1 MotifSize=5 SeqCount=3 SeqSize=158 GC=45.57% Ns=0 MaxPossibleHits=292
##seq_name	start	end	strand	motif	pvalue	score	score_pct	match
1	30	34	+	1-motifA	0.0078125	4.874	73.4	CTCGC
1	31	35	-	1-motifA	0.0341796875	2.482	37.4	TCGCG
2	4	8	+	1-motifA	0.0166015625	3.860	58.2	GTCGA
2	16	20	-	1-motifA	0.0078125	4.874	73.4	GCGAG
2	42	46	+	1-motifA	0.01953125	3.725	56.1	GTCTA
2	43	47	-	1-motifA	0.015625	3.867	58.3	TCTAG
3	28	32	+	1-motifA	0.01953125	3.725	56.1	GTCTA
```

One can also use minimotif to get basic information about motifs and sequences.
By only using minimotif with one of these at a time, the following is output:

```sh
$ bin/minimotif -m test/motif.jaspar
----------------------------------------
Motif: 1-motifA (N1 L1)
MaxScore=6.64	Threshold=[exceeds max]
Motif PWM:
	A	C	G	T
1:	-3.74	1.66	-1.12	-1.73
2:	-7.23	0.38	-2.81	1.35
3:	-1.43	1.34	-2.59	-0.09
4:	-3.06	-7.23	1.02	0.88
5:	1.27	-0.49	-0.36	-3.36
Score=-24.14	-->     p=1
Score=-12.07	-->     p=0.82
Score=0.00	-->     p=0.11
Score=3.32	-->     p=0.022
Score=6.64	-->     p=0.00098
----------------------------------------

$ bin/minimotif -s test/dna.fa
##seq_num	seq_name	size	gc_pct	n_count
1	1	55	49.09	0
2	2	70	45.71	0
3	3	33	39.39	0
```

This mode shows the internal PWM representation of motifs, as well P-values
for the min and max possible scores (with some in-between scores). Basic
information about sequences is output, including size, GC percent, and the
number of non-DNA/RNA letters found.

Finally, scanning can be restricted to only parts of the input sequences as
specified in a BED file. This will, of course, result in nearly linear
speed-ups to the runtime proportional to the fraction of the input sequences
being scanned. Example output:

```
##minimotif v1.2 [ -t 0.04 -m test/motif.jaspar -s test/dna.fa -x test/dna.bed ]
##MotifCount=1 MotifSize=5 BedCount=2 BedSize=73 SeqCount=3 SeqSize=158 GC=45.57% Ns=0
##bed_range	bed_name	seq_name	start	end	strand	motif	pvalue	score	score_pct	match
1:1-35(+)	A	1	30	34	+	1-motifA	0.0078125	4.874	73.4	CTCGC
2:11-48(-)	B	2	16	20	-	1-motifA	0.0078125	4.874	73.4	GCGAG
2:11-48(-)	B	2	43	47	-	1-motifA	0.015625	3.867	58.3	TCTAG
```

## Benchmarking

Using GNU Time on my MacbookPro M1 and the following equivalent commands to
record time elapsed and peak memory usage.

### Default minimotif settings (low-mem mode active):
```sh
minimotif -v -t 0.0001 -m motifs.txt -s seqs.fa > res.txt
```
### minimotif with multi-threading (and low-mem mode implicitly disabled):
```sh
minimotif -j4 -v -t 0.0001 -m motifs.txt -s seqs.fa > res.txt
```
### fimo with Q-values disabled and immediate printing of results:
```sh
fimo --verbosity 1 --thresh 0.0001 --text motifs.txt seqs.fa > res.txt
```

|                                |     `minimotif`    |   `minimotif -j4`  |      `fimo`      |
|:------------------------------:|:------------------:|:------------------:|:----------------:|
| 100x1Kbp (100Kbp) +  10 motifs |    0.02s,   4.30MB |    0.02s,   9.91MB |    0.23s, 3.92MB |
| 100x1Kbp (100Kbp) + 100 motifs |    0.20s,   5.89MB |    0.07s,  12.31MB |    1.96s, 4.44MB |
| 100x10Kbp (1Mbp)  +  10 motifs |    0.10s,   4.28MB |    0.06s,  11.44MB |    2.44s, 4.20MB |
| 100x10Kbp (1Mbp)  + 100 motifs |    0.70s,   6.02MB |    0.23s,  15.80MB |   23.24s, 4.77MB |
|   TAIR10 (120Mbp) +  10 motifs |    6.89s,  41.08MB |    2.36s, 153.10MB | 4m41.99s, 4.01MB |
|   TAIR10 (120Mbp) + 100 motifs | 1m06.29s,  41.59MB |   19.14s, 152.10MB |     (not run)    |
|   GRCh38 (3.2Gbp) +  10 motifs | 3m04.80s, 249.50MB | 1m02.30s,   3.09GB |     (not run)    |

### It's still not fast enough!

If you are unfortunate enough to be working with genomes sized in the billions
and find this is not fast enough, then if you are willing to lose out on the
dependency-free convenience of minimotif I recommend trying out
[MOODS](https://github.com/jhkorhonen/MOODS). This library makes use of several
filtering algorithms to significantly speed up scanning (whereas minimotif
dumbly scores every possible match for all motifs across all sequences). I have
found after some brief testing that when scanning hundreds of motifs across
sequences in the Mbp-Gbp range several-fold speed-ups can be achieved.
(Alternatively, if CPU time is meaningless to you and you have access to a large
number of cores you can surpass even these impressive scanning times by
making liberal use of minimotif's `-j` flag.) See the latest MOODS
[paper](https://academic.oup.com/bioinformatics/article/33/4/514/2726114) for
details.

## Extra scripts

A few extra utilities are included in the `scripts/` folder. These take the
minimotif results via `stdin` and output their results to `stdout`.

- `add_qvals.sh`: Calculate Benjamini-Hochberg adjusted P-values (or Q-values)
  and add them as a tenth column. (Note: Do *not* deduplicate/remove
  overlapping hits before calculating Q-values.) Currently not compatible with
  minimotif run using the `-x` flag.
- `dedup_hits.sh`: Remove lower-scoring overlapping hits (of the same motif).
- `sort_coord.sh`: Sort the results by coordinate. Currently not compatible
  with minimotif run using the `-x` flag.
- `sort_motif.sh`: Sort the results by motif name.
- `sort_pval.sh`: Sort the results by P-value.
- `to_bed.sh`: Convert the output to a BED6+4 format.

All of these scripts use the `sort` program. Extra arguments (e.g.
`--buffer-size`) can be used by setting a `SORT_ARGS` variable. (Be careful
not to change the sorting behaviour.)

In this example, the output of minimotif is piped to `add_qvals.sh` to first
calculate Q-values, then to `dedup_hits.sh` to remove overlapping hits,
and finally coordinate sorted with `sort_coord.sh`:

```sh
bin/minimotif -t 0.1 -m test/motif.homer -s test/dna.fa \
  | scripts/add_qvals.sh \
  | scripts/dedup_hits.sh \
  | scripts/sort_coord.sh \
  > res.clean.txt
```

## Compatible motif formats

The format will be auto-detected by minimotif. A brief overview of the
requirements for each format follows.

### [HOCOMOCO](https://hocomoco11.autosome.ru)

Example HOCOMOCO motif:

```
>1-motifA
144.000000003	180.000000003	79.000000001	97.000000001
179.000000003	57.000000001	181.000000003	83.000000001
176.000000003	53.000000001	255.000000003	16.0000000004
9.0000000002	3.0000000001	16.0000000004	472.000000006
10.0000000002	5.0000000001	476.000000006	9.0000000002
473.000000006	5.0000000001	8.0000000002	14.0000000002
16.0000000004	415.000000006	27.0000000004	42.000000001
12.0000000002	8.0000000002	438.000000006	42.000000001
39.000000001	159.000000003	8.0000000002	294.000000006
144.000000003	233.000000003	72.000000001	51.000000001
268.000000006	57.000000001	93.000000001	82.000000001
```

HOCOMOCO motifs have a header line starting with `>` followed by the motif
name. Only mononucleotide count matrices (PCM) can be used. The counts are
split into four columns (A,C,G,T/U). These counts need not be integers. The
headers cannot contain the tab character.

### [JASPAR](https://jaspar.genereg.net/docs/)

Example JASPAR motif:

```
>1-motifA
A [   3   0  16   5 106 ]
C [ 139  57 111   0  31 ]
G [  20   6   7  89  34 ]
T [  13 112  41  81   4 ]
```

JASPAR motifs each have a header line starting with `>` followed by the motif
name. Counts for each DNA letter at each position follow this line. Each row
starts the character for each DNA letter, and is followed by counts for each
position in the motif enclosed by square brackets. These counts must be
integers.

### [HOMER](http://homer.ucsd.edu/homer/motif/creatingCustomMotifs.html)

Example HOMER motif:

```
>CYCKA	1-motifA	5.30478607528482
0.015	0.798	0.113	0.074
0.000	0.325	0.033	0.641
0.094	0.630	0.040	0.236
0.028	0.000	0.510	0.463
0.607	0.177	0.192	0.024
```

HOMER motif header lines have (at least) three tab-separated essential elements,
in addition to needing to begin with the `>` character:

1. Consensus sequence (minimotif ignores this but requires something be here)

2. One-word name

3. Minimum logodds score (minimotif ignores this but will warn if missing and
   running with `-v`/`-w`)

This header line is immediately followed by the motif probability values (first
column for A, second for C, third for G, and fourth for T/U).

### [MEME](https://meme-suite.org/meme/doc/meme-format.html)

Example MEME motif:

```
MEME version 5

ALPHABET= ACGT

strands: + -

Background letter frequencies
A 0.249493 C 0.257606 G 0.249493 T 0.243408

MOTIF 1-motifA
letter-probability matrix: alength= 4 w= 5 nsites= 175 E= 0
 0.015363  0.797841  0.113101  0.073695 
 0.000245  0.325285  0.033308  0.641162 
 0.093911  0.629765  0.039881  0.236444 
 0.027719  0.000057  0.509626  0.462598 
 0.606981  0.176988  0.191848  0.024182 
```

There are three essential elements to a MEME motif file:

1. MEME version information

2. For each motif, a line that starts with `MOTIF ` and contains the one-word
   motif name (extra text is ignored)

3. For each motif, a line that starts with `letter-probability matrix` and is
   immediately followed by the motif probability values (first column for A,
   second for C, third for G, fourth for T/U)

Other optional lines:

- The `ALPHABET` line is checked by minimotif; if using a custom alphabet
  definition then it is ignored (it is up to the user to ensure the custom
  alphabet represents DNA/RNA)

- The `strands:` line is briefly examined, and if it does not match the minimotif
  settings for which strands to scan a warning will be emitted if using `-v`/`-w`

- Background values will be used if a line starting with `Background letter
  frequencies` and immediately followed by probabilities for A,C,G,T/U is found

Both full and minimal MEME motif files can be used (alongside its derivatives
DREME and STREME).

