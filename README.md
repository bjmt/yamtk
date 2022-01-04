# minimotif: A small DNA/RNA scanner for HOMER/JASPAR/MEME motifs

## Motivation

[fimo](https://meme-suite.org/meme/tools/fimo) is a great motif scanning
program, with complex functionality for scanning any kind of motifs/sequences
such as calculating P and Q values. It also has a tiny memory footprint, even
when scanning large numbers of motifs and sequences billions long. Unfortunately
all of this would appear to have come at a cost, as fimo can be excruciatingly
slow when attempting to scan a large collection of motifs against a moderately
sized genome. minimotif was created as a result.

This was achieved by sacrificing two things: complexity and memory. minimotif
has fewer options; chiefly among missing options are Q-values. This is because
calculating Q-values requires juggling all results within memory, which can be
rather problematic when results are in the millions. (fimo has the option to
output results immediately and skip Q-value calculation, but even so fimo is
quite slow.) Secondly, minimotif loads input sequences completely into memory.
This means that you are scanning the human genome, be prepared that minimotif
will need over 3GB of memory. (Alternatively, use the `-l` flag to limit this
to about 250MB and sacrifice a bit of performance in some cases.)

## Installation

```sh
git clone https://github.com/bjmt/minimotif
cd minimotif
make
```

This will create the final binary as `bin/minimotif` within the project folder.

## Usage

```
minimotif v1.0  Copyright (C) 2022  Benjamin Jean-Marie Tremblay

Usage:  minimotif [options] [ -m motifs.txt | -1 CONSENSUS ] -s sequences.fa

 -m <str>   Filename of text file containing motifs. Acceptable formats: MEME,
            JASPAR, HOMER. Must be 1-50 bases wide.
 -1 <str>   Instead of -m, scan a single consensus sequence. Ambiguity letters
            are allowed. Must be 1-50 bases wide. The -b, -t, -p and -n flags
            are unused.
 -s <str>   Filename of fasta-formatted file containing DNA/RNA sequences to
            scan. Use '-' for stdin. Omitting -s will cause minimotif to print
            the parsed motifs instead of scanning. Alternatively, solely
            providing -s and not -m/-1 will cause minimotif to return sequence
            stats. Any spaces found are not read into the final scanned
            sequence. Non-standard characters (i.e. other than ACGTU) will be
            read but are treated as gaps during scanning.
 -o <str>   Filename to output results. By default output goes to stdout.
 -b <dbl>   Comma-separated background probabilities for A,C,G,T. By default
            the background probability values from the motif file (MEME only)
            are used, or a uniform background is assumed. Used in PWM
            generation.
 -f         Only scan the forward strand.
 -t <dbl>   Threshold P-value. Default: 1e-05.
 -p <int>   Pseudocount for PWM generation. Default: 1. Must be a positive
            integer.
 -n <int>   Number of motif sites used in PWM generation. Default: 1000.
 -d         Deduplicate motif/sequence names. Default: abort. Duplicates will
            have the motif/sequence and line numbers appended.
 -r         Trim motif (JASPAR only) and sequence names to the first word.
 -l         Low memory mode. Only allows a single sequence in memory at a
            time. Reading sequences from stdin is disabled. If scanning many
            smaller sequences with large numbers of motifs, the impact on
            performance may be significant.
 -g         Print a progress bar during scanning. This turns off some of the
            messages printed by -w. Note that it's only useful if there is
            more than one input motif.
 -v         Verbose mode. Recommended when using for the first time with new
            motifs/sequences, as warnings about potential issues will only be
            printed when -v/-w are set.
 -w         Very verbose mode. Only recommended for debugging purposes.
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
##minimotif v1.0 [ -r -t 0.04 -m test/motif.jaspar -s test/dna.fa ]
##MotifCount=1 MotifSize=5 SeqCount=3 SeqSize=158 GC=45.57% Ns=0
##seqname	start	end	strand	motif	pvalue	score	score_pct	match
1  	30	34	+	1-motifA	0.0078125	4.874	73.4	CTCGC
1  	31	35	-	1-motifA	0.0341796875	2.482	37.4	TCGCG
2	4	8	+	1-motifA	0.0166015625	3.860	58.2	GTCGA
2	42	46	+	1-motifA	0.01953125	3.725	56.1	GTCTA
2	16	20	-	1-motifA	0.0078125	4.874	73.4	GCGAG
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

$ bin/minimotif -r -s test/dna.fa
##seqnum	line_num	seqname	size	gc_pct	n_count
1	1	1  	55	49.09	0
2	5	2	70	45.71	0
3	10	3	33	39.39	0
```

This mode shows the internal PWM representation of motifs, as well P-values
for the min and max possible scores (with some in-between scores). Basic
information about sequences is output, including size, GC percent, and the
number of non-DNA/RNA letters found.

## Benchmarking

Using GNU Time on my MacbookPro M1 and the following equivalent commands to
record time elapsed and peak memory usage (Q-values turned off in fimo):

```sh
time -v minimotif -v -t 0.0001 -m motifs.txt -s seqs.fa > res.txt
time -v minimotif -l -v -t 0.0001 -m motifs.txt -s seqs.fa > res.txt
time -v fimo --verbosity 1 --thresh 0.0001 --text motifs.txt seqs.fa > res.txt
```
|                                          |     `minimotif`    |   `minimotif -l`  |      `fimo`      |
|:----------------------------------------:|:------------------:|:-----------------:|:----------------:|
|         100 x 1Kbp (100Kbp) +  10 motifs |    0.05s,  4.42MB  |    0.06s,  4.28MB |    0.23s, 3.92MB |
|         100 x 1Kbp (100Kbp) + 100 motifs |    0.22s,  6.05MB  |    0.27s,  7.00MB |    1.96s, 4.44MB |
|          100 x 10Kbp (1Mbp) +  10 motifs |    0.15s,  5.38MB  |    0.16s,  4.34MB |    2.44s, 4.20MB |
|          100 x 10Kbp (1Mbp) + 100 motifs |    1.13s,  6.89MB  |    1.25s,  6.98MB |   23.24s, 4.77MB |
| Arabidopsis genome (120Mbp) +  10 motifs |   12.12s, 150.50MB |   13.75s, 33.39MB | 4m41.99s, 4.01MB |
| Arabidopsis genome (120Mbp) + 100 motifs | 1m58.48s, 157.90MB | 2m17.56s, 33.58MB |     (not run)    |

From the benchmarks the speed advantage of minimotif over fimo is clear. Also
obvious however, is the associated high memory usage costs. By sacrificing a
bit of speed, this can be alleviated somewhat when scanning very large
sequences by using the low-memory option, `-l`. When using this option only one
sequence is ever kept in memory, meaning that the max memory usage will be tied
to the size of the largest single sequence within the fasta file, instead of
the size of all sequences added together. In the example benchmark above when
scanning the Arabidopsis genome, minimotif only requests enough memory to hold
the largest chromosome (~30Mbp) rather then the entire genome (~120Mbp).

## Compatible motif formats

The format will be auto-detected by minimotif. A brief overview of the
requirements for each format follows.

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
position in the motif enclosed by square brackets.

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

HOMER motif header lines have three tab-separated essential elements, in
addition to needing to begin with the `>` character:

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

Both full and minimal MEME motif files can be used.

