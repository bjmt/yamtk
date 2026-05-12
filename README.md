# yamtk: Yet Another Motif ToolKit

* Motif scanning: [yamtk scan](#yamscan)
* Deduplicate overlapping motif hits: [yamtk dedup](#yamdedup)
* Motif enrichment: [yamtk enr](#yamenr)
* Higher-order sequence shuffling: [yamtk shuf](#yamshuf)
* Miscellaneous utility scripts: [Extra scripts](#extra-scripts)

## Installation

Requires a C compiler (tested with gcc/clang), GNU Make, and [Zlib](https://zlib.net).

```sh
git clone https://github.com/bjmt/yamtk  # Or download a recent release
cd yamtk
make
```

This will create the final binary within the project folder. Use `make install`
to move it to `/usr/local/bin`.

## Motivation

I occasionally find myself needing to scan motifs against the Arabidopsis
genome from the command line. Usually having the MEME suite installed locally,
I resort to using [fimo](https://meme-suite.org/meme/tools/fimo). Inevitably
the wait time exceeds my limited patience. Eventually I decided to perform my
own re-write of fimo, only this time I would only include features I need in
addition to making sure it could run fast enough to satisfy me. This effort led
to creating a faster replacement, yamscan, and then later a smattering of
additional related programs/utilities as my needs demanded.

I would also like to mention that the fast performance of these programs was
made possible by making use of several functions from Dr. Heng Li's amazing
[klib](https://github.com/attractivechaos/klib) C library.

## yamscan

A regular DNA/RNA scanner with a focus on simplicity and speed.

### Usage

```
yamtk v2.1.0  Copyright (C) 2026  Benjamin Jean-Marie Tremblay
Usage:  yamtk scan [options] [ -m motifs.txt | -1 CONSENSUS ] -s sequences.fa

 -m <str>   Motif file (MEME, JASPAR, HOMER, HOCOMOCO PCM). 1-50 bases wide.
 -1 <str>   Scan a single consensus sequence (ambiguity letters allowed).
            1-50 bases wide. Ignores -b, -t, -0, -p, -n.
 -s <str>   Input FASTA/FASTQ. Can be gzipped. Use '-' for stdin.
            Omit -s to print parsed motifs; omit -m/-1 to print sequence stats.
            Non-ACGTU characters are treated as gaps during scanning.
 -x <str>   BED file of ranges to restrict scanning to. Col 4 = range name,
            col 6 = strand. Can be gzipped. Disables -f.
 -o <str>   Output file. Default: stdout.
 -b <dbl,dbl,dbl,dbl>  Background probabilities for A,C,G,T. Default: from
            motif file (MEME only) or uniform.
 -f         Only scan the forward strand.
 -t <dbl>   P-value threshold. Default: 0.0001.
 -0         Report all hits with score >= 0 (no p-value threshold).
 -p <int>   Pseudocount for PWM generation. Default: 1.
 -n <int>   Motif sites for PPM->PCM conversion. Default: 1000.
 -M         Mask lower-case bases (skip scanning at those positions).
 -d         Deduplicate motif/sequence names (default: abort on duplicates).
            Incompatible with -x.
 -r         Do not trim motif (HOCOMOCO/JASPAR) and sequence names to the
            first word.
 -l         Load all sequences into memory (default: one at a time). Auto-set
            with stdin or -j > 1.
 -j <int>   Threads. Default: 1. Capped at number of input motifs.
 -g         Show progress bar (only useful with multiple motifs).
 -v / -w / -h   Verbose / very-verbose / help.
```

### Output

yamscan reports basic information about matches, including coordinates,
scores, P-values, percent of the scores from the max, and the actual match.
Additional information about the motifs and sequences is included in the
header, which can be used for calculating Q-values after the fact. The
coordinates are 1-based.

Example output:

```
##yamscan v2.0.0 [ -t 0.04 -m test/motif.jaspar -s test/dna.fa ]
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

One can also use yamscan to get basic information about motifs and sequences.
By only using yamscan with one of these at a time, the following is output:

```sh
$ yamtk scan -m test/motif.jaspar
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

$ yamtk scan -s test/dna.fa
##seq_num	seq_name	size	gc_pct	n_count
1	1	55	49.09	0
2	2	70	45.71	0
3	3	33	39.39	0

$ yamtk scan -s test/dna.fa -x test/dna.bed
##bed_range	bed_name	seq_num	seq_name	size	gc_pct	n_count
1:1-35(+)	A	1	1	35	51.43	0
2:11-48(-)	B	2	2	38	50.00	0
```

This mode shows the internal PWM representation of motifs, as well P-values
for the min and max possible scores (with some in-between scores). Basic
information about sequences is output, including size, GC percent, and the
number of non-DNA/RNA letters found. If a BED file is also provided then the
sequence stats are restricted to those ranges.

Finally, scanning can be restricted to only parts of the input sequences as
specified in a BED file. This will, of course, result in nearly linear
speed-ups to the runtime proportional to the fraction of the input sequences
being scanned. Example output:

```
##yamscan v2.0.0 [ -t 0.04 -m test/motif.jaspar -s test/dna.fa -x test/dna.bed ]
##MotifCount=1 MotifSize=5 BedCount=2 BedSize=73 SeqCount=3 SeqSize=158 GC=45.57% Ns=0
##bed_range	bed_name	seq_name	start	end	strand	motif	pvalue	score	score_pct	match
1:1-35(+)	A	1	30	34	+	1-motifA	0.0078125	4.874	73.4	CTCGC
2:11-48(-)	B	2	16	20	-	1-motifA	0.0078125	4.874	73.4	GCGAG
2:11-48(-)	B	2	43	47	-	1-motifA	0.015625	3.867	58.3	TCTAG
```

### Comparing yamscan and fimo

The two programs have slightly different defaults, so right out of the box they
will not give identical values for scores and P-values. However with some slight
tweaking to even things out we can see the outputs are nearly identical.

fimo:
```sh
$ fimo --bfile --uniform-- --motif-pseudo 1 --no-qvalue --text --thresh 0.02 test/motif.meme test/dna.fa
Using motif +motif of width 5.
Using motif -motif of width 5.
motif_id	motif_alt_id	sequence_name	start	stop	strand	score	p-value	q-value	matched_sequence
motif		1	30	34	+	4.87379	0.00781		CTCGC
motif		2	4	8	+	3.83495	0.0166		GTCGA
motif		2	16	20	-	4.87379	0.00781		CTCGC
motif		2	42	46	+	3.69903	0.0195		GTCTA
motif		2	43	47	-	3.91262	0.0137		CTAGA
motif		3	28	32	+	3.69903	0.0195		GTCTA
```

yamscan (also manually setting the nsites value found in the motif file):
```sh
$ yamtk scan -b 0.25,0.25,0.25,0.25 -p 1 -n 175 -s test/dna.fa -t 0.02 -m test/motif.meme
##yamscan v2.0.0 [ -b 0.25,0.25,0.25,0.25 -p 1 -n 175 -s test/dna.fa -t 0.02 -m test/motif.meme ]
##MotifCount=1 MotifSize=5 SeqCount=3 SeqSize=158 GC=45.57% Ns=0 MaxPossibleHits=292
##seq_name	start	end	strand	motif	pvalue	score	score_pct	match
1	30	34	+	motif	0.0078125	4.877	73.4	CTCGC
2	4	8	+	motif	0.0166015625	3.843	57.8	GTCGA
2	16	20	-	motif	0.0078125	4.877	73.4	GCGAG
2	42	46	+	motif	0.01953125	3.704	55.7	GTCTA
2	43	47	-	motif	0.013671875	3.919	59.0	TCTAG
3	28	32	+	motif	0.01953125	3.704	55.7	GTCTA
```

(Actually, you might notice one small difference: yamscan always returns the
sequence on the forward strand for each hit, whereas fimo will return reverse
complement sequences when hits are on the reverse strand. If you prefer the
fimo behaviour you can pipe the yamscan output to `scripts/flip_rc.sh`)

### Benchmarking

Using GNU Time on my MacbookPro M1 and the following equivalent commands to
record time elapsed and peak memory usage. fimo is from MEME v5.4.1.

Default yamscan settings (low-mem mode active):
```sh
yamtk scan -v -t 0.0001 -m motifs.txt -s seqs.fa > res.txt
```
yamtk scan with multi-threading (and low-mem mode implicitly disabled):
```sh
yamtk scan -j4 -v -t 0.0001 -m motifs.txt -s seqs.fa > res.txt
```
fimo with Q-values disabled and immediate printing of results:
```sh
fimo --verbosity 1 --thresh 0.0001 --text motifs.txt seqs.fa > res.txt
```
homer using motifs with logodds scores matching P = 0.0001:
```sh
scanMotifGenomeWide.pl motifs.txt seqs.fa > res.txt
```

|                                |      `yamscan`     |    `yamscan -j4`   |      `fimo`      |      `homer`      |
|:------------------------------:|:------------------:|:------------------:|:----------------:|:-----------------:|
| 100x1Kbp (100Kbp) +  10 motifs |    0.02s,   4.30MB |    0.02s,   9.91MB |    0.23s, 3.92MB |    0.48s,  6.64MB |
| 100x1Kbp (100Kbp) + 100 motifs |    0.20s,   5.89MB |    0.07s,  12.31MB |    1.96s, 4.44MB |    1.43s,  6.63MB |
| 100x10Kbp (1Mbp)  +  10 motifs |    0.10s,   4.28MB |    0.06s,  11.44MB |    2.44s, 4.20MB |    1.45s,  6.67MB |
| 100x10Kbp (1Mbp)  + 100 motifs |    0.70s,   6.02MB |    0.23s,  15.80MB |   23.24s, 4.77MB |   10.32s,  6.66MB |
|   TAIR10 (120Mbp) +  10 motifs |    6.89s,  41.08MB |    2.36s, 153.10MB | 4m41.99s, 4.01MB | 1m55.18s, 34.11MB |
|   TAIR10 (120Mbp) + 100 motifs | 1m06.29s,  41.59MB |   19.14s, 152.10MB |     (not run)    |     (not run)     |
|   GRCh38 (3.2Gbp) +  10 motifs | 3m04.80s, 249.50MB | 1m02.30s,   3.09GB |     (not run)    |     (not run)     |

### It's still not fast enough!

If you are unfortunate enough to be working with genomes sized in the billions
and find this is not fast enough, then if you are willing to lose out on the
dependency-free convenience of yamscan I recommend trying out
[MOODS](https://github.com/jhkorhonen/MOODS). This library makes use of several
filtering algorithms to significantly speed up scanning (whereas yamscan
dumbly scores every possible match for all motifs across all sequences). I have
found after some brief testing that when scanning hundreds of motifs across
sequences in the Mbp-Gbp range several-fold speed-ups can be achieved. (In my
limited testing I also noticed that for smaller scanning jobs MOODS can itself
be several times slower due to the it's large startup cost -- an additional
tradeoff to keep in mind.) Alternatively, if CPU time is meaningless to you and
you have access to a large number of cores you can surpass even these
impressive scanning times by making liberal use of yamscan's `-j` flag. See
the latest MOODS
[paper](https://academic.oup.com/bioinformatics/article/33/4/514/2726114) for
details.

## yamdedup

Remove overlapping motif hits (or any type of sequence range).

This program is meant to work with direct yamscan output, and can work with a
piped `stdin`. (However it can also work with regular BED files!) This means it 
expects sequence/motif/strand combinations to be in ascending order by their
start coordinates. Separate sequence/motif/strand combinations can be
interleaved (which is the case when yamscan is run using multiple threads),
but within each unique combination the entries must be coordinate sorted.
It tries its utmost to use the least amount of memory required to faithfully
remove overlapping hits across the entire file, but for complex inputs it
may still end up using several dozen MBs. As an example, running yamdedup
with a BED file containing ~30,000,000 unique hits duplicated three times
(final range count of ~120,000,000, or <2 GB gzip-compressed) for ~500 motifs
(interleaved) across the Arabidopsis genome runs in about 48 seconds and
reaches 80 MB peak memory usage. Of course this will vary wildly depending on
how much of the input is overlapping; if all ranges in a file are overlapping
in one big chain, then yamdedup will be forced to store the entire input in
memory at a time.

### Usage

```sh
yamtk v2.0.0  Copyright (C) 2022  Benjamin Jean-Marie Tremblay               
Usage:  yamtk dedup [options] -i [ results.txt | ranges.bed ]                    
                                                                              
 -i <str>   Filename of yamscan results file or a tab-delimited BED file      
            with at least six columns: (1) sequence name, (2) 0-based start,  
            (3) end, (4) motif name, (5) motif score, and (6) the strand.     
            If a yamscan results file is provided the P-value is used for     
            deciding which overlapping hit(s) to discard, otherwise for BED   
            the score column is used and lower scores are discarded. The input
            is assumed to be partially sorted: different sequence/motif/strand
            combinations can be interleaved, but the individual combinations  
            themselves must be sorted by their start coordinates. The         
            yamscan program outputs its results this way, so no additional    
            sorting is needed. Can be gzipped. Use '-' for stdin.             
 -o <str>   Filename to output results. By default output goes to stdout.     
 -s         Ignore strand when finding overlapping ranges.                    
 -m         Ignore motif name when finding overlapping ranges.                
 -0         Ignore scores when removing overlapping ranges, causing yamdedup  
            to simply remove overlapping ranges in the order they appear.     
 -S         Sort on range size instead of score/p-value (keeping larger ones).
 -r         Sort in opposite order (i.e., keep lower scores, higher p-values, 
            or smaller ranges).                                               
 -y         Force yamdedup to treat the input as a yamscan output file.       
 -b         Force yamdedup to treat the input as a BED file.                  
 -v         Verbose mode.                                                     
 -w         Very verbose mode.                                                
 -h         Print this help message. 
```

### Examples

Let us consider the following basic scenario:

```sh
$ yamtk scan -t 0.2 -m test/motif.jaspar -s test/dna.fa | head -n6
##yamscan v2.0.0 [ -t 0.2 -m test/motif.jaspar -s test/dna.fa ]
##MotifCount=1 MotifSize=5 SeqCount=3 SeqSize=158 GC=45.57% Ns=0 MaxPossibleHits=292
##seq_name	start	end	strand	motif	pvalue	score	score_pct	match
1	3	7	+	1-motifA	0.060546875	1.459	22.0	GCTGA
1	7	11	+	1-motifA	0.1640625	-1.165	-17.6	ACTGA
1	11	15	+	1-motifA	0.0654296875	1.236	18.6	ATCGA
```

In this example we can see that all three hits overlap, but the middle hit has
a much worse score. yamdedup will recognize this and simply remove only that hit:

```sh
$ yamtk scan -t 0.2 -m test/motif.jaspar -s test/dna.fa | head -n6 | yamtk dedup -i-
##yamscan v2.0.0 [ -t 0.2 -m test/motif.jaspar -s test/dna.fa ]
##MotifCount=1 MotifSize=5 SeqCount=3 SeqSize=158 GC=45.57% Ns=0 MaxPossibleHits=292
##seq_name	start	end	strand	motif	pvalue	score	score_pct	match
1	3	7	+	1-motifA	0.060546875	1.459	22.0	GCTGA
1	11	15	+	1-motifA	0.0654296875	1.236	18.6	ATCGA
```

Again, yamdedup won't mind if the input is a properly formatted BED file:

```sh
$ yamtk scan -t 0.2 -m test/motif.jaspar -s test/dna.fa | head -n6 | scripts/to_bed.sh | yamtk dedup -i-
1	2	7	1-motifA	12	+	1.459	22.0	0.060546875	.
1	10	15	1-motifA	11	+	1.236	18.6	0.0654296875	.
```

There are a few options available for controlling the behaviour of yamdedup.
For example, the default behaviour for BED files is to prioritize keeping
higher scores, but this can be reversed using `-r`:

```sh
$ yamtk scan -t 0.2 -m test/motif.jaspar -s test/dna.fa | head -n6 | scripts/to_bed.sh | yamtk dedup -i- -r
1	6	11	1-motifA	7	+	-1.165	-17.6	0.1640625	.
```

Now we can see that the middle lower-scoring hit was kept instead.

Other options include ignoring the strand of motif hits (`-s`) if you wish to
consider hits on opposite strands to be seen as overlapping, or even ignore the
motif names (`-m`) themselves to remove any overlapping range period.
Overlapping hits can be removed in the order they appear, irrespective of their
hit score, using `-0`, or even consider the size of the ranges (`-S`) instead
of the scores (perhaps more useful for non-motif BED inputs). (In the future I
may consider adding additional options such as requiring specific amounts of
overlap.)

## yamenr

Compute motif enrichment between a positive and negative (or shuffled-positive)
sequence set. Outputs a TSV with per-motif enrichment statistics and
Benjamini-Hochberg FDR-corrected q-values.

### Usage

```
yamtk v2.1.0  Copyright (C) 2026  Benjamin Jean-Marie Tremblay
Usage:  yamtk enr [options] -i positives.fa -m motifs.txt

 -i <str>   Positives FASTA. Can be gzipped. Use '-' for stdin (requires -n).
 -n <str>   Negatives FASTA. If omitted, positives are shuffled (see -k, -s).
 -m <str>   Motif file (MEME, JASPAR, HOMER, HOCOMOCO PCM).
 -o <str>   Output TSV file. Default: stdout.
 -b <dbl,dbl,dbl,dbl>  Background probabilities for A,C,G,T. Default: from
            motif file (MEME only) or uniform.
 -p <int>   Pseudocount for PWM generation. Default: 1.
 -N <int>   Motif sites for PPM->PCM conversion. Default: 1000.
 -d         Deduplicate motif/sequence names (default: abort on duplicates).
 -r         Do not trim motif (HOCOMOCO/JASPAR) and sequence names to the
            first word.
 -t <dbl>   P-value threshold for a hit. Default: 0.0001.
 -T <str>   Test mode: 'seqs' (default), 'sites', or 'ranksum'.
            seqs    = Fisher's exact on per-sequence hit presence.
            sites   = Fisher's exact on per-position hit rate.
            ranksum = Threshold-free Mann-Whitney U on max PWM score.
 -f         Only scan the forward strand.
 -k <int>   Shuffle k-mer size when -n is absent. Default: 2.
 -s <uint>  RNG seed for shuffling. Default: time-seeded.
 -q <dbl>   Only report rows with q-value <= this. Default: 1.0.
 -j <int>   Threads. Default: 1.
 -v / -w / -h   Verbose / very-verbose / help.
```

### Examples

Test enrichment against an external negative set:

```sh
$ yamtk enr -i positives.fa -n negatives.fa -m motifs.meme
```

Use a shuffled-positive null (k=2, deterministic seed) and report only
significant motifs (q ≤ 0.05) with 4 threads:

```sh
$ yamtk enr -i positives.fa -m motifs.meme -k 2 -s 42 -q 0.05 -j 4
```

Use threshold-free ranksum mode (AUC effect size):

```sh
$ yamtk enr -i positives.fa -n negatives.fa -m motifs.meme -T ranksum
```

### Output

Tab-separated with three comment header lines followed by one row per motif,
sorted ascending by q-value:

```
##yamenr v2.1.0 [ ... ]
##MotifCount=N PosSeqs=N NegSeqs=N NegSource=... TestMode=seqs Effect=seq_fold
##motif  motif_id  pos_n  pos_seq_hits  pos_site_hits  neg_n  neg_seq_hits  neg_site_hits  effect  log2_effect  pvalue  qvalue
```

`effect` and `log2_effect` meaning depends on `-T`:
- `seqs` / `sites` — fold enrichment (positive/negative hit rate ratio)
- `ranksum` — AUC (0.5 = no signal; range 0–1) / log-odds of AUC

## yamshuf

A regular DNA/RNA sequence shuffler with a focus on simplicity and speed.

### Usage

```
yamtk v2.2.0  Copyright (C) 2026  Benjamin Jean-Marie Tremblay
Usage:  yamtk shuf [options] -i sequences.fa

 -i <str>   Input FASTA/FASTQ. Can be gzipped. Use '-' for stdin.
            Non-ACGTU characters are treated as N (except with -l or -k 1).
 -k <int>   k-mer size for shuffling. Default: 3. k=1 uses Fisher-Yates;
            max k for Euler/Markov: 9.
 -o <str>   Output file. Default: stdout.
 -s <int>   RNG seed. Default: 4.
 -m         Markov shuffling: generates sequences with similar k-mer
            frequencies. Best for large sequences.
 -l         Linear k-mer shuffle (fast Fisher-Yates over k-mer blocks).
 -r <int>   Repeat shuffle N times per sequence; index appended to name.
 -R         Reset RNG to seed before each sequence instead of just once.
 -x <str>   BED file of ranges to restrict shuffling to. Bases inside
            ranges are shuffled in place; bases outside pass through
            unchanged. Incompatible with -p.
 -n         Output RNA instead of DNA. Only applies when k > 1 and -l
            is not used.
 -p         Print k-mer counts instead of shuffling (-i, -k, -o only).
 -v / -w / -h   Verbose / very-verbose / help.
```

### Overview of shuffling algorithms

yamshuf uses a few different shuffling implementations. When k = 1, it performs
a Fisher-Yates shuffle of the input sequences. For higher values of k, there
are three available methods. The default method (Euler) involves constructing
a k-mer edge graph from the k-mer counts of the input sequences and finding a
random Eulerian walk through all of the available k-mers (thus ending up with
the exact same number of k-mers in the final shuffled sequences). The `-m`
flag triggers the use of the Markov method, where a Markov chain is constructed
from the k-mer counts in the sequences and the resulting probabilities are used
for generating new sequences of equal length (thus ending up with new sequences
with similar, but not exactly the same, k-mer counts). Finally, the `-l` flag
results in the input sequences being split linearly into chunks (with length k)
which are then shuffled via the Fisher-Yates method.

### Benchmarking

Using GNU Time on my MacbookPro M1 and the following equivalent commands to
record time elapsed and peak memory usage. A fasta file containing 10 sequences
(each 10 Mbp long) with the alphabet ACGTN is used as input.

Default yamshuf settings (Euler shuffling, or regular Fisher-Yates for k = 1):
```sh
yamtk shuf -k $K -i 100Mbp.fa > shuffled.fa
```
yamtk shuf with Markov shuffling:
```sh
yamtk shuf -m -k $K -i 100Mbp.fa > shuffled.fa
```
yamtk shuf with linear shuffling:
```sh
yamtk shuf -l -k $K -i 100Mbp.fa > shuffled.fa
```
MEME v5.4.1 `fasta-shuffle-letters` tool (using the uShuffle library):
```sh
fasta-shuffle-letters -k $K 100Mbp.fa shuffled.fa
```

|                 |   `yamshuf`    |  `yamshuf -m`  |  `yamshuf -l`  |`fasta-shuffle-letters`|
|:---------------:|:--------------:|:--------------:|:--------------:|:---------------------:|
| 10x10Mbp, k = 1 | 0.40s, 19.46MB |      (n/a)     |      (n/a)     |     0.63s,  52.82MB   |
| 10x10Mbp, k = 2 | 1.80s, 19.49MB | 1.66s, 19.49MB | 0.36s, 19.52MB |     2.62s, 416.46MB   |
| 10x10Mbp, k = 3 | 1.89s, 19.42MB | 1.70s, 19.42MB | 0.32s, 19.42MB |     3.40s, 416.50MB   |
| 10x10Mbp, k = 4 | 1.93s, 19.49MB | 1.71s, 19.42MB | 0.30s, 19.42MB |     4.36s, 416.50MB   |
| 10x10Mbp, k = 5 | 2.01s, 19.46MB | 1.87s, 19.47MB | 0.29s, 19.42MB |     6.86s, 425.68MB   |
| 10x10Mbp, k = 6 | 2.06s, 19.66MB | 1.77s, 19.55MB | 0.34s, 19.42MB |    13.42s, 416.98MB   |
| 10x10Mbp, k = 7 | 2.53s, 20.35MB | 2.25s, 20.05MB | 0.30s, 19.42MB |    33.55s, 418.54MB   |
| 10x10Mbp, k = 8 | 2.63s, 23.87MB | 2.32s, 22.54MB | 0.28s, 19.46MB |       (not run)       |
| 10x10Mbp, k = 9 | 5.51s, 41.26MB | 4.56s, 34.72MB | 0.30s, 19.42MB |       (not run)       |

For some reason the `fasta-shuffle-letters` program is very memory hungry, much
more so than the original uShuffle program. However the standalone uShuffle
program did not include a fasta reader, meaning it could only take sequences
as command line arguments, severely limiting the maximum sequence size (and
thus I cannot benchmark it with the above sequence sizes).

### Limitations of the Euler and Markov methods

yamshuf is a rather limited program in that it only recognizes a
five letter alphabet (ACGTN or ACGUN; all other characters are recognized as
N). This allows the program to use hard-coded constants and graph structures.
This is in contrast to other tools such as uShuffle (and my own programs
[universalmotif](https://bioconductor.org/packages/universalmotif) and
[sequence-utils](https://github.com/bjmt/sequence-utils)) which are generically
coded to allow for any alphabet, but in turn (likely) allow for fewer compiler
optimizations. Additionally, yamshuf constructs a k-mer edge graph for the
complete set of possible k-mers for any k, even if some k-mers are absent in
the input sequences. This means that for increasing values of k, the graph
structure increases exponentially. As a result shuffling with high values of k
is impractical. uShuffle gets around this by building a k-mer edge graph using
only available k-mers, thus allowing for much higher values of k when shuffling
short sequences (for longer sequences, so many k-mers will likely be present
that it becomes affected by the same issue of impractically large edge graphs).

If there is not a requirement that the shuffled sequences have the exact same
k-mer counts as the input sequences, then none of these limitations apply when
using the linear method (`-l`). This is because yamshuf merely moves around
chunks of the input sequences without needing to count k-mers or build any
edge graph. (In fact, the higher the value of k, the fewer chunks there are
to move around, thus increasing the speed of the shuffling.)

## Extra scripts

A few extra utilities are included in the `scripts/` folder. These take the
yamscan results via `stdin` and output their results to `stdout`.

### Utilities requiring the `sort` program

- `add_qvals.sh`: Calculate Benjamini-Hochberg adjusted P-values (or Q-values)
  and add them as a tenth column. (Note: Do *not* deduplicate/remove
  overlapping hits before calculating Q-values.) Currently not compatible with
  yamscan run using the `-x` flag. (Generally I find it doesn't make much
  sense to calculate Q-values when scanning for motifs.) Be warned that this can
  be quite slow for big inputs, since it has to sort everything twice.
- `dedup_hits.sh`: Remove lower-scoring overlapping hits (of the same motif).
  Currently not compatible with yamscan run using the `-x` flag. This should
  only be used for very small inputs (eg <100,000 rows) unless you are willing
  to wait a while. *Note: As of yamscan v1.4 this scripts has been deprecated
  in favour of the yamdedup program.*
- `sort_coord.sh`: Sort the results by coordinate.
- `sort_motif.sh`: Sort the results by motif name.
- `sort_pval.sh`: Sort the results by P-value.

Extra arguments can be provided by setting a `SORT_ARGS` environmental variable
in the shell. For example:

```sh
$ SORT_ARGS="--buffer-size=1G" scripts/add_qvals.sh < results.txt
```

### Simple operations/format conversions

- `flip_rc.sh`: Reverse complement sequence matches on the reverse strand. This
  can be useful if you wish to see matches from the strand the motif was
  matched from, as the default is to always return the sequence from the forward
  strand.
- `std_kmers.sh`: Filter the output of `yamtk shuf -p` to only include k-mers
  containing standard letters (ACGT or ACGU).
- `to_bed.sh`: Convert the results to a BED6+4 format.
- `to_gff3.sh`: Convert the results to GFF3.
- `to_gtf.sh`: Convert the results to GTF/GFF2.

See the scripts for a description of the output formats.

### Example

In this example, the output of yamscan is first piped to `flip_rc.sh` to
reverse complement the reverse strand motif matches, then to `add_qvals.sh`
to calculate Q-values, to `dedup_hits.sh` to remove overlapping hits,
coordinate sorted with `sort_coord.sh`, and finally converted to BED with
`to_bed.sh`:

```sh
yamscan -t 0.05 -m test/motif.homer -s test/dna.fa \
  | scripts/flip_rc.sh \
  | scripts/add_qvals.sh \
  | scripts/dedup_hits.sh \
  | scripts/sort_coord.sh \
  | scripts/to_bed.sh \
  > res.clean.txt
```

## Compatible motif formats

The format will be auto-detected by yamscan. A brief overview of the
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

1. Consensus sequence (yamscan ignores this but requires something be here)

2. One-word name

3. Minimum logodds score (yamscan ignores this but will warn if missing and
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

- The `ALPHABET` line is checked by yamscan; if using a custom alphabet
  definition then it is ignored (it is up to the user to ensure the custom
  alphabet represents DNA/RNA)

- The `strands:` line is briefly examined, and if it does not match the yamscan
  settings for which strands to scan a warning will be emitted if using `-v`/`-w`

- Background values will be used if a line starting with `Background letter
  frequencies` and immediately followed by probabilities for A,C,G,T/U is found

Both full and minimal MEME motif files can be used (alongside its derivatives
DREME and STREME).

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

