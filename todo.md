## Todo

- yamseq: sequence manipulation (subsetting, masking, insertion, duplication,
  rna/dna conversion, reverse complement, seq sizes)
- yamshuf: DNA/RNA shuffler
  + streme paper says it preserves positions of separator characters; I should do this too
- yame: motif elicitation
  + (future) Seed-width trimming during refinement: after each refinement pass,
    check IC of the two outermost columns and drop them if below threshold; continue
    with the narrowed PPM. Would improve p-value sensitivity (fewer uninformative
    PWM degrees of freedom) in addition to the output-time IC trimming already
    implemented (MIN_IC_BITS).
- yamconv: convert between motif formats
- yamseed: seed sequences with motifs of interest

- yamscan: Make only using -x and -s output the subset sequences instead of just
  info about the ranges? --> actually make it a separate program to manipulate seqs
- get rid of infinite for loops, always use the `#define`'d bounds

- multithread low-mem mode?

