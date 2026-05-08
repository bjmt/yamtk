## Todo

- yamenr:
  + sequence name duplication shouldn't matter, only for motif names
  + peak memory is printed twice
  + some -v text is printed to stdout?
  + print alt motif name too?
  + is the number of motif hits properly calculated for seqs mode?
  + NegSource: don't print filepath

- yamsort: use file re-reading to perform nicer sorting
- yamcompare: tomtom clone
- yamrefine: optimize motif for maximum enrichemnt in a +ve set vs -ve
  + optimize for addn specificity (extend edges) or less (trim low IC edges/letters)
- yamseq: sequence manipulation (subsetting, masking, insertion, formating, shuffling?)
- yamshuf: DNA/RNA shuffler
  + streme paper says it preserves positions of separator characters; I should do this too
- yame: motif elicitation
- yamconv: convert between motif formats
- yamseed: seed sequences with motifs of interest

- yamscan: Make only using -x and -s output the subset sequences instead of just
  info about the ranges? --> actually make it a separate program to manipulate seqs
- get rid of infinite for loops, always use the `#define`'d bounds

- multithread low-mem mode?

