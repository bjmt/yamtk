## Todo

- rename repo to yam-toolkit (yet another motif toolkit)
  + yamscan
  + yamdedup: remove overlapping hits from yamscan
  + yamsort: use file re-reading to perform nicer sorting
  + yamcompare: tomtom clone
  + yamrefine: optimize motif for maximum enrichemnt in a +ve set vs -ve
    - optimize for addn specificity (extend edges) or less (trim low IC edges/letters)
    - maximum enrichment in target sequences vs maximize under-enrichment in background
  + yamseq: sequence manipulation (subsetting, masking, insertion, formating, shuffling?)
  + yamshuf: DNA/RNA shuffler
    - streme paper says it preserves positions of separator characters; I should do this too
  + yame: motif elicitation
  + yamenr: motif enrichment
  + yamconv: convert between motif formats
  + yamseed: seed sequences with motifs of interest

- yamscan: Make only using -x and -s output the subset sequences instead of just
  info about the ranges? --> actually make it a separate program to manipulate seqs
- minimotif: add free(line) when doing a badexit() in read functions
- minimotif: in `peek_through_seqs`, do I really need to do `kseq_rewind`?
  + probably doesn't matter either way, the function doesn't actually do anything
    computationally intensive, it only resets a few values in the stream struct
- minimotif: maybe add an option to taking in stdin that doesn't require reading
  in the whole the whole sequence input before scanning; this would be giving up
  making a nice header of course... --> could be triggered by a -L flag or some
  such to force low-mem mode when stdin is used
- minimotif: should I really be ignoring the nsites value in meme motif?

- calculate MaxPossibleHits when `-x` is used

- get rid of infinite for loops, always use the `#define`'d bounds

- consensus element scanning: set # of mismatches
  + note that it's not so simple, since ambiguity codes
    mess up simple score decreases to control exact #
    of mismatches

- multithread low-mem mode?

- make a program to compare a set of +ve and
  -ve sequences and find the best score thresh

