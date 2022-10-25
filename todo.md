## Todo

- consensus element scanning: set # of mismatches
  + note that it's not so simple, since ambiguity codes
    mess up simple score decreases to control exact #
    of mismatches

- multithread low-mem mode?

- make c program for dedupping
  + iterate through each motif/sequence/strand combination;
    read in just the matching data, sort it in mem, then print
    dedupped

- make a program to compare a set of +ve and
  -ve sequences and find the best score thresh

- make a tomtom clone

- add an option to only scan in specified bed regions
  + instead of fixing `add_qvals.sh` and `dedup_hits.sh`, make c versions

- bug?
  + `echo "> A\nATGTGCAC\n" | bin/minimotif -s-`
  + `echo "> A\nATGTGCAC\n" | bin/minimotif -w -s-`

- benchmark different default values for:
  + `MAX_MOTIF_SIZE`
  + `ALLOC_CHUNK_SIZE`
  + `SEQ_REALLOC_SIZE`

