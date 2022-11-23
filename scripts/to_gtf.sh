#!/bin/bash

# Output format: GTF/GFF2
# Score col: min(1000, -10*log10(P-value))
# Name attribute: [motif name]_[sequence name][strand]

awk '
  BEGIN{
    IFS = "\t"
    OFS = "\t"
  }
  function score(x) {
    y = -10 * log(x) / log(10)
    y = int(y * 10) / 10.0
    y = y > 1000 ? 1000 : y
    return y
  }
  /^[^##]/ {
    if (NF == 9) {
      print $1,"yamscan","nucleotide_motif",$2,$3,score($6),$4,".","name \"" $5 "_" $1 $4 "\"; motif_id \"" $5 "\"; p_value \"" $6 "\"; score \"" $7 "\"; match \"" $9 "\";"
    } else if (NF == 10) {
      print $1,"yamscan","nucleotide_motif",$2,$3,score($6),$4,".","name \"" $5 "_" $1 $4 "\"; motif_id \"" $5 "\"; p_value \"" $6 "\"; score \"" $7 "\"; match \"" $9 "\"; q_value \"" $10 "\";"
    } else if (NF == 11) {
      print $3,"yamscan","nucleotide_motif",$4,$5,score($8),$6,".","name \"" $7 "_" $3 $6 "\"; motif_id \"" $7 "\"; p_value \"" $8 "\"; score \"" $9 "\"; match \"" $11 "\"; bed_range \"" $1 "\"; bed_id \"" $2 "\";"
    } else if (NF == 12) {
      print $3,"yamscan","nucleotide_motif",$4,$5,score($8),$6,".","name \"" $7 "_" $3 $6 "\"; motif_id \"" $7 "\"; p_value \"" $8 "\"; score \"" $9 "\"; match \"" $11 "\"; bed_range \"" $1 "\"; bed_id \"" $2 "\"; q_value \"" $12 "\";"
    } else {
      print "Error: Input is malformed; expected 9-12 fields, found " NF > "/dev/stderr"
      exit 1
    }
  }' < /dev/stdin

