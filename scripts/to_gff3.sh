#!/bin/bash

awk '
  BEGIN{
    IFS = "\t"
    OFS = "\t"
    print "##gff-version 3"
  }
  function score(x) {
    y = -10 * log(x) / log(10)
    y = int(y * 10) / 10.0
    y = y > 1000 ? 1000 : y
    return y
  }
  /^[^##]/ {
    if (NF == 9) {
      print $1,"minimotif","nucleotide_motif",$2,$3,score($6),$4,".","ID=" $5 "_" $1 $4 ";Name=" $5 ";P_Value=" $6 ";Score=" $7 ";Match=" $9 ";"
    } else if (NF == 10) {
      print $1,"minimotif","nucleotide_motif",$2,$3,score($6),$4,".","ID=" $5 "_" $1 $4 ";Name=" $5 ";P_Value=" $6 ";Score=" $7 ";Match=" $9 ";Q_value=" $10 ";"
    } else if (NF == 11) {
      print $3,"minimotif","nucleotide_motif",$4,$5,score($8),$6,".","ID=" $7 "_" $3 $6 ";Name=" $7 ";P_Value=" $8 ";Score=" $9 ";Match=" $11 ";Bed_Range=" $1 ";Bed_ID=" $2 ";"
    } else if (NF == 12) {
      print $3,"minimotif","nucleotide_motif",$4,$5,score($8),$6,".","ID=" $7 "_" $3 $6 ";Name=" $7 ";P_Value=" $8 ";Score=" $9 ";Match=" $11 ";Bed_Range=" $1 ";Bed_ID=" $2 ";Q_value=" $12 ";"
    } else {
      print "Error: Input is malformed; expected 9-12 fields, found " NF > "/dev/stderr"
      exit 1
    }
  }' < /dev/stdin

