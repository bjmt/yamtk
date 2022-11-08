#!/bin/bash

awk '
  BEGIN {
    IFS = "\t"
    OFS = "\t"
  }
  function flip ( m )
  {
    RC = ""
    if (m != "") {
      split(m, x, "")
      for (i=1;i<=length(m);i++) {
        if (x[i] == "A") {
          RC = "T" RC
        } else if (x[i] == "a") {
          RC = "t" RC
        } else if (x[i] == "C") {
          RC = "G" RC
        } else if (x[i] == "c") {
          RC = "g" RC
        } else if (x[i] == "G") {
          RC = "C" RC
        } else if (x[i] == "g") {
          RC = "c" RC
        } else if (x[i] == "T") {
          RC = "A" RC
        } else if (x[i] == "t") {
          RC = "a" RC
        } else {
          RC = x[i] RC
        }
      }
    }
    return RC
  }
  /^##/ {
    print $0
  }
  /^[^##]/ {
    if (NF == 9 || NF == 10) {
      if ($4 == "-") {
        RC = flip($9)
        if (NF == 9) {
          print $1,$2,$3,$4,$5,$6,$7,$8,RC
        } else {
          print $1,$2,$3,$4,$5,$6,$7,$8,RC,$10
        }
      } else {
        print $0
      }
    } else if (NF == 11 || NF == 12) {
      if ($6 == "-") {
        RC = flip($11)
        if (NF == 11) {
          print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,RC
        } else {
          print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,RC,$12
        }
      } else {
        print $0
      }
    } else {
      print "Error: Input is malformed; expected 9-12 fields, found " NF > "/dev/stderr"
      exit 1
    }
  }
  ' < /dev/stdin

