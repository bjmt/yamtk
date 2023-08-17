#!/bin/bash

# Filter the k-mer count output from yamshuf to only include ACGT/ACGU.

awk -F "\t" '
  NR == 1 {
    printf "%s", $1
    for (i = 2; i < NF; i++) {
      keepCol[i] = gsub(/N/, ".", $i)
      if (keepCol[i] == 0) {
        printf "\t%s", $i
      }
    }
    printf "\n"
  }
  NR > 1 {
    printf "%s", $1
    for (i = 2; i < NF; i++) {
      if (keepCol[i] == 0) {
        printf "\t%s", $i
      }
    }
    printf "\n"
  }
  ' < /dev/stdin

