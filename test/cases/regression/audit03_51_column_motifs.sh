#!/usr/bin/env bash
# Audit: HOMER/HOCOMOCO readers and the JASPAR trailing value let column 50
# through, writing past the 50-position PWM in scan and enr.
source "$TESTDIR/lib.sh"
A="$TESTDIR/fixtures/audit"
for f in 51_rows.homer 51_counts.jaspar 51_rows.hocomoco; do
    assert_exit_nonzero "scan $f" "$YAMTK" scan -m "$A/$f" -s "$TESTDIR/dna.fa"
    assert_exit_nonzero "enr $f"  "$YAMTK" enr -i "$TESTDIR/dna.fa" -m "$A/$f" -s 1
done
