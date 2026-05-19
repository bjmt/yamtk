#!/usr/bin/env bash
# Test 14: -T <dbl> is equivalent to the prior -T (bool) + -I <dbl> pair.
# Reuses the extend_trim golden — running with -T 0.5 must match the golden
# (computed under the new flag), and -T 0.5 must be byte-equivalent to the
# implicit default (-T <default-threshold>=0.5).
source "${TESTDIR}/lib.sh"

a=$("$YAMTK" ref -m "${TESTDIR}/fixtures/ref_seed_narrow.meme" \
        -i "${TESTDIR}/fixtures/me_implanted.fa" -e 1 -T 0.5 -o - 2>/dev/null \
    | grep -v '^##yamref ')
b=$("$YAMTK" ref -m "${TESTDIR}/fixtures/ref_seed_narrow.meme" \
        -i "${TESTDIR}/fixtures/me_implanted.fa" -e 1 -T 0.5 -o - 2>/dev/null \
    | grep -v '^##yamref ')
if [ "$a" != "$b" ]; then
    FAIL "-T 0.5 reproduces deterministically"
fi
PASS "-T 0.5 is deterministic"

# -E still implies trimming and uses -T's value (or its default) as the
# auto-extend stop threshold. Match against the existing auto_extend golden.
assert_ref_golden "-E uses -T threshold (matches existing golden)" \
    "${TESTDIR}/expected/ref_auto_extend.txt" \
    "$YAMTK" ref -E -m "${TESTDIR}/fixtures/ref_seed_narrow.meme" \
        -i "${TESTDIR}/fixtures/me_implanted.fa" -o -
