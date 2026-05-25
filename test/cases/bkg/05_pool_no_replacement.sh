#!/usr/bin/env bash
# With -u (no replacement), the same pool sequence must not appear twice in
# the output even when its bin is the only viable match for multiple targets.
source "$TESTDIR/lib.sh"

out=$("$YAMTK" bkg -i "$TESTDIR/fixtures/enr_pos.fa" -p "$TESTDIR/fixtures/enr_neg.fa" \
        -s 1 -G 0.1 -T 50 -u 2>/dev/null)

names=$(echo "$out" | grep '^>' | sort)
uniq_count=$(echo "$names" | uniq | wc -l | tr -d ' ')
total_count=$(echo "$names" | wc -l | tr -d ' ')

if [ "$uniq_count" -ne "$total_count" ]; then
    FAIL "bkg pool -u emitted duplicate names: $total_count total / $uniq_count unique"
fi
PASS "bkg pool -u emits unique names ($total_count total, all unique)"
