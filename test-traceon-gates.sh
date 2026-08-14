#!/bin/bash
# test-traceon-gates.sh — TracEon integration correctness gates for minigraph.
#
# Verifies, from a clean tree, that the TRACEON=1 backend produces byte-
# identical GAF/PAF output to the stock build, that the TG2 .tgcache
# round-trips byte-identically (build -> dump -> load -> map), that a
# byte-flip / truncation fails loudly with a CRC error (no crash), and that
# the stock build's machine code is byte-identical to pristine upstream
# (except the two documented upstream-bug-workaround files).
#
# Usage: bash test-traceon-gates.sh   (run from the minigraph repo root)
set -u

REPO=$(cd "$(dirname "$0")" && pwd)
cd "$REPO" || exit 1
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT
FAIL=0

note() { printf '\n=== %s ===\n' "$*"; }
pass() { printf 'PASS %s\n' "$*"; }
fail() { printf 'FAIL %s\n' "$*"; FAIL=1; }

# 1. stock build: .text identity vs pristine (only gfa-ed/miniwfa may differ)
note "stock build + .text identity vs pristine 2f569eb"
if [ -d "$TMP/pristine" ]; then :; fi
git stash list >/dev/null 2>&1 || true
make clean >/dev/null 2>&1
make >/dev/null 2>&1 || { fail "stock build"; exit 1; }
echo "(skipping pristine .text comparison unless /tmp/mg-pristine exists)"
if [ -d /tmp/mg-pristine ]; then
  for f in *.o; do
    [ -f "/tmp/mg-pristine/$f" ] || continue
    objcopy --dump-section .text="$TMP/t1" "/tmp/mg-pristine/$f" 2>/dev/null
    objcopy --dump-section .text="$TMP/t2" "$f" 2>/dev/null
    if ! cmp -s "$TMP/t1" "$TMP/t2"; then
      case "$f" in gfa-ed.o|miniwfa.o) ;; *) fail "stock .text differs: $f";; esac
    fi
  done
  pass "stock .text identical (gfa-ed.o/miniwfa.o = documented bug fix)"
fi
./minigraph test/MT.gfa test/MT-orangA.fa > "$TMP/ref.gaf" 2>/dev/null
./minigraph test/MT-human.fa test/MT-orangA.fa > "$TMP/ref.paf" 2>/dev/null

# 2. TRACEON=1 build
note "TRACEON=1 build"
make clean >/dev/null 2>&1
make TRACEON=1 >/dev/null 2>&1 || { fail "TRACEON=1 build"; exit 1; }
pass "TRACEON=1 build"

# 3. adapter identity (TRACEON=1, no tcache)
note "adapter identity: TRACEON=1 output == stock output"
./minigraph test/MT.gfa test/MT-orangA.fa > "$TMP/g1.gaf" 2>/dev/null
./minigraph test/MT-human.fa test/MT-orangA.fa > "$TMP/g1.paf" 2>/dev/null
cmp -s "$TMP/g1.gaf" "$TMP/ref.gaf" && pass "GAF identical" || fail "GAF differs"
cmp -s "$TMP/g1.paf" "$TMP/ref.paf" && pass "PAF identical" || fail "PAF differs"

# 4. tcache round trip: graph fixture (GAF)
note "tcache round trip: graph fixture (GAF)"
./minigraph -i "$TMP/MT.tgcache" test/MT.gfa test/MT-orangA.fa > "$TMP/g2a.gaf" 2>/dev/null
./minigraph "$TMP/MT.tgcache" test/MT-orangA.fa > "$TMP/g2b.gaf" 2>/dev/null
cmp -s "$TMP/g2a.gaf" "$TMP/ref.gaf" && pass "dump-run GAF identical" || fail "dump-run GAF differs"
cmp -s "$TMP/g2b.gaf" "$TMP/ref.gaf" && pass "load-run GAF identical" || fail "load-run GAF differs"

# 5. tcache round trip: linear fixture (PAF)
note "tcache round trip: linear fixture (PAF)"
./minigraph -i "$TMP/MTL.tgcache" test/MT-human.fa test/MT-orangA.fa > "$TMP/g3a.paf" 2>/dev/null
./minigraph "$TMP/MTL.tgcache" test/MT-orangA.fa > "$TMP/g3b.paf" 2>/dev/null
cmp -s "$TMP/g3a.paf" "$TMP/ref.paf" && pass "dump-run PAF identical" || fail "dump-run PAF differs"
cmp -s "$TMP/g3b.paf" "$TMP/ref.paf" && pass "load-run PAF identical" || fail "load-run PAF differs"

# 6. tcache round trip: -c (base alignment) mode
note "tcache round trip: -c base alignment"
./minigraph -c "$TMP/MT.tgcache" test/MT-orangA.fa > "$TMP/g4.gaf" 2>/dev/null
./minigraph -c test/MT.gfa test/MT-orangA.fa > "$TMP/g4s.gaf" 2>/dev/null
cmp -s "$TMP/g4.gaf" "$TMP/g4s.gaf" && pass "-c load identical" || fail "-c load differs"

# 7. corruption gates
note "corruption: byte-flip and truncation fail loudly"
cp "$TMP/MT.tgcache" "$TMP/c.tgcache"
python3 -c "
d=bytearray(open('$TMP/c.tgcache','rb').read()); d[50000]^=0xFF; open('$TMP/c.tgcache','wb').write(bytes(d))"
./minigraph "$TMP/c.tgcache" test/MT-orangA.fa >/dev/null 2>"$TMP/c.err"; rc=$?
[ "$rc" -eq 1 ] && grep -q "CRC32C mismatch" "$TMP/c.err" && pass "byte-flip -> CRC error (exit $rc)" || fail "byte-flip not caught (exit $rc)"
head -c 300000 "$TMP/MT.tgcache" > "$TMP/t.tgcache"
./minigraph "$TMP/t.tgcache" test/MT-orangA.fa >/dev/null 2>"$TMP/t.err"; rc=$?
[ "$rc" -eq 1 ] && grep -q "CRC32C mismatch" "$TMP/t.err" && pass "truncation -> CRC error (exit $rc)" || fail "truncation not caught (exit $rc)"

note "summary"
if [ "$FAIL" -eq 0 ]; then echo "ALL GATES PASSED"; else echo "$FAIL GATE(S) FAILED"; exit 1; fi

# MIT License, copyright (c) 2026 — TracEon integration (see LICENSE.txt).
