#!/usr/bin/env bash
# compile_and_run_tests.sh
# Compiles all shared units and the console self-test runner, then runs it.
# Usage: bash fpc/compile_and_run_tests.sh [--verbose]

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SHARED="$SCRIPT_DIR/_shared"
OUT="$SCRIPT_DIR/_build"
RUNNER_SRC="$SCRIPT_DIR/_unittests/run_self_tests_console.pas"
RUNNER_BIN="$OUT/run_self_tests_console"
FPC_FLAGS="-Mobjfpc -Fu$SHARED -FU$OUT"
VERBOSE=1

for arg in "$@"; do
  [[ "$arg" == "--verbose" ]] && VERBOSE=1
done

run() {
  if [[ $VERBOSE -eq 1 ]]; then
    "$@"
  else
    "$@" > /dev/null 2>&1
  fi
}

echo "=== jpmath test suite ==="
echo "Shared units : $SHARED"
echo "Build output : $OUT"
echo

mkdir -p "$OUT"

# ---------------------------------------------------------------------------
# 1. Compile shared units
# ---------------------------------------------------------------------------
echo "Compiling shared units..."

UNITS=(
  jpmtypes
  jpmanneal jpmarith jpmbessel jpmchebyshev jpmcomplex
  jpmcontinued jpmderivative jpmdiffeq jpmgeometry jpmintegration
  jpmlstsqr jpmmatrices jpmoptimize jpmpolynomials jpmroots
  jpmsignal jpmsimplex jpmsort jpmspecial jpmstats
  jpmappointment jpminterpolation
)

ERRORS=0
for u in "${UNITS[@]}"; do
  SRC="$SHARED/$u.pas"
  if [[ ! -f "$SRC" ]]; then
    echo "  SKIP (not found): $u"
    continue
  fi
  if fpc $FPC_FLAGS "$SRC" > "$OUT/${u}_compile.log" 2>&1; then
    [[ $VERBOSE -eq 1 ]] && echo "  OK: $u"
  else
    echo "  FAIL (compile): $u"
    cat "$OUT/${u}_compile.log"
    ERRORS=$((ERRORS + 1))
  fi
done

if [[ $ERRORS -gt 0 ]]; then
  echo
  echo "ERROR: $ERRORS unit(s) failed to compile. Aborting."
  exit 1
fi
echo "  All shared units compiled."
echo

# ---------------------------------------------------------------------------
# 2. Compile the console self-test runner
#    Source: _unittests/run_self_tests_console.pas
# ---------------------------------------------------------------------------
echo "Compiling test runner..."
if ! fpc $FPC_FLAGS -o"$RUNNER_BIN" "$RUNNER_SRC" > "$OUT/runner_compile.log" 2>&1; then
  echo "  FAIL (compile test runner)"
  cat "$OUT/runner_compile.log"
  exit 1
fi
echo "  Test runner compiled."
echo

# ---------------------------------------------------------------------------
# 3. Run the tests
# ---------------------------------------------------------------------------
echo "Running self_tests..."
echo
"$RUNNER_BIN"
STATUS=$?
echo
if [[ $STATUS -eq 0 ]]; then
  echo "=== All tests passed ==="
else
  echo "=== Some tests FAILED (exit code $STATUS) ==="
fi
exit $STATUS
