#!/usr/bin/env bash
# Runs all finite verifications (C1)--(C9), stopping on the first failure.
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p .build logs
log=logs/verification.log
: > "$log"
run() {
  label="$1"; shift
  echo "==== $label ====" | tee -a "$log"
  "$@" 2>&1 | tee -a "$log"
}
compiler="${CXX:-g++}"
"$compiler" -std=c++17 -O2 -Wall -Wextra verify_R_finite.cpp -o .build/verify_R_finite
run 'C1 finite quintinomial enumeration' ./.build/verify_R_finite
run 'C2, C3 affine collisions and integral lines (primary)' python3 verify_collision_lattices.py
run 'C2, C3 affine collisions and integral lines (independent exterior minors)' python3 verify_collision_lattices_alt.py
run 'C4 fifth-power supports, independent algorithms and CSV' python3 verify_fifth_norm.py
run 'C5, C6 binary prime exceptions and Phi_13 over F4' python3 verify_prime_exceptions.py
run 'C7 binary mixing obstructions' python3 verify_binary_mixing.py
run 'C8 quaternary mixing obstructions' python3 verify_quaternary_mixing.py
run 'C9 small cores' python3 verify_binary_cores.py
echo 'ALL CHECKS C1--C9 PASSED' | tee -a "$log"
