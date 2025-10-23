#!/bin/bash
# Smoke test for Poisson CLI test runner
# Verifies basic execution with one manufactured solution pair

set -e

RUNNER="$1"
if [ -z "$RUNNER" ]; then
    echo "Usage: $0 <path_to_poisson_tests_cli>"
    exit 1
fi

echo "=== Poisson CLI Smoke Test ==="
echo "Runner: $RUNNER"

# Test 1: Single pair poly2
echo "Test 1: poly2 on 8x8 grid"
$RUNNER --grid 8 8 --pair poly2

# Test 2: Single pair sinpi
echo "Test 2: sinpi on 16x16 grid"
$RUNNER --grid 16 16 --pair sinpi

# Test 3: Multiple pairs
echo "Test 3: Multiple pairs"
$RUNNER --grid 8 8 --pair poly2 --pair sinpi

echo "=== All smoke tests passed ==="
exit 0
