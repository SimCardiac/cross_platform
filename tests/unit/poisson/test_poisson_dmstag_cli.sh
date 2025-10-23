#!/bin/bash
# DMStag CLI smoke test script
# Usage: test_poisson_dmstag_cli.sh <path-to-poisson_dmstag_tests_cli>

set -e

if [ $# -ne 1 ]; then
    echo "Usage: $0 <path-to-poisson_dmstag_tests_cli>"
    exit 1
fi

RUNNER="$1"

if [ ! -x "$RUNNER" ]; then
    echo "Error: $RUNNER is not executable or does not exist"
    exit 1
fi

echo "=== DMStag CLI Smoke Test ==="
echo "Testing: $RUNNER"
echo

# Test 1: Single pair (poly2)
echo "Test 1: Single pair poly2 on 8x8 grid"
$RUNNER --grid 8 8 --pair poly2
if [ $? -ne 0 ]; then
    echo "FAIL: poly2 test failed"
    exit 1
fi
echo "PASS"
echo

# Test 2: Single pair (sinpi)
echo "Test 2: Single pair sinpi on 16x16 grid"
$RUNNER --grid 16 16 --pair sinpi
if [ $? -ne 0 ]; then
    echo "FAIL: sinpi test failed"
    exit 1
fi
echo "PASS"
echo

# Test 3: Multiple pairs
echo "Test 3: Multiple pairs on 8x8 grid"
$RUNNER --grid 8 8 --pair poly2 --pair sinpi
if [ $? -ne 0 ]; then
    echo "FAIL: multiple pairs test failed"
    exit 1
fi
echo "PASS"
echo

echo "=== All DMStag smoke tests passed ==="
exit 0
