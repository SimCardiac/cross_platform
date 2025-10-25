#!/bin/bash
# ==============================================================================
# Quick Parallel Efficiency Test for Unified PDE Solver
# Optimized version with reduced iterations for CI/CD
# ==============================================================================

set -e  # Exit on error

if [ $# -lt 1 ]; then
    echo "Usage: $0 <path_to_unified_pde_solver_2D> [problem_type] [nx] [ny]"
    exit 1
fi

SOLVER=$1
PROBLEM_TYPE=${2:-poisson}
NX=${3:-64}
NY=${4:-64}

echo "========================================"
echo "Quick Parallel Efficiency Test"
echo "========================================"
echo "Solver: $SOLVER"
echo "Problem type: $PROBLEM_TYPE"
echo "Grid size: ${NX}x${NY}"
echo "========================================"

# Check if solver exists
if [ ! -f "$SOLVER" ]; then
    echo "Error: Solver not found at $SOLVER"
    exit 1
fi

# Check if mpirun is available
if ! command -v mpirun &> /dev/null; then
    echo "Warning: mpirun not found. Running serial test only."
    # Run a simple serial test
    if [ "$PROBLEM_TYPE" == "heat" ]; then
        $SOLVER -problem_type heat -nx $NX -ny $NY -dt 0.001 -T 0.001 > /dev/null 2>&1
    else
        $SOLVER -problem_type poisson -nx $NX -ny $NY > /dev/null 2>&1
    fi
    echo "✓ PASS: Serial execution successful"
    exit 0
fi

# Array to store results
declare -a NUM_PROCS=(1 2 4)
declare -a TIMES=()

echo ""
echo "Running tests..."
echo "========================================"

# Run tests with different numbers of processes
for NP in "${NUM_PROCS[@]}"; do
    echo "Testing with $NP process(es)..."
    
    START_TIME=$(date +%s.%N)
    if [ "$PROBLEM_TYPE" == "heat" ]; then
        # Use very short simulation for heat equation
        mpirun -np $NP $SOLVER -problem_type heat -nx $NX -ny $NY -dt 0.001 -T 0.005 > /dev/null 2>&1
    else
        mpirun -np $NP $SOLVER -problem_type poisson -nx $NX -ny $NY > /dev/null 2>&1
    fi
    END_TIME=$(date +%s.%N)
    
    ELAPSED=$(echo "$END_TIME - $START_TIME" | bc)
    TIMES+=($ELAPSED)
    
    echo "  Time: ${ELAPSED}s"
done

echo ""
echo "========================================"
echo "Results Summary"
echo "========================================"
printf "%-12s %-12s %-12s %-12s\n" "Processes" "Time(s)" "Speedup" "Efficiency(%)"
echo "----------------------------------------"

# Calculate speedup and efficiency
SERIAL_TIME=${TIMES[0]}
for i in "${!NUM_PROCS[@]}"; do
    NP=${NUM_PROCS[$i]}
    TIME=${TIMES[$i]}
    
    SPEEDUP=$(echo "scale=3; $SERIAL_TIME / $TIME" | bc)
    EFFICIENCY=$(echo "scale=2; ($SPEEDUP / $NP) * 100" | bc)
    
    printf "%-12d %-12.3f %-12.3f %-12.2f\n" $NP $TIME $SPEEDUP $EFFICIENCY
done

echo "========================================"
echo ""

# Check if efficiency is reasonable
EFFICIENCY_2=$(echo "scale=2; (${TIMES[0]} / ${TIMES[1]} / 2) * 100" | bc)
EFFICIENCY_4=$(echo "scale=2; (${TIMES[0]} / ${TIMES[2]} / 4) * 100" | bc)

echo "Parallel efficiency check:"
echo "  2 processes: ${EFFICIENCY_2}% (expected >60%)"
echo "  4 processes: ${EFFICIENCY_4}% (expected >35%)"

# More relaxed thresholds for quick testing
PASS_2=$(echo "$EFFICIENCY_2 > 60" | bc)
PASS_4=$(echo "$EFFICIENCY_4 > 35" | bc)

if [ "$PASS_2" -eq 1 ] && [ "$PASS_4" -eq 1 ]; then
    echo ""
    echo "✓ PASS: Parallel efficiency is acceptable"
    exit 0
else
    echo ""
    echo "✗ FAIL: Parallel efficiency is below expected thresholds"
    echo "Note: This may be due to small problem size or system load"
    exit 1
fi
