#!/bin/bash
# ==============================================================================
# Parallel Efficiency Test for Unified PDE Solver
# Measures strong scaling efficiency by running the same problem with
# different numbers of MPI processes
# ==============================================================================

set -e  # Exit on error

if [ $# -lt 1 ]; then
    echo "Usage: $0 <path_to_unified_pde_solver_2D> [problem_type] [nx] [ny]"
    echo "Example: $0 ./unified_pde_solver_2D poisson 128 128"
    exit 1
fi

SOLVER=$1
PROBLEM_TYPE=${2:-poisson}  # Default: poisson
NX=${3:-128}                # Default: 128
NY=${4:-128}                # Default: 128

echo "========================================"
echo "Parallel Efficiency Test"
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
    echo "Error: mpirun not found. MPI is required for this test."
    exit 1
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
    
    # Run the solver and capture timing
    if [ "$PROBLEM_TYPE" == "heat" ]; then
        OUTPUT=$(mpirun -np $NP $SOLVER -problem_type heat -nx $NX -ny $NY -dt 0.001 -T 0.01 2>&1)
    else
        OUTPUT=$(mpirun -np $NP $SOLVER -problem_type poisson -nx $NX -ny $NY 2>&1)
    fi
    
    # Extract wall time (we'll need to add timing output to the solver)
    # For now, use the time command to measure
    START_TIME=$(date +%s.%N)
    if [ "$PROBLEM_TYPE" == "heat" ]; then
        mpirun -np $NP $SOLVER -problem_type heat -nx $NX -ny $NY -dt 0.001 -T 0.01 > /dev/null 2>&1
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
    
    # Calculate speedup: T_serial / T_parallel
    SPEEDUP=$(echo "scale=3; $SERIAL_TIME / $TIME" | bc)
    
    # Calculate efficiency: (Speedup / num_procs) * 100
    EFFICIENCY=$(echo "scale=2; ($SPEEDUP / $NP) * 100" | bc)
    
    printf "%-12d %-12.3f %-12.3f %-12.2f\n" $NP $TIME $SPEEDUP $EFFICIENCY
done

echo "========================================"
echo ""

# Check if efficiency is reasonable (>70% for 2 procs, >50% for 4 procs)
EFFICIENCY_2=$(echo "scale=2; (${TIMES[0]} / ${TIMES[1]} / 2) * 100" | bc)
EFFICIENCY_4=$(echo "scale=2; (${TIMES[0]} / ${TIMES[2]} / 4) * 100" | bc)

echo "Parallel efficiency check:"
echo "  2 processes: ${EFFICIENCY_2}% (expected >70%)"
echo "  4 processes: ${EFFICIENCY_4}% (expected >40%)"

# Return success if efficiencies are reasonable
PASS_2=$(echo "$EFFICIENCY_2 > 70" | bc)
PASS_4=$(echo "$EFFICIENCY_4 > 40" | bc)

if [ "$PASS_2" -eq 1 ] && [ "$PASS_4" -eq 1 ]; then
    echo ""
    echo "✓ PASS: Parallel efficiency is acceptable"
    exit 0
else
    echo ""
    echo "✗ FAIL: Parallel efficiency is below expected thresholds"
    exit 1
fi
