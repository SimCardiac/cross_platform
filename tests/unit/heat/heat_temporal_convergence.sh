#!/bin/bash

# Heat Equation Temporal Convergence Test Script
# Tests convergence order with respect to time step dt
# Usage: heat_temporal_convergence.sh [path_to_heat_DMStag_2D_executable]

# Get executable path from argument or use default
HEAT_EXEC="${1:-./heat_DMStag_2D}"

# Check if executable exists
if [ ! -f "$HEAT_EXEC" ]; then
  echo "Error: Executable not found: $HEAT_EXEC"
  exit 1
fi

echo "=== Heat Equation Temporal Convergence Test ==="
echo "Testing temporal (time) convergence with implicit Euler"
echo "Executable: $HEAT_EXEC"
echo ""
echo "Fixed grid: 128x128 (spatial error negligible)"
echo "Variable time step: dt is halved each test"
echo ""
printf "%-15s %-15s %-15s %-12s\n" "dt" "L2_Error" "Conv_Rate" "Steps"
echo "-----------------------------------------------------------------------"

declare -a errors
declare -a dts=(0.01 0.005 0.0025 0.00125 0.000625)

# Fixed grid size (fine enough so spatial error is negligible)
n=128
T_final=0.1
for i in "${!dts[@]}"; do
  dt=${dts[$i]}
  
  # Calculate number of time steps
  steps=$(echo "scale=0; $T_final/$dt" | bc)
  
  # Run the test and capture output
  output=$(mpirun -np 4 "$HEAT_EXEC" -nx $n -ny $n -dt $dt -T $T_final -heat_check_error 2>&1)
  
  # Parse the L2 error line
  error_line=$(echo "$output" | grep "||u(T) - u_exact(T)||_L2")
  if [ -n "$error_line" ]; then
    error=$(echo "$error_line" | awk '{print $NF}')
    errors[$i]=$error
    
    # Calculate convergence rate if not first iteration
    if [ $i -gt 0 ]; then
      prev_error=${errors[$((i-1))]}
      # For temporal convergence, dt is halved, so rate = log(error_old/error_new)/log(2)
      rate=$(echo "scale=4; l($prev_error/$error)/l(2)" | bc -l)
      printf "%-15s %-15s %-15s %-12s\n" "$dt" "$error" "$rate" "$steps"
    else
      printf "%-15s %-15s %-15s %-12s\n" "$dt" "$error" "---" "$steps"
    fi
  else
    echo "Error: Could not parse output for dt=$dt"
    echo "Output was:"
    echo "$output"
  fi
done

echo ""
echo "Expected temporal convergence rate: ~1.0 (first-order in time)"
echo "Formula: rate = log(error_dt / error_{dt/2}) / log(2)"
echo "Note: Implicit Euler is first-order accurate in time"
echo "      Grid is fixed at 128x128 to minimize spatial error"
