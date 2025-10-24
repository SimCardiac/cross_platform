#!/bin/bash

# Heat Equation Convergence Test Script (Edge-Centered)
# Usage: heat_edge_convergence.sh [path_to_heat_DMStag_edge_2D_executable]

# Get executable path from argument or use default
HEAT_EXEC="${1:-./heat_DMStag_edge_2D}"

# Check if executable exists
if [ ! -f "$HEAT_EXEC" ]; then
  echo "Error: Executable not found: $HEAT_EXEC"
  exit 1
fi

echo "=== Heat Equation Convergence Test (Edge-Centered) ==="
echo "Testing spatial convergence with edge-centered DMStag"
echo "Executable: $HEAT_EXEC"
echo ""
printf "%-10s %-12s %-15s %-15s %-12s\n" "Grid" "h" "L2_Error" "Conv_Rate" "T_final"
echo "-----------------------------------------------------------------------"

declare -a errors
declare -a sizes=(16 32 64 128)

for i in "${!sizes[@]}"; do
  n=${sizes[$i]}
  
  # For heat equation, dt should scale with h^2 for stability
  h=$(echo "scale=8; 1.0/$n" | bc)
  dt=$(echo "scale=8; 0.2*$h*$h" | bc)
  
  # Run the test and capture output
  output=$(mpirun -np 4 "$HEAT_EXEC" -nx $n -ny $n -dt $dt -T 0.01 -heat_check_error 2>&1)
  
  # Parse the L2 error line
  error_line=$(echo "$output" | grep "||u(T) - u_exact(T)||_L2")
  if [ -n "$error_line" ]; then
    error=$(echo "$error_line" | awk '{print $NF}')
    errors[$i]=$error
    
    # Calculate convergence rate if not first iteration
    if [ $i -gt 0 ]; then
      prev_error=${errors[$((i-1))]}
      rate=$(echo "scale=4; l($prev_error/$error)/l(2)" | bc -l)
      printf "%-10s %-12s %-15s %-15s %-12s\n" "${n}x${n}" "$h" "$error" "$rate" "0.01"
    else
      printf "%-10s %-12s %-15s %-15s %-12s\n" "${n}x${n}" "$h" "$error" "---" "0.01"
    fi
  else
    echo "Error: Could not parse output for n=$n"
    echo "Output was:"
    echo "$output"
  fi
done

echo ""
echo "Expected spatial convergence rate: ~2.0 (second-order accurate)"
echo "Formula: rate = log(error_h / error_{h/2}) / log(2)"
echo "Note: dt scales as h^2 to maintain stability (dt = 0.2*h^2)"
echo "      Edge-centered DOF on DOWN edges (horizontal cell boundaries)"
