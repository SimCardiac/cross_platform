#!/bin/bash

# Vertex-Centered Poisson Convergence Test Script
# Usage: poisson_vertex_convergence.sh [path_to_poisson_vertex_2D_executable]

# Get executable path from argument or use default
POISSON_EXEC="${1:-./poisson_vertex_2D}"

# Check if executable exists
if [ ! -f "$POISSON_EXEC" ]; then
  echo "Error: Executable not found: $POISSON_EXEC"
  exit 1
fi

echo "=== Vertex-Centered Poisson Convergence Test ==="
echo "Testing spatial convergence with DMDA vertex-centered grid"
echo "Executable: $POISSON_EXEC"
echo ""
printf "%-10s %-12s %-15s %-15s %-12s\n" "Grid" "h" "L2_Error" "Conv_Rate" "Iterations"
echo "-----------------------------------------------------------------------"

declare -a errors
declare -a sizes=(16 32 64 128)

for i in "${!sizes[@]}"; do
  n=${sizes[$i]}
  
  # Run the test and capture output
  output=$(mpirun -np 4 "$POISSON_EXEC" \
    -nx $n -ny $n \
    -poisson_check_error \
    -convergence_test \
    -ksp_type gmres \
    -pc_type bjacobi \
    -sub_pc_type ilu \
    -ksp_rtol 1e-12 \
    2>&1)
  
  # Parse the CONVERGENCE line: CONVERGENCE: N h error iterations
  conv_line=$(echo "$output" | grep "CONVERGENCE:")
  if [ -n "$conv_line" ]; then
    h=$(echo "$conv_line" | awk '{print $3}')
    error=$(echo "$conv_line" | awk '{print $4}')
    iters=$(echo "$conv_line" | awk '{print $5}')
    errors[$i]=$error
    
    # Calculate convergence rate if not first iteration
    if [ $i -gt 0 ]; then
      prev_error=${errors[$((i-1))]}
      rate=$(echo "scale=4; l($prev_error/$error)/l(2)" | bc -l)
      printf "%-10s %-12s %-15s %-15s %-12s\n" "${n}x${n}" "$h" "$error" "$rate" "$iters"
    else
      printf "%-10s %-12s %-15s %-15s %-12s\n" "${n}x${n}" "$h" "$error" "---" "$iters"
    fi
  else
    echo "Error: Could not parse output for n=$n"
  fi
done

echo ""
echo "Expected spatial convergence rate: ~2.0 (second-order accurate)"
echo "Formula: rate = log(error_h / error_{h/2}) / log(2)"
