#!/bin/bash
# Staggered (Mixed) Heat Convergence Test (Item 6)
HEAT_EXEC="${1:-./heat_staggered_2D}"
if [ ! -f "$HEAT_EXEC" ]; then echo "Error: $HEAT_EXEC not found"; exit 1; fi
echo "=== Staggered (Mixed) Heat Convergence Test ==="
echo "Executable: $HEAT_EXEC"
echo ""
printf "%-10s %-12s %-15s %-15s %-12s\n" "Grid" "h" "L2_Error" "Conv_Rate" "Steps"
echo "-----------------------------------------------------------------------"
declare -a errors
declare -a sizes=(16 32 64 128)
for i in "${!sizes[@]}"; do
  n=${sizes[$i]}
  output=$(mpirun -np 4 "$HEAT_EXEC" -nx $n -ny $n -T 0.05 -heat_check_error -convergence_test 2>&1)
  conv_line=$(echo "$output" | grep "CONVERGENCE:")
  if [ -n "$conv_line" ]; then
    h=$(echo "$conv_line" | awk '{print $3}')
    error=$(echo "$conv_line" | awk '{print $4}')
    steps=$(echo "$conv_line" | awk '{print $5}')
    errors[$i]=$error
    if [ $i -gt 0 ]; then
      prev_error=${errors[$((i-1))]}
      rate=$(echo "scale=4; l($prev_error/$error)/l(2)" | bc -l)
      printf "%-10s %-12s %-15s %-15s %-12s\n" "${n}x${n}" "$h" "$error" "$rate" "$steps"
    else
      printf "%-10s %-12s %-15s %-15s %-12s\n" "${n}x${n}" "$h" "$error" "---" "$steps"
    fi
  fi
done
echo ""
echo "Expected spatial convergence rate: ~2.0 (second-order accurate)"
echo "dt = h^2 (coupled space-time refinement)"
