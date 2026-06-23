#!/bin/bash
EXEC="${1:-./stokes_projection_2D}"
[ ! -f "$EXEC" ] && echo "Error: not found" && exit 1
echo "=== Stokes Projection Convergence ==="
printf "%-8s %-12s %-14s %-10s\n" "Grid" "h" "L2_u" "Rate"
echo "--------------------------------------------"
declare -a eu
for n in 16 32 64 128; do
  out=$("$EXEC" -nx $n -ny $n -T 0.01 -stokes_check_error -convergence_test 2>&1)
  h=$(echo "$out" | awk '/CONVERGENCE/{print $2}')
  e=$(echo "$out" | awk '/CONVERGENCE/{print $3}')
  eu+=("$e"); i=$((${#eu[@]}-1))
  if [ $i -gt 0 ]; then
    r=$(echo "scale=3; l(${eu[$((i-1))]}/$e)/l(2)" | bc -l)
    printf "%-8s %-12s %-14s %-10s\n" "${n}x${n}" "$h" "$e" "$r"
  else
    printf "%-8s %-12s %-14s %-10s\n" "${n}x${n}" "$h" "$e" "---"
  fi
done
echo "Expected: ~2.0, dt = h^2"
