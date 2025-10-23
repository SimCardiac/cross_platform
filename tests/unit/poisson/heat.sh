# cd /Users/pengfei/Documents/GitHub/cross_platform/build && cat > spatial_convergence_test.sh << 'EOF'
#!/bin/bash
echo "=== Spatial Convergence Test (Fixed dt) ==="
echo "Using fixed small dt=0.0001 to isolate spatial error"
echo ""
printf "%-10s %-12s %-15s %-15s\n" "Grid" "h" "L2_Error" "Conv_Rate"
echo "------------------------------------------------------------"

declare -a errors
declare -a sizes=(16 32 64 128)
dt=0.0001

for i in "${!sizes[@]}"; do
  n=${sizes[$i]}
  h=$(echo "scale=8; 1.0/$n" | bc)
  
  # Run the test and capture error
  output=$(mpirun -np 4 ./ex_heat_stagger -nx $n -ny $n -dt $dt -T 0.01 -heat_check_error 2>&1)
  error=$(echo "$output" | grep "L2" | awk '{print $NF}')
  errors[$i]=$error
  
  # Calculate convergence rate if not first iteration
  if [ $i -gt 0 ]; then
    prev_error=${errors[$((i-1))]}
    rate=$(echo "scale=4; l($prev_error/$error)/l(2)" | bc -l)
    printf "%-10s %-12s %-15s %-15s\n" "${n}x${n}" "$h" "$error" "$rate"
  else
    printf "%-10s %-12s %-15s %-15s\n" "${n}x${n}" "$h" "$error" "---"
  fi
done

echo ""
echo "Expected spatial convergence rate: ~2.0"
# EOF
# chmod +x spatial_convergence_test.sh && bash spatial_convergence_test.sh