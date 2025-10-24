#!/bin/bash

# Stokes Equation Convergence Test Script
# Tests spatial convergence for projection method Stokes solver

EXECUTABLE="$1"

if [ -z "$EXECUTABLE" ]; then
    echo "Usage: $0 <path_to_executable>"
    exit 1
fi

if [ ! -f "$EXECUTABLE" ]; then
    echo "Error: Executable $EXECUTABLE not found"
    exit 1
fi

echo "=== Stokes Equation Convergence Test (Projection Method) ==="
echo "Testing spatial convergence with staggered grid"
echo "Executable: $EXECUTABLE"
echo ""
echo "Grid       h            U_Error         V_Error         P_Error         U_Rate     V_Rate     P_Rate"
echo "------------------------------------------------------------------------------------------------------------"

# Grid sizes to test
GRIDS=(16 32 64 128)

# Arrays to store results
declare -a errors_u
declare -a errors_v
declare -a errors_p
declare -a hs

# Run simulations
for N in "${GRIDS[@]}"; do
    h=$(echo "scale=8; sqrt(2.0) / $N" | bc -l)
    hs+=($h)
    
    output=$(mpirun -np 4 "$EXECUTABLE" -nx $N -ny $N -convergence_test 2>&1)
    
    # Parse CONVERGENCE line: N h u_err v_err p_err
    conv_line=$(echo "$output" | grep "^CONVERGENCE:")
    
    if [ -z "$conv_line" ]; then
        echo "Error: Could not parse output for N=$N"
        echo "Output was:"
        echo "$output"
        exit 1
    fi
    
    # Extract errors
    u_err=$(echo "$conv_line" | awk '{print $4}')
    v_err=$(echo "$conv_line" | awk '{print $5}')
    p_err=$(echo "$conv_line" | awk '{print $6}')
    
    errors_u+=($u_err)
    errors_v+=($v_err)
    errors_p+=($p_err)
done

# Compute and print convergence rates
for i in "${!GRIDS[@]}"; do
    N=${GRIDS[$i]}
    h=${hs[$i]}
    u_err=${errors_u[$i]}
    v_err=${errors_v[$i]}
    p_err=${errors_p[$i]}
    
    # Format h
    h_fmt=$(printf "%.8f" $h)
    
    if [ $i -eq 0 ]; then
        # First grid, no rate
        printf "%-10s %-12s %-15s %-15s %-15s %-10s %-10s %-10s\n" \
               "${N}x${N}" "$h_fmt" "$u_err" "$v_err" "$p_err" "---" "---" "---"
    else
        # Compute convergence rates
        u_err_prev=${errors_u[$((i-1))]}
        v_err_prev=${errors_v[$((i-1))]}
        p_err_prev=${errors_p[$((i-1))]}
        
        u_rate=$(echo "scale=4; l($u_err_prev / $u_err) / l(2.0)" | bc -l)
        v_rate=$(echo "scale=4; l($v_err_prev / $v_err) / l(2.0)" | bc -l)
        p_rate=$(echo "scale=4; l($p_err_prev / $p_err) / l(2.0)" | bc -l)
        
        printf "%-10s %-12s %-15s %-15s %-15s %-10s %-10s %-10s\n" \
               "${N}x${N}" "$h_fmt" "$u_err" "$v_err" "$p_err" "$u_rate" "$v_rate" "$p_rate"
    fi
done

echo ""
echo "Expected convergence rates:"
echo "  Velocity (u, v): ~2.0 (second-order accurate)"
echo "  Pressure (p):    ~2.0 (second-order accurate)"
echo "Formula: rate = log(error_h / error_{h/2}) / log(2)"
echo "Note: Projection method with staggered grid (MAC grid)"
