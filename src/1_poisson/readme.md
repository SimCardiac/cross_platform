

## cell-centered 有限差分格式求解 poisson 方程

```bash
source .bashrc && cmake -S . -B build && cmake --build build --target poisson_DMStag_2D 
./tests/unit/1_poisson/poisson_convergence.sh build/poisson_DMStag_2D
```

```
=== Poisson Equation Convergence Test ===
Testing spatial convergence with element-centered DMStag
Executable: build/poisson_DMStag_2D

Grid       h            L2_Error        Conv_Rate       Iterations  
-----------------------------------------------------------------------
16x16      0.0625       0.00160948      ---             629         
32x32      0.03125      0.000401789     2.0021          2602        
64x64      0.015625     0.000100411     2.0005          10707       
128x128    0.0078125    2.51004e-05     2.0001          40001       
256x256    0.00390625   0.00218438      -6.4552         40001       

Expected spatial convergence rate: ~2.0 (second-order accurate)
Formula: rate = log(error_h / error_{h/2}) / log(2)
```