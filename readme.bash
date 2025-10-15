




gcc      14.2.0 
openmpi  5.0.8
petsc    3.23.4
openblas 0.3.27
kokkos   4.6.01

目前这些步骤可能有些多余，先通过 homebrew安装 gcc，再通过spack编译一个版本的gcc，后面 spack 编译软件都用这个gcc来编译。

brew install gcc
spack compiler find
spack install gcc@14.2.0 %gcc@15.2.0

然后用 gcc 14.2.0 编译其他软件。


spack load gcc@14.2.0  %gcc@15.2.0

spack env create cross && spack env activate cross
spack add openmpi@5.0.8 %gcc@14.2.0
spack add openblas@0.3.27 %gcc@14.2.0
spack add kokkos@4.6.01 %gcc@14.2.0
spack add petsc@3.23.4 %gcc@14.2.0
spack concretize -f