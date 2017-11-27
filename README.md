# 栅格地理计算并行算子 （PaRGO）

# 1. 简介

Parallel Raster-based Geocomputation Operators

# 2. 编译安装

## 2.1. 编译环境

PaRGO依赖GDAL，根据不同并行版本依赖MPI（MSMPI，OpenMPI，MPICH，或Intel MPI）、CUDA。

+ Windows 10 with Visual Studio 2013, MSMPI-v8, GDAL-1.11.4
+ Windows 10 with mingw64 (GCC-4.9.3), MSMPI-v8, GDAL-1.11.5
+ CentOS 6.2 (cluster) with GCC-4.8.4, MPICH-3.1.4, GDAL-1.9.0
+ Red Hat Server 6.2 (cluster) with ICC-12.1.0, Intel MPI 4.0, GDAL-1.11.5
+ macOS 10.12.6 with Clang-8.0 (or GCC-4.9.3), OpenMPI-1.10.4, GDAL-1.11.4 (Framework)



+ Windows 10-64bit with Visual Studio 2012

  ```shell
  cmake -G "Visual Studio 12 2013 Win64" C:/z_code/DTA/PaRGO -DMPI=1 -DOPENMP=1 -DINSTALL_PREFIX=D:/compile/bin/pargo
  ```

  ​

