# 栅格地理计算并行算子 （PaRGO）

# 1. 简介

Parallel Raster-based Geocomputation Operators

# 2. 编译安装

## 2.1. 编译环境

PaRGO依赖GDAL，根据不同并行版本依赖MPI（MSMPI，OpenMPI，MPICH，或Intel MPI）、OpenMP、CUDA。

+ Windows 10 with Visual Studio 2010/2013/2015, MSMPI-v8.1, GDAL-1.11.4
+ Windows 10 with msys2/mingw64 (GCC-9.1.0), MSMPI-v8.1, GDAL-3.0
+ CentOS 6.2 (cluster) with GCC-4.8.4, MPICH-3.1.4, GDAL-1.9.0
+ Red Hat Server 6.2 (cluster) with ICC-12.1.0, Intel MPI 4.0, GDAL-1.11.5
+ macOS 10.14.5 with Clang-10.0 with Xcode, OpenMPI-4.0.1, GDAL-2.4.2 (brew installed)


+ Windows 10-64bit with Visual Studio 2013

  ```shell
  cmake -G "Visual Studio 12 2013 Win64" C:/z_code/DTA/PaRGO ..
  ```
  By default, the install directory is /path/to/source/bin,
  and can also be specified by add `-DINSTALL_PREFIX` argument, 
  e.g., `-DINSTALL_PREFIX=D:/compile/bin/pargo`
  ​

