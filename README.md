# Parallel Raster-based Geocomputation Operators (PaRGO) Version 2
## 1. introduction

Coming soon...

Usage on Windows is listed below.

## 2. Prerequisites
**See user manual in /doc for more details**
### 2.1. MSVC, CMake
### 2.2. GDAL
### 2.3. MS-MPI
### 2.4. Depending on specific computing platforms, PaRGO has dependency on MPI, MPI+OpenMP, or CUDA.


  ​
## 3. Installation

```shell
git clone .......
cd PaRGO
mkdir build
cd build
cmake -G "Visual Studio 10 2010 Win64" ../PaRGO -DUSE_MPI_DEBUGGER=1 ..
#(open VS2010, run ALL_BUILD and ALL_ISNTALL)

or:
cmake -G "Visual Studio 14 2015 Win64" ../PaRGO -DUSE_MPI_DEBUGGER=1 ..
#(open VS2015, run ALL_BUILD and ALL_ISNTALL)


```
  By default, the install directory is /path/to/source/bin,
  and can also be specified by add `-DINSTALL_PREFIX` argument, 
  e.g., `-DINSTALL_PREFIX=D:/compile/bin/pargo`

## 4. Run

examples (replace `C:\lib\mpi\v8.1\Bin\mpiexec.exe` to your installed path):

```shell
MultiScaleLE:
C:\lib\mpi\v8.1\Bin\mpiexec.exe -n 4 D:\src\PaRGO\vs2010\build\apps\morphology\Debug\multiScaleLE.exe D:\arcgis-data\pargo\ywzdem5m.tif D:\arcgis-data\pargo\le\moore20.nbr D:\arcgis-data\pargo\le\output.tif 1

Slope:
C:\lib\mpi\v8.1\Bin\mpiexec.exe -n 4 D:\src\PaRGO\vs2010\build_gdal1\apps\morphology\Debug\slope.exe D:\arcgis-data\pargo\fcm\dem_nenjiang_10.tif D:\arcgis-data\pargo\moore.nbr D:\arcgis-data\pargo\fcm\slope_nenjiang_10.tif

IDW:
C:\lib\mpi\8.1\Bin\mpiexec.exe -n 8 C:\src\PaRGO\vs2010\build\apps\spatial\Release\idw.exe -sample D:\data-arcgis\pargo\idw\sc_p_3w.shp -mask D:\data-arcgis\pargo\idw\idw_bounding_mask_10m.tif -dataNbr D:\data-arcgis\pargo\neigh.nbr -computeNbr D:\data-arcgis\pargo\neigh.nbr -out D:\data-arcgis\pargo\idw\out\temp.tif -resolution 10 -fieldIndex 3 -idwExp 2 -searchPointNum 12 -searchRange 0 -blockSize 5000 -dcmp space

FCM:
C:\lib\mpi\8.1\Bin\mpiexec.exe -n 8 C:\src\PaRGO\vs2010\build\apps\spatial\Release\fcm.exe -inputs D:\data-arcgis\pargo\fcm\normalized\twi_1.tif,D:\data-arcgis\pargo\fcm\normalized\plan_1.tif,D:\data-arcgis\pargo\fcm\normalized\prof_1.tif,D:\data-arcgis\pargo\fcm\normalized\slp_1.tif, -dataNbr D:\data-arcgis\pargo\moore.nbr -computeNbr D:\data-arcgis\pargo\moore.nbr -out D:\data-arcgis\pargo\fcm\out\temp.tif -clusterNum 5 -maxIter 5 -tolerance 0.00001 -weight 2 -dcmp space
```

## 5. Develop based on PaRGO

  1. create program files (`.cpp` or/and `.h`) in a proper folder (e.g., `PaRGO\apps\spatial\myAlgorithm`), referring to the path of existing programs.
  2. modify the `CMkaeLists.txt` of the created folder (e.g., `PaRGO\apps\spatial\CMkaeLists.txt`)
    > FILE(GLOB MYALGOFILES ./myAlgorithm/*.cpp)
      SET(MYALGOFILES ${MYALGOFILES} ${GPRO_SRCS})
      ADD_EXECUTABLE(myalgorithm ${MYALGOFILES})
      SET(SPATIAL_TARGETS idw
                    fcm
                    myalgorithm
                    )

## 5. Tested platforms
+ Windows 10 with Visual Studio 2010/2013/2015, MSMPI-v8.1, GDAL-1.11.4
+ Windows 10 with msys2/mingw64 (GCC-9.1.0), MSMPI-v8.1, GDAL-3.0
+ CentOS 6.2 (cluster) with GCC-4.8.4, MPICH-3.1.4, GDAL-1.9.0
+ Red Hat Server 6.2 (cluster) with ICC-12.1.0, Intel MPI 4.0, GDAL-1.11.5
+ macOS 10.14.5 with Clang-10.0 with Xcode, OpenMPI-4.0.1, GDAL-2.4.2 (brew installed)
+ Windows 10-64bit with Visual Studio 2013