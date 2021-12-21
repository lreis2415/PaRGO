Parallel Raster-based Geocomputation Operators (PaRGO) Version 2 - Introduction
==========

[toc]

# 1. Introduction

PaRGO is a C++ parallel programming framework for raster-based geocomputation, featuring support for easy implementation of various geocomputation algorithms in a serial style but with parallel performance on different parallel platforms. In PaRGO Version 2, we have improved its load-balancing performance using the idea of the spatial computational domain. The main features of PaRGO (or PaRGO V2) are:

1. support for one parallel program running on different parallel platforms/models: MPI and MPI+OpenMP for CPU in Beowulf and SMP clusters, and CUDA for GPU.
2. support for implementation of raster-based geocomputation algorithms with different characteristics: local, focal, zonal, and global algorithms.
3. enables flexible, practical, and effective load-balancing strategy in multiple modes: the intensity ratio mode, the estimate function mode, and the preliminary experiment mode, for uniform or nonuniform data and computation spatial distribution.

Usage on Windows is documented below.

# 2. Installation for Windows

## 2.1 MSVC

PaRGO is a C++ project. Please install `Visual Studio 2010 (VS2010) or later versions`. 

**To choose a suitable version**: any newly released version of `VS` that compatible with GDAL is fine if you only want to run existing algorithms in PaRGO.  `Visual Studio 2010` is specifically recommended if you intend to develop new algorithms with PaRGO, since it is the last version that supports the MPI Cluster Debugger.

## 2.2 CMake

Any newly released version of CMake is fine. You can download the installer from [CMake Official Website](https://cmake.org/download/).

## 2.3 GDAL

`GDAL 1.x` and `GDAL 2.x` are both supported. 

Take GDAL 2.4.4 as an example:

1. From http://download.gisinternals.com/release.php download the 

   - [release-1900-x64-gdal-2-4-4-mapserver-7-4-3.zip](http://download.gisinternals.com/sdk/downloads/release-1900-x64-gdal-2-4-4-mapserver-7-4-3.zip) 

   - [release-1900-x64-gdal-2-4-4-mapserver-7-4-3-libs.zip](http://download.gisinternals.com/sdk/downloads/release-1900-x64-gdal-2-4-4-mapserver-7-4-3-libs.zip)


2. Unzip them into **the same folder** (e.g., D:\lib\gdal\2-4-4-vs2015x6)

3. Modify system environmental variables: (replace “D:\lib\gdal\2-4-4-vs2015x64” to your path)

   - `GDAL_ROOT`= `C:\lib\gdal\2-4-4-vs2015x64`

   - `GDAL_DATA`= `C:\lib\gdal\2-4-4-vs2015x64\bin\gdal-data`

   - `GDAL_PATHS`= `C:\lib\gdal\2-4-4-vs2015x64\bin; C:\lib\gdal\2-4-4-vs2015x64\bin\gdal\apps; C:\lib\gdal\2-4-4-vs2015x64\bin\proj\apps; C:\lib\gdal\2-4-4-vs2015x64\bin\curl;`

   - Add `%GDAL_PATHS%` to `PATH`

To test if the GDAL is installed successfully, open the command line (CMD) and type:

```
gdalinfo
```

The CMD would print the usages of `gdalinfo` if success.

## 2.4 MS-MPI

Download `MS-MPI v6 or later versions` from [GitHub (v10 or later)](https://github.com/microsoft/Microsoft-MPI/releases) or [Microsoft Archived Websites (v8.1)](https://www.microsoft.com/en-us/download/details.aspx?id=55494). 

**To choose a suitable version**: note that `VS2010` supports up to `MS-MPI v8.x`. So please install `v8.1` if you want to develop new algorithms with PaRGO.

Then:

1. Install the `msmpisdk.msi` and `MSMpiSetup.exe` to a designated path. Take the default path as an example.
2. Add those paths to system environment variables:
   - `MSMPI_BIN`=`C:\Program Files\Microsoft MPI\Bin\`
   - `MSMPI_INC`=`C:\Program Files (x86)\Microsoft SDKs\MPI\Include\`
   - `MSMPI_LIB32`=`C:\Program Files (x86)\Microsoft SDKs\MPI\Lib\x86\`
   - `MSMPI_LIB64`=`C:\Program Files (x86)\Microsoft SDKs\MPI\Lib\x64\`

To test if the MS-MPI is installed successfully, open the command line and type:

```shell
mpiexec
```

The CMD would print the usages of `mpiexec` if success.

## 2.5 (Optional) OpenMP and CUDA

Depending on specific computing platforms, PaRGO has dependency on OpenMP or CUDA.

# 3. Build

After installation, you can download the PaRGO project from GitHub using `git clone` or just download the zip file.

Enter the root folder of the PaRGO project, and build it from the command line in the following steps.

1. Create the folder

    ```shell
    cd PaRGO
    mkdir build
    cd build
    ```
    
2. Compile the project

    If you are using VS2010, use

    ```
    cmake -G "Visual Studio 10 2010 Win64" ../PaRGO -DUSE_MPI_DEBUGGER=1 ..
    ```

    If you are using VS2015, use

    ```
    cmake -G "Visual Studio 14 2015 Win64" ../PaRGO -DUSE_MPI_DEBUGGER=1 ..
    ```

    By default, the install directory is /path/to/source/bin, which can also be specified by adding `-DINSTALL_PREFIX` argument, e.g.:

    ```
    cmake -G "Visual Studio 10 2010 Win64" ../PaRGO -DUSE_MPI_DEBUGGER=1 -DINSTALL_PREFIX=D:/compile/bin/pargo ..
    ```

4. Build the project

   1. using VS GUI:

      Run `ALL_BUILD` and `INSTALL`  in the "Solution explorer -> CMakePredefinedTargets"

   2. using command line:

      Build in the visual studio command line prompt. You can find it in the Windows startup menu, in the Visual Studio folder. Note to choose the x64 version if your Windows is x64.

      For VS2010, the folder is named "Microsoft Visual studio 2010". It also can be found in the path like `"C:\ProgramData\Microsoft\Windows\Start Menu\Programs\Microsoft Visual Studio 2010\Visual Studio Tools\Visual Studio x64 Win64 command prompt (2010).lnk"`

      For VS2015 (and later), the folder is named "Visual Studio 2010".

      In the command prompt, `cd` to the "build" directory and use

      ```
      msbuild ALL_BUILD.vcxproj /p:Configuration=Release
      msbuild INSTALL.vcxproj /p:Configuration=Release
      
      ```

Reference: [CMake and Visual Studio | Cognitive Waves (wordpress.com)](https://cognitivewaves.wordpress.com/cmake-and-visual-studio/)

# 4. Run

examples:

```shell
MultiScaleLE:
mpiexec -n 4 D:\src\PaRGO\vs2010\build\apps\morphology\Debug\multiScaleLE.exe D:\arcgis-data\pargo\ywzdem5m.tif D:\arcgis-data\pargo\le\moore20.nbr D:\arcgis-data\pargo\le\output.tif 1

Slope:
mpiexec -n 4 D:\src\PaRGO\vs2010\build_gdal1\apps\morphology\Debug\slope.exe D:\arcgis-data\pargo\fcm\dem_nenjiang_10.tif D:\arcgis-data\pargo\moore.nbr D:\arcgis-data\pargo\fcm\slope_nenjiang_10.tif

IDW:
mpiexec -n 8 C:\src\PaRGO\vs2010\build\apps\spatial\Release\idw.exe -sample D:\data-arcgis\pargo\idw\sc_p_3w.shp -mask D:\data-arcgis\pargo\idw\idw_bounding_mask_10m.tif -dataNbr D:\data-arcgis\pargo\neigh.nbr -computeNbr D:\data-arcgis\pargo\neigh.nbr -out D:\data-arcgis\pargo\idw\out\temp.tif -resolution 10 -fieldIndex 3 -idwExp 2 -searchPointNum 12 -searchRange 0 -blockSize 5000 -dcmp space

FCM:
mpiexec -n 8 C:\src\PaRGO\vs2010\build\apps\spatial\Release\fcm.exe -inputs D:\data-arcgis\pargo\fcm\normalized\twi_1.tif,D:\data-arcgis\pargo\fcm\normalized\plan_1.tif,D:\data-arcgis\pargo\fcm\normalized\prof_1.tif,D:\data-arcgis\pargo\fcm\normalized\slp_1.tif, -dataNbr D:\data-arcgis\pargo\moore.nbr -computeNbr D:\data-arcgis\pargo\moore.nbr -out D:\data-arcgis\pargo\fcm\out\temp.tif -clusterNum 5 -maxIter 5 -tolerance 0.00001 -weight 2 -dcmp space
```

Note the "\PaRGO\vs2010\build\apps\spatial\\`Release`\fcm.exe" has significant better performance than the `Debug` path.

# 5. Tested Platforms

+ Windows 10 with Visual Studio 2010/2013/2015, MSMPI-v8.1, GDAL-1.11.4, GDAL-2.4.4
+ Windows 10 with msys2/mingw64 (GCC-9.1.0), MSMPI-v8.1, GDAL-3.0
+ CentOS 6.2 (cluster) with GCC-4.8.4, MPICH-3.1.4, GDAL-1.9.0
+ Red Hat Server 6.2 (cluster) with ICC-12.1.0, Intel MPI 4.0, GDAL-1.11.5
+ macOS 10.14.5 with Clang-10.0 with Xcode, OpenMPI-4.0.1, GDAL-2.4.2 (brew installed)
+ Windows 10-64bit with Visual Studio 2013

# 6. Development Based on PaRGO

## 6.1 Create a new program

  1. create program files (`.cpp` and `.h`) in a proper folder (e.g., `PaRGO\apps\spatial\myAlgorithm`), referring to the path of existing programs.
  2. modify the `CMkaeLists.txt` of the created folder (e.g., `PaRGO\apps\spatial\CMkaeLists.txt`). For a new algorithm, there are four places need to be replaced. Replace the `myAlgorithm` with the directory of your app, and `MYALGOFILES` with a name you like:

    FILE(GLOB MYALGOFILES ./myAlgorithm/*.cpp)
      
    SET(MYALGOFILES ${MYALGOFILES} ${GPRO_SRCS})
      
    ADD_EXECUTABLE(myAlgorithm ${MYALGOFILES})
      
    SET(SPATIAL_TARGETS idw
                        fcm
                        myalgorithm
                        )

Now you can open Visual Studio and start to program in the PaRGO way!

## 6.2 Programming with PaRGO & MPI

PaRGO encapsulates lots of parallel details (MPI functions) to provide serial programming experience. If you are unfamiliar with MPI and parallel programming,  you will find it relatively easy to use PaRGO. If you are familiar with MPI, you will also find it convenient to make parallel program.

### Write a local algorithm

To write a simple parallel local algorithm, you don't have to know any parallel programming knowledge. Please refer to the `reclassify` algorithm for a simple local algorithm.

### MPI basics

MPI is the Message Passing Interface. By using MPI functions, multiple processes are assigned with each others' computation tasks, and communicate with each other to exchange intermediate results. Algorithms in PaRGO almost all include the following mostly used functions. While it's fine to copy some of the functions from an existing algorithm in PaRGO (e.g., the `reclassify` algorithm) to write a simple local geocomputation algorithm, it's better to understand some basic functions before coding. 

1. `MPI_Init(int argc, char* argv[])` 
2. `MPI_Finalize()`
3. `MPI_Comm_size(MPI_Comm comm, int *size)`
4. `MPI_Comm_rank(MPI_Comm comm, int *rank)`
5. `MPI_Send(void* buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm)`
6. `MPI_Recv(void* buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status *status) `

Full API document please see [MPI Reference - Message Passing Interface | Microsoft Docs](https://docs.microsoft.com/en-us/message-passing-interface/mpi-reference).



### Write a focal algorithm



### Write a global algorithm



## 6.3 Debug

Set your algorithm as `startup project` in VS, and use the `Debug` mode to try out your program in **serial** for the first time. If it goes without error, you can start to try it in **parallel**. Take the `reclassify` algorithm running with 4 processes as an example, the property settings should be like the following.

### MPI Cluster Debugger

It is recommended to debug using MPI Cluster Debugger for better error locating, in which you can see the outputs of every processes in separate windows. To configure the MPI Cluster Debugger, right click your program in the solution explore and open **properties -> Configuration Properties -> Debug -> MPI Cluster Debugger**. 

- Run Environment: `localhost/4` 
- Application  Parameter: `-input D:\data\dem.tif -output D:\data\output.tif`

Two optional properties can specify the path and version of MPI and GDAL:

- MPIExec Command: `"C:\Program Files\Microsoft MPI\Bin\mpiexec.exe"`
- MPIExec Arguments: `-env PATH C:\lib\gdal\2-4-4-vs2015x64\bin`

### Local Windows Debugger

If you only want to use one window for outputting, you can choose the `Local Windows Debugger`, and set:

- Command: `"C:\Program Files\Microsoft MPI\Bin\mpiexec.exe"`
- Command Arguments: `-n 4 "$(TargetPath)"`

## 6.4 Run

Programs compiled in the `Release` mode would have significantly better performance than the `Debug` mode. Switch to the `Release` mode in your VS, right click your project and `Build` it, and executable files would be generated at path like `C:\src\PaRGO\vs2010\build\apps\spatial\Release\fcm.exe`. You can run it through command line or just in VS.



# 7. Common Problems

1. I cannot compile my MPI project in Visual Studio.

   Errors like `cannot open source file "mpi.h"` and `error LNK2019: unresolve external symbol...` may be due to the VS is not configured correctly. The VS needs manual configuration for a new MPI-based project.

   Right click the project in the "Solution Explorer" and click "properties":

   - In **Configuration Properties -> C/C++ -> General -> Additional Include Directories**, append `C:\Program Files (x86)\Microsoft SDKs\MPI\Include`. 
   - In Configuration Properties -> Linker ->
     -  **General -> Additional Library Directories**, append `C:\Program Files (x86)\Microsoft SDKs\MPI\Lib\x86`
     -  **Input -> Additional Dependencies**, append `msmpi.lib`

