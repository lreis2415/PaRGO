# ----------------------------------------------------------------------------
#  Root CMake file for PaRGO.
#
#     Please from the off-tree build directory, invoke:
#       $ cmake -G "<Generator>" <path/to/PaRGO>
#
#       Optional cmake options:
#         -DINSTALL_PREFIX=<path/to/install>
#         -DUSE_MPI_DEBUGGER When set to 1, apps will not included in folders and MPI debugger in VS2010 works!
#         -DLLVM_ROOT_DIR Specific the root directory of brew installed LLVM, e.g., /usr/local/opt/llvm
#
#  Routine testing platforms and compilers include:
#     1. Windows 10 with Visual Studio 2010/2013/2015, MSMPI-v8.1, GDAL-1.11.4
#     2. Windows 10 with msys2/mingw64 (GCC-9.1.0), MSMPI-v8.1, GDAL-3.0
#     3. CentOS 6.2 (cluster) with GCC-4.8.4, MPICH-3.1.4, GDAL-1.9.0
#     4. Red Hat Server 6.2 (cluster) with ICC-12.1.0, Intel MPI 4.0, GDAL-1.11.5
#     5. macOS 10.14.5 with Clang-10.0 with Xcode, OpenMPI-4.0.1, GDAL-2.4.2 (brew installed)
#     TODO: Add CUDA test enviroments.
#
#  Created and maintained by Liang-Jun Zhu (zlj@lreis.ac.cn)
#  Latest updated: Jul. 21, 2019
#  Copyright (C) 2018-2019 Lreis, IGSNRR, CAS
# ----------------------------------------------------------------------------
### Disable in-source builds to prevent source tree corruption.
IF(" ${CMAKE_SOURCE_DIR}" STREQUAL " ${CMAKE_BINARY_DIR}")
  MESSAGE(FATAL_ERROR "
FATAL: In-source builds are not allowed. 
       You should create a separate directory for build files."
)
ENDIF()
CMAKE_MINIMUM_REQUIRED(VERSION 3.1 FATAL_ERROR)
if(POLICY CMP0015)
  cmake_policy(SET CMP0015 NEW)
endif()
if(POLICY CMP0054)
  cmake_policy(SET CMP0054 NEW)
endif()
if(POLICY CMP0074)
  cmake_policy(SET CMP0074 NEW)
endif()

### Set coding language
PROJECT(PaRGO CXX C)

### PaRGO Project name with configuration.
SET(PRJNAME "PaRGO")

### Set cmake find path
set(CMAKE_MODULE_PATH
        ${CMAKE_CURRENT_SOURCE_DIR}/cmake
        ${CMAKE_MODULE_PATH}
        )

### Set default or specified installation directory.
IF (INSTALL_PREFIX)
  SET(CMAKE_INSTALL_PREFIX ${INSTALL_PREFIX} CACHE PATH "Set specified install path" FORCE)
ELSE()
  SET(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR}/bin CACHE PATH "Set default install path" FORCE)
ENDIF ()
SET(INSTALL_DIR ${CMAKE_INSTALL_PREFIX})

### In case of Makefiles if the user does not setup CMAKE_BUILD_TYPE, assume it's Release.
IF(CMAKE_GENERATOR MATCHES "Makefiles|Ninja" AND "${CMAKE_BUILD_TYPE}" STREQUAL "")
  SET(CMAKE_BUILD_TYPE Release)
ENDIF()

### Allow project folders in MSVC or XCODE.
IF(MSVC OR XCODE)
  SET_PROPERTY(GLOBAL PROPERTY USE_FOLDERS ON)
ENDIF()

### OS platform.
if (WIN32)
  add_definitions(-DWINDOWS)
  if (MSVC)
    add_definitions(-DMSVC)
  endif ()
elseif (APPLE)
  add_definitions(-DMACOSX)
  if (XCODE)
    add_definitions(-DXCODE)
  endif ()
else ()
  add_definitions(-DLINUX)
endif ()

### Assume built-in pthreads on macOS X
IF(APPLE)
    set(CMAKE_THREAD_LIBS_INIT "-lpthread")
    set(CMAKE_HAVE_THREADS_LIBRARY 1)
    set(CMAKE_USE_WIN32_THREADS_INIT 0)
    set(CMAKE_USE_PTHREADS_INIT 1)
    set(THREADS_PREFER_PTHREAD_FLAG ON)
ENDIF()

### Break in case of popular CMake configuration mistakes.
IF(NOT CMAKE_SIZEOF_VOID_P GREATER 0)
  MESSAGE(FATAL_ERROR "CMake fails to determine the bitness of the target platform.
  Please check your CMake and compiler installation."
)
ENDIF()

### Set compiling version, MPI and OpenMP are required!
SET(WITH_OPENMP 1)
SET(WITH_MPI 1)
IF (CUDA)  # TODO
  SET(WITH_CUDA 1)
ENDIF()

IF (NOT USE_MPI_DEBUGGER)
  SET(USE_MPI_DEBUGGER 0)
ENDIF()

### Detect compiler and target platform architecture.
INCLUDE(Utils)
SET(ENABLE_CXX11 1)
INCLUDE(DetectCXXCompiler)

### Add standard paths or specified paths for Find libraries and headers.
INCLUDE(AddFindPaths)

### Find packages.
### Find MPI which is required
FIND_PACKAGE(MPI REQUIRED)
IF(MPI_FOUND)
  SET(WITH_MPI 1)
  SET(PRJNAME ${PRJNAME}_MPI)
  MESSAGE(STATUS "Compiling with MPI...")
ELSE()
  MESSAGE(FATAL_ERROR "FATAL: Could not find MPI, which is required!")
ENDIF()

### Find OMP which is required
IF(XCODE)
  # To use OpenMP in Xcode built by CMake, the recommended way is to 
  #   install llvm by 'brew install llvm', which has builtin 'libomp'
  IF (NOT LLVM_ROOT_DIR)
    SET(LLVM_ROOT_DIR "/usr/local/opt/llvm")
  ENDIF ()
  FIND_PACKAGE(LLVM)
  IF(LLVM_FOUND)
    geo_include_directories(${LLVM_LIBRARY_DIRS}/clang/${LLVM_VERSION_BASE_STRING}/include)
	SET(OpenMP_CXX_FLAGS "-fopenmp")
    SET(OpenMP_C_FLAGS "-fopenmp")
	FIND_LIBRARY(OpenMP_LIBRARY
		NAMES "omp"
		HINTS ${LLVM_LIBRARY_DIRS}) 
    SET(OPENMP_FOUND 1)
  ENDIF ()
ELSE ()  
  FIND_PACKAGE(OpenMP)
ENDIF ()
IF(OPENMP_FOUND)
  SET(WITH_OPENMP 1)
  SET(PRJNAME ${PRJNAME}_OMP)
  ADD_DEFINITIONS(-DSUPPORT_OMP)
  INCLUDE(SetFlagOpenMP)
  MESSAGE(STATUS "Compiling with OpenMP...")
ENDIF()

FIND_PACKAGE(GDAL REQUIRED)
IF (GDAL_FOUND)
  SET(WITH_GDAL 1)
  ADD_DEFINITIONS(-DUSE_GDAL)
  INCLUDE(SetFlagGDAL)
  geo_include_directories(${GDAL_INCLUDE_DIR})
  MESSAGE(STATUS "Compiling with GDAL")
ELSE()
  MESSAGE(FATAL_ERROR "FATAL: Could not find GDAL, which is required!")
ENDIF ()

### Set final project name
PROJECT(${PRJNAME})

### Use statically or dynamically linked CRT? Default: dynamic
IF(MSVC)
  INCLUDE(CRTLinkage)
ENDIF(MSVC)

SET(GPRO_DIR ${CMAKE_CURRENT_SOURCE_DIR}/gpro)
SET(APPS_DIR ${CMAKE_CURRENT_SOURCE_DIR}/apps)

# GPRO source code
FILE(GLOB GPRO_SRCS ${GPRO_DIR}/*.cpp ${GPRO_DIR}/*.h)
# Adds include directories in such a way that directories from the TauDEM source tree go first.
geo_include_directories(${GPRO_DIR} ${GDAL_INCLUDE_DIR} ${MPI_INCLUDE_PATH})
ADD_SUBDIRECTORY(${APPS_DIR})

unset(USE_MPI_DEBUGGER CACHE)
### Build platform.
STATUS("")
STATUS("  Platform:")
IF(NOT CMAKE_VERSION VERSION_LESS 2.8.11 AND NOT BUILD_INFO_SKIP_TIMESTAMP)
  STRING(TIMESTAMP TIMESTAMP "" UTC)
  IF(TIMESTAMP)
    STATUS("    Timestamp:"    ${TIMESTAMP})
  ENDIF()
ENDIF()
STATUS("    Host:"             ${CMAKE_HOST_SYSTEM_NAME} ${CMAKE_HOST_SYSTEM_VERSION} ${CMAKE_HOST_SYSTEM_PROCESSOR})
IF(CMAKE_CROSSCOMPILING)
  STATUS("    Target:"         ${CMAKE_SYSTEM_NAME} ${CMAKE_SYSTEM_VERSION} ${CMAKE_SYSTEM_PROCESSOR})
ENDIF()
STATUS("    CMake:"            ${CMAKE_VERSION})
STATUS("    CMake generator:"  ${CMAKE_GENERATOR})
STATUS("    CMake build tool:" ${CMAKE_BUILD_TOOL})
IF(MSVC)
  STATUS("    MSVC:"           ${MSVC_VERSION})
ENDIF()
IF(CMAKE_GENERATOR MATCHES Xcode)
  STATUS("    Xcode:"          ${XCODE_VERSION})
ENDIF()
IF(NOT CMAKE_GENERATOR MATCHES "Xcode|Visual Studio")
  STATUS("    Configuration:"  ${CMAKE_BUILD_TYPE})
ENDIF()

### C/C++ options.
IF(CMAKE_CXX_COMPILER_VERSION)
  SET(GEO_COMPILER_STR "${CMAKE_CXX_COMPILER} ${CMAKE_CXX_COMPILER_ARG1} (ver ${CMAKE_CXX_COMPILER_VERSION})")
ELSE()
  SET(GEO_COMPILER_STR "${CMAKE_CXX_COMPILER} ${CMAKE_CXX_COMPILER_ARG1}")
ENDIF()
STRING(STRIP "${GEO_COMPILER_STR}" GEO_COMPILER_STR)

STATUS("")
STATUS("  C/C++:")
IF(ENABLE_CXX11 OR HAVE_CXX11)
STATUS("    C++11:" HAVE_CXX11 THEN YES ELSE NO)
ENDIF()
STATUS("    C++ Compiler:"           ${GEO_COMPILER_STR})
STATUS("    C++ flags (Release):"    ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE})
STATUS("    C++ flags (Debug):"      ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG})
STATUS("    C Compiler:"             ${CMAKE_C_COMPILER} ${CMAKE_C_COMPILER_ARG1})
STATUS("    C flags (Release):"      ${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_RELEASE})
STATUS("    C flags (Debug):"        ${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_DEBUG})
IF(WIN32)
  STATUS("    Linker flags (Release):" ${CMAKE_EXE_LINKER_FLAGS} ${CMAKE_EXE_LINKER_FLAGS_RELEASE})
  STATUS("    Linker flags (Debug):"   ${CMAKE_EXE_LINKER_FLAGS} ${CMAKE_EXE_LINKER_FLAGS_DEBUG})
ELSE()
  STATUS("    Linker flags (Release):" ${CMAKE_SHARED_LINKER_FLAGS} ${CMAKE_SHARED_LINKER_FLAGS_RELEASE})
  STATUS("    Linker flags (Debug):"   ${CMAKE_SHARED_LINKER_FLAGS} ${CMAKE_SHARED_LINKER_FLAGS_DEBUG})
ENDIF()

### Dependencies.
STATUS("    GDAL: ${GDAL_LIBRARIES} ${GDAL_INCLUDE_DIR}")
STATUS("    MPI: ${MPI_LIBRARIES} ${MPI_INCLUDE_PATH}")

### Auxiliary.
STATUS("")
STATUS("  Install path:" ${CMAKE_INSTALL_PREFIX})
STATUS("")