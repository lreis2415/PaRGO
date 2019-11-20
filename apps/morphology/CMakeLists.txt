### Set each program separatly.
# relief
FILE(GLOB RELIEFFILES ./relief/*.cpp)
SET(RELIEFFILES ${RELIEFFILES} ${GPRO_SRCS})
# slope 
FILE(GLOB SLPFILES ./slope/*.cpp)
SET(SLPFILES ${SLPFILES} ${GPRO_SRCS})
## slopeFFD
#SET(SLPFFDDIR ./slopeFFD)
#SET(SLPFFDFILES ${SLPFFDDIR}/slopeOperator.cpp ${SLPFFDDIR}/pSlope.cpp ${GPRO_SRCS})
## slopeMD
#SET(SLPMDDIR ./slopeMD)
#SET(SLPMDFILES ${SLPMDDIR}/slopeOperator.cpp ${SLPMDDIR}/pSlope.cpp ${GPRO_SRCS})
## slopeSD
#SET(SLPSDDIR ./slopeSD)
#SET(SLPSDFILES ${SLPSDDIR}/slopeOperator.cpp ${SLPSDDIR}/pSlope.cpp ${GPRO_SRCS})
## slopeSFD
#SET(SLPSFDDIR ./slopeSFD)
#SET(SLPSFDFILES ${SLPSFDDIR}/slopeOperator.cpp ${SLPSFDDIR}/pSlope.cpp ${GPRO_SRCS})
## slopeTFD
#SET(SLPTFDDIR ./slopeTFD)
#SET(SLPTFDFILES ${SLPTFDDIR}/slopeOperator.cpp ${SLPTFDDIR}/pSlope.cpp ${GPRO_SRCS})
## slopeTFDW
#SET(SLPTFDWDIR ./slopeTFDW)
#SET(SLPTFDWFILES ${SLPTFDWDIR}/slopeOperator.cpp ${SLPTFDWDIR}/pSlope.cpp ${GPRO_SRCS})
# multiScaleLE
FILE(GLOB LEFILES ./MultiScaleLE/*.cpp)
SET(LEFILES ${LEFILES} ${GPRO_SRCS})

ADD_EXECUTABLE(relief ${RELIEFFILES})
ADD_EXECUTABLE(slope ${SLPFILES})
#ADD_EXECUTABLE(slopeFFD ${SLPFFDFILES})
#ADD_EXECUTABLE(slopeMD ${SLPMDFILES})
#ADD_EXECUTABLE(slopeSD ${SLPSDFILES})
#ADD_EXECUTABLE(slopeSFD ${SLPSFDFILES})
#ADD_EXECUTABLE(slopeTFD ${SLPTFDFILES})
#ADD_EXECUTABLE(slopeTFDW ${SLPTFDWFILES})
ADD_EXECUTABLE(multiScaleLE ${LEFILES})

SET(MORPH_TARGETS relief
                  slope
#                  slopeFFD
#                  slopeMD
#                  slopeSD
#                  slopeSFD
#                  slopeTFD
#                  slopeTFDW
                  multiScaleLE
                  )

foreach (c_target ${MORPH_TARGETS})
    target_link_libraries(${c_target} ${MPI_LIBRARIES} ${GDAL_LIBRARIES})
    ### For Xcode generator, specific LLVM-Clang installed by brew, and add link library of OpenMP explicitly.
    IF(XCODE AND LLVM_FOUND)
        TARGET_LINK_LIBRARIES(${c_target} ${OpenMP_LIBRARY})
        SET_TARGET_PROPERTIES(${c_target} PROPERTIES XCODE_ATTRIBUTE_CC ${LLVM_ROOT_DIR}/bin/clang)
        SET_TARGET_PROPERTIES(${c_target} PROPERTIES XCODE_ATTRIBUTE_CXX ${LLVM_ROOT_DIR}/bin/clang++)
        SET_TARGET_PROPERTIES(${c_target} PROPERTIES XCODE_ATTRIBUTE_COMPILER_INDEX_STORE_ENABLE "No")
    ENDIF()
    install(TARGETS ${c_target} DESTINATION ${INSTALL_DIR})
    LIST(APPEND APP_TARGETS ${c_target})
    IF (MSVC_VERSION GREATER 1600 OR (MSVC_VERSION EQUAL 1600 AND NOT USE_MPI_DEBUGGER) OR XCODE)
        SET_PROPERTY(TARGET ${c_target} PROPERTY FOLDER "Morphology")
    ENDIF ()
endforeach (c_target ${MORPH_TARGETS})

SET(APP_TARGETS "${APP_TARGETS}" PARENT_SCOPE)
