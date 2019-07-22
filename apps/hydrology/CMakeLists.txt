### Set each program separatly.
# Remove pit and flat
FILE(GLOB PITREMFILES ./pitremove/*.cpp)
SET(PITREMFILES ${PITREMFILES} ${GPRO_SRCS})

# Single flow direction D8
FILE(GLOB D8FILES ./flowdird8/*.cpp)
SET(D8FILES ${D8FILES} ${GPRO_SRCS})

# Multi-flow direction MFD-md
FILE(GLOB MFDMDFILES ./flowdirmfdmd/*.cpp)
SET(MFDMDFILES ${MFDMDFILES} ${GPRO_SRCS})

# SCA D8
FILE(GLOB SCAD8FILES ./scad8/*.cpp)
SET(SCAD8FILES ${SCAD8FILES} ${GPRO_SRCS})

# SCA Differential
FILE(GLOB SCADIFFERFILES ./scadiffer/*.cpp)
SET(SCADIFFERFILES ${SCADIFFERFILES} ${GPRO_SRCS})

# SCA MFD-md
FILE(GLOB SCAMFDMDFILES ./scamfdmd/*.cpp)
SET(SCAMFDMDFILES ${SCAMFDMDFILES} ${GPRO_SRCS})

# TWI algorithm
FILE(GLOB TWIFILES ./twi/*.cpp)
SET(TWIFILES ${TWIFILES} ${GPRO_SRCS})

ADD_EXECUTABLE(pitremove ${PITREMFILES})
ADD_EXECUTABLE(flowdird8 ${D8FILES})
ADD_EXECUTABLE(flowdirmfdmd ${MFDMDFILES})
ADD_EXECUTABLE(scad8 ${SCAD8FILES})
ADD_EXECUTABLE(scadiffer ${SCADIFFERFILES})
ADD_EXECUTABLE(scamfdmd ${SCAMFDMDFILES})
ADD_EXECUTABLE(twi ${TWIFILES})

SET(HYDRO_TARGETS pitremove
                  flowdird8
                  flowdirmfdmd
                  scad8
                  scadiffer
                  scamfdmd
                  twi
                  )

foreach (c_target ${HYDRO_TARGETS})
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
    IF (MSVC OR XCODE)
        SET_PROPERTY(TARGET ${c_target} PROPERTY FOLDER "Hydrology")
    ENDIF ()
endforeach (c_target ${HYDRO_TARGETS})

SET(APP_TARGETS "${APP_TARGETS}" PARENT_SCOPE)
