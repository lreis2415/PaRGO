### Set each program separatly.
# Inverse Distance Weight for interpolation.
FILE(GLOB IDWFILES ./idw/*.cpp)
SET(IDWFILES ${IDWFILES} ${GPRO_SRCS})
ADD_EXECUTABLE(idw ${IDWFILES})

FILE(GLOB FCMFILES ./fcm/*.cpp)
SET(FCMFILES ${FCMFILES} ${GPRO_SRCS})
ADD_EXECUTABLE(fcm ${FCMFILES})

FILE(GLOB RECLASSIFYFILES ./reclassify/*.cpp)
SET(RECLASSIFYFILES ${RECLASSIFYFILES} ${GPRO_SRCS})
ADD_EXECUTABLE(reclassify ${RECLASSIFYFILES})

SET(SPATIAL_TARGETS idw
                    fcm
                    reclassify
                    )

foreach (c_target ${SPATIAL_TARGETS})
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
        SET_PROPERTY(TARGET ${c_target} PROPERTY FOLDER "Spatial")
    ENDIF ()
endforeach (c_target ${SPATIAL_TARGETS})

SET(APP_TARGETS "${APP_TARGETS}" PARENT_SCOPE)
