### Set each program separatly.

FILE(GLOB DEMO1FILES ./demo1/*.cpp)
SET(DEMO1FILES ${DEMO1FILES} ${GPRO_SRCS})
ADD_EXECUTABLE(demo1_reclassify ${DEMO1FILES})


FILE(GLOB DEMO2FILES ./demo2/*.cpp)
SET(DEMO2FILES ${DEMO2FILES} ${GPRO_SRCS})
ADD_EXECUTABLE(demo2_slope ${DEMO2FILES})

######## $1: variable, name it whatever
######## $2: directory name of your app
######## Add your app here, replace $1 and $2 ########
# FILE(GLOB $1 ./$2/*.cpp)
# SET($1 ${$1} ${GPRO_SRCS})
# ADD_EXECUTABLE($2 ${$1})
######## Add your app here, replace $1 and $2 ########

SET(DEMO_TARGETS demo1_reclassify
                 demo2_slope
                ######## Add your app here, replace $2 ########
                #$2
                ######## Add your app here, replace $2 ########
                )

foreach (c_target ${DEMO_TARGETS})
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
        SET_PROPERTY(TARGET ${c_target} PROPERTY FOLDER "Demo")
    ENDIF ()
endforeach (c_target ${DEMO_TARGETS})

SET(APP_TARGETS "${APP_TARGETS}" PARENT_SCOPE)
