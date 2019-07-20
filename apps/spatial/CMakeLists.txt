### Set each program separatly.
# Inverse Distance Weight for interpolation.
FILE(GLOB IDWFILES ./idw/*.cpp)
SET(IDWFILES ${IDWFILES} ${GPRO_SRCS})

ADD_EXECUTABLE(idw ${IDWFILES})

SET(SPATIAL_TARGETS idw
                    )

foreach (c_target ${SPATIAL_TARGETS})
    target_link_libraries(${c_target} ${MPI_LIBRARIES} ${GDAL_LIBRARIES})
    install(TARGETS ${c_target} DESTINATION ${INSTALL_DIR})
    LIST(APPEND APP_TARGETS ${c_target})
    IF (MSVC)
        SET_PROPERTY(TARGET ${c_target} PROPERTY FOLDER "Spatial")
    ENDIF ()
endforeach (c_target ${SPATIAL_TARGETS})

SET(APP_TARGETS "${APP_TARGETS}" PARENT_SCOPE)
