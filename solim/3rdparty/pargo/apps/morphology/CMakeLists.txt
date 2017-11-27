### Set each program separatly.
# relief
SET(RELIEFDIR ./relief)
SET(RELIEFFILES ${RELIEFDIR}/reliefOperator.cpp ${RELIEFDIR}/pRelief.cpp ${GPRO_SRCS})
# slope using ?? algorithm
SET(SLPDIR ./slope)
SET(SLPFILES ${SLPDIR}/slopeOperator.cpp ${SLPDIR}/pSlope.cpp ${GPRO_SRCS})

ADD_EXECUTABLE(relief ${RELIEFFILES})
ADD_EXECUTABLE(slope ${SLPFILES})

SET(MORPH_TARGETS relief
                  slope
                  )

foreach (c_target ${MORPH_TARGETS})
    target_link_libraries(${c_target} ${MPI_LIBRARIES} ${GDAL_LIBRARIES})
    install(TARGETS ${c_target} DESTINATION ${INSTALL_DIR})
    LIST(APPEND APP_TARGETS ${c_target})
    IF (MSVC)
        SET_PROPERTY(TARGET ${c_target} PROPERTY FOLDER "Morphology")
    ENDIF ()
endforeach (c_target ${MORPH_TARGETS})

SET(APP_TARGETS "${APP_TARGETS}" PARENT_SCOPE)
