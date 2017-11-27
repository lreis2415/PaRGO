### Set each program separatly.
# Remove pit and falt
SET(PITREMDIR ./pitremove)
SET(PITREMFILES ${PITREMDIR}/pitRemoveOperator.cpp ${PITREMDIR}/pPitRemove.cpp ${GPRO_SRCS})
# Single flow direction D8
SET(D8DIR ./flowdird8)
SET(D8FILES ${D8DIR}/d8Operator.cpp ${D8DIR}/pD8.cpp ${GPRO_SRCS})
# Multi-flow direction MFD-md
SET(MFDMDDIR ./flowdirmfdmd)
SET(MFDMDFILES ${MFDMDDIR}/mfdOperator.cpp ${MFDMDDIR}/pMFDmd.cpp ${GPRO_SRCS})
# SCA D8
SET(SCAD8DIR ./scad8)
SET(SCAD8FILES ${SCAD8DIR}/scaOperator.cpp ${SCAD8DIR}/pSCA.cpp ${GPRO_SRCS})
# SCA Differential
SET(SCADIFFERDIR ./scadiffer)
SET(SCADIFFERFILES ${SCADIFFERDIR}/SCAOperator.cpp ${SCADIFFERDIR}/pSCA.cpp ${GPRO_SRCS})
# SCA MFD-md
SET(SCAMFDMDDIR ./scamfdmd)
SET(SCAMFDMDFILES ${SCAMFDMDDIR}/scaOperator.cpp ${SCAMFDMDDIR}/pSCA.cpp ${GPRO_SRCS})
# TWI ?? algorithm
SET(TWIDIR ./twi)
SET(TWIFILES ${TWIDIR}/TWIOperator.cpp ${TWIDIR}/pTWI.cpp ${GPRO_SRCS})

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
    install(TARGETS ${c_target} DESTINATION ${INSTALL_DIR})
    LIST(APPEND APP_TARGETS ${c_target})
    IF (MSVC)
        SET_PROPERTY(TARGET ${c_target} PROPERTY FOLDER "Hydrology")
    ENDIF ()
endforeach (c_target ${HYDRO_TARGETS})

SET(APP_TARGETS "${APP_TARGETS}" PARENT_SCOPE)
