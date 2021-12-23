# GPRO source code
SET(GPRO_SRCS ${GPRO_SRCS} PARENT_SCOPE)
# Adds include directories in such a way that directories from the TauDEM source tree go first.
geo_include_directories(${GPRO_DIR} ${GDAL_INCLUDE_DIR} ${MPI_INCLUDE_PATH})

ADD_SUBDIRECTORY(./hydrology)
ADD_SUBDIRECTORY(./morphology)
ADD_SUBDIRECTORY(./spatial)
ADD_SUBDIRECTORY(./demo)

### For CLion to implement the "make install" command
ADD_CUSTOM_TARGET(install_${PROJECT_NAME}
        $(MAKE) install
        DEPENDS ${APP_TARGETS}
        COMMENT "Installing ${PROJECT_NAME}")
