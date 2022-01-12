# Define collections of links that have to be made for morphologica
set(MORPH_LIBS_CORE ${ARMADILLO_LIBRARY} ${ARMADILLO_LIBRARIES} ${HDF5_C_LIBRARIES} ${LAPACK_LIBRARIES})
set(MORPH_LIBS_GL ${OpenCV_LIBS} OpenGL::GL Freetype::Freetype glfw)

# Agent-based modelling
add_subdirectory (agent)

# Reaction-diffusion style modelling
add_subdirectory (rd)
