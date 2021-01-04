# Define collections of links that have to be made for morphologica
set(MORPH_LIBS_CORE ${ARMADILLO_LIBRARY} ${ARMADILLO_LIBRARIES} ${HDF5_C_LIBRARIES} ${LAPACK_LIBRARIES} jsoncpp_lib)
set(MORPH_LIBS_GL ${OpenCV_LIBS} OpenGL::GL Freetype::Freetype glfw)

# James model 0. 2-D RD Karbowski with variable numbers of guidance
# molecules and TC populations and runtime setting of parameters to
# examine behaviour.
# add_executable(james0 james1.cpp)
# target_compile_definitions(james0 PUBLIC FLT=float COMPILE_PLOTTING)
# target_link_libraries(james0 ${MORPH_LIBS_CORE} ${MORPH_LIBS_GL})

# Competition based model based on the model in the BarrelEmerge paper
add_executable(rettec rettec.cpp)
target_compile_definitions(rettec PUBLIC FLT=float COMPILE_PLOTTING AXONCOMP)
target_link_libraries(rettec ${MORPH_LIBS_CORE} ${MORPH_LIBS_GL})

# Competition model with no graphics
add_executable(rettec_c rettec.cpp)
target_compile_definitions(rettec_c PUBLIC FLT=float AXONCOMP)
target_link_libraries(rettec_c ${MORPH_LIBS_CORE})

# Difference here is that AXONCOMP is NOT defined and so a
# model is compiled from rd_rettec_nocomp.h, which implements a
# non-competition based stopping criterion, intended to be a
# 'Sperry-like' model.
add_executable(rettec_nc rettec.cpp)
target_compile_definitions(rettec_nc PUBLIC FLT=float COMPILE_PLOTTING DEBUG)
target_link_libraries(rettec_nc ${MORPH_LIBS_CORE} ${MORPH_LIBS_GL})

# No-graphics Sperry-like model
add_executable(rettec_nc_c rettec.cpp)
target_compile_definitions(rettec_nc_c PUBLIC FLT=float)
target_link_libraries(rettec_nc_c ${MORPH_LIBS_CORE})
