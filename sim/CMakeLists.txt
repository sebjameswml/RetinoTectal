# James model 0. 2-D RD Karbowski with variable numbers of guidance
# molecules and TC populations and runtime setting of parameters to
# examine behaviour.
# add_executable(james0 james1.cpp)
# target_compile_definitions(james0 PUBLIC FLT=float COMPILE_PLOTTING)
# target_link_libraries(james0 ${MORPH_LDFLAGS} ${HDF5LIBLINK} ${JSONLIBLINK})

# Divisive normalization of a_i along with competition (linear)
add_executable(rettec rettec.cpp)
target_compile_definitions(rettec PUBLIC FLT=float COMPILE_PLOTTING DNCOMP)
target_link_libraries(rettec ${MORPH_LDFLAGS} ${HDF5LIBLINK} ${JSONLIBLINK} ${OpenCV_LIBRARIES} glfw)

add_executable(rettec_c rettec.cpp)
target_compile_definitions(rettec_c PUBLIC FLT=float DNCOMP)
target_link_libraries(rettec_c ${MORPH_LDFLAGS} ${HDF5LIBLINK} ${JSONLIBLINK} ${OpenCV_LIBRARIES})
