# Competition agent model
add_executable(agent1 agent1.cpp)
target_compile_definitions(agent1 PUBLIC FLT=float VISUALISE=1)
target_link_libraries(agent1 ${MORPH_LIBS_CORE} ${MORPH_LIBS_GL})

#add_executable(agent1d agent1.cpp)
#target_compile_definitions(agent1d PUBLIC FLT=double VISUALISE=1)
#target_link_libraries(agent1d ${MORPH_LIBS_CORE} ${MORPH_LIBS_GL})

add_executable(agent1c agent1.cpp)
target_compile_definitions(agent1c PUBLIC FLT=float)
target_link_libraries(agent1c ${MORPH_LIBS_CORE} ${MORPH_LIBS_GL})

add_executable(agent1_eval agent1_eval.cpp)
target_compile_definitions(agent1_eval PUBLIC FLT=float VISUALISE=1)
target_link_libraries(agent1_eval ${MORPH_LIBS_CORE} ${MORPH_LIBS_GL})

add_executable(agent1_evalc agent1_eval.cpp)
target_compile_definitions(agent1_evalc PUBLIC FLT=float)
target_link_libraries(agent1_evalc ${MORPH_LIBS_CORE} ${MORPH_LIBS_GL})

add_executable(tissuevis tissuevis.cpp)
target_compile_definitions(tissuevis PUBLIC FLT=float VISUALISE=1)
target_link_libraries(tissuevis ${MORPH_LIBS_CORE} ${MORPH_LIBS_GL})

add_executable(testnet testnet.cpp)
target_compile_definitions(testnet PUBLIC FLT=float VISUALISE=1)
target_link_libraries(testnet ${MORPH_LIBS_CORE} ${MORPH_LIBS_GL})

# Run anneal1, watching each individual simulation
add_executable(anneal1 anneal1.cpp)
target_compile_definitions(anneal1 PUBLIC FLT=float VISUALISE=1)
target_link_libraries(anneal1 ${MORPH_LIBS_CORE} ${MORPH_LIBS_GL})

# Run anneal1, watching the optimisation varables
add_executable(anneal1c anneal1.cpp)
target_compile_definitions(anneal1c PUBLIC FLT=float OPTVIS=1)
target_link_libraries(anneal1c ${MORPH_LIBS_CORE} ${MORPH_LIBS_GL})

# Run anneal1 with no visuals at all
add_executable(anneal1c0 anneal1.cpp)
target_compile_definitions(anneal1c0 PUBLIC FLT=float)
target_link_libraries(anneal1c0 ${MORPH_LIBS_CORE} ${MORPH_LIBS_GL})

add_executable(anneal1read anneal1read.cpp)
target_compile_definitions(anneal1read PUBLIC FLT=float)
target_link_libraries(anneal1read ${MORPH_LIBS_CORE} ${MORPH_LIBS_GL})

# For 4D searches. anneal1read.cpp would need quite a bit of hacking
# to make this work. I DID make anneal.cpp capable of handling >3
# search dimensions though, and it would be nice to do the GJ
# optimization with s, m_g, r_j and m_j and then vis the results
# nicely. However, I can just look at the results in Octave, using the
# HDF5.
#
#add_executable(anneal1read4 anneal1read.cpp)
#target_compile_definitions(anneal1read4 PUBLIC FLT=float D=4)
#target_link_libraries(anneal1read4 ${MORPH_LIBS_CORE} ${MORPH_LIBS_GL})
