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

add_executable(tissuevis tissuevis.cpp)
target_compile_definitions(tissuevis PUBLIC FLT=float VISUALISE=1)
target_link_libraries(tissuevis ${MORPH_LIBS_CORE} ${MORPH_LIBS_GL})

add_executable(testnet testnet.cpp)
target_compile_definitions(testnet PUBLIC FLT=float VISUALISE=1)
target_link_libraries(testnet ${MORPH_LIBS_CORE} ${MORPH_LIBS_GL})

add_executable(search1 search1.cpp)
target_compile_definitions(search1 PUBLIC FLT=float VISUALISE=1)
target_link_libraries(search1 ${MORPH_LIBS_CORE} ${MORPH_LIBS_GL})

add_executable(search1c search1.cpp)
target_compile_definitions(search1c PUBLIC FLT=float)
target_link_libraries(search1c ${MORPH_LIBS_CORE} ${MORPH_LIBS_GL})

add_executable(search2 search2.cpp)
target_compile_definitions(search2 PUBLIC FLT=float VISUALISE=1)
target_link_libraries(search2 ${MORPH_LIBS_CORE} ${MORPH_LIBS_GL})

add_executable(search2c search2.cpp)
target_compile_definitions(search2c PUBLIC FLT=float)
target_link_libraries(search2c ${MORPH_LIBS_CORE} ${MORPH_LIBS_GL})
