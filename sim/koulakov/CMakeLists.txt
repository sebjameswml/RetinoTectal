# Koulakov and Tsigankov model(s)

# Models from the paper 'A stochastic model for retinocollicular map
# development'

add_executable(koulakov1d_fig4 koulakov1d_fig4.cpp)
target_link_libraries(koulakov1d_fig4 ${MORPH_LIBS_CORE} ${MORPH_LIBS_GL})

add_executable(koulakov1d_fig5 koulakov1d_fig5.cpp)
target_link_libraries(koulakov1d_fig5 ${MORPH_LIBS_CORE} ${MORPH_LIBS_GL})

add_executable(koulakov1d_fig6 koulakov1d_fig6.cpp)
target_link_libraries(koulakov1d_fig6 ${MORPH_LIBS_CORE} ${MORPH_LIBS_GL})
