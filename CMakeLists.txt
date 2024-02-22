# Min version 3.9 for OpenMP find_package
cmake_minimum_required(VERSION 3.9)
project(RetinoTectal)

message(STATUS "Install prefix: ${CMAKE_INSTALL_PREFIX}")
message(STATUS "  (This can be changed with `cmake -DCMAKE_INSTALL_PREFIX=/some/place`")

# From CMAKE_SYSTEM work out which of __OSX__, __GLN__, __NIX__ or __WIN__ are required
message(STATUS "Operating system: " ${CMAKE_SYSTEM})
if(CMAKE_SYSTEM MATCHES Linux.*)
  set(EXTRA_HOST_DEFINITION "-D__GLN__")
elseif(CMAKE_SYSTEM MATCHES BSD.*)
  set(EXTRA_HOST_DEFINITION "-D__NIX__")
elseif(APPLE)
  set(EXTRA_HOST_DEFINITION "-D__OSX__")
else()
  message(ERROR "Operating system not supported: " ${CMAKE_SYSTEM})
endif()

# morphologica uses c++-17 language features
set(CMAKE_CXX_STANDARD 17)

# Add the host definition to CXXFLAGS along with other switches,
# depending on OS/Compiler and your needs/preferences
if (APPLE)
  set(CMAKE_CXX_FLAGS "${EXTRA_HOST_DEFINITION} -Wall -Wfatal-errors -g -O3")
else() # assume g++ (or a gcc/g++ mimic like Clang)
  set(CMAKE_CXX_FLAGS "${EXTRA_HOST_DEFINITION} -Wall -Wfatal-errors -g -Wno-unused-result -Wno-unknown-pragmas -march=native -O0")
endif()

# Tell clang to be quiet about brace initialisers:
if(CMAKE_CXX_COMPILER_ID MATCHES Clang)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-missing-braces")
endif()

# Add OpenMP flags here, if necessary
find_package(OpenMP)
if(OpenMP_FOUND)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

# Additional GL compiler flags.
set(OpenGL_GL_PREFERENCE "GLVND")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DGL3_PROTOTYPES -DGL_GLEXT_PROTOTYPES")
if(APPLE)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DGL_SILENCE_DEPRECATION")
endif()

# Tell the program where the morph fonts are, to compile them into the binary
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DMORPH_FONTS_DIR=\"\\\"${PROJECT_SOURCE_DIR}/morphologica/fonts\\\"\"")

# Lib finding - we have to do this, so that the correct include paths
# are present for OpenGL, etc; they're all included by morphologica
# header files.
find_package(HDF5 REQUIRED)
find_package(Armadillo REQUIRED)
find_package(OpenGL REQUIRED)
find_package(glfw3 3.3 REQUIRED)
find_package(Freetype REQUIRED)

# Define collections of includes for the dependencies
set(MORPH_INC_CORE ${ARMADILLO_INCLUDE_DIR} ${ARMADILLO_INCLUDE_DIRS} ${HDF5_INCLUDE_DIR})
set(MORPH_INC_GL ${OPENGL_INCLUDE_DIR} ${GLFW3_INCLUDE_DIR} ${FREETYPE_INCLUDE_DIRS})
include_directories(${MORPH_INC_CORE} ${MORPH_INC_GL})

# Assuming that you installed morphologica in-tree (i.e. 'next to' schnakenberg.cpp).
set(MORPH_INCLUDE_PATH "${PROJECT_SOURCE_DIR}/morphologica" CACHE PATH "The path to morphologica")
include_directories(BEFORE ${MORPH_INCLUDE_PATH}/include) # Allows GL3/gl3.h to be found
include_directories(BEFORE ${MORPH_INCLUDE_PATH})         # Allows morph/Header.h to be found

if(OpenMP_FOUND)
  link_libraries(OpenMP::OpenMP_CXX)
endif()

if(NOT EXISTS "${MORPH_INCLUDE_PATH}/morph")
  message(FATAL_ERROR "Make sure you have installed morphologica at the location: ${MORPH_INCLUDE_PATH}. The suggested way to do this is to git clone morphologica inside the ${PROJECT_NAME} base directory")
endif()

add_subdirectory(sim)

# For debugging of variables:
option(DEBUG_VARIABLES OFF)
if(DEBUG_VARIABLES)
  get_cmake_property(_variableNames VARIABLES)
  foreach (_variableName ${_variableNames})
    message(STATUS "${_variableName}=${${_variableName}}")
  endforeach()
endif(DEBUG_VARIABLES)
