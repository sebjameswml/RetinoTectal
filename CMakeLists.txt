# Min version 3.1 for the modern way to select c++-11/14/17/etc
cmake_minimum_required(VERSION 3.1)
project(RetinoTectal)

message(STATUS "Install prefix: ${CMAKE_INSTALL_PREFIX}")
message(STATUS "  (This can be changed with `cmake -DCMAKE_INSTALL_PREFIX=/some/place`")

set (CMAKE_CXX_STANDARD 17)

# From CMAKE_SYSTEM work out which of __OSX__, __GLN__, __NIX__ or
# __WIN__ are required
message(STATUS "Operating system: " ${CMAKE_SYSTEM})
if(CMAKE_SYSTEM MATCHES Linux.*)
  set(MORPH_HOST_DEFINITION "-D__GLN__")
elseif(CMAKE_SYSTEM MATCHES BSD.*)
  set(MORPH_HOST_DEFINITION "-D__NIX__")
elseif(APPLE)
  set(MORPH_HOST_DEFINITION "-D__OSX__")
elseif(CMAKE_SYSTEM MATCHES Win.*)
  set(MORPH_HOST_DEFINITION "-D__WIN__")
else()
  message(ERROR "Operating system not supported: " ${CMAKE_SYSTEM})
endif()

set(CMAKE_CXX_STANDARD 17)

# Add the host definition to CXXFLAGS along with other switches
#
#OpenMP. With cmake 3.9+, this is the right way to do it:
find_package(OpenMP)

# Set up compiler flags
if (APPLE)
  set(CMAKE_CXX_FLAGS "${MORPH_HOST_DEFINITION} -Wall -g -O3 -Wno-missing-braces")
else()
  # To use Intel compiler, call cmake as: cmake -DCMAKE_CXX_COMPILER=icpc ..
  if (CMAKE_CXX_COMPILER_ID MATCHES Intel)
    set(CMAKE_CXX_FLAGS "${MORPH_HOST_DEFINITION} -Wall -g -std=c++17 -xHOST -O3 -D__ICC__")
  else() # GCC or Clang
    set(CMAKE_CXX_FLAGS "${MORPH_HOST_DEFINITION} -Wall -g -Wfatal-errors -Wno-unused-result -Wno-unknown-pragmas -march=native -O3")
  endif()
endif()

if(OpenMP_FOUND)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

# Additional GL compiler flags.
#
# Following `cmake --help-policy CMP0072`
set(OpenGL_GL_PREFERENCE "GLVND")
# maybe -DGL_GLEXT_PROTOTYPES
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DGL3_PROTOTYPES -DGL_GLEXT_PROTOTYPES")

if(APPLE)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DGL_SILENCE_DEPRECATION")
endif()

# Lib finding - we have to do this, so that the correct include paths
# are present for OpenGL, OpenCV, etc; they're all included by
# morphologica header files.
find_package(OpenCV REQUIRED)
find_package(OpenGL REQUIRED)
find_package(HDF5 REQUIRED)
find_package(Armadillo REQUIRED)
find_package(jsoncpp REQUIRED)
find_package(glfw3 REQUIRED)
if(OpenMP_FOUND)
  link_libraries(OpenMP::OpenMP_CXX)
endif()

include_directories(${OpenCV_INCLUDE_DIRS} ${OPENGL_INCLUDE_DIR} ${GLUT_INCLUDE_DIR} ${X11_INC_SEARCH_PATH} ${HDF5_INCLUDE_DIR})

# The include path for morphologica defaults to in-tree, but user can change it if necessary
set(MORPH_INCLUDE_PATH "${PROJECT_SOURCE_DIR}/morphologica" CACHE PATH "The path to the morphologica headers (default: here, in the RetinoTectal directory)")

# Show message saying where morphologica is
message(STATUS "morphologica path: ${MORPH_INCLUDE_PATH}")
message(STATUS "  (This can be changed with, e.g.: `cmake -DMORPH_INCLUDE_PATH=/usr/local/include`")

# Use our include path for the morph directories
include_directories(${MORPH_INCLUDE_PATH})

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
