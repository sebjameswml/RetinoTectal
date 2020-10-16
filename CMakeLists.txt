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
#find_package(OpenMP)
#if(OpenMP_CXX_FOUND)
#    target_link_libraries(MyTarget PUBLIC OpenMP::OpenMP_CXX)
#endif()
#
if (APPLE)
  set(CMAKE_CXX_FLAGS "${MORPH_HOST_DEFINITION} -Wall -g -Xpreprocessor -fopenmp")
else()
  if (CMAKE_CXX_COMPILER_ID MATCHES Intel) # Intel compiler
    set(CMAKE_CXX_FLAGS "${MORPH_HOST_DEFINITION} -Wall -g -Wfatal-errors -xHOST -O3 -qopenmp -D__ICC__")
  else() # GCC
    set(CMAKE_CXX_FLAGS "${MORPH_HOST_DEFINITION} -Wall -g -Wfatal-errors -Wno-unused-result -Wno-unknown-pragmas -march=native -O3 -fopenmp")
  endif()
endif()

# Additional GL compiler flags.
#
# Following `cmake --help-policy CMP0072`
set(OpenGL_GL_PREFERENCE "GLVND")
# maybe -DGL_GLEXT_PROTOTYPES
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DGL3_PROTOTYPES -DGL_GLEXT_PROTOTYPES")

# Lib finding - we have to do this, so that the correct include paths
# are present for OpenGL, OpenCV, etc; they're all included by
# morphologica header files.
find_package(OpenCV REQUIRED)
find_package(OpenGL REQUIRED)
#find_package(GLUT REQUIRED)
#find_package(X11 REQUIRED)
find_package(HDF5 REQUIRED)
find_package(Armadillo REQUIRED)
find_package(jsoncpp REQUIRED)

# glfw3 library finding. However, it's possible to compile from
# source, and there are instructions on the glfw3 website on how to do
# this with cmake.
#
# glfw3 provides a cmake find_package script, though it doesn't set
# lib paths etc. Then you just link with "glfw3" if you happened to
# compile the static libglfw3.a or with "glfw" if you compiled the
# dynamic libglfw.so
find_package(glfw3 REQUIRED)

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
