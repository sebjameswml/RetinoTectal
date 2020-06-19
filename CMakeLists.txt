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
find_package(GLUT REQUIRED)
find_package(X11 REQUIRED)
find_package(HDF5 REQUIRED)

# glfw3 library finding. However, it's possible to compile from
# source, and there are instructions on the glfw3 website on how to do
# this with cmake.
#
# glfw3 provides a cmake find_package script, though it doesn't set
# lib paths etc. Then you just link with "glfw3" if you happened to
# compile the static libglfw3.a or with "glfw" if you compiled the
# dynamic libglfw.so
find_package(glfw3 REQUIRED)
# You can find GLFW library with pkgconfig, but it's switched off for now
if(0)
  find_package(PkgConfig)
  if(PKG_CONFIG_FOUND)
    # Use pkg-config to check GLFW3 is present:
    pkg_check_modules(GLFW3 REQUIRED glfw3)
    if (GLFW3_FOUND)
      # Good....
      message(INFO, "pkg-config found glfw3. Excellent.")
    else(GLFW3_FOUND)
      # Bad....
      message(ERROR, "glfw3 was NOT found by pkg-config. Uh oh.")
    endif(GLFW3_FOUND)
  else()
    message(WARNING "There's no pkg-config on this system to check for libglfw3. You may need to try `sudo apt-get install pkg-config`")
  endif()
endif()

include_directories(${OpenCV_INCLUDE_DIRS} ${OPENGL_INCLUDE_DIR} ${GLUT_INCLUDE_DIR} ${X11_INC_SEARCH_PATH} ${HDF5_INCLUDE_DIR})

# Find libmorphologica using pkgconfig
find_package(PkgConfig)
if(PKG_CONFIG_FOUND)
  message(STATUS "Have pkg_config, searching for libmorphologica...")

  if(APPLE)
    # Manually set PKG_CONFIG_PATH on Mac, as by default, mac ports will search only in /opt/local
    set(ENV{PKG_CONFIG_PATH}  "ENV{PKG_CONFIG_PATH}:${CMAKE_INSTALL_PREFIX}/lib/pkgconfig")
  endif(APPLE)

  pkg_check_modules(MORPH REQUIRED libmorphologica)
  if (MORPH_FOUND)
    message(STATUS "pkg_config found libmorphologica")
    include_directories(${MORPH_INCLUDEDIR})
    include_directories(${MORPH_INCLUDEDIR}/morph)
    message(STATUS "Included morph includedir:  ${MORPH_INCLUDEDIR} and ${MORPH_INCLUDEDIR}/morph")
  else()
    message(FATAL_ERROR "You need libmorphologica from github.com/ABRG_Models/morphologica")
  endif(MORPH_FOUND)
else()
  message(WARNING "There's no pkg-config on this system to check for libmorphologica. "
    "You can set MORPH_LIB with `cmake -DMORPH_LIB=/pathto/lib` ...")
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
