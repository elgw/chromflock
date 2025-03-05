cmake_minimum_required(VERSION 3.9)

# Might be of interest:
# cmake -DCMAKE_BUILD_TYPE=Debug ..
# todo:
# - see that the install works
# - create a deb file
# - build lua from here

project(chromflock
  VERSION 0.4.4
  LANGUAGES C)

include(ExternalProject)

set (CMAKE_C_STANDARD 11)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# Don't install to /usr/local/ use /usr/
set( CMAKE_INSTALL_PREFIX
  /usr )

# Replace variables by their values
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/src/cf_version.h.in ${CMAKE_CURRENT_SOURCE_DIR}/src/cf_version.h)

set(mflock_SOURCES
  src/mflock_cli.c
  src/mflock.c
  src/functional.c
  src/cmmwrite.c
  src/wio.c
  src/hsvrgb.c
  src/liveview.c
  src/ellipsoid.c
  src/contact_pairs_io.c
  src/cf_util.c)

set(aflock_SOURCES
  src/aflock.c
  src/wio.c
  src/contact_pairs_io.c
  src/cf_util.c
  src/ellipsoid.c
  src/oscp.c
)

set(chromflock_SOURCES,
  src/chromflock.c)

include_directories(src/lua-5.3.5/src/)

add_executable( mflock ${mflock_SOURCES})
add_executable( aflock ${aflock_SOURCES})
add_executable( chromflock ${aflock_SOURCES})

# TODO add_custom_command for building lua

#
# OpenMP for mflock (TODO: This is optional)
#

find_package(OpenMP)
if(OpenMP_C_FOUND)
  message("Found OpenMP")
  target_link_libraries(mflock OpenMP::OpenMP_C)
  target_link_libraries(aflock OpenMP::OpenMP_C)
else()
  message("OpenMP not found automatically, trying manual")

  set(
    CMAKE_C_FLAGS
    "${CMAKE_C_FLAGS} -fopenmp=omp /openmp"
  )
  SET(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -fopenmp=omp -lomp")

endif()

#
# Math library, if needed
#

find_library(MATH_LIBRARY m)
if(MATH_LIBRARY)
  target_link_libraries(mflock ${MATH_LIBRARY})
  target_link_libraries(aflock ${MATH_LIBRARY})
  target_link_libraries(chromflock ${MATH_LIBRARY})
endif()

#
# libz
#
find_package(ZLIB REQUIRED)
target_link_libraries(mflock ZLIB::ZLIB)
target_link_libraries(aflock ZLIB::ZLIB)
target_link_libraries(chromflock ZLIB::ZLIB)

find_package(SDL2 REQUIRED)
target_link_libraries(
  mflock
  SDL2::SDL2
  SDL2::SDL2main
)

target_link_libraries(mflock
  ${CMAKE_CURRENT_SOURCE_DIR}/src/lua-5.3.5/src/liblua.a
  )

# Note the value of CMAKE_INSTALL_PREFIX
include(GNUInstallDirs)
install(TARGETS mflock)
