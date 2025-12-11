cmake_minimum_required(VERSION 3.9)

# Might be of interest:
# - Build in debug mode:
#   $ cmake -DCMAKE_BUILD_TYPE=Debug ..
#
# - Built a .deb file
#   $ cpack -G DEB
#
# todo:
# - see that the install works
# - create a deb file
# - build lua from here

project(chromflock
  VERSION 0.4.5
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

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/docs/version.rst.in ${CMAKE_CURRENT_SOURCE_DIR}/docs/version.rst)

set(mflock_SOURCES
  src/mflock_cli.c
  src/mflock.c
  src/functional.c
  src/cmmwrite.c
  src/wio.c
  src/hsvrgb.c
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

set(chromflock_SOURCES
  src/chromflock.c
  src/cc2cpm.c
    src/string2any.c
    src/any2string.c
    src/oscp.c
    src/sprite2cmap.c
    src/chromflock_init.c
    src/balance.c
    src/cf_util.c
    src/contact_pairs_io.c
  )

add_executable( mflock ${mflock_SOURCES})
add_executable( aflock ${aflock_SOURCES})
add_executable( chromflock ${chromflock_SOURCES})

#
# Handle large files
#

add_definitions(-D_FILE_OFFSET_BITS=64)

#
# OpenMP for mflock
#

find_package(OpenMP)
if(OpenMP_C_FOUND)
  message("Found OpenMP")
  target_link_libraries(mflock OpenMP::OpenMP_C)
  target_link_libraries(aflock OpenMP::OpenMP_C)
else()
  message(STATUS "OpenMP not found automatically, trying manual")
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

#
# SDL2
#
find_package(SDL2)
if(SDL2_FOUND)
 message(STATUS "SDL2 found")
target_link_libraries(
  mflock
  SDL2::SDL2
  SDL2::SDL2main
)
add_definitions(-DSDL)
target_sources(mflock PRIVATE src/liveview.c)
else()
  message(STATUS "SDL2 not found, some features will be disabled")
endif()

#
# LUA 5.3.5 -- note: has to be built separately
#
include_directories(src/lua-5.3.5/src/)
target_link_libraries(mflock
  ${CMAKE_CURRENT_SOURCE_DIR}/src/lua-5.3.5/src/liblua.a
)


# Enable -flto (GCC) and similar if available, see
# https://cmake.org/cmake/help/latest/policy/CMP0069.html
# https://cmake.org/cmake/help/latest/module/CheckIPOSupported.html
include(CheckIPOSupported)
check_ipo_supported(RESULT result OUTPUT output)
if(result)
  message("Enabling IPO")
  set_property(TARGET aflock PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
  set_property(TARGET mflock PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
  set_property(TARGET chromflock PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
  set(CMAKE_INTERPROCEDURAL_OPTIMIZATION_JOBS 8)
else()
  message(WARNING "IPO is not supported: ${output}")
endif()



# Note the value of CMAKE_INSTALL_PREFIX
include(GNUInstallDirs)
install(TARGETS mflock aflock chromflock)

#
# .deb package
#
# Ref:
# https://cmake.org/cmake/help/book/mastering-cmake/chapter/Packaging%20With%20CPack.html

SET(CPACK_DEBIAN_PACKAGE_MAINTAINER "Erik Wernersson")
# resolve dependencies automatically
set(CPACK_DEBIAN_PACKAGE_SHLIBDEPS ON)
set(CPACK_DEBIAN_PACKAGE_HOMEPAGE "https://www.github.com/elgw/deconwolf/")
SET(CPACK_DEBIAN_PACKAGE_SECTION "science")
INCLUDE(CPack)
