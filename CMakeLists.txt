#
#  File        : CMakeLists.txt
#                ( cmake-based build system configuration file )
#
#  Basic usage :
#
#                 - cmake .
#                 - make
#
#                ( https://gmic.eu )
#
#  Author      : Sander Knopper
#                ( https://github.com/saknopper )
#
#                Based on work by Ștefan Talpalaru
#                ( https://github.com/stefantalpalaru )
#
#  Licenses    : This file is 'dual-licensed', you have to choose one
#                 of the two licenses below to apply.
#
#                CeCILL-C
#                The CeCILL-C license is close to the GNU LGPL.
#                ( http://cecill.info/licences/Licence_CeCILL-C_V1-en.html )
#
#            or  CeCILL v2.1
#                The CeCILL license is compatible with the GNU GPL.
#                ( http://cecill.info/licences/Licence_CeCILL_V2.1-en.html )
#
#  This software is governed either by the CeCILL or the CeCILL-C license
#  under French law and abiding by the rules of distribution of free software.
#  You can  use, modify and or redistribute the software under the terms of
#  the CeCILL or CeCILL-C licenses as circulated by CEA, CNRS and INRIA
#  at the following URL: "http://cecill.info".
#
#  As a counterpart to the access to the source code and  rights to copy,
#  modify and redistribute granted by the license, users are provided only
#  with a limited warranty  and the software's author,  the holder of the
#  economic rights,  and the successive licensors  have only  limited
#  liability.
#
#  In this respect, the user's attention is drawn to the risks associated
#  with loading,  using,  modifying and/or developing or reproducing the
#  software by the user in light of its specific status of free software,
#  that may mean  that it is complicated to manipulate,  and  that  also
#  therefore means  that it is reserved for developers  and  experienced
#  professionals having in-depth computer knowledge. Users are therefore
#  encouraged to load and test the software's suitability as regards their
#  requirements in conditions enabling the security of their systems and/or
#  data to be ensured and,  more generally, to use and operate it in the
#  same conditions as regards security.
#
#  The fact that you are presently reading this means that you have had
#  knowledge of the CeCILL and CeCILL-C licenses and that you accept its terms.
#

cmake_minimum_required(VERSION 3.14.0)

project(gmic CXX C)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
  message("Build directory is equal to source directory. Binaries will be put in the src directory.")
  message("")

  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/src/" CACHE FILEPATH "Output directory for all targets." FORCE)
  set(GMIC_BINARIES_PATH ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
else()
  set(GMIC_BINARIES_PATH ${CMAKE_BINARY_DIR})
endif()

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

include(FeatureSummary)
include(GNUInstallDirs)

find_package(CImg REQUIRED)
find_package(GMicStdlib REQUIRED)
find_package(ZLIB REQUIRED)
find_package(PkgConfig REQUIRED)

find_package(PNG)
find_package(JPEG)
find_package(TIFF)
pkg_check_modules(FFTW fftw3)
pkg_check_modules(CURL libcurl)

# X11 uniquement sur les systèmes non-Apple/non-Windows (miroir du Makefile)
if(NOT APPLE AND NOT WIN32)
  pkg_check_modules(X11 x11)
endif()

# options controlling the build process
option(BUILD_LIB "Build the GMIC shared library" ON)
option(BUILD_LIB_STATIC "Build the GMIC static library" ON)
option(BUILD_CLI "Build the CLI interface" ON)
option(BUILD_MAN "Build the manpage" ON)
option(BUILD_BASH_COMPLETION "Build Bash completion" ON)
option(CUSTOM_CFLAGS "Override default compiler optimization flags" OFF)
option(ENABLE_DYNAMIC_LINKING "Dynamically link the binaries to the GMIC shared library" OFF)
option(ENABLE_LTO "Enable -flto (Link Time Optimizer) on gcc and clang" OFF)
option(ENABLE_CHECK_SHARED_IMAGES "Enable extra checking of shared images (adds -Dgmic_check_shared_images, slows down code)" ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# compile flags
set(CMAKE_POSITION_INDEPENDENT_CODE True)

list(APPEND GMIC_CXX_COMPILE_FLAGS
  -Dgmic_is_parallel
  -Dcimg_use_abort
  "-Dcimg_date=\"\""
  "-Dcimg_time=\"\""
  -Dcimg_use_zlib
  -pedantic
)

# -Dcimg_use_vt100: enabled on Unix only (mirrors Makefile behaviour)
if(NOT APPLE AND NOT WIN32)
  list(APPEND GMIC_CXX_COMPILE_FLAGS -Dcimg_use_vt100)
endif()

# -Dgmic_check_shared_images: extra shared-image checking, enabled by default (mirrors Makefile CHECK_SHARED_IMAGES_CFLAGS)
if(ENABLE_CHECK_SHARED_IMAGES)
  list(APPEND GMIC_CXX_COMPILE_FLAGS -Dgmic_check_shared_images)
endif()

if(APPLE)
    list(APPEND GMIC_CXX_COMPILE_FLAGS -mmacosx-version-min=10.8 -stdlib=libc++ -Wno-error=c++11-narrowing -Wc++11-extensions -fpermissive)
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    list(APPEND GMIC_CXX_COMPILE_FLAGS -Wno-error=c++11-narrowing -fpermissive)
else()
    list(APPEND GMIC_CXX_COMPILE_FLAGS -Wno-error=narrowing -fpermissive)
endif()

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  list(APPEND GMIC_CXX_COMPILE_FLAGS
    -Wall
    -Wextra
    -Wfatal-errors
    -Werror=unknown-pragmas
    -Werror=unused-label
    -Wshadow
  )
endif()

if(NOT "${PRERELEASE_TAG}" STREQUAL "")
    list(APPEND GMIC_CXX_COMPILE_FLAGS "-Dgmic_prerelease=\"${PRERELEASE_TAG}\"")
endif()

if (ENABLE_LTO)
  # https://stackoverflow.com/a/47370726/160386
  include(CheckIPOSupported)
  check_ipo_supported(RESULT supported OUTPUT error)
  if(supported)
    message(STATUS "IPO / LTO enabled")
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
  else()
    message(STATUS "IPO / LTO not supported: <${error}>")
  endif()
endif()

if(ENABLE_DYNAMIC_LINKING)
  if(NOT BUILD_LIB)
    message(FATAL_ERROR "ENABLE_DYNAMIC_LINKING needs BUILD_LIB")
  endif()
  set(CMAKE_SKIP_RPATH TRUE)
endif()

string(TOLOWER "${CMAKE_BUILD_TYPE}" CMAKE_BUILD_TYPE_LOWER)
if (${CMAKE_BUILD_TYPE_LOWER} STREQUAL "debug")
    list(PREPEND GMIC_CXX_COMPILE_FLAGS -g -Wall -Wextra -pedantic -Dcimg_verbosity=3)
endif()
if (${CMAKE_BUILD_TYPE_LOWER} STREQUAL "relwithdebinfo")
    list(PREPEND GMIC_CXX_COMPILE_FLAGS -g)
endif()

if (NOT CUSTOM_CFLAGS)
    if (${CMAKE_BUILD_TYPE_LOWER} STREQUAL "debug")
        list(PREPEND GMIC_CXX_COMPILE_FLAGS -Og)
    else()
        list(PREPEND GMIC_CXX_COMPILE_FLAGS -O3 -mtune=generic)
    endif()
endif()

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
  message(STATUS "OpenMP found — enabling parallel processing (-Dcimg_use_openmp)")
  list(APPEND GMIC_CXX_COMPILE_FLAGS -Dcimg_use_openmp)
else()
  message(STATUS "OpenMP not found — parallel processing disabled")
endif()

set(GMIC_OPTIONAL_COMPILE_FLAGS "")
set(GMIC_OPTIONAL_LIBS "")

if(PNG_FOUND)
  message(STATUS "libpng found — enabling PNG support (-Dcimg_use_png)")
  list(APPEND GMIC_OPTIONAL_COMPILE_FLAGS -Dcimg_use_png)
  list(APPEND GMIC_OPTIONAL_LIBS PNG::PNG)
endif()

if(JPEG_FOUND)
  message(STATUS "libjpeg found — enabling JPEG support (-Dcimg_use_jpeg)")
  list(APPEND GMIC_OPTIONAL_COMPILE_FLAGS -Dcimg_use_jpeg)
  list(APPEND GMIC_OPTIONAL_LIBS JPEG::JPEG)
endif()

if(TIFF_FOUND)
  message(STATUS "libtiff found — enabling TIFF support (-Dcimg_use_tiff)")
  list(APPEND GMIC_OPTIONAL_COMPILE_FLAGS -Dcimg_use_tiff)
  list(APPEND GMIC_OPTIONAL_LIBS TIFF::TIFF)
endif()

if(FFTW_FOUND)
  message(STATUS "fftw3 found — enabling FFT support (-Dcimg_use_fftw3)")
  list(APPEND GMIC_OPTIONAL_COMPILE_FLAGS -Dcimg_use_fftw3)
  find_library(FFTW3_THREADS_LIB fftw3_threads)
  if(FFTW3_THREADS_LIB)
    list(APPEND GMIC_OPTIONAL_LIBS ${FFTW_LIBRARIES} ${FFTW3_THREADS_LIB})
  else()
    list(APPEND GMIC_OPTIONAL_LIBS ${FFTW_LIBRARIES})
  endif()
  list(APPEND GMIC_OPTIONAL_INCLUDE_DIRS ${FFTW_INCLUDE_DIRS})
endif()

if(CURL_FOUND)
  message(STATUS "libcurl found — enabling network download support (-Dcimg_use_curl)")
  list(APPEND GMIC_OPTIONAL_COMPILE_FLAGS -Dcimg_use_curl)
  list(APPEND GMIC_OPTIONAL_LIBS ${CURL_LIBRARIES})
  list(APPEND GMIC_OPTIONAL_INCLUDE_DIRS ${CURL_INCLUDE_DIRS})
endif()

if(X11_FOUND)
  message(STATUS "X11 found — enabling display support (-Dcimg_display=1)")
  list(APPEND GMIC_OPTIONAL_COMPILE_FLAGS -Dcimg_display=1)
  list(APPEND GMIC_OPTIONAL_LIBS ${X11_LIBRARIES})
  list(APPEND GMIC_OPTIONAL_INCLUDE_DIRS ${X11_INCLUDE_DIRS})
endif()

# source files
set(GMIC_LIB_SOURCES src/gmic.cpp)

function(gmic_apply_common_deps target_name)
  target_link_libraries(${target_name} PRIVATE ZLIB::ZLIB)

  if(OpenMP_CXX_FOUND)
    target_link_libraries(${target_name} PRIVATE OpenMP::OpenMP_CXX)
  endif()

  if(GMIC_OPTIONAL_LIBS)
    target_link_libraries(${target_name} PRIVATE ${GMIC_OPTIONAL_LIBS})
  endif()

  if(GMIC_OPTIONAL_COMPILE_FLAGS)
    target_compile_options(${target_name} PRIVATE ${GMIC_OPTIONAL_COMPILE_FLAGS})
  endif()

  if(GMIC_OPTIONAL_INCLUDE_DIRS)
    target_include_directories(${target_name} PRIVATE ${GMIC_OPTIONAL_INCLUDE_DIRS})
  endif()

endfunction()

if(BUILD_LIB)
  add_library(libgmic SHARED ${GMIC_LIB_SOURCES})
  target_compile_options(libgmic PRIVATE ${GMIC_CXX_COMPILE_FLAGS} -Dgmic_core)
  set_target_properties(libgmic PROPERTIES SOVERSION "1" OUTPUT_NAME "gmic")
  target_link_libraries(libgmic
    PRIVATE
      CImg::CImg
      GMicStdlib::Stdlib
  )
  target_include_directories(libgmic
    PUBLIC
      $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
  )

  gmic_apply_common_deps(libgmic)

  install(TARGETS libgmic EXPORT GmicTargets
    PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
    LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
    ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  )
endif()


if(BUILD_LIB_STATIC)
  add_library(libgmicstatic STATIC ${GMIC_LIB_SOURCES})
  target_compile_options(libgmicstatic PRIVATE ${GMIC_CXX_COMPILE_FLAGS} -Dgmic_core)
  set_target_properties(libgmicstatic PROPERTIES OUTPUT_NAME "gmic")
  target_link_libraries(libgmicstatic
    PRIVATE
      CImg::CImg
      GMicStdlib::Stdlib
  )
  target_include_directories(libgmicstatic
    PUBLIC
      $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
  )

  gmic_apply_common_deps(libgmicstatic)

  install(TARGETS libgmicstatic EXPORT GmicTargets ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}")
endif()

if(BUILD_LIB OR BUILD_LIB_STATIC)
  install(FILES src/gmic.h DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}")
endif()

if(BUILD_CLI)
  if(NOT ENABLE_DYNAMIC_LINKING AND NOT BUILD_LIB_STATIC)
    message(FATAL_ERROR "BUILD_CLI requires either BUILD_LIB_STATIC=ON or ENABLE_DYNAMIC_LINKING=ON")
  endif()
  add_executable(gmic src/gmic_cli.cpp)
  target_compile_options(gmic PRIVATE ${GMIC_CXX_COMPILE_FLAGS})
  if(ENABLE_DYNAMIC_LINKING)
    target_link_libraries(gmic
      PUBLIC
        libgmic
      PRIVATE
        CImg::CImg
    )
  else()
    target_link_libraries(gmic PRIVATE libgmicstatic)
  endif()

  gmic_apply_common_deps(gmic)

  if(NOT ${CMAKE_BUILD_TYPE_LOWER} STREQUAL "debug")
    find_program(STRIP_TOOL strip)
    if(STRIP_TOOL)
      add_custom_command(TARGET gmic POST_BUILD
        COMMAND ${STRIP_TOOL} $<TARGET_FILE:gmic>
        COMMENT "Stripping gmic binary (like 'strip gmic' in the Makefile)"
      )
    endif()
  endif()

  install(TARGETS gmic
    RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
    LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  )
endif()

if(BUILD_MAN OR BUILD_BASH_COMPLETION)
  if(NOT BUILD_CLI)
    message(FATAL_ERROR "BUILD_MAN and BUILD_BASH_COMPLETION require BUILD_CLI=ON")
  endif()
endif()

if(BUILD_MAN)
  file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/man)
  add_custom_command(
    OUTPUT ${CMAKE_BINARY_DIR}/man/gmic.1
    DEPENDS gmic
    COMMAND ${CMAKE_COMMAND} -E env
      "LD_LIBRARY_PATH=${GMIC_BINARIES_PATH}"
      "DYLD_LIBRARY_PATH=${GMIC_BINARIES_PATH}"
      ${GMIC_BINARIES_PATH}/gmic
        ${CMAKE_SOURCE_DIR}/src/gmic_stdlib.gmic
        it ${CMAKE_SOURCE_DIR}/src/gmic_stdlib.gmic
        reference man > ${CMAKE_BINARY_DIR}/man/gmic.1
  )
  add_custom_target(man ALL DEPENDS ${CMAKE_BINARY_DIR}/man/gmic.1)
  install(FILES ${CMAKE_BINARY_DIR}/man/gmic.1 DESTINATION ${CMAKE_INSTALL_MANDIR}/man1/)
endif()

if(BUILD_BASH_COMPLETION)
  file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/resources)
  add_custom_command(
    OUTPUT ${CMAKE_BINARY_DIR}/resources/gmic_bashcompletion.sh
    DEPENDS gmic
    COMMAND ${CMAKE_COMMAND} -E env
      "LD_LIBRARY_PATH=${GMIC_BINARIES_PATH}"
      "DYLD_LIBRARY_PATH=${GMIC_BINARIES_PATH}"
      ${GMIC_BINARIES_PATH}/gmic
        ${CMAKE_SOURCE_DIR}/src/gmic_stdlib.gmic
        it ${CMAKE_SOURCE_DIR}/src/gmic_stdlib.gmic
        parse_cli bashcompletion
        > ${CMAKE_BINARY_DIR}/resources/gmic_bashcompletion.sh
  )
  add_custom_target(bashcompletion ALL DEPENDS ${CMAKE_BINARY_DIR}/resources/gmic_bashcompletion.sh)
  install(FILES ${CMAKE_BINARY_DIR}/resources/gmic_bashcompletion.sh
    DESTINATION ${CMAKE_INSTALL_DATADIR}/bash-completion/completions
    RENAME gmic
  )
endif()

include(CMakePackageConfigHelpers)
file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/GmicConfig.cmake.in" "@PACKAGE_INIT@\ninclude(\${CMAKE_CURRENT_LIST_DIR}/GmicTargets.cmake)\n")
configure_package_config_file(
  ${CMAKE_CURRENT_BINARY_DIR}/GmicConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/GmicConfig.cmake
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/gmic
)
install(EXPORT GmicTargets DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/gmic)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/GmicConfig.cmake DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/gmic)
feature_summary(WHAT ALL)
