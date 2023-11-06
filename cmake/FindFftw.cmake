find_package(PkgConfig)

pkg_check_modules(FFTW3 fftw3>=3.0)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Fftw
  REQUIRED_VARS FFTW3_LIBRARIES
)

add_library(Fftw::Fftw INTERFACE IMPORTED)
set_target_properties(Fftw::Fftw PROPERTIES
  INTERFACE_INCLUDE_DIRECTORIES "${FFTW3_INCLUDE_DIRS}"
  INTERFACE_LINK_LIBRARIES "${FFTW3_LIBRARIES}"
)

find_library(FFTW3_THREADS_LIB fftw3_threads PATHS ${FFTW3_LIBRARY_DIRS})
if(NOT FFTW3_THREADS_LIB STREQUAL "FFTW3_THREADS_LIB-NOTFOUND")
  add_library(Fftw::Threads INTERFACE IMPORTED)
  set_target_properties(Fftw::Threads PROPERTIES
    INTERFACE_LINK_LIBRARIES "-lfftw3_threads"
  )
endif()
