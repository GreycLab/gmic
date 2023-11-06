# The CMake modules shipped with OpenCV do not contain targets.
find_package(PkgConfig)

pkg_check_modules(OPENCV opencv4)

if(NOT OPENCV_FOUND)
  pkg_check_modules(OPENCV opencv)
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(OpenCV
  REQUIRED_VARS OPENCV_LIBRARIES
)

# OpenCV 4 can contain broken path in .pc file so we need to filter it.
# https://github.com/opencv/opencv/pull/17377
set(OPENCV_VALID_INCLUDE_DIRS)
foreach(dir ${OPENCV_INCLUDE_DIRS})
  if(EXISTS ${dir})
    list(APPEND OPENCV_VALID_INCLUDE_DIRS ${dir})
  endif()
endforeach()

add_library(OpenCV::OpenCV INTERFACE IMPORTED)
set_target_properties(OpenCV::OpenCV PROPERTIES
  INTERFACE_INCLUDE_DIRECTORIES "${OPENCV_VALID_INCLUDE_DIRS}"
  INTERFACE_LINK_LIBRARIES "${OPENCV_LIBRARIES}"
)
