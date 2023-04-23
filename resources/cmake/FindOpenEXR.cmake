find_package(PkgConfig)

pkg_check_modules(OPENEXR OpenEXR)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(OpenEXR
  REQUIRED_VARS OPENEXR_LIBRARIES
)

add_library(OpenEXR::OpenEXR INTERFACE IMPORTED)
set_target_properties(OpenEXR::OpenEXR PROPERTIES
  INTERFACE_INCLUDE_DIRECTORIES "${OPENEXR_INCLUDE_DIRS}"
  INTERFACE_LINK_LIBRARIES "${OPENEXR_LIBRARIES}"
)
