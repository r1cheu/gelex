# Locate the LAPACKE C interface (header + linkable symbols) and expose it as
# the imported target LAPACKE::LAPACKE.
#
# Handles the common conda/distro layouts where lapacke.h lives under a
# backend-specific subdirectory, and the OpenBLAS case where the LAPACKE symbols
# are bundled into the LAPACK library instead of a standalone liblapacke.
# Requires LAPACK::LAPACK (find_package(LAPACK)) for that fallback.

include(FindPackageHandleStandardArgs)
include(CheckCXXSymbolExists)

find_path(
  LAPACKE_INCLUDE_DIR
  NAMES lapacke.h
  PATH_SUFFIXES "" lapacke openblas include/lapacke)

find_library(LAPACKE_LIBRARY NAMES lapacke)

if(NOT LAPACKE_LIBRARY
   AND TARGET LAPACK::LAPACK
   AND LAPACKE_INCLUDE_DIR)
  set(CMAKE_REQUIRED_LIBRARIES LAPACK::LAPACK)
  set(CMAKE_REQUIRED_INCLUDES ${LAPACKE_INCLUDE_DIR})
  check_cxx_symbol_exists(LAPACKE_dpotrf lapacke.h LAPACKE_EMBEDDED_IN_LAPACK)
  unset(CMAKE_REQUIRED_LIBRARIES)
  unset(CMAKE_REQUIRED_INCLUDES)
endif()

find_package_handle_standard_args(
  LAPACKE
  REQUIRED_VARS LAPACKE_INCLUDE_DIR
  HANDLE_COMPONENTS)

if(LAPACKE_FOUND AND NOT TARGET LAPACKE::LAPACKE)
  add_library(LAPACKE::LAPACKE INTERFACE IMPORTED)
  set_property(TARGET LAPACKE::LAPACKE PROPERTY INTERFACE_INCLUDE_DIRECTORIES
                                                "${LAPACKE_INCLUDE_DIR}")
  if(LAPACKE_LIBRARY)
    set_property(TARGET LAPACKE::LAPACKE PROPERTY INTERFACE_LINK_LIBRARIES
                                                  "${LAPACKE_LIBRARY}")
  elseif(TARGET LAPACK::LAPACK)
    set_property(TARGET LAPACKE::LAPACKE PROPERTY INTERFACE_LINK_LIBRARIES
                                                  LAPACK::LAPACK)
  endif()
endif()

mark_as_advanced(LAPACKE_INCLUDE_DIR LAPACKE_LIBRARY)
