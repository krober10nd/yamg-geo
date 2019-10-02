# - Check for the presence of VTK
#
# The following variables are set when VTK is found:
#  HAVE_VTK       = Set to true, if all components of VTK
#                          have been found.
#  VTK_INCLUDES   = Include path for the header files of VTK
#  VTK_LIBRARIES  = Link these to use VTK

## -----------------------------------------------------------------------------
## Check for the header files

find_path (VTK_INCLUDES 
    NAMES "vtkCommon.h"
    PATHS /usr/local/include /usr/include $ENV{VTK_ROOT}/include
  )

## -----------------------------------------------------------------------------
## Check for the library

find_library (VTK_LIBRARIES 
    NAMES "vtkCommon"
    PATHS /usr/local/lib /usr/lib /lib $ENV{VTK_ROOT}/lib
  )

## -----------------------------------------------------------------------------
## Actions taken when all components have been found

if (VTK_INCLUDES AND VTK_LIBRARIES)
  set (HAVE_VTK TRUE)
else (VTK_INCLUDES AND VTK_LIBRARIES)
  if (NOT VTK_FIND_QUIETLY)
    if (NOT VTK_INCLUDES)
      message (STATUS "Unable to find VTK header files!")
    endif (NOT VTK_INCLUDES)
    if (NOT VTK_LIBRARIES)
      message (STATUS "Unable to find VTK library files!")
    endif (NOT VTK_LIBRARIES)
  endif (NOT VTK_FIND_QUIETLY)
endif (VTK_INCLUDES AND VTK_LIBRARIES)

if (HAVE_VTK)
  if (NOT VTK_FIND_QUIETLY)
    message (STATUS "Found components for VTK")
    message (STATUS "VTK_INCLUDES = ${VTK_INCLUDES}")
    message (STATUS "VTK_LIBRARIES     = ${VTK_LIBRARIES}")
  endif (NOT VTK_FIND_QUIETLY)
else (HAVE_VTK)
  if (VTK_FIND_REQUIRED)
    message (FATAL_ERROR "Could not find VTK!")
  endif (VTK_FIND_REQUIRED)
endif (HAVE_VTK)

mark_as_advanced (
  HAVE_VTK
  VTK_LIBRARIES
  VTK_INCLUDES
  )


