# - Check for the presence of P4EST
#
# The following variables are set when P4EST is found:
#  HAVE_P4EST       = Set to true, if all components of P4EST
#                          have been found.
#  P4EST_INCLUDES   = Include path for the header files of P4EST
#  P4EST_LIBRARIES  = Link these to use P4EST

## -----------------------------------------------------------------------------
## Check for the header files

find_path (P4EST_INCLUDES 
    NAMES "p4est.h"
    PATHS /usr/local/include /usr/include $ENV{P4EST_ROOT}/local/include/
  )

## -----------------------------------------------------------------------------
## Check for the library

find_library (P4EST_LIBRARIES 
    NAMES "p4est"
    PATHS /usr/local/lib /usr/lib /lib $ENV{P4EST_ROOT}/local/lib/
  )

## -----------------------------------------------------------------------------
## Actions taken when all components have been found

if (P4EST_INCLUDES AND P4EST_LIBRARIES)
  set (HAVE_P4EST TRUE)
else (P4EST_INCLUDES AND P4EST_LIBRARIES)
  if (NOT P4EST_FIND_QUIETLY)
    if (NOT P4EST_INCLUDES)
      message (STATUS "Unable to find P4EST header files!")
    endif (NOT P4EST_INCLUDES)
    if (NOT P4EST_LIBRARIES)
      message (STATUS "Unable to find P4EST library files!")
    endif (NOT P4EST_LIBRARIES)
  endif (NOT P4EST_FIND_QUIETLY)
endif (P4EST_INCLUDES AND P4EST_LIBRARIES)

if (HAVE_P4EST)
  if (NOT P4EST_FIND_QUIETLY)
    message (STATUS "Found components for P4EST")
    message (STATUS "P4EST_INCLUDES = ${P4EST_INCLUDES}")
    message (STATUS "P4EST_LIBRARIES     = ${P4EST_LIBRARIES}")
  endif (NOT P4EST_FIND_QUIETLY)
else (HAVE_P4EST)
  if (P4EST_FIND_REQUIRED)
    message (FATAL_ERROR "Could not find P4EST!")
  endif (P4EST_FIND_REQUIRED)
endif (HAVE_P4EST)

mark_as_advanced (
  HAVE_P4EST
  P4EST_LIBRARIES
  P4EST_INCLUDES
  )


