# - Check for the presence of SC
#
# The following variables are set when SC is found:
#  HAVE_SC       = Set to true, if all components of SC
#                          have been found.
#  SC_INCLUDES   = Include path for the header files of SC
#  SC_LIBRARIES  = Link these to use SC

## -----------------------------------------------------------------------------
## Check for the header files

find_path (SC_INCLUDES 
    NAMES "sc.h"
    PATHS /usr/local/include /usr/include $ENV{SC_ROOT}/local/include
  )

## -----------------------------------------------------------------------------
## Check for the library

find_library (SC_LIBRARIES 
    NAMES "sc"
    PATHS /usr/local/lib /usr/lib /lib $ENV{SC_ROOT}/local/lib
  )

## -----------------------------------------------------------------------------
## Actions taken when all components have been found

if (SC_INCLUDES AND SC_LIBRARIES)
  set (HAVE_SC TRUE)
else (SC_INCLUDES AND SC_LIBRARIES)
  if (NOT SC_FIND_QUIETLY)
    if (NOT SC_INCLUDES)
      message (STATUS "Unable to find SC header files!")
    endif (NOT SC_INCLUDES)
    if (NOT SC_LIBRARIES)
      message (STATUS "Unable to find SC library files!")
    endif (NOT SC_LIBRARIES)
  endif (NOT SC_FIND_QUIETLY)
endif (SC_INCLUDES AND SC_LIBRARIES)

if (HAVE_SC)
  if (NOT SC_FIND_QUIETLY)
    message (STATUS "Found components for SC")
    message (STATUS "SC_INCLUDES = ${SC_INCLUDES}")
    message (STATUS "SC_LIBRARIES     = ${SC_LIBRARIES}")
  endif (NOT SC_FIND_QUIETLY)
else (HAVE_SC)
  if (SC_FIND_REQUIRED)
    message (FATAL_ERROR "Could not find SC!")
  endif (SC_FIND_REQUIRED)
endif (HAVE_SC)

mark_as_advanced (
  HAVE_SC
  SC_LIBRARIES
  SC_INCLUDES
  )


