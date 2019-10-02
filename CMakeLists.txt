CMAKE_MINIMUM_REQUIRED(VERSION 3.11)
# accept the silly policy to enable the use of _root vars
CMAKE_POLICY(SET CMP0074 NEW)
PROJECT(YAMG)

# The version number.
set (YAMG_VERSION_MAJOR 0)
set (YAMG_VERSION_MINOR 1)

set (CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMake/Modules ${CMAKE_MODULE_PATH})

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

enable_language(CXX)

set (YAMG_LIBRARIES)

find_package(P4EST REQUIRED)
if(HAVE_P4EST)
    message(STATUS "Found P4EST --- ${P4EST_LIBRARIES}")

    include_directories(${P4EST_INCLUDES})

    set (YAMG_LIBRARIES ${P4EST_LIBRARIES} ${YAMG_LIBRARIES})
    add_definitions(-DP4EST_FOUND)
else()
    message(SEND_ERROR, "Cannot find p4est - required.")
endif()

find_package(SC REQUIRED)
if(HAVE_SC)
    message(STATUS "Found SC --- ${SC_LIBRARIES}")

    include_directories(${SC_INCLUDES})

    set (YAMG_LIBRARIES ${SC_LIBRARIES} ${YAMG_LIBRARIES})
    add_definitions(-DSC_FOUND)
else()
    message(SEND_ERROR, "Cannot find SC - required.")
endif()

FIND_PACKAGE(VTK NO_MODULE)
if(VTK_FOUND)
  message(STATUS "Found VTK: ${VTK_DIR} (found version \"${VTK_VERSION}\")")

  include(${VTK_USE_FILE})
  set (YAMG_LIBRARIES ${VTK_LIBRARIES} ${YAMG_LIBRARIES})
  add_definitions(-DVTK_FOUND)
else()
  message(SEND_ERROR, "Cannot find VTK - required.")
endif()

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
else()
        message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()

if(CMAKE_COMPILER_IS_GNUCXX)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-literal-suffix -Wno-deprecated")
endif()
if(${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
  message(STATUS "Adding magic sauce for Intel compiler.")
  SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -xHOST -g")
endif()

find_package(MPI REQUIRED)
if(MPI_FOUND)
  include_directories(${MPI_INCLUDE_PATH})
  set (YAMG_LIBRARIES ${MPI_LIBRARIES} ${YAMG_LIBRARIES})
  add_definitions(-DMPI_FOUND)
endif()

include_directories(include ${CMAKE_BINARY_DIR}/include)

set(CXX_SOURCES ${CMAKE_BINARY_DIR}/Yamg.cpp)

# target_link_library
add_library(yamg SHARED ${CXX_SOURCES})
target_link_libraries(yamg ${YAMG_LIBRARIES})

# add tests 
ADD_SUBDIRECTORY(tests EXCLUDE_FROM_ALL)

ADD_EXECUTABLE(yamg-vp ${CMAKE_BINARY_DIR}/yamg-vp.cpp)
TARGET_LINK_LIBRARIES(yamg-vp ${YAMG_LIBRARIES} yamg)
