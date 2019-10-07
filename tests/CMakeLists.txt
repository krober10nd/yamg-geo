PROJECT(YAMG_TESTS)
cmake_minimum_required(VERSION 3.0)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/tests/bin)

set(YAMG_TESTS_SRC ${CMAKE_SOURCE_DIR}/tests/src)

enable_language(CXX)

FIND_PACKAGE(VTK NO_MODULE)
if(VTK_FOUND)
  message(STATUS "Found VTK: ${VTK_DIR} (found version \"${VTK_VERSION}\")")

  include(${VTK_USE_FILE})
  set (YAMG_LIBRARIES ${VTK_LIBRARIES} ${YAMG_LIBRARIES})
  add_definitions(-DVTK_FOUND)
else()
  message(SEND_ERROR, "Cannot find VTK - required.")
endif()

include_directories(include ${CMAKE_SOURCE_DIR}/include)

# test volume 
ADD_EXECUTABLE(test_Volume ${YAMG_TESTS_SRC}/test_Volume.cpp)

TARGET_LINK_LIBRARIES(test_Volume ${YAMG_LIBRARIES} yamg)

