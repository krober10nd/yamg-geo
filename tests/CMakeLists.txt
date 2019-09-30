set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ../bin)
set(YAMG_TEST_SRC ${CMAKE_SOURCE_DIR}/src)

#include_directories(${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/tests/include)

ADD_EXECUTABLE(test_Area ${YAMG_TEST_SRC}/test_Area.cpp)
TARGET_LINK_LIBRARIES(test_Area ${YAMG_LIBRARIES})

