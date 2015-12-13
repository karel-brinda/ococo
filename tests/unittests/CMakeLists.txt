include_directories(../..)

list( APPEND CMAKE_CXX_FLAGS "-Wall -std=c++11")

add_executable(tests tests.cpp)
target_link_libraries(tests ${googletest_LIB} ococo_lib)
