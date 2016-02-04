include_directories(../..)

list( APPEND CMAKE_CXX_FLAGS "-Wall -std=c++11")

add_executable(unittests unittests.cpp)

find_package (Threads REQUIRED)

target_link_libraries(unittests ${googletest_LIB} ${CMAKE_THREAD_LIBS_INIT} ${Boost_LIBRARIES} ${htslib_LIB} ${zlib_LIB})
