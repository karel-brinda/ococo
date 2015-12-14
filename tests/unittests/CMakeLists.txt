include_directories(../..)

list( APPEND CMAKE_CXX_FLAGS "-Wall -std=c++11")

add_executable(tests tests.cpp)

#set(PThreadLib -pthread)
find_package (Threads REQUIRED)

#target_link_libraries(tests ${googletest_LIB} ococo_lib ${PThreadLib} )
target_link_libraries(tests ${googletest_LIB} ococo_lib ${CMAKE_THREAD_LIBS_INIT} )
