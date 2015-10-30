cmake_minimum_required (VERSION 2.6)
project (ococo)


#FIND_PACKAGE( Boost 1.46 COMPONENTS program_options system log log_setup thread REQUIRED )
FIND_PACKAGE( Boost 1.46 COMPONENTS program_options REQUIRED )
INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} )

list( APPEND CMAKE_CXX_FLAGS "-std=gnu++11  -O3 -DNDEBUG ${CMAKE_CXX_FLAGS}")
add_definitions( -Wall -D_FILE_OFFSET_BITS=64 )

add_executable( ococo ococo.cpp)
target_link_libraries( ococo ${Boost_LIBRARIES})

