cmake_minimum_required (VERSION 2.6)
project (ococo)

include(ExternalProject)

FIND_PACKAGE( Boost 1.46 COMPONENTS program_options system log log_setup thread REQUIRED )
#FIND_PACKAGE( Boost 1.46 COMPONENTS program_options REQUIRED )
INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} )

include("contrib/htslib.cmake")
include("contrib/zlib.cmake")
include("contrib/googletest.cmake")

add_subdirectory(tests)

#list( APPEND CMAKE_CXX_FLAGS "-std=gnu++11  -O3 -DNDEBUG ${CMAKE_CXX_FLAGS}")
list( APPEND CMAKE_CXX_FLAGS "-std=c++11")
add_definitions( -Wall -D_FILE_OFFSET_BITS=64 -g)

add_library (ococo_lib ococo.cpp)

add_executable( ococo main.cpp)
add_dependencies( ococo htslib zlib ococo_lib)
target_link_libraries( ococo ococo_lib ${Boost_LIBRARIES}  ${htslib_LIB} ${zlib_LIB} )
