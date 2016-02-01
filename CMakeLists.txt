cmake_minimum_required (VERSION 2.6)
project (ococo)

include(ExternalProject)

option(DEBUGGING_MODE "Build in the debugging mode" OFF)
set(DEBUGGING_SEVERITY "trace" CACHE STRING  "Verbosity of debugging mode")
#
# FIX!!!
#

###############
# COMPILATION #
###############

if(DEBUGGING_MODE)
	find_package( Boost 1.46 COMPONENTS program_options system log log_setup thread REQUIRED )
	add_definitions(-DDEBUGGING_MODE)
else(DEBUGGING_MODE)
	find_package( Boost 1.46 COMPONENTS program_options REQUIRED )
endif(DEBUGGING_MODE)

include_directories( ${Boost_INCLUDE_DIR} )

include("contrib/htslib.cmake")
include("contrib/zlib.cmake")
include("contrib/googletest.cmake")

add_subdirectory(tests/unit_tests)

#list(APPEND CMAKE_CXX_FLAGS "-std=c++11")


#list( APPEND CMAKE_CXX_FLAGS "-std=gnu++11  -O3 -DNDEBUG ${CMAKE_CXX_FLAGS}")
#list(APPEND CMAKE_CXX_FLAGS "-Werrror -Wshadow -std=c++11")
#add_definitions(-D_FILE_OFFSET_BITS=64 -g)


if(DEBUGGING_MODE)
	list(APPEND CMAKE_CXX_FLAGS "-std=c++11 -Wall -Wextra -pedantic -Wshadow")	
else(DEBUGGING_MODE)
	list(APPEND CMAKE_CXX_FLAGS "-std=c++11 -DNDEBUG -march=native -O3")
endif(DEBUGGING_MODE)


add_executable(ococo main.cpp)
add_dependencies(ococo htslib zlib)
target_link_libraries(ococo ${Boost_LIBRARIES} ${htslib_LIB} ${zlib_LIB})


################
# INSTALLATION #
################

install(TARGETS ococo DESTINATION bin)


#########
# TESTS #
#########

enable_testing()
add_test(NAME unit_tests COMMAND unittests)
#add_test(NAME integration_tests COMMAND "cd tests/integration_test/ && ./integration_tests.bats")
