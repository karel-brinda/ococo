cmake_minimum_required (VERSION 2.6)
project (ococo)


###################
# MODES SWITCHERS #
###################

option(DEBUG "Build in the debugging mode" OFF)
option(OCOCO32 "Increase size of stats per position from 16b to 32b." OFF)
option(VERBOSE_VCF "Print also VCF records without updates." OFF)
option(INSTALL_DEBUG_SCRIPTS "Install debugging scripts." OFF)

# FIX!!!
# (seems not to work correctly)
set(DEBUGGING_SEVERITY "trace" CACHE STRING  "Verbosity of debugging mode")

if(DEBUG)
	find_package( Boost 1.46 COMPONENTS program_options system log log_setup thread REQUIRED )
	add_definitions(-DDEBUGGING_MODE)
	#list(APPEND CMAKE_CXX_FLAGS "-std=c++11 -Wall -Wextra -pedantic -Wshadow")	
else(DEBUG)
	find_package( Boost 1.46 COMPONENTS program_options REQUIRED )
	#list(APPEND CMAKE_CXX_FLAGS "-std=c++11 -DNDEBUG -march=native -O3")
endif(DEBUG)

if(OCOCO32)
	add_definitions(-DOCOCO32)
endif(OCOCO32)

if(VERBOSE_VCF)
	add_definitions(-DVERBOSE_VCF)
endif(VERBOSE_VCF)


######################
# EXTERNAL LIBRARIES #
######################

include(ExternalProject)

include_directories( ${Boost_INCLUDE_DIR} )

include("contrib/htslib.cmake")
include("contrib/zlib.cmake")
include("contrib/googletest.cmake")


################
# MAIN PROGRAM #
################

add_subdirectory(tests/unit_tests)

if(DEBUG)
	list( APPEND CMAKE_CXX_FLAGS "-std=c++11 -Wall -Wextra -pedantic -Wshadow" )
else(DEBUG)
	list( APPEND CMAKE_CXX_FLAGS "-std=c++11 -DNDEBUG -march=native -O3" )
endif(DEBUG)

add_executable(ococo main.cpp)
add_dependencies(ococo htslib zlib)
target_link_libraries(ococo ${Boost_LIBRARIES} ${htslib_LIB} ${zlib_LIB})


################
# INSTALLATION #
################

install(TARGETS ococo DESTINATION bin)

if(INSTALL_DEBUG_SCRIPTS)
	install(PROGRAMS misc/ococo_test_bam.sh DESTINATION bin)
	install(PROGRAMS misc/ococo_test_bam.sh DESTINATION bin)
endif(INSTALL_DEBUG_SCRIPTS)


#########
# TESTS #
#########

enable_testing()
add_test(NAME unit_tests COMMAND unittests)
#add_test(NAME integration_tests COMMAND "cd tests/integration_test/ && ./integration_tests.bats")
