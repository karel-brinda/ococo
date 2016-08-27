cmake_minimum_required (VERSION 2.6)
project (ococo C CXX)


##############################
# COMPILATION MODE SWITCHERS #
##############################

if(NOT CMAKE_BUILD_TYPE)
	message(STATUS "Setting build type to 'RELEASE' as none was specified.")
	set(CMAKE_BUILD_TYPE RELEASE CACHE STRING "Choose the type of build." FORCE)
	set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "DEBUG" "RELEASE")
endif()

if(CMAKE_BUILD_TYPE MATCHES DEBUG)
	message("DEBUG mode")
endif(CMAKE_BUILD_TYPE MATCHES DEBUG)

if(CMAKE_BUILD_TYPE MATCHES RELEASE)
	message("RELEASE mode")
endif(CMAKE_BUILD_TYPE MATCHES RELEASE)

option(INSTALL_DEBUG_SCRIPTS "Install debugging scripts." OFF)
set(DEBUGGING_SEVERITY "trace" CACHE STRING  "Verbosity of debugging mode.")

option(BUILD_TESTS "Build tests." OFF)


##########
# HTSLIB #
##########

if (CMAKE_GENERATOR STREQUAL "Unix Makefiles")
	set(MAKE_COMMAND "$(MAKE)")
else()
	find_program(MAKE_COMMAND NAMES make gmake)
endif()

include(ExternalProject)
ExternalProject_Add(htslib
	PREFIX ${CMAKE_BINARY_DIR}/ext/htslib.tmp
	SOURCE_DIR "${CMAKE_BINARY_DIR}/ext/htslib"
	BUILD_IN_SOURCE 1
	UPDATE_COMMAND ""
	CONFIGURE_COMMAND ""
	BUILD_COMMAND ${MAKE_COMMAND} lib-static
	INSTALL_COMMAND ""
	)

include_directories(${CMAKE_BINARY_DIR}/ext/htslib)


########
# ZLIB #
########

find_package(ZLIB REQUIRED)
if (ZLIB_FOUND)
	include_directories(${ZLIB_INCLUDE_DIRS})
else()
	message (FATAL_ERROR "zlib not found.")
endif(ZLIB_FOUND)

include_directories(${htslib_INSTALL}/include)

################
# MAIN PROGRAM #
################

include_directories( "src" )
add_subdirectory(tests/unit_tests)

set(CMAKE_CXX_FLAGS         "${CMAKE_CXX_FLAGS}         -std=c++11 -Wall -Wextra -Wno-missing-field-initializers -Wshadow -g ")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -march=native -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG   "${CMAKE_CXX_FLAGS_DEBUG}   -O0 -pedantic")

add_executable(ococo src/main.cpp)
add_dependencies(ococo htslib ococo_core)
add_library(ococo_core
	src/caller.h
	src/params.cpp
	src/params.h
	src/misc.cpp
	src/misc.h
	src/ococo.h
	src/types.h
	src/stats.h
	src/version.h
	)
target_link_libraries(ococo_core ${CMAKE_BINARY_DIR}/ext/htslib/libhts.a ${ZLIB_LIBRARIES})
target_link_libraries(ococo ${ZLIB_LIBRARIES} ococo_core)


################
# INSTALLATION #
################

install(TARGETS ococo DESTINATION bin)

if(INSTALL_DEBUG_SCRIPTS)
	install(PROGRAMS scripts/ococo_test_bam.sh DESTINATION bin)
endif(INSTALL_DEBUG_SCRIPTS)


#########
# TESTS #
#########

add_test(NAME unit_tests COMMAND unittests)
#add_test(NAME integration_tests COMMAND "cd tests/integration_test/ && ./integration_tests.bats")
