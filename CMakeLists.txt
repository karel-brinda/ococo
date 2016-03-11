cmake_minimum_required (VERSION 2.6)
project (ococo C CXX)


##############################
# COMPILATION MODE SWITCHERS #
##############################

if(NOT CMAKE_BUILD_TYPE)
    message(STATUS "Setting build type to 'DEBUG' as none was specified.")
    set(CMAKE_BUILD_TYPE DEBUG CACHE STRING "Choose the type of build." FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "DEBUG" "RELEASE" "RELWITHDEBINFO")
endif()

if(CMAKE_BUILD_TYPE MATCHES DEBUG)
    message("DEBUG mode")
endif(CMAKE_BUILD_TYPE MATCHES DEBUG)

if(CMAKE_BUILD_TYPE MATCHES RELEASE)
    message("RELEASE mode")
endif(CMAKE_BUILD_TYPE MATCHES RELEASE)

if(CMAKE_BUILD_TYPE MATCHES RELWITHDEBINFO)
    message("RELWITHDEBINFO mode")
endif(CMAKE_BUILD_TYPE MATCHES RELWITHDEBINFO)

#option(DEBUG "Build in the debugging mode" OFF)
option(OCOCO32 "Increase size of stats per position from 16b to 32b." OFF)
option(VERBOSE_VCF "Print also VCF records without updates." OFF)
option(INSTALL_DEBUG_SCRIPTS "Install debugging scripts." OFF)
set(DEBUGGING_SEVERITY "trace" CACHE STRING  "Verbosity of debugging mode.")

if(CMAKE_BUILD_TYPE MATCHES DEBUG)
    find_package( Boost 1.46 COMPONENTS program_options system log log_setup thread REQUIRED )
    add_definitions(-DDEBUGGING_MODE)
    add_definitions(-DDEBUGGING_SEVERITY=${DEBUGGING_SEVERITY})
else(CMAKE_BUILD_TYPE MATCHES DEBUG)
    find_package( Boost 1.46 COMPONENTS program_options REQUIRED )
endif(CMAKE_BUILD_TYPE MATCHES DEBUG)

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

include("cmake-ext/htslib.cmake")
include("cmake-ext/zlib.cmake")
include("cmake-ext/googletest.cmake")


################
# MAIN PROGRAM #
################

add_subdirectory(tests/unit_tests)

set(CMAKE_CXX_FLAGS                "${CMAKE_CXX_FLAGS}                -std=c++11 -Wall -Wextra -Wno-missing-field-initializers -Wshadow ")
set(CMAKE_CXX_FLAGS_RELEASE        "${CMAKE_CXX_FLAGS_RELEASE}        -O3 -march=native -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -O3 -march=native -DNDEBUG -g")
set(CMAKE_CXX_FLAGS_DEBUG          "${CMAKE_CXX_FLAGS_DEBUG}          -O0 -g -pedantic")

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
