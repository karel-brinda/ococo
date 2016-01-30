cmake_minimum_required (VERSION 2.6)
project (ococo)

include(ExternalProject)


###############
# COMPILATION #
###############

find_package( Boost 1.46 COMPONENTS program_options system log log_setup thread REQUIRED )
include_directories( ${Boost_INCLUDE_DIR} )

include("contrib/htslib.cmake")
include("contrib/zlib.cmake")
include("contrib/googletest.cmake")

add_subdirectory(tests/unit_tests)

#list( APPEND CMAKE_CXX_FLAGS "-std=gnu++11  -O3 -DNDEBUG ${CMAKE_CXX_FLAGS}")
#list(APPEND CMAKE_CXX_FLAGS "-Werrror -Wshadow -std=c++11")
list(APPEND CMAKE_CXX_FLAGS "-Wall -Wextra -pedantic -Wshadow -std=c++11")
#add_definitions(-D_FILE_OFFSET_BITS=64 -g)

#add_library(ococo_lib ococo.cpp)
#target_link_libraries(ococo_lib ${htslib_LIB} ${zlib_LIB} )

add_executable(ococo main.cpp)
#add_dependencies(ococo htslib zlib ococo_lib)
add_dependencies(ococo htslib zlib)
#target_link_libraries(ococo ococo_lib ${Boost_LIBRARIES} ${htslib_LIB} ${zlib_LIB})
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
