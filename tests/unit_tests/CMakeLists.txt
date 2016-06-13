include_directories(../..)

if (BUILD_TESTS)

	list(APPEND CMAKE_CXX_FLAGS "-Wall -std=c++11")

	add_executable(unittests unittests.cpp)

	find_package (Threads REQUIRED)

	enable_testing()
	find_package(GTest REQUIRED)
	include_directories(${GTEST_INCLUDE_DIRS})

	target_link_libraries(unittests ${GTEST_BOTH_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} ${Boost_LIBRARIES} ${htslib_LIB} ${zlib_LIB} ococo_core)

endif (BUILD_TESTS)
