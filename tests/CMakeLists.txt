include_directories(..)

add_executable(tests tests.cpp)
target_link_libraries(tests ${googletest_LIB} ococo_lib)
