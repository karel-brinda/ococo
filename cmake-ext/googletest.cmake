include(ExternalProject)

set(googletest_PREFIX ${CMAKE_BINARY_DIR}/cmake-ext/googletest-prefix)
set(googletest_INSTALL ${CMAKE_BINARY_DIR}/cmake-ext/googletest-install)

ExternalProject_Add(googletest
    PREFIX ${googletest_PREFIX}
    GIT_REPOSITORY "https://github.com/google/googletest.git"
    GIT_TAG "ff5ffd457e0"
    INSTALL_DIR ${googletest_INSTALL}
    UPDATE_COMMAND ""
    CMAKE_ARGS
        -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
        -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
        -DCMAKE_INSTALL_PREFIX=${googletest_INSTALL}
)

include_directories(${googletest_INSTALL}/include)
set(googletest_LIB ${googletest_INSTALL}/lib/libgtest.a ${googletest_INSTALL}/lib/libgtest_main.a)

