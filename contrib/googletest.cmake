set(googletest_PREFIX ${CMAKE_BINARY_DIR}/contrib/googletest-prefix)
set(googletest_INSTALL ${CMAKE_BINARY_DIR}/contrib/googletest-install)

ExternalProject_Add(googletest
    PREFIX ${googletest_PREFIX}
    GIT_REPOSITORY "https://github.com/google/googletest.git"
    #GIT_TAG "fix-signed-one-bit-bitfield"
    UPDATE_COMMAND ""
    BUILD_IN_SOURCE 1
    INSTALL_DIR ${googletest_INSTALL}
    #CMAKE_ARGS -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
    #           -DCMAKE_INSTALL_PREFIX=${googletest_INSTALL}
    #           -DAMD64=ON
    LOG_DOWNLOAD ON
    LOG_CONFIGURE ON
    LOG_BUILD ON 
    )

include_directories(${googletest_INSTALL}/include)
set(googletest_LIB ${googletest_INSTALL}/lib/libgtest.a ${googletest_INSTALL}/lib/libgtest_main.a)
