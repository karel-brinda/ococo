# build zlib
set(zlib_PREFIX ${CMAKE_BINARY_DIR}/contrib/zlib-prefix)
set(zlib_INSTALL ${CMAKE_BINARY_DIR}/contrib/zlib-install)

ExternalProject_Add(zlib
    PREFIX ${zlib_PREFIX}
#    GIT_REPOSITORY "https://github.com/jtkukunas/zlib.git"
#    GIT_TAG "e176b3c23ace88d5ded5b8f8371bbab6d7b02ba8"
    GIT_REPOSITORY "https://github.com/madler/zlib.git"
    GIT_TAG "v1.2.8"
    UPDATE_COMMAND ""
    BUILD_IN_SOURCE 1
    CONFIGURE_COMMAND ${zlib_PREFIX}/src/zlib/configure --prefix=${zlib_INSTALL} --static
    INSTALL_DIR ${zlib_INSTALL}
#    CMAKE_ARGS -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
#               -DCMAKE_INSTALL_PREFIX=${zlib_INSTALL}
#               -DAMD64=ON
    LOG_DOWNLOAD 1
    LOG_INSTALL 1
    )

include_directories(${zlib_INSTALL}/include)
set(zlib_LIB ${zlib_INSTALL}/lib/libz.a)
