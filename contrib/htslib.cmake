set(htslib_PREFIX ${CMAKE_BINARY_DIR}/contrib/htslib-prefix)
set(htslib_INSTALL ${CMAKE_BINARY_DIR}/contrib/htslib-install)

if (CMAKE_GENERATOR STREQUAL "Unix Makefiles")
    # when using the makefile generator, use the special variable $(MAKE) to invoke make
    # this enables the jobserver to work correctly
    set(MAKE_COMMAND "$(MAKE)")
else()
    # invoke make explicitly
    # in this case, we assume the parent build system is running in parallel already so no -j flag is added
    find_program(MAKE_COMMAND NAMES make gmake)
endif()

ExternalProject_Add(htslib
    PREFIX ${htslib_PREFIX}
    GIT_REPOSITORY "https://github.com/nsubtil/htslib.git"
    GIT_TAG "fix-signed-one-bit-bitfield"
    UPDATE_COMMAND ""
    BUILD_IN_SOURCE 1
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ${MAKE_COMMAND} lib-static
    INSTALL_COMMAND ${MAKE_COMMAND} install prefix=${htslib_INSTALL}
    LOG_DOWNLOAD 1
    )

add_dependencies(htslib zlib)
include_directories(${htslib_INSTALL}/include)
set(htslib_LIB ${htslib_INSTALL}/lib/libhts.a)
