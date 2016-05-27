#pragma once

#include <cstdarg>
#include <cstdio>
#include <iostream>
#include <cstdlib>
#include <string>

#include "version.h"

namespace ococo {

void print_version();

void fatal_error(const char *format, ...);

void error(const char *format, ...);

void warning(const char *format, ...);

void info(const char *format, ...);

bool file_exists(const std::string &fn);

/*
 * Get a right full mask (right n bits set to 1)
 *
 * T - type
 * size - number of 1's
 */
template <typename T, int size>
constexpr T right_full_mask() {
    static_assert(size <= 8 * sizeof(T), "Exceeding data type borders.");
    return (size == 0) ? 0
                       : (((static_cast<T>(0x1) << (size - 1)) - 1) << 1) | 1;
}

}
