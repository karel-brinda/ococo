#pragma once

#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <string>

namespace ococo {

void fatal_error(const char *format, ...) {
    va_list args;
    va_start(args, format);
    fprintf(stderr, "[ococo:fatal-error]: ");
    vfprintf(stderr, format, args);
    va_end(args);
}

void error(const char *format, ...) {
    va_list args;
    va_start(args, format);
    fprintf(stderr, "[ococo:error]: ");
    vfprintf(stderr, format, args);
    va_end(args);
}

void warning(const char *format, ...) {
    va_list args;
    va_start(args, format);
    fprintf(stderr, "[ococo:warning]: ");
    vfprintf(stderr, format, args);
    va_end(args);
}

void info(const char *format, ...) {
    va_list args;
    va_start(args, format);
    fprintf(stderr, "[ococo]: ");
    vfprintf(stderr, format, args);
    va_end(args);
}

bool file_exists(const std::string &fn) {
    FILE *file;

    file = fopen(fn.c_str(), "r");
    if (file) {
        fclose(file);
        return true;
    }
    return false;
}

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
