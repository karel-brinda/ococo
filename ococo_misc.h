#ifndef _OCOCO_MISC_H_
#define _OCOCO_MISC_H_

#pragma once

#include <string>
#include <cstdlib>
#include <cstdarg>
#include <cstdio>

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
template <typename T, int size> constexpr T right_full_mask() {
    static_assert(size <= 8 * sizeof(T), "Exceeding data type borders.");
    return (size == 0) ? 0
                       : (((static_cast<T>(0x1) << (size - 1)) - 1) << 1) | 1;
}

/*
 * Get a left full mask (left n bits set to 1)
 *
 * T - type
 * size - number of 1's
 */
template <typename T, int size> constexpr T left_full_mask() {
    return right_full_mask<T, size>() << (8 * sizeof(T) - size);
}

/*
 * Get bits ending at a right coordinate
 *
 * T - type
 * right - the rightmost bit (from right)
 * size - number of bits
 */
template <typename T, int size, int rightmost_bit>
inline T get_right_bits(T pattern) {
    static_assert(rightmost_bit + size <= 8 * sizeof(T),
                  "Exceeding data type borders.");
    return (pattern >> rightmost_bit) & right_full_mask<T, size>();
}

/*
 * Get bits starting at a left coordinate
 *
 * T - type
 * left - the leftmost bit (from left)
 * size - number of bits
 */
template <typename T, int size, int leftmost_bit>
inline T get_left_bits(T pattern) {
    return get_right_bits<T, size, 8 * sizeof(T) - leftmost_bit - size>(
        pattern);
}

#ifdef not_finished_remove
/*
 * Set bits end at a right coordinate
 *
 * T - type
 * right - the rightmost bit (from right)
 * size - number of bits
 */
template <typename T, int size, int rightmost_bit>
inline T set_right_bits(T pattern) {
    static_assert(rightmost_bit + size <= 8 * sizeof(T),
                  "Exceeding data type borders.");
    return ((pattern & right_full_mask<T, size>()) << rightmost_bit);
}

/*
 * Set bits starting at a left coordinate
 *
 * T - type
 * left - the leftmost bit (from left)
 * size - number of bits
 */
template <typename T, int size, int leftmost_bit>
inline T set_left_bits(T pattern) {
    return set_right_bits<T, 8 * sizeof(T) - leftmost_bit - size, size>(
        pattern);
}
#endif
};

#endif