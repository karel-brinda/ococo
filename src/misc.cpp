#include "misc.h"

void ococo::print_version() {
    // clang-format off
    std::cerr <<
           "\n"
           "Program: ococo (Online consensus caller, call cons. from an unsorted SAM/BAM stream)\n"
           "Version: " << OCOCO_VERSION  << "\n"
           "Contact: Karel Brinda <kbrinda@hsph.harvard.edu>\n";
    // clang-format on
    std::cerr << std::endl;
}

void ococo::fatal_error(const char *format, ...) {
    va_list args;
    va_start(args, format);
    fprintf(stderr, "[ococo:fatal-error]: ");
    vfprintf(stderr, format, args);
    va_end(args);
}

void ococo::error(const char *format, ...) {
    va_list args;
    va_start(args, format);
    fprintf(stderr, "[ococo:error]: ");
    vfprintf(stderr, format, args);
    va_end(args);
}

void ococo::warning(const char *format, ...) {
    va_list args;
    va_start(args, format);
    fprintf(stderr, "[ococo:warning]: ");
    vfprintf(stderr, format, args);
    va_end(args);
}

void ococo::info(const char *format, ...) {
    va_list args;
    va_start(args, format);
    fprintf(stderr, "[ococo]: ");
    vfprintf(stderr, format, args);
    va_end(args);
}

bool ococo::file_exists(const std::string &fn) {
    FILE *file;

    file = fopen(fn.c_str(), "r");
    if (file) {
        fclose(file);
        return true;
    }
    return false;
}

double ococo::realtime()
{
        struct timeval tp;
        //struct timezone tzp;
        //gettimeofday(&tp, &tzp);
        gettimeofday(&tp, nullptr);
        return tp.tv_sec + tp.tv_usec * 1e-6;
}

double ococo::cputime()
{
        struct rusage r;
        getrusage(RUSAGE_SELF, &r);

        //todo: check also memory
        //std::cerr << r.ru_maxrss << std::endl;
        return r.ru_utime.tv_sec + r.ru_stime.tv_sec + 1e-6 * (r.ru_utime.tv_usec + r.ru_stime.tv_usec);
}
