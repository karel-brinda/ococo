#include "ococo.h"

#ifdef DEBUGGING_MODE
#define BOOST_LOG_DYN_LINK
#endif

//#include <boost/format.hpp>

#include <climits>
#include <cstdio>
#include <cstdlib>

#ifdef DEBUGGING_MODE

#ifndef DEBUGGING_SEVERITY
#define DEBUGGING_SEVERITY trace
#endif

#include <boost/log/core.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/trivial.hpp>

namespace logging = boost::log;

void boost_logging_init() {
    logging::core::get()->set_filter(logging::trivial::severity >=
                                     logging::trivial::DEBUGGING_SEVERITY);
}

#endif

/*
    CONFIGURATION - STATISTICS
    --------------------------
*/

#ifdef OCOCO32
typedef uint32_t OCOCO_BASIC_TYPE;
#else
typedef uint16_t OCOCO_BASIC_TYPE;
#endif

constexpr uint32_t BITS_PER_COUNTER = (sizeof(OCOCO_BASIC_TYPE) * 8 - 4) / 4;
typedef ococo::stats_t<OCOCO_BASIC_TYPE, BITS_PER_COUNTER, 4> STATS_T;

/*
    --------------------------
*/


int main(int argc, const char *argv[]) {
    int main_return_code = 0;

#ifdef DEBUGGING_MODE
    boost_logging_init();

    /*
     BOOST_LOG_TRIVIAL(trace) << "A trace severity message";
     BOOST_LOG_TRIVIAL(debug) << "A debug severity message";
     BOOST_LOG_TRIVIAL(info) << "An informational severity message";
     BOOST_LOG_TRIVIAL(warning) << "A warning severity message";
     BOOST_LOG_TRIVIAL(error) << "An error severity message";
     BOOST_LOG_TRIVIAL(fatal) << "A fatal severity message";
     */

    BOOST_LOG_TRIVIAL(info) << "Ococo started.";
#endif

    /*
     * Default configuration.
     */

#ifdef DEBUGGING_MODE
    BOOST_LOG_TRIVIAL(info) << "Parsing command-line parameters.";
#endif
    ococo::consensus_params_t tmp_params = ococo::consensus_params_t(argc, argv);
    if (!tmp_params.correctly_initialized){
        return -1;
    }

    return main_return_code;
}