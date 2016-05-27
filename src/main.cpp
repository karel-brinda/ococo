#include "ococo.h"

//#ifdef DEBUGGING_MODE
#define BOOST_LOG_DYN_LINK
//#endif

//#include <boost/format.hpp>

#include <climits>
#include <cstdio>
#include <cstdlib>

#include <boost/log/core.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/trivial.hpp>

#ifndef DEBUGGING_SEVERITY
#define DEBUGGING_SEVERITY trace
#endif

namespace logging = boost::log;

#ifdef DEBUGGING_MODE
const bool debugging = true;
#else
const bool debugging = false;
#endif

/*
 --------------------------
 */


int main(int argc, const char *argv[]) {
    if(debugging){
        
        logging::core::get()->set_filter(logging::trivial::severity >=
                                         logging::trivial::DEBUGGING_SEVERITY);
        
        /*
         BOOST_LOG_TRIVIAL(trace) << "A trace severity message";
         BOOST_LOG_TRIVIAL(debug) << "A debug severity message";
         BOOST_LOG_TRIVIAL(info) << "An informational severity message";
         BOOST_LOG_TRIVIAL(warning) << "A warning severity message";
         BOOST_LOG_TRIVIAL(error) << "An error severity message";
         BOOST_LOG_TRIVIAL(fatal) << "A fatal severity message";
         */
        
        BOOST_LOG_TRIVIAL(info) << "Ococo started.";
    }
    
    /*
     * Default configuration.
     */    
    ococo::params_t params = ococo::params_t(argc, argv);
    if (!params.correctly_initialized){
        return EXIT_FAILURE;
    }
    
    switch (params.counter_configuration){
            
        case ococo::OCOCO16:
        {
            ococo::caller_t<uint16_t, 3, 4, debugging> caller(&params);
            if (!caller.correctly_initialized){
                return EXIT_FAILURE ;
            }
            caller.run();
            return caller.return_code;
        }
            
        case ococo::OCOCO32:
        {
            ococo::caller_t<uint32_t, 7, 4, debugging> caller(&params);
            if (!caller.correctly_initialized){
                return EXIT_FAILURE ;
            }
            caller.run();
            return caller.return_code;
        }
            
        case ococo::OCOCO64:
        {
            ococo::caller_t<uint64_t, 15, 4, debugging> caller(&params);
            if (!caller.correctly_initialized){
                return EXIT_FAILURE ;
            }
            caller.run();
            return caller.return_code;
        }
    }
    
    return EXIT_FAILURE;
}
