#include "ococo.h"

#define BOOST_LOG_DYN_LINK

#include <climits>
#include <cstdio>
#include <cstdlib>

int main(int argc, const char *argv[]) {
    /* Use default configuration */
    ococo::params_t params = ococo::params_t(argc, argv);
    if (!params.correctly_initialized) {
        return EXIT_FAILURE;
    }

    switch (params.counter_configuration) {
        case ococo::OCOCO16: {
            ococo::caller_t<uint16_t, 3, 4> caller(&params);
            if (!caller.correctly_initialized) {
                return EXIT_FAILURE;
            }
            caller.run();
            return caller.return_code;
        }

        case ococo::OCOCO32: {
            ococo::caller_t<uint32_t, 7, 4> caller(&params);
            if (!caller.correctly_initialized) {
                return EXIT_FAILURE;
            }
            caller.run();
            return caller.return_code;
        }

        case ococo::OCOCO64: {
            ococo::caller_t<uint64_t, 15, 4> caller(&params);
            if (!caller.correctly_initialized) {
                return EXIT_FAILURE;
            }
            caller.run();
            return caller.return_code;
        }
    }

    return EXIT_FAILURE;
}
