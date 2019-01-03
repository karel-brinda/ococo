/* The MIT License

   Copyright (c) 2016-2019 Karel Brinda (kbrinda@hsph.harvard.edu)

   Permission is hereby granted, free of charge, to any person obtaining a copy
   of this software and associated documentation files (the "Software"), to deal
   in the Software without restriction, including without limitation the rights
   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
   copies of the Software, and to permit persons to whom the Software is
   furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included in
   all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

#include "ococo.h"

#include <climits>
#include <cstdio>
#include <cstdlib>

using namespace ococo;

int main(int argc, const char **argv) {
    /* Use the default configuration */
    params_t params = params_t(argc, argv);
    if (!params.correctly_initialized) {
        return EXIT_FAILURE;
    }

    switch (params.counter_configuration) {
        case OCOCO16: {
            ococo_t<uint16_t, 3, 4> ococo(&params);
            if (!ococo.correctly_initialized) {
                return EXIT_FAILURE;
            }
            ococo.run();
            return ococo.return_code;
        }

        case OCOCO32: {
            ococo_t<uint32_t, 7, 4> ococo(&params);
            if (!ococo.correctly_initialized) {
                return EXIT_FAILURE;
            }
            ococo.run();
            return ococo.return_code;
        }

        case OCOCO64: {
            ococo_t<uint64_t, 15, 4> ococo(&params);
            if (!ococo.correctly_initialized) {
                return EXIT_FAILURE;
            }
            ococo.run();
            return ococo.return_code;
        }
    }

    return EXIT_FAILURE;
}
