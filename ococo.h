#pragma once

#include "ococo_misc.h"
#include "consensus_params.h"
#include "consensus_functions.h"
#include "ococo_types.h"
#include "stats.h"

#include "htslib/kseq.h"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <exception>
#include <iostream>
#include <sstream>
#include <string>
#include <zlib.h>

namespace ococo {

KSEQ_INIT(gzFile, gzread)

/*********************
 *** Configuration ***
 *********************/

const int fasta_line_l = 50;
const int stats_delim_l = 10;


/***********************************
 *                                 *
 *   Consensus calling functions   *
 *                                 *
 ***********************************/

}
