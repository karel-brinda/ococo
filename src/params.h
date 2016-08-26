#pragma once

#include "types.h"
#include "misc.h"

#include <cstdio>
#include <iostream>
#include <string>

#include <boost/program_options.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include <htslib/faidx.h>
#include <htslib/khash.h>
#include <htslib/kseq.h>
#include <htslib/kstring.h>
#include <htslib/sam.h>
#include <htslib/kseq.h>

/****************************
 *** Consensus parameters ***
 ****************************/

namespace ococo {

		  enum mode_t { BATCH, REALTIME };

		  enum strategy_t {
					 NO_UPDATES,
					 STOCHASTIC,
					 STOCHASTIC_AMB,
					 MAJORITY,
					 MAJORITY_AMB,
					 count
		  };

		  enum counter_configuration_t {
					 OCOCO16,
					 OCOCO32,
					 OCOCO64,
		  };

		  struct params_t {
					 bool correctly_initialized;
					 int return_code;

					 std::string command;

					 /*
					  * Counter parameters
					  */
					 counter_configuration_t counter_configuration;
					 std::string counters_str;
					 std::string counters_str_descr;
					 int32_t stats_bits_per_position;
					 int32_t stats_bits_per_nucleotide;

					 /*
					  * Input parameters
					  */
					 std::string sam_fn;
					 std::string fasta_in_fn;
					 std::string stats_in_fn;


					 /*
					  * Output parameters
					  */
					 bool verbose;

					 std::string vcf_fn;
					 std::string fasta_out_fn;
					 std::string stats_out_fn;
					 std::string pileup_fn;
					 std::string log_fn;

					 /*
					  * Files
					  */

					 FILE *vcf_file;
					 FILE *pileup_file;
					 FILE *fasta_out_file;
					 samFile *sam_file;
					 FILE *log_file;


					 /*
					  * Consensus calling parameters
					  */

					 mode_t mode;
					 strategy_t strategy;

					 /* minimum mapping quality for update */
					 int32_t min_mapq;

					 /* minimum base quality for update */
					 int32_t min_baseq;

					 /* initial values for counters corresponding to ref */
					 int32_t init_ref_weight;

					 /* minimum coverage for update (does not include init_ref_weight */
					 int32_t min_coverage;

					 /* threshold for having majority */
					 double majority_threshold;

					 /* auxiliary */
					 std::string strategy_str;
					 std::string mode_str;
					 int64_t n_upd;

					 /*
					  * Array of consensus calling functions
					  */
					 char (*cons_alg[strategy_t::count])(const pos_stats_uncompr_t &psu,
										  const params_t &params);

					 params_t();

					 params_t(int argc, const char *argv[]);

					 ~params_t();

					 void parse_commandline(int argc, const char *argv[]);

					 void init_default_values();
		  };

}
