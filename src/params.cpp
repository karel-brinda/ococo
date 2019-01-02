#include "params.h"
#include "consensus.h"

#include <getopt.h>

/****************************
 *** Consensus parameters ***
 ****************************/

void ococo::params_t::init_default_values() {
    verbose               = false;
    counters_str          = "ococo32";
    counter_configuration = OCOCO32;
    mode                  = BATCH;
    mode_str              = "batch";
    strategy              = MAJORITY;
    strategy_str          = "majority";
    min_mapq              = default_q;
    min_baseq             = default_Q;
    init_ref_weight       = default_w;
    min_coverage          = default_c;
    majority_threshold    = default_M;
    coverage_filter       = default_C;

    cons_alg[strategy_t::NO_UPDATES]     = &cons_call_no_updates;
    cons_alg[strategy_t::STOCHASTIC]     = &cons_call_stoch;
    cons_alg[strategy_t::STOCHASTIC_AMB] = &cons_call_stoch_amb;
    cons_alg[strategy_t::MAJORITY]       = &cons_call_maj;
    cons_alg[strategy_t::MAJORITY_AMB]   = &cons_call_maj_amb;

    in_sam_file = nullptr;

    out_vcf_file    = nullptr;
    out_pileup_file = nullptr;
    out_fasta_file  = nullptr;
    out_sam_file    = nullptr;
    out_log_file    = nullptr;

    n_upd = 0;

    correctly_initialized = true;
    return_code           = 0;
}

ococo::params_t::params_t() { init_default_values(); }

ococo::params_t::params_t(int argc, const char **argv) {
    init_default_values();
    parse_commandline(argc, argv);
}

ococo::params_t::~params_t() {
    /*
     * Close files.
     */

    if (in_sam_file != nullptr) {
        int error_code = sam_close(in_sam_file);
        if (error_code != 0) {
            ococo::error("Input SAM file could not be closed.\n");
            return_code = -1;
        }
    }

    if (out_sam_file != nullptr) {
        int error_code = sam_close(out_sam_file);
        if (error_code != 0) {
            ococo::error("Output SAM file could not be closed.\n");
            return_code = -1;
        }
    }

    if (out_vcf_file != nullptr) {
        int error_code = fclose(out_vcf_file);
        if (error_code != 0) {
            ococo::error("Output VCF file could not be closed.\n");
            return_code = -1;
        }
    }

    if (out_pileup_file != nullptr) {
        int error_code = fclose(out_pileup_file);
        if (error_code != 0) {
            return_code = error_code;
            ococo::error("Output pileup file could not be closed.\n");
            return_code = -1;
        }
    }

    if (out_fasta_file != nullptr) {
        int error_code = fclose(out_fasta_file);
        if (error_code != 0) {
            ococo::error("Output FASTA consensus file could not be closed.\n");
            return_code = -1;
        }
    }

    if (out_log_file != nullptr) {
        int error_code = fclose(out_log_file);
        if (error_code != 0) {
            ococo::warning("Log file could not be closed.\n");
        }
    }
}

void ococo::params_t::print_help() {
    print_version();

    std::cerr
        <<
        // clang-format off

        // "---------------------------------------------------------------------------------"
           "Usage:   ococo -i <SAM/BAM file> [other options]\n\n"
        // "---------------------------------------------------------------------------------"
           "Input:\n"
           "  -i, --input FILE      input SAM/BAM file (- for standard input)\n"
           "  -f, --fasta-ref FILE  initial FASTA reference (otherwise seq of N's is used)\n"
           "  -s, --stats-in FILE   input statistics\n\n"
        // "---------------------------------------------------------------------------------"
           "Output:\n"
           "  -O, --output FILE     output SAM/BAM file (- for standard output)\n"
           "  -F, --fasta-cons FILE FASTA file with the consensus\n"
           "  -S, --stats-out FILE  output statistics\n"
           "  -V, --vcf-cons FILE   VCF file with updates of consensus (- for standard output)\n"
           "  -P, --pileup FILE     truncated pileup (- for standard output)\n\n"
        // "---------------------------------------------------------------------------------"
           "Variant and consensus calling:\n"
           "  -x, --counters STR    counter configuration: [ococo32]\n"
           "                           - ococo16 (3b/counter, 16b/position)\n"
           "                           - ococo32 (7b/counter, 32b/position)\n"
           "                           - ococo64 (15b/counter, 64b/position)\n"
           "  -m, --mode STR        mode: [batch]\n"
           "                           - real-time (updates reported immediately)\n"
           "                           - batch (updates reported after end of algn stream)\n"
        //    "  -t, --strategy STR    strategy for updates: [majority]\n"
        //    "                           - majority (update to majority base)\n"
        //    "                           - stochastic (update to stochastically drawn base)\n"
        //    "                           - no-updates (no updates, only counters updated)\n"
           "  -q, --min-MQ INT      skip alignments with mapping quality smaller than INT [" << default_q << "]\n"
           "  -Q, --min-BQ INT      skip bases with base quality smaller than INT [" << default_Q <<"]\n"
           "  -c, --min-cov INT     minimum coverage required for an update [" << default_c <<"]\n"
           "  -w, --ref-weight INT  initial counter value for nucleotides from the ref ["<< default_w <<"]\n"
           "  -M, --maj-thres FLOAT majority threshold [" << default_M << "]\n\n"
       // "---------------------------------------------------------------------------------"
           "Coverage normalization:\n"
           "  -C, --cov-filt INT    skip alignments when coverage greater than INT [" << default_C << "]\n\n"
       // "---------------------------------------------------------------------------------"
           "Others:\n"
           "  -v, --version         print version information and exit\n"
           "  -h, --help            print this message and exit\n"
           //"  --log FILE            auxiliary log file\n"
           "  --verbose             verbose mode (report every update of a counter)\n\n"
        // "---------------------------------------------------------------------------------"
           "Examples:\n"
           "   ococo -i test.bam -f test.fa -m real-time -V -\n"
           "   ococo -x ococo64 -i test.bam -f test.fa -P - -V variants.vcf\n\n"
        // "---------------------------------------------------------------------------------"
           "Note:\n"
           "   For more details, see the manual page 'man ./ococo.1'.\n"
        // "---------------------------------------------------------------------------------"
        // clang-format on
        << std::endl;
}

void ococo::params_t::parse_commandline(int argc, const char **argv) {
    /* Parse cmd parameters */
    std::stringstream cmd;
    for (int32_t i = 0; i < argc; i++) {
        cmd << argv[i];
        if (i != argc - 1) {
            cmd << " ";
        }
    }
    command = cmd.str();

    if (argc == 1) {
        print_help();
        exit(1);
    }

    const struct option lopts[] = {
        {"version", no_argument, nullptr, 'v'},
        {"help", no_argument, nullptr, 'h'},
        //
        {"input", required_argument, nullptr, 'i'},
        {"fasta-ref", required_argument, nullptr, 'f'},
        {"stats-in", required_argument, nullptr, 's'},
        //
        {"fasta-cons", required_argument, nullptr, 'F'},
        {"stats-out", required_argument, nullptr, 'S'},
        {"vcf-cons", required_argument, nullptr, 'V'},
        {"pileup", required_argument, nullptr, 'P'},
        {"output", required_argument, nullptr, 'O'},
        {"cov-filt", required_argument, nullptr, 'C'},
        {"log", required_argument, nullptr, 'L'},
        {"verbose", no_argument, nullptr, 'W'},
        //
        {"counters", required_argument, nullptr, 'x'},
        {"mode", required_argument, nullptr, 'm'},
        {"strategy", required_argument, nullptr, 't'},
        {"min-MQ", required_argument, nullptr, 'q'},
        {"min-BQ", required_argument, nullptr, 'Q'},
        {"ref-weight", required_argument, nullptr, 'w'},
        {"min-cov", required_argument, nullptr, 'c'},
        {"min-coverage", required_argument, nullptr, 'c'},  // deprec
        {"maj-thres", required_argument, nullptr, 'M'},
        {"majority-threshold", required_argument, nullptr, 'M'},  // deprec
        //
        {nullptr, 0, nullptr, 0}};

    int c;
    using std::string;
    while ((c = getopt_long(argc, (char *const *)argv,
                            "vhi:f:s:F:S:V:P:L:W:x:m:t:q:Q:w:c:M:O:C:", lopts,
                            nullptr)) >= 0) {
        switch (c) {
            case 'v': {
                print_version();
                exit(0);
                break;
            }
            case 'h': {
                print_help();
                exit(0);
                break;
            }
            case 'i': {
                in_sam_fn = optarg;
                break;
            }
            case 'f': {
                in_fasta_fn = optarg;
                break;
            }
            case 's': {
                in_stats_fn = optarg;
                break;
            }
            case 'F': {
                out_fasta_fn = optarg;
                break;
            }
            case 'S': {
                out_stats_fn = optarg;
                break;
            }
            case 'V': {
                out_vcf_fn = optarg;
                break;
            }
            case 'P': {
                out_pileup_fn = optarg;
                break;
            }
            case 'O': {
                out_sam_fn = optarg;
                break;
            }
            case 'C': {
                coverage_filter = atoi(optarg);
                break;
            }
            case 'L': {
                out_log_fn = optarg;
                break;
            }
            case 'W': {
                verbose = true;
                break;
            }
            case 'x': {
                counters_str = optarg;

                if (counters_str.compare("ococo16") == 0) {
                    counter_configuration = OCOCO16;
                    counters_str_descr =
                        "ococo16 (16 bits per position, 3bits per nucleotide "
                        "counter)";
                } else if (counters_str.compare("ococo32") == 0) {
                    counter_configuration = OCOCO32;
                    counters_str_descr =
                        "ococo32 (32 bits per position, 7bits per nucleotide "
                        "counter)";
                } else if (counters_str.compare("ococo64") == 0) {
                    counter_configuration = OCOCO64;
                    counters_str_descr =
                        "ococo64 (64 bits per position, 15bits per nucleotide "
                        "counter)";
                } else {
                    ococo::error(
                        "Unknown counter configuration '%s'. Possible modes "
                        "are 'ococo16', 'ococo32', and 'ococo64'.\n",
                        counters_str.c_str());
                    exit(1);
                }

                break;
            }
            case 'm': {
                mode_str = optarg;

                if (mode_str.compare("batch") == 0) {
                    mode = ococo::mode_t::BATCH;
                } else if (mode_str.compare("real-time") == 0) {
                    mode = ococo::mode_t::REALTIME;
                } else {
                    ococo::error(
                        "Unknown mode '%s'. Possible modes are 'batch' and "
                        "'real-time'.\n",
                        mode_str.c_str());
                    exit(1);
                }

                break;
            }
            case 't': {
                strategy_str = optarg;

                if (strategy_str.compare("stochastic") == 0) {
                    strategy = ococo::strategy_t::STOCHASTIC;
                } else if (strategy_str.compare("no-updates") == 0) {
                    strategy = ococo::strategy_t::NO_UPDATES;
                } else if (strategy_str.compare("majority") == 0) {
                    strategy = ococo::strategy_t::MAJORITY;
                } else {
                    ococo::error(
                        "Unknown strategy '%s'. Possible strategies are "
                        "'majority', 'stochastic', and 'no-updates'.\n",
                        strategy_str.c_str());
                    exit(1);
                }

                break;
            }
            case 'q': {
                min_mapq = atoi(optarg);
                break;
            }
            case 'Q': {
                min_baseq = atoi(optarg);
                break;
            }
            case 'w': {
                init_ref_weight = atoi(optarg);
                break;
            }
            case 'c': {
                min_coverage = atoi(optarg);
                break;
            }
            case 'M': {
                majority_threshold = atof(optarg);
                break;
            }
            case '?': {
                ococo::error("Unknown error");
                exit(1);
                break;
            }
        }
    }
    if (in_sam_fn.size() == 0) {
        ococo::error("SAM/BAM file must be specified (option '-i').\n");
        exit(1);
    }

    if (out_pileup_fn.size() != 0 && out_pileup_fn.compare(out_vcf_fn) == 0) {
        ococo::error(
            "Pileup and VCF files cannot be the same (both currently '%s').\n",
            out_pileup_fn.c_str());
        exit(1);
    }

    ococo::info("Ococo starting: %s\n", counters_str_descr.c_str());
}
