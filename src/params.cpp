#include "params.h"
#include "consensus.h"

#include <getopt.h>

/****************************
 *** Consensus parameters ***
 ****************************/

void ococo::params_t::init_default_values() {
    verbose               = false;
    counters_str          = "ococo16";
    counter_configuration = OCOCO16;
    mode                  = BATCH;
    mode_str              = "batch";
    strategy              = MAJORITY;
    strategy_str          = "majority";
    min_mapq              = 1;
    min_baseq             = 13;
    init_ref_weight       = 0;
    min_coverage          = 2;
    majority_threshold    = 0.60;

    cons_alg[strategy_t::NO_UPDATES]     = &cons_call_no_updates;
    cons_alg[strategy_t::STOCHASTIC]     = &cons_call_stoch;
    cons_alg[strategy_t::STOCHASTIC_AMB] = &cons_call_stoch_amb;
    cons_alg[strategy_t::MAJORITY]       = &cons_call_maj;
    cons_alg[strategy_t::MAJORITY_AMB]   = &cons_call_maj_amb;

    vcf_file       = nullptr;
    pileup_file    = nullptr;
    fasta_out_file = nullptr;
    sam_file       = nullptr;
    log_file       = nullptr;

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

    if (sam_file != nullptr) {
        int error_code = sam_close(sam_file);
        if (error_code != 0) {
            ococo::error("SAM file could not be closed.\n");
            return_code = -1;
        }
    }

    if (vcf_file != nullptr) {
        int error_code = fclose(vcf_file);
        if (error_code != 0) {
            ococo::error("VCF file could not be closed.\n");
            return_code = -1;
        }
    }

    if (pileup_file != nullptr) {
        int error_code = fclose(pileup_file);
        if (error_code != 0) {
            return_code = error_code;
            ococo::error("Pileup file could not be closed.\n");
            return_code = -1;
        }
    }

    if (fasta_out_file != nullptr) {
        int error_code = fclose(fasta_out_file);
        if (error_code != 0) {
            ococo::error("FASTA consensus file could not be closed.\n");
            return_code = -1;
        }
    }

    if (log_file != nullptr) {
        int error_code = fclose(log_file);
    }
}

void ococo::params_t::print_help() {
    std::cerr
        <<

        // clang-format off
        // "---------------------------------------------------------------------------------"
           "Generic options:\n"
           "  -v, --version         print version and exit\n"
           "  -h, --help            print this message and exit\n\n"           
        // "---------------------------------------------------------------------------------"
           "Input options:\n"
           "  -i, --input FILE      input SAM/BAM file (- for standard input)\n"
           "  -f, --fasta-ref FILE  initial FASTA reference (otherwise sequence of N's \n"
           "                           considered as the reference)\n"
           "  -s, --stats-in arg    input statistics.\n\n"
        // "---------------------------------------------------------------------------------"
           "Output options:\n"
           "  -F, --fasta-cons FILE FASTA file with consensus\n"
           "  -S, --stats-out FILE  outputs statistics\n"
           "  -V, --vcf-cons FILE   VCF file with updates of consensus (- for standard output)\n"
           "  -P, --pileup FILE     truncated pileup (- for standard output)\n"
           "  --log FILE            auxiliary log file\n"
           "  --verbose             verbose mode\n\n"
        // "---------------------------------------------------------------------------------"
           "Parameters of consensus calling:\n"
           "  -x, --counters STR    counters configuration: [ococo16]\n"
           "                           - ococo16 (3b/counter, 16b/position)\n"
           "                           - ococo32 (7b/counter, 32b/position)\n"
           "                           - ococo64 (15b/counter, 64b/position)\n"
           "  -m, --mode STR        mode:  [batch]\n"
           "                           - real-time / batch\n"
           "  -t, --strategy STR    strategy for updates: [majority]\n"
           "                           - no-updates / majority / stochastic\n"
           //"  -a [ --allow-amb ]                    Allow updates to ambiguous "
           //"nucleotides.\n"
           "  -q, --min-MQ INT      skip alignments with mapping quality smaller than INT [1]\n"
           "  -Q, --min-BQ INT      skip bases with base quality smaller than INT [13]\n"
           "  -w, --ref-weight INT  initial counter value for nucleotides from ref [0]\n"
           "  -c, --min-cov INT     minimum coverage required for update [2]\n"
           "  -M, --maj-thres FLOAT majority threshold [0.6]"
        // "---------------------------------------------------------------------------------"
        // clang-format on
        << std::endl;
}

void ococo::params_t::parse_commandline(int argc, const char **argv) {
    /* Save cmd parameters */

    std::stringstream cmd;
    for (int32_t i = 0; i < argc; i++) {
        cmd << argv[i];
        if (i != argc - 1) {
            cmd << " ";
        }
    }
    command = cmd.str();

    /* Parse cmd parameters */

    const struct option lopts[] = {
        {"version", no_argument, NULL, 'v'},
        {"help", no_argument, NULL, 'h'},
        //
        {"input", required_argument, NULL, 'i'},
        {"fasta-ref", required_argument, NULL, 'f'},
        {"stats-in", required_argument, NULL, 's'},
        //
        {"fasta-cons", required_argument, NULL, 'F'},
        {"stats-out", required_argument, NULL, 'S'},
        {"vcf-cons", required_argument, NULL, 'V'},
        {"pileup", required_argument, NULL, 'P'},
        {"log", required_argument, NULL, 'L'},
        {"verbose", required_argument, NULL, 'W'},
        //
        {"counters", required_argument, NULL, 'x'},
        {"mode", required_argument, NULL, 'm'},
        {"strategy", required_argument, NULL, 's'},
        {"min-MQ", required_argument, NULL, 'q'},
        {"min-BQ", required_argument, NULL, 'Q'},
        {"ref-weight", required_argument, NULL, 'w'},
        {"min-cov", required_argument, NULL, 'c'},
        {"min-coverage", required_argument, NULL, 'c'},  // deprec
        {"maj-thres", required_argument, NULL, 'M'},
        {"majority-threshold", required_argument, NULL, 'M'},  // deprec
        //
        {NULL, 0, NULL, 0}};

    int c;
    while ((c = getopt_long(argc, (char *const *)argv,
                            "vhi:f:s:F:S:V:P:L:W:x:m:s:q:Q:w:c:M:", lopts,
                            NULL)) >= 0) {
        switch (c) {
            case 'h':
                // print_help();
                break;
            case 'v':
                // print_help();
                break;
            case 1:
                break;
            case '?':
                correctly_initialized = false;
                return_code           = -1;
                // return -1;
        }
    }

    try {
        namespace po = boost::program_options;

        po::options_description options_generic("Generic options");
        options_generic.add_options()
            //
            ("version,v", "Print version and exit.")
            //
            ("help,h", "Print this message and exit.")
            //
            ;

        po::options_description options_input("Input options");
        options_input.add_options()
            //
            ("input,i", po::value<std::string>(&sam_fn)->required(),
             "Input SAM/BAM file (- for standard input).")
            //
            ("fasta-ref,f", po::value<std::string>(&fasta_in_fn),
             "Initial FASTA reference (if not provided, sequence of N's is "
             "considered as the reference).")
            //
            ("stats-in,s", po::value<std::string>(&stats_in_fn),
             "Input statistics.")
            //
            ;

        po::options_description options_output("Output options");
        options_output.add_options()
            //
            ("fasta-cons,F", po::value<std::string>(&fasta_out_fn),
             "FASTA file with consensus.")
            //
            ("stats-out,S", po::value<std::string>(&stats_out_fn),
             "Outputs statistics.")
            //
            ("vcf-cons,V", po::value<std::string>(&vcf_fn),
             "VCF file with updates of consensus (- for standard output).")
            //
            ("pileup,P", po::value<std::string>(&pileup_fn),
             "Truncated pileup (- for standard output).")
            //
            ("log", po::value<std::string>(&log_fn), "Auxiliary log file.")
            //
            ("verbose", "Verbose mode.")
            //
            ;

        po::options_description options_consensus(
            "Parameters of consensus calling");
        options_consensus.add_options()
            //
            ("counters,x",
             po::value<std::string>(&counters_str)->default_value(counters_str),
             "Counters configuration: \n - ococo16 (3b/counter, "
             "16b/position)\n - ococo32 (7b/counter, 32b/position)\n - ococo64 "
             "(15b/counter, 64b/position)")
            //
            ("mode,m",
             po::value<std::string>(&mode_str)->default_value(mode_str),
             "Mode: real-time / batch.")
            //
            ("strategy,t",
             po::value<std::string>(&strategy_str)->default_value(strategy_str),
             "Strategy for updates: no-updates / majority / stochastic.")
            //
            ("allow-amb,a", "Allow updates to ambiguous nucleotides.")
            //
            ("min-MQ,q", po::value<int32_t>(&min_mapq)->default_value(min_mapq),
             "Skip alignments with mapping quality smaller than INT.")
            //
            ("min-BQ,Q",
             po::value<int32_t>(&min_baseq)->default_value(min_baseq),
             "Skip bases with base quality smaller than INT.")
            //
            ("ref-weight,w", po::value<int32_t>(&init_ref_weight)
                                 ->default_value(init_ref_weight),
             "Initial counter value for nucleotides from the reference.")
            //
            ("min-coverage,c",
             po::value<int32_t>(&min_coverage)->default_value(min_coverage),
             "Minimum coverage required for update.")
            //
            ("majority-threshold,M", po::value<double>(&majority_threshold)
                                         ->default_value(majority_threshold),
             "Majority threshold.")
            //
            ;

        po::options_description options_all;
        options_all.add(options_generic)
            .add(options_input)
            .add(options_output)
            .add(options_consensus);

        po::variables_map vm;
        try {
            po::store(
                po::command_line_parser(argc, argv).options(options_all).run(),
                vm);  // can throw

            if (vm.count("version")) {
                std::cout << std::endl;
                print_version();
                std::cout << std::endl;
                exit(0);
            }

            if (vm.count("help")) {
                std::cout << options_all << "\n";
                exit(0);
            }

            po::notify(vm);
            // throws on error, so do after help in case there
            // are any problems
            if (vm.count("strategy")) {
                if (strategy_str.compare("no-updates") == 0) {
                    strategy = ococo::strategy_t::NO_UPDATES;
                } else if (strategy_str.compare("majority") == 0) {
                    if (vm.count("allow-amb") == 0) {
                        strategy = ococo::strategy_t::MAJORITY;
                    } else {
                        strategy = ococo::strategy_t::MAJORITY_AMB;
                    }
                } else if (strategy_str.compare("stochastic") == 0) {
                    if (vm.count("allow-amb") == 0) {
                        strategy = ococo::strategy_t::STOCHASTIC;
                    } else {
                        strategy = ococo::strategy_t::STOCHASTIC_AMB;
                    }
                } else {
                    ococo::error(
                        "Unknown strategy '%s'. Possible strategies "
                        "are 'majority' and 'stochastic'.\n",
                        strategy_str.c_str());
                    correctly_initialized = false;
                    return_code           = -1;
                    return;
                }
            }

            if (vm.count("mode")) {
                if (mode_str.compare("batch") == 0) {
                    mode = ococo::mode_t::BATCH;
                } else if (mode_str.compare("real-time") == 0) {
                    mode = ococo::mode_t::REALTIME;
                } else {
                    ococo::error(
                        "Unknown mode '%s'. Possible modes are "
                        "'batch' and 'real-time'.\n",
                        mode_str.c_str());
                    correctly_initialized = false;
                    return_code           = -1;
                    return;
                }
            }

            if (vm.count("verbose")) {
                verbose = true;
            }

            if (vm.count("counters")) {
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
                    correctly_initialized = false;
                    return_code           = -1;
                    return;
                }
            }
            ococo::info("Ococo starting: %s\n", counters_str_descr.c_str());

        } catch (po::error &e) {
            std::cout << options_all << "\n";
            ococo::error("%s.\n", e.what());
            correctly_initialized = false;
            return_code           = -1;
            return;
        }

    } catch (std::exception &e) {
        ococo::error("Unhandled Exception: %s.\n", e.what());
        correctly_initialized = false;
        return_code           = -1;
        return;
    }
}
