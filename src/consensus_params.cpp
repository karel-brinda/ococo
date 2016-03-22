#include "consensus_params.h"
#include "consensus_functions.h"

/****************************
 *** Consensus parameters ***
 ****************************/

void ococo::consensus_params_t::init_default_values() {
    mode=BATCH;
    strategy=MAJORITY;
    min_mapq=1;
    min_baseq=13;
    init_ref_weight=0;
    min_coverage=2;
    majority_threshold=0.60;

    cons_alg[strategy_t::NO_UPDATES] = &cons_call_no_updates;
    cons_alg[strategy_t::STOCHASTIC] = &cons_call_stoch;
    cons_alg[strategy_t::STOCHASTIC_AMB] = &cons_call_stoch_amb;
    cons_alg[strategy_t::MAJORITY] = &cons_call_maj;
    cons_alg[strategy_t::MAJORITY_AMB] = &cons_call_maj_amb;

    vcf_file = nullptr;
    pileup_file = nullptr;
    fasta_out_file = nullptr;
    sam_file = nullptr;

    correctly_initialized=true;
}

ococo::consensus_params_t::consensus_params_t(){
	init_default_values();
}

ococo::consensus_params_t::consensus_params_t(int argc, const char *argv[]){
	init_default_values();
	parse_commandline(argc, argv);
}

ococo::consensus_params_t::~consensus_params_t(){
    /*
     * Close files.
     */

    if (sam_file != nullptr) {
        int error_code = sam_close(sam_file);
        if (error_code != 0) {
            ococo::error("SAM file could not be closed.\n");
            //main_return_code = -1;
        }
    }

    if (vcf_file != nullptr) {
        int error_code = fclose(vcf_file);
        if (error_code != 0) {
            ococo::error("VCF file could not be closed.\n");
            //main_return_code = -1;
        }
    }

    if (pileup_file != nullptr) {
        int error_code = fclose(pileup_file);
        if (error_code != 0) {
            ococo::error("Pileup file could not be closed.\n");
            //main_return_code = -1;
        }
    }

    if (fasta_out_file != nullptr) {
        int error_code = fclose(fasta_out_file);
        if (error_code != 0) {
            ococo::error("FASTA consensus file could not be closed.\n");
            //main_return_code = -1;
        }
    }
}

void ococo::consensus_params_t::parse_commandline(int argc, const char *argv[]){

    try {

        std::stringstream min_mq_message;
        min_mq_message
            << "Skip alignments with mapping quality smaller than INT. ["
            << min_mapq << "]";

        std::stringstream min_bq_message;
        min_bq_message << "Skip bases with base quality smaller than INT. ["
                       << min_baseq << "]";

        std::stringstream ref_weight_message;
        ref_weight_message
            << "Initial counter value for nucleotides from the reference. ["
            << init_ref_weight << "]";

        std::stringstream min_coverage_message;
        min_coverage_message << "Minimum coverage required for update. ["
                             << min_coverage << "]";

        std::stringstream majority_threshold_message;
        majority_threshold_message << "Majority threshold. ["
                                   << majority_threshold << "]";

        namespace po = boost::program_options;

        po::options_description options_generic("Generic options");
        options_generic.add_options()
            //
            ("version,v",
            "Print version and exit.");


        po::options_description options_input("Input options");
        options_input.add_options()
            //
            ("input,i",
            po::value<std::string>(&sam_fn)->required(),
            "Input SAM/BAM file (- for standard input).")
            //
            (
            "fasta-ref,f", po::value<std::string>(&fasta_in_fn),
            "Initial FASTA reference (if not provided, sequence of N's is "
            "considered as the reference).")
            //
            ("stats-in,s",po::value<std::string>(&stats_in_fn),
            "Input statistics.")
            //
            ;

        po::options_description options_output("Output options");
        options_output.add_options()
            //
            (
            "fasta-cons,F", po::value<std::string>(&fasta_out_fn),
            "FASTA file with consensus.")
            //
            (
            "stats-out,S", po::value<std::string>(&stats_out_fn),
            "Outputs statistics.")
            //
            (
            "vcf-cons,V", po::value<std::string>(&vcf_fn),
            "VCF file with updates of consensus (- for standard output)."
            )
            //
            (
            "pileup,P", po::value<std::string>(&pileup_fn),
            "Truncated pileup (- for standard output).")
            //
            ("verbose",
            "Verbose mode.")
            //
            ;


        po::options_description options_consensus("Parameters of consensus calling");
        options_consensus.add_options()
            //
            (
            "mode,m", po::value<std::string>(&mode_str),
            "Mode: real-time / batch. [batch]")
            //
            (
            "strategy,t", po::value<std::string>(&strategy_str),
            "Strategy for updates: no-updates / majority / stochastic. "
            "[majority]")
            //
            ("allow-amb,a", "Allow updates to ambiguous nucleotides.")
            //
            (
            "min-MQ,q", po::value<int32_t>(&min_mapq),
            min_mq_message.str().c_str())
            //
            (
            "min-BQ,Q", po::value<int32_t>(&min_baseq),
            min_bq_message.str().c_str())
            //
            (
            "ref-weight,w", po::value<int32_t>(&init_ref_weight),
            ref_weight_message.str().c_str())
            //
            (
            "min-coverage,c", po::value<int32_t>(&min_coverage),
            min_coverage_message.str().c_str())
            //
            (
            "majority-threshold,M",
            po::value<double>(&majority_threshold),
            majority_threshold_message.str().c_str()
            )
            //
            ;

        po::options_description options_all;
        options_all.add(options_generic).add(options_input).add(options_output).add(options_consensus);

        po::variables_map vm;
        try {

            po::store(po::command_line_parser(argc, argv)
                          .options(options_all)
                          .run(),
                      vm);  // can throw

            if (vm.count("version")) {
                std::cout<<std::endl;
                print_version();
                std::cout<<std::endl;
                exit(0);
            }

            po::notify(vm); // throws on error, so do after help in case there
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
                    ococo::error("Unknown strategy '%s'. Possible strategies "
                                 "are 'majority' and 'stochastic'.\n",
                                 strategy_str.c_str());
                    correctly_initialized=false;
                    return;
                }
            }

            if (vm.count("mode")) {
                if (mode_str.compare("batch") == 0) {
                    mode = ococo::mode_t::BATCH;
                } else if (mode_str.compare("real-time") == 0) {
                    mode = ococo::mode_t::REALTIME;
                } else {
                    ococo::error("Unknown mode '%s'. Possible modes are "
                                 "'batch' and 'real-time'.\n",
                                 mode_str.c_str());
                    correctly_initialized=false;
                    return;
                }
            }

        } catch (po::error &e) {
            std::cout << options_all << "\n";
            ococo::error("%s.\n", e.what());
            correctly_initialized=false;
            return;
        }

    } catch (std::exception &e) {
        ococo::error("Unhandled Exception: %s.\n", e.what());
        correctly_initialized=false;
        return;
    }
}