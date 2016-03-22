#include "ococo.h"

const char *OCOCO_VERSION = "0.1.0";

#ifdef DEBUGGING_MODE
#define BOOST_LOG_DYN_LINK
#endif

//#include <boost/format.hpp>

#include <boost/program_options.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

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

void print_version(){
    std::cerr << "Program: Ococo (online consensus caller, "
        << "http://github.com/karel-brinda/ococo)." << std::endl
        << "Version: " << OCOCO_VERSION << " ("
        << 8 * sizeof(OCOCO_BASIC_TYPE) << "bit variant"
        << ", counter size " << BITS_PER_COUNTER << "bits";

#ifdef DEBUGGING_MODE
    std::cerr << ", debugging mode";
#endif

#ifdef VERBOSE_VCF
    std::cerr << ", verbose VCF";
#endif

    std::cerr << ")." << std::endl;
}

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

    ococo::consensus_params_t tmp_params = ococo::consensus_params_t();

/*
 * Parse command-line parameters.
 */

#ifdef DEBUGGING_MODE
    BOOST_LOG_TRIVIAL(info) << "Parsing command-line parameters.";
#endif

    try {

        std::stringstream min_mq_message;
        min_mq_message
            << "Skip alignments with mapping quality smaller than INT. ["
            << tmp_params.min_mapq << "]";

        std::stringstream min_bq_message;
        min_bq_message << "Skip bases with base quality smaller than INT. ["
                       << tmp_params.min_baseq << "]";

        std::stringstream ref_weight_message;
        ref_weight_message
            << "Initial counter value for nucleotides from the reference. ["
            << tmp_params.init_ref_weight << "]";

        std::stringstream min_coverage_message;
        min_coverage_message << "Minimum coverage required for update. ["
                             << tmp_params.min_coverage << "]";

        std::stringstream majority_threshold_message;
        majority_threshold_message << "Majority threshold. ["
                                   << tmp_params.majority_threshold << "]";

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
            po::value<std::string>(&tmp_params.sam_fn)->required(),
            "Input SAM/BAM file (- for standard input).")
            //
            (
            "fasta-ref,f", po::value<std::string>(&tmp_params.fasta_in_fn),
            "Initial FASTA reference (if not provided, sequence of N's is "
            "considered as the reference).")
            //
            ("stats-in,s",po::value<std::string>(&tmp_params.stats_in_fn),
            "Input statistics.")
            //
            ;

        po::options_description options_output("Output options");
        options_output.add_options()
            //
            (
            "fasta-cons,F", po::value<std::string>(&tmp_params.fasta_out_fn),
            "FASTA file with consensus.")
            //
            (
            "stats-out,S", po::value<std::string>(&tmp_params.stats_out_fn),
            "Outputs statistics.")
            //
            (
            "vcf-cons,V", po::value<std::string>(&tmp_params.vcf_fn),
            "VCF file with updates of consensus (- for standard output)."
            )
            //
            (
            "pileup,P", po::value<std::string>(&tmp_params.pileup_fn),
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
            "mode,m", po::value<std::string>(&tmp_params.mode_str),
            "Mode: real-time / batch. [batch]")
            //
            (
            "strategy,t", po::value<std::string>(&tmp_params.strategy_str),
            "Strategy for updates: no-updates / majority / stochastic. "
            "[majority]")
            //
            ("allow-amb,a", "Allow updates to ambiguous nucleotides.")
            //
            (
            "min-MQ,q", po::value<int32_t>(&tmp_params.min_mapq),
            min_mq_message.str().c_str())
            //
            (
            "min-BQ,Q", po::value<int32_t>(&tmp_params.min_baseq),
            min_bq_message.str().c_str())
            //
            (
            "ref-weight,w", po::value<int32_t>(&tmp_params.init_ref_weight),
            ref_weight_message.str().c_str())
            //
            (
            "min-coverage,c", po::value<int32_t>(&tmp_params.min_coverage),
            min_coverage_message.str().c_str())
            //
            (
            "majority-threshold,M",
            po::value<double>(&tmp_params.majority_threshold),
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
                if (tmp_params.strategy_str.compare("no-updates") == 0) {
                    tmp_params.strategy = ococo::strategy_t::NO_UPDATES;
                } else if (tmp_params.strategy_str.compare("majority") == 0) {
                    if (vm.count("allow-amb") == 0) {
                        tmp_params.strategy = ococo::strategy_t::MAJORITY;
                    } else {
                        tmp_params.strategy = ococo::strategy_t::MAJORITY_AMB;
                    }
                } else if (tmp_params.strategy_str.compare("stochastic") == 0) {
                    if (vm.count("allow-amb") == 0) {
                        tmp_params.strategy = ococo::strategy_t::STOCHASTIC;
                    } else {
                        tmp_params.strategy = ococo::strategy_t::STOCHASTIC_AMB;
                    }
                } else {
                    ococo::error("Unknown strategy '%s'. Possible strategies "
                                 "are 'majority' and 'stochastic'.\n",
                                 tmp_params.strategy_str.c_str());
                    return EXIT_FAILURE;
                }
            }

            if (vm.count("mode")) {
                if (tmp_params.mode_str.compare("batch") == 0) {
                    tmp_params.mode = ococo::mode_t::BATCH;
                } else if (tmp_params.mode_str.compare("real-time") == 0) {
                    tmp_params.mode = ococo::mode_t::REALTIME;
                } else {
                    ococo::error("Unknown mode '%s'. Possible modes are "
                                 "'batch' and 'real-time'.\n",
                                 tmp_params.mode_str.c_str());
                    return EXIT_FAILURE;
                }
            }

        } catch (po::error &e) {
            std::cout << options_all << "\n";
            ococo::error("%s.\n", e.what());
            return EXIT_FAILURE;
        }

    } catch (std::exception &e) {
        ococo::error("Unhandled Exception: %s.\n", e.what());
        return EXIT_FAILURE;
    }

    /*
     * Read SAM headers.
     */

    ococo::info("Initialing SAM/BAM reader.\n");

    hts_itr_t *iter = nullptr;

    samFile *in = nullptr;
    bam1_t *b = nullptr;
    bam_hdr_t *header = nullptr;

    STATS_T *stats = nullptr;

#ifdef DEBUGGING_MODE
    BOOST_LOG_TRIVIAL(info) << "SAM/BAM reader initialization: reading '"
                            << tmp_params.sam_fn.c_str() << "'.";
#endif

    in = sam_open(tmp_params.sam_fn.c_str(), "r");
    if (in == nullptr) {
        ococo::fatal_error("Problem with opening SAM/BAM file ('%s').\n",
                           tmp_params.sam_fn.c_str());
        main_return_code = -1;
        goto cleaning;
    }
    if ((header = sam_hdr_read(in)) == 0) {
        ococo::fatal_error("SAM/BAM headers are missing or corrupted.\n");
        main_return_code = -1;
        goto cleaning;
    }

    stats = new (std::nothrow) STATS_T(tmp_params, *header);
    if (stats == nullptr || !stats->check_allocation()) {
        ococo::fatal_error("Allocation of the main structure failed.\n");
        main_return_code = -1;
        goto cleaning;
    }

    /*
     * Load FASTA and stats.
     */

    if (!tmp_params.stats_in_fn.empty() && !tmp_params.fasta_in_fn.empty()) {
        ococo::fatal_error("Initial FASTA reference and input statistics "
                           "cannot be used at the same time.\n");
        main_return_code = -1;
        goto cleaning;
    }

    if (!tmp_params.stats_in_fn.empty()) {
        ococo::info("Loading statistics ('%s').\n", tmp_params.stats_in_fn.c_str());
#ifdef DEBUGGING_MODE
        BOOST_LOG_TRIVIAL(info) << "Importing statistics: '" << tmp_params.stats_in_fn
                                << "'.";
#endif
        int error_code = stats->import_stats(tmp_params.stats_in_fn);
        if (error_code != 0) {
            ococo::fatal_error("Import of statistics failed (file '%s').\n",
                               tmp_params.stats_in_fn.c_str());
            main_return_code = -1;
            goto cleaning;
        }
    } else {
#ifdef DEBUGGING_MODE
        BOOST_LOG_TRIVIAL(info) << "No file with statistics provided.";
#endif

        if (!tmp_params.fasta_in_fn.empty()) {
            ococo::info("Loading reference ('%s').\n", tmp_params.fasta_in_fn.c_str());

#ifdef DEBUGGING_MODE
            BOOST_LOG_TRIVIAL(info) << "Loading FASTA: '" << tmp_params.fasta_in_fn
                                    << "'.";
#endif

            int error_code = stats->load_fasta(tmp_params.fasta_in_fn);
            if (error_code != 0) {
                ococo::fatal_error("Loading of FASTA failed (file '%s').\n",
                                   tmp_params.fasta_in_fn.c_str());
                main_return_code = -1;
                goto cleaning;
            }
        }

        else {
            ococo::info("Neither reference, nor statistics provided. Going to "
                        "consider sequence of N's as a reference.\n");
        }
    }

    /*
     * Open VCF file.
     */

    if (tmp_params.vcf_fn.size() > 0) {
        ococo::info("Opening VCF stream ('%s').\n", tmp_params.vcf_fn.c_str());

#ifdef DEBUGGING_MODE
        BOOST_LOG_TRIVIAL(info) << "Open VCF: '" << tmp_params.vcf_fn << "'.";
#endif

        if (tmp_params.vcf_fn == std::string("-")) {
            tmp_params.vcf_file = stdout;
        } else {
            tmp_params.vcf_file = fopen(tmp_params.vcf_fn.c_str(), "w+");
            if (tmp_params.vcf_file == nullptr) {
                ococo::fatal_error("Problem with opening VCF file '%s'.\n",
                                   tmp_params.vcf_fn.c_str());
                main_return_code = -1;
                goto cleaning;
            }
        }

        std::stringstream cmd;
        for (int32_t i = 0; i < argc; i++) {
            cmd << argv[i];
            if (i != argc - 1) {
                cmd << " ";
            }
        }

        char buf[PATH_MAX + 1];
        char *res = realpath(tmp_params.fasta_in_fn.c_str(), buf);
        std::string fasta_full_path;
        if (res) {
            fasta_full_path = std::string(buf);
        } else {
            fasta_full_path = tmp_params.fasta_in_fn;
        }

        stats->print_vcf_header(tmp_params.vcf_file, cmd.str(), fasta_full_path);
    } else {
#ifdef DEBUGGING_MODE
        BOOST_LOG_TRIVIAL(info) << "No VCF file required.";
#endif
    }

    /*
     * Open pileup file.
     */

    if (tmp_params.pileup_fn.size() > 0) {
        ococo::info("Opening pileup stream ('%s').\n", tmp_params.pileup_fn.c_str());

#ifdef DEBUGGING_MODE
        BOOST_LOG_TRIVIAL(info) << "Open pileup: '" << tmp_params.pileup_fn << "'.";
#endif

        if (tmp_params.pileup_fn == std::string("-")) {
            tmp_params.pileup_file = stdout;
        } else {
            tmp_params.pileup_file = fopen(tmp_params.pileup_fn.c_str(), "w+");
            if (tmp_params.pileup_file == nullptr) {
                ococo::fatal_error("Problem with opening pileup file '%s'.\n",
                                   tmp_params.pileup_fn.c_str());
                main_return_code = -1;
                goto cleaning;
            }
        }

    } else {
#ifdef DEBUGGING_MODE
        BOOST_LOG_TRIVIAL(info) << "No pileup file required.";
#endif
    }

    /*
     * Open consensus FASTA file.
     */

    if (tmp_params.fasta_out_fn.size() > 0) {
        tmp_params.fasta_out_file = fopen(tmp_params.fasta_out_fn.c_str(), "w+");

        ococo::info("Opening consensus file ('%s').\n", tmp_params.fasta_out_fn.c_str());

#ifdef DEBUGGING_MODE
        BOOST_LOG_TRIVIAL(info) << "Open FASTA for consensus: '" << tmp_params.fasta_out_fn
                                << "'.";
#endif

        tmp_params.fasta_out_file = fopen(tmp_params.fasta_out_fn.c_str(), "w+");

        if (tmp_params.fasta_out_file == nullptr) {
            ococo::fatal_error(
                "Problem with opening FASTA for consensus: '%s'.\n",
                tmp_params.fasta_out_fn.c_str());
            main_return_code = -1;
            goto cleaning;
        }
    } else {
#ifdef DEBUGGING_MODE
        BOOST_LOG_TRIVIAL(info) << "No FASTA file for consensus required.";
#endif
    }

    /*
     * Process alignments.
     */
    ococo::info("Starting the main loop.\n");

#ifdef DEBUGGING_MODE
    BOOST_LOG_TRIVIAL(info) << "Starting the main loop.";
#endif

    int32_t r;
    b = bam_init1();
    while ((r = sam_read1(in, header, b)) >=
           0) { // read one alignment from `in'
        const char *rname = bam_get_qname(b);
        const uint8_t *seq = bam_get_seq(b);
        const uint8_t *qual = bam_get_qual(b);
        const uint32_t *cigar = bam_get_cigar(b);
        const int32_t n_cigar = b->core.n_cigar;
        //+b->core.l_qname
        const int32_t seqid = b->core.tid;
        const int64_t mappping_pos = b->core.pos;
        const int32_t mapq = b->core.qual;
        const int32_t flags = b->core.flag;

#ifdef DEBUGGING_MODE
        BOOST_LOG_TRIVIAL(debug)
            << "Reading alignment: rname='" << rname << ", chrom=" << seqid
            << ", pos=" << mappping_pos << ", mapq=" << mapq
            << ", flags=" << flags;
#endif

        if ((flags & BAM_FUNMAP) != 0) {
#ifdef DEBUGGING_MODE
            BOOST_LOG_TRIVIAL(debug) << "Discarded: read is not aligned.";
#endif
            continue;
        }

        if (!stats->seq_active[seqid]) {
#ifdef DEBUGGING_MODE
            BOOST_LOG_TRIVIAL(debug)
                << "Discarded: consensus calling is off for this chromosome.";
#endif
            continue;
        }

        if (mapq < stats->params.min_mapq) {
#ifdef DEBUGGING_MODE
            BOOST_LOG_TRIVIAL(debug)
                << "Discarded: mapping quality is too low.";
#endif
            continue;
        }

        int32_t ref_pos = mappping_pos;
        for (int32_t cigar_grp = 0, read_pos = 0; cigar_grp < n_cigar;
             cigar_grp++) {
            const int32_t op = bam_cigar_op(cigar[cigar_grp]);
            const int32_t ol = bam_cigar_oplen(cigar[cigar_grp]);

            const int32_t next_read_pos = read_pos + ol;
            switch (op) {
            case BAM_CMATCH:
            case BAM_CDIFF:
            case BAM_CEQUAL:

                for (; read_pos < next_read_pos; read_pos++, ref_pos++) {
                    const uint8_t nt16 = bam_seqi(seq, read_pos);
                    const uint8_t nt4 = ococo::nt16_nt4[nt16];
                    const char nt256 = ococo::nt16_nt256[nt16];
                    const int32_t bq = qual[read_pos];
                    
                    if (bq != 0xff && bq < (stats->params.min_baseq)) {
#ifdef DEBUGGING_MODE
                        BOOST_LOG_TRIVIAL(trace)
                            << "Omitting base (too low base quality): chrom="
                            << seqid << ", pos=" << ref_pos
                            << ", nucl=" << nt256 << ", quality=" << bq << ".";
#endif
                        continue;
                    }

                    if (nt4 == 0x4) {
#ifdef DEBUGGING_MODE
                        BOOST_LOG_TRIVIAL(trace)
                            << "Omitting base (ambiguous nucleotide): chrom="
                            << seqid << ", pos=" << ref_pos
                            << ", nucl=" << nt256 << ", quality=" << bq << ".";
#endif
                        continue;
                    }

#ifdef DEBUGGING_MODE
                    BOOST_LOG_TRIVIAL(trace)
                        << "Incrementing counter: chrom=" << seqid
                        << ", pos=" << ref_pos << ", nucl=" << nt256
                        << ", quality=" << bq << ". Old state: "
                        << stats->debug_str_counters(seqid, ref_pos) << ",";
#endif

                    stats->seq_stats[seqid][ref_pos] =
                        stats->increment(stats->seq_stats[seqid][ref_pos], nt4);

#ifdef DEBUGGING_MODE
                    BOOST_LOG_TRIVIAL(trace)
                        << "           ...new state: "
                        << stats->debug_str_counters(seqid, ref_pos) << ".";
#endif

                    if (stats->params.mode == ococo::mode_t::REALTIME) {
                        stats->call_consensus_position(tmp_params.vcf_file, tmp_params.pileup_file,
                                                       seqid, ref_pos);
#ifdef DEBUGGING_MODE
                        BOOST_LOG_TRIVIAL(trace)
                            << "Consensus called. New state: "
                            << stats->debug_str_counters(seqid, ref_pos) << ".";
#endif
                    }
                }

                break;

            case BAM_CDEL:
            case BAM_CREF_SKIP:
                ref_pos += ol;
                break;

            case BAM_CSOFT_CLIP:
                read_pos += ol;
                break;

            case BAM_CBACK:
                ref_pos -= ol;
                break;

            case BAM_CINS:
                read_pos += ol;
                break;

            case BAM_CPAD:
            case BAM_CHARD_CLIP:
                break;
            }
        }

#ifdef DEBUGGING_MODE
        BOOST_LOG_TRIVIAL(debug) << "Alignment of '" << rname
                                 << "' incorporated into statistics.";
#endif
    }

    /*
     * Calling final consensus and export stats.
     */

    if (stats->params.mode == ococo::mode_t::BATCH) {
#ifdef DEBUGGING_MODE
        BOOST_LOG_TRIVIAL(info) << "Calling consensus for the entire reference "
                                   "sequence (batch mode).";
#endif
        stats->call_consensus(tmp_params.vcf_file, tmp_params.pileup_file);

        if (tmp_params.fasta_out_fn.size() > 0) {
#ifdef DEBUGGING_MODE
            BOOST_LOG_TRIVIAL(info) << "Saving FASTA: '" << tmp_params.fasta_out_fn
                                    << "'.";
#endif
            int error_code = stats->save_fasta(tmp_params.fasta_out_fn);
            if (error_code != 0) {
                ococo::error("FASTA '%s' could not be saved.\n",
                             tmp_params.fasta_out_fn.c_str());
                main_return_code = -1;
            }
        } else {
#ifdef DEBUGGING_MODE
            BOOST_LOG_TRIVIAL(info) << "FASTA not saved.";
#endif
        }
    }

    if (tmp_params.stats_out_fn.size() > 0) {
#ifdef DEBUGGING_MODE
        BOOST_LOG_TRIVIAL(info) << "Saving statistics: '" << tmp_params.stats_out_fn
                                << "'.";
#endif

        ococo::info("Saving statistics ('%s').\n", tmp_params.stats_out_fn.c_str());

        int error_code = stats->export_stats(tmp_params.stats_out_fn);
        if (error_code != 0) {
            ococo::error("Statistics could not be saved ('%s').\n",
                         tmp_params.stats_out_fn.c_str());
            main_return_code = -1;
        }
    } else {
#ifdef DEBUGGING_MODE
        BOOST_LOG_TRIVIAL(info) << "Statistics not saved.";
#endif
    }

cleaning:

/*
 * Free memory.
 */

#ifdef DEBUGGING_MODE
    BOOST_LOG_TRIVIAL(info) << "Freeing memory.";
#endif
    hts_itr_destroy(iter);
    bam_destroy1(b);
    bam_hdr_destroy(header);

    /*
     * Close files.
     */

    if (in != nullptr) {
        int error_code = sam_close(in);
        if (error_code != 0) {
            ococo::error("SAM file could not be closed.\n");
            main_return_code = -1;
        }
    }

    if (tmp_params.vcf_file != nullptr) {
        int error_code = fclose(tmp_params.vcf_file);
        if (error_code != 0) {
            ococo::error("VCF file could not be closed.\n");
            main_return_code = -1;
        }
    }

    if (tmp_params.pileup_file != nullptr) {
        int error_code = fclose(tmp_params.pileup_file);
        if (error_code != 0) {
            ococo::error("Pileup file could not be closed.\n");
            main_return_code = -1;
        }
    }

    if (tmp_params.fasta_out_file != nullptr) {
        int error_code = fclose(tmp_params.fasta_out_file);
        if (error_code != 0) {
            ococo::error("FASTA consensus file could not be closed.\n");
            main_return_code = -1;
        }
    }

    if (stats != nullptr) {
        delete stats;
    }

    if (main_return_code == 0) {
        ococo::info("Ococo successfully finished. Bye.\n");
    }

#ifdef DEBUGGING_MODE
    BOOST_LOG_TRIVIAL(info) << "Ococo finished.";
#endif

    return main_return_code;
}