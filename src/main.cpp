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

    ococo::consensus_params_t tmp_params = ococo::consensus_params_t(argc, argv);

/*
 * Parse command-line parameters.
 */

#ifdef DEBUGGING_MODE
    BOOST_LOG_TRIVIAL(info) << "Parsing command-line parameters.";
#endif


    /*
     * Read SAM headers.
     */

    ococo::info("Initialing SAM/BAM reader.\n");

    hts_itr_t *iter = nullptr;

    bam1_t *b = nullptr;
    bam_hdr_t *header = nullptr;

    STATS_T *stats = nullptr;

#ifdef DEBUGGING_MODE
    BOOST_LOG_TRIVIAL(info) << "SAM/BAM reader initialization: reading '"
                            << tmp_params.sam_fn.c_str() << "'.";
#endif

    tmp_params.sam_file = sam_open(tmp_params.sam_fn.c_str(), "r");
    if (tmp_params.sam_file == nullptr) {
        ococo::fatal_error("Problem with opening SAM/BAM file ('%s').\n",
                           tmp_params.sam_fn.c_str());
        main_return_code = -1;
        goto cleaning;
    }
    if ((header = sam_hdr_read(tmp_params.sam_file)) == 0) {
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
    while ((r = sam_read1(tmp_params.sam_file, header, b)) >=
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