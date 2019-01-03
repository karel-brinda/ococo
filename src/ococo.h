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

#pragma once

#include <cstdlib>
#include "consensus.h"
#include "io.h"
#include "params.h"
#include "stats.h"

namespace ococo {

/*! @struct
    @abstract                      Structure for metadata for 1 sequence.
    @field correctly_initialized   Flag for a correct initialization from the
                                   input files.
    @field return_code             Return code of the caller. 0 if everything
                                   ok.
    @field b                       Structure for one alignment.
    @field header                  SAM/BAM header.
    @field stats                   Pileup statistics.
    @field params                  Program arguments.
    @field t_real                  Initial timestamp.
*/
template <typename T, int counter_size, int refbase_size>
struct ococo_t {
    bool correctly_initialized;
    int return_code;

    bam1_t *b;
    bam_hdr_t *header;

    stats_t<T, counter_size, refbase_size> *stats;

    params_t *params;
    double t_real;

    /*! @func
        @abstract  Open all files and load headers.
    */
    ococo_t(params_t *params_);
    ~ococo_t();

    /*! @func
        @abstract  Check whether the alignment passes filters.
    */
    bool check_read(int32_t seqid, int32_t flags, int32_t mapq);

    /*! @func
        @abstract  Run calling.
    */
    void run();
};

template <typename T, int counter_size, int refbase_size>
ococo_t<T, counter_size, refbase_size>::ococo_t(params_t *params_)
    : params(params_) {
    /*
     * Read SAM headers.
     */

    t_real = realtime();
    ococo::info("Initializing the SAM/BAM reader.\n");

    correctly_initialized = true;
    return_code           = EXIT_SUCCESS;

    b      = nullptr;
    header = nullptr;
    stats  = nullptr;

    params->in_sam_file = sam_open(params->in_sam_fn.c_str(), "r");
    if (params->in_sam_file == nullptr) {
        ococo::fatal_error("Problem with opening the SAM/BAM file ('%s').\n",
                           params->in_sam_fn.c_str());
        correctly_initialized = false;
        return;
    }

    if ((header = sam_hdr_read(params->in_sam_file)) == 0) {
        ococo::fatal_error("SAM/BAM headers are missing or corrupted.\n");
        correctly_initialized = false;
        return;
    }

    stats = new (std::nothrow)
        stats_t<T, counter_size, refbase_size>(params, *header);
    if (stats == nullptr || !stats->check_allocation()) {
        ococo::fatal_error(
            "Allocation of the table for nucleotide statistics failed.\n");
        correctly_initialized = false;
        return;
    }

    /*
     * Open output SAM stream
     */

    if (params->out_sam_fn.size() > 0) {
        params->out_sam_file = sam_open(params->out_sam_fn.c_str(), "w");
        if (params->out_sam_file == nullptr) {
            ococo::fatal_error(
                "Problem with opening the SAM/BAM file ('%s').\n",
                params->out_sam_fn.c_str());
            correctly_initialized = false;
            return;
        }
    }

    /*
     * Load FASTA and stats.
     */

    if (!params->in_stats_fn.empty() && !params->in_fasta_fn.empty()) {
        ococo::fatal_error(
            "Initial FASTA reference and input statistics "
            "cannot be used at the same time.\n");
        correctly_initialized = false;
        return;
    }

    if (!params->in_stats_fn.empty()) {
        ococo::info("Loading previously saved statistics ('%s').\n",
                    params->in_stats_fn.c_str());

        int error_code = stats->import_stats(params->in_stats_fn);
        if (error_code != 0) {
            ococo::fatal_error("Import of the statistics failed (file '%s').\n",
                               params->in_stats_fn.c_str());
            correctly_initialized = false;
            return;
        }
    } else {
        if (!params->in_fasta_fn.empty()) {
            ococo::info("Loading the reference ('%s').\n",
                        params->in_fasta_fn.c_str());

            int error_code = stats->load_fasta(params->in_fasta_fn);
            if (error_code != 0) {
                ococo::fatal_error("Loading of the FASTA failed (file '%s').\n",
                                   params->in_fasta_fn.c_str());
                correctly_initialized = false;
                return;
            }
        }

        else {
            ococo::info(
                "Neither reference, nor statistics provided. Going "
                "to use a sequence of N's as the reference.\n");
        }
    }

    /*
     * Open VCF file.
     */

    if (params->out_vcf_fn.size() > 0) {
        ococo::info("Opening the VCF stream ('%s').\n",
                    params->out_vcf_fn.c_str());

        if (params->out_vcf_fn == std::string("-")) {
            params->out_vcf_file = stdout;
        } else {
            params->out_vcf_file = fopen(params->out_vcf_fn.c_str(), "w+");
            if (params->out_vcf_file == nullptr) {
                ococo::fatal_error("Problem with opening the VCF file '%s'.\n",
                                   params->out_vcf_fn.c_str());
                correctly_initialized = false;
                return;
            }
        }

        char buf[PATH_MAX + 1];
        char *res = realpath(params->in_fasta_fn.c_str(), buf);
        std::string fasta_full_path;
        if (res) {
            fasta_full_path = std::string(buf);
        } else {
            fasta_full_path = params->in_fasta_fn;
        }

        print_vcf_header(params->out_vcf_file);
    }

    /*
     * Open pileup file.
     */

    if (params->out_pileup_fn.size() > 0) {
        ococo::info("Opening the pileup stream ('%s').\n",
                    params->out_pileup_fn.c_str());

        if (params->out_pileup_fn == std::string("-")) {
            params->out_pileup_file = stdout;
        } else {
            params->out_pileup_file =
                fopen(params->out_pileup_fn.c_str(), "w+");
            if (params->out_pileup_file == nullptr) {
                ococo::fatal_error(
                    "Problem with opening the pileup file '%s'.\n",
                    params->out_pileup_fn.c_str());
                correctly_initialized = false;
                return;
            }
        }
    }

    /*
     * Open consensus FASTA file.
     */

    if (params->out_fasta_fn.size() > 0) {
        ococo::info("Opening the consensus FASTA file ('%s').\n",
                    params->out_fasta_fn.c_str());
        params->out_fasta_file = fopen(params->out_fasta_fn.c_str(), "w+");

        if (params->out_fasta_file == nullptr) {
            ococo::fatal_error(
                "Problem with opening the consensus FASTA file: '%s'.\n",
                params->out_fasta_fn.c_str());
            correctly_initialized = false;
            return;
        }
    }

    /*
     * Open log file.
     */

    if (params->out_log_fn.size() > 0) {
        ococo::info("Opening the log file ('%s').\n",
                    params->out_log_fn.c_str());
        params->out_log_file = fopen(params->out_log_fn.c_str(), "w+");

        if (params->out_log_file == nullptr) {
            ococo::fatal_error("Problem with opening the log file: '%s'.\n",
                               params->out_log_fn.c_str());
            correctly_initialized = false;
            return;
        }
    }
}

/*
//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
*/

template <typename T, int counter_size, int refbase_size>
bool ococo_t<T, counter_size, refbase_size>::check_read(int32_t seqid,
                                                        int32_t flags,
                                                        int32_t mapq) {
    if ((flags & BAM_FUNMAP) != 0) {
        return false;
    }

    if (!stats->seq_active[seqid]) {
        return false;
    }

    if (mapq < stats->params->min_mapq) {
        return false;
    }

    return true;
}

template <typename T, int counter_size, int refbase_size>
void ococo_t<T, counter_size, refbase_size>::run() {
    /*
     * Process alignments.
     *
     * Notes
     * - if reading bam fails - a critical error
     * - if writing bam header fails - a critical error
     * - if writing bam fails - a non-critical error
     */
    ococo::info("Starting the main loop.\n");

    int32_t return_value;
    b              = bam_init1();
    int64_t n_upd0 = 0;
    int64_t i_read = 0;

    if (stats->params->out_sam_file != nullptr) {
        int error_code = sam_hdr_write(stats->params->out_sam_file, header);
        if (error_code != 0) {
            return_code = EXIT_FAILURE;
            ococo::error("Construction of the SAM header failed (error %d)",
                         error_code);
            return;
        }
    }

    while ((return_value = sam_read1(params->in_sam_file, header, b)) >= 0) {
        /* SAM/BAM variables */
        const char *rname         = bam_get_qname(b);
        const uint8_t *seq        = bam_get_seq(b);
        const uint8_t *qual       = bam_get_qual(b);
        const uint32_t *cigar     = bam_get_cigar(b);
        const int32_t n_cigar     = b->core.n_cigar;
        const int32_t seqid       = b->core.tid;
        const int64_t mapping_pos = b->core.pos;
        const int32_t mapq        = b->core.qual;
        const int32_t flags       = b->core.flag;

        /* filtration on the alignment level */
        bool read_ok = check_read(seqid, flags, mapq);
        if (!read_ok) {
            continue;
        }

        /* coverage statistics for the region of the alignment */
        int32_t low_cov_thres = stats->params->coverage_filter;
        int32_t npos_low_cov  = 0;
        int32_t npos_high_cov = 0;
        int32_t pseudo_rlen   = 0;

        if (low_cov_thres < 0) {
            low_cov_thres = 424242;
        }

        /* iteration over individual bases */
        int32_t ref_pos = mapping_pos;
        for (int32_t cigar_grp = 0, read_pos = 0; cigar_grp < n_cigar;
             cigar_grp++) {
            const int32_t op = bam_cigar_op(cigar[cigar_grp]);
            const int32_t ol = bam_cigar_oplen(cigar[cigar_grp]);

            const int32_t next_read_pos = read_pos + ol;
            switch (op) {
                case BAM_CMATCH:
                case BAM_CDIFF:
                case BAM_CEQUAL:

                    for (; read_pos < next_read_pos; ++read_pos, ++ref_pos) {
                        /* filtration on the level of base */
                        const uint8_t nt16 = bam_seqi(seq, read_pos);
                        const uint8_t nt4  = ococo::nt16_nt4[nt16];
                        const int32_t bq   = qual[read_pos];

                        if (bq != 0xff && bq < (stats->params->min_baseq)) {
                            continue;
                        }

                        if (nt4 == 0x4) {
                            continue;
                        }

                        /* updating counters */
                        pos_stats_uncompr_t psu;
                        psu.decompress<T, counter_size, refbase_size>(
                            stats->seq_stats[seqid][ref_pos]);
                        psu.increment(nt4);

                        /* updating coverage statistics */
                        if (psu.sum - 1 < low_cov_thres) {
                            ++npos_low_cov;
                        } else {
                            if (2 * (psu.sum - 1) > 3 * low_cov_thres) {
                                ++npos_high_cov;
                            }
                        }
                        ++pseudo_rlen;

                        /* consensus calling for the current position */
                        if (stats->params->mode == ococo::mode_t::REALTIME) {
                            stats->call_consensus_position_uncompressed(
                                params->out_vcf_file, params->out_pileup_file,
                                seqid, ref_pos, psu);
                        }

                        /* compressing the counters a putting them back to the
                         * statistics */
                        stats->seq_stats[seqid][ref_pos] =
                            psu.compress<T, counter_size, refbase_size>();
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
            }  // switch (op)

        }  // for (int32_t cigar_grp

        /* Filtering the alignment based on coverage. */
        if (stats->params->out_sam_file != nullptr) {
            // at least 5% pos. low coverage => print read
            // if(npos_low_cov * 20 >= pseudo_rlen){
            if (npos_low_cov > 0 && npos_high_cov < 0.5 * pseudo_rlen) {
                int error_code =
                    sam_write1(stats->params->out_sam_file, header, b);
                if (error_code != 0) {
                    return_code = EXIT_FAILURE;
                    ococo::error("Writing SAM failed (error %d)", error_code);
                    break;
                }
            }
        }

        /* Logging the number of updates from this alignment. */
        if (stats->params->out_log_file != nullptr) {
            fprintf(stats->params->out_log_file,
                    "%" PRIu64 "\t%s\t%" PRIu64 "\n", i_read, rname,
                    stats->params->n_upd - n_upd0);
            n_upd0 = stats->params->n_upd;
        }

        i_read += 1;
    }  // while ((r = sam_read1

    if (return_value < -1) {
        ococo::error("Truncated BAM stream (error %" PRId32
                     "). Ococo will still try to print results.\n",
                     return_value);
        return_code = EXIT_FAILURE;
    }

    /*
     * Call final consensus and export stats.
     */

    if (stats->params->mode == ococo::mode_t::BATCH) {
        stats->call_consensus(params->out_vcf_file, params->out_pileup_file);

        if (params->out_fasta_fn.size() > 0) {
            int error_code = stats->save_fasta(params->out_fasta_fn);
            if (error_code != 0) {
                ococo::error("FASTA '%s' could not be saved.\n",
                             params->out_fasta_fn.c_str());
                return_code = EXIT_FAILURE;
            }
        }
    }

    if (params->out_stats_fn.size() > 0) {
        ococo::info("Saving the obtained statistics ('%s').\n",
                    params->out_stats_fn.c_str());

        int error_code = stats->export_stats(params->out_stats_fn);
        if (error_code != 0) {
            ococo::error("THe statistics could not be saved ('%s').\n",
                         params->out_stats_fn.c_str());
            return_code = EXIT_FAILURE;
        }
    }
}

/*
//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
*/

template <typename T, int counter_size, int refbase_size>
ococo_t<T, counter_size, refbase_size>::~ococo_t() {
    bam_destroy1(b);
    bam_hdr_destroy(header);

    if (params->in_sam_file != nullptr) {
        int error_code = sam_close(params->in_sam_file);
        if (error_code != 0) {
            ococo::error("Input SAM file could not be closed.\n");
            return_code = -1;
        }
    }

    if (params->out_sam_file != nullptr) {
        int error_code = sam_close(params->out_sam_file);
        if (error_code != 0) {
            ococo::error("Output SAM file could not be closed.\n");
            return_code = -1;
        }
    }

    if (params->out_vcf_file != nullptr) {
        int error_code = fclose(params->out_vcf_file);
        if (error_code != 0) {
            ococo::error("Output VCF file could not be closed.\n");
            return_code = -1;
        }
    }

    if (params->out_pileup_file != nullptr) {
        int error_code = fclose(params->out_pileup_file);
        if (error_code != 0) {
            return_code = error_code;
            ococo::error("Output pileup file could not be closed.\n");
            return_code = -1;
        }
    }

    if (params->out_fasta_file != nullptr) {
        int error_code = fclose(params->out_fasta_file);
        if (error_code != 0) {
            ococo::error("Output FASTA consensus file could not be closed.\n");
            return_code = -1;
        }
    }

    if (params->out_log_file != nullptr) {
        int error_code = fclose(params->out_log_file);
        if (error_code != 0) {
            ococo::warning("Log file could not be closed.\n");
        }
    }

    if (stats != nullptr) {
        delete stats;
    }

    if (return_code == EXIT_SUCCESS && correctly_initialized == true) {
        ococo::info("Ococo successfully finished. Bye.\n");
        ococo::info("%.3f sec; CPU: %.3f sec\n", realtime() - t_real,
                    cputime());
    }
}
}  // namespace ococo
