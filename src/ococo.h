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

#include "bamfiles.h"
#include "debugging.h"
#include "io.h"
#include "logfile.h"
#include "params.h"
#include "pileupfile.h"
#include "stats.h"
#include "vcffile.h"

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
template <typename T>
struct Ococo {
    bool correctly_initialized;
    int return_code;

    stats_t<T> *stats;

    params_t *params;
    double t_real;

    /*! @func
        @abstract  Open all files and load headers.
    */
    Ococo(params_t *params_) : params(params_) {
        /*
         * Read SAM headers.
         */

        t_real = realtime();
        info("Initializing the SAM/BAM reader.\n");

        correctly_initialized = true;
        return_code           = EXIT_SUCCESS;

        stats = nullptr;

        /*
         * Load FASTA and stats.
         */

        if (!params->in_stats_fn.empty() && !params->in_fasta_fn.empty()) {
            fatal_error(
                "Initial FASTA reference and input statistics "
                "cannot be used at the same time.\n");
            correctly_initialized = false;
            return;
        }

        if (!params->in_stats_fn.empty()) {
            info("Loading previously saved statistics ('%s').\n",
                 params->in_stats_fn.c_str());

            int error_code = stats->import_stats(params->in_stats_fn);
            if (error_code != 0) {
                fatal_error("Import of the statistics failed (file '%s').\n",
                            params->in_stats_fn.c_str());
                correctly_initialized = false;
                return;
            }
        } else {
            if (!params->in_fasta_fn.empty()) {
                info("Loading the reference ('%s').\n",
                     params->in_fasta_fn.c_str());

                int error_code = stats->load_fasta(params->in_fasta_fn);
                if (error_code != 0) {
                    fatal_error("Loading of the FASTA failed (file '%s').\n",
                                params->in_fasta_fn.c_str());
                    correctly_initialized = false;
                    return;
                }
            }

            else {
                info(
                    "Neither reference, nor statistics provided. Going "
                    "to use a sequence of N's as the reference.\n");
            }
        }

        if (params->out_vcf_fn.size() > 0) {
            char buf[PATH_MAX + 1];
            char *res = realpath(params->in_fasta_fn.c_str(), buf);
            std::string fasta_full_path;
            if (res) {
                fasta_full_path = std::string(buf);
            } else {
                fasta_full_path = params->in_fasta_fn;
            }
        }

        /*
         * Open consensus FASTA file.
         */

        if (params->out_fasta_fn.size() > 0) {
            info("Opening the consensus FASTA file ('%s').\n",
                 params->out_fasta_fn.c_str());
            params->out_fasta_file = fopen(params->out_fasta_fn.c_str(), "w+");

            if (params->out_fasta_file == nullptr) {
                fatal_error(
                    "Problem with opening the consensus FASTA file: '%s'.\n",
                    params->out_fasta_fn.c_str());
                correctly_initialized = false;
                return;
            }
        }
    }

    ~Ococo() {
        if (params->out_fasta_file != nullptr) {
            int error_code = fclose(params->out_fasta_file);
            if (error_code != 0) {
                error("Output FASTA consensus file could not be closed.\n");
                return_code = -1;
            }
        }

        if (stats != nullptr) {
            delete stats;
        }

        if (return_code == EXIT_SUCCESS && correctly_initialized == true) {
            info("Ococo successfully finished. Bye.\n");
            info("%.3f sec; CPU: %.3f sec\n", realtime() - t_real, cputime());
        }
    }

    /*! @func
        @abstract  Check whether the alignment passes filters.
    */
    bool check_read(int32_t seqid, int32_t flags, int32_t mapq) {
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

    /*! @func
        @abstract  Run calling.
    */
    void run() {
        /*
         * Process alignments.
         *
         * Notes
         * - if reading bam fails - a critical error
         * - if writing bam header fails - a critical error
         * - if writing bam fails - a non-critical error
         */

        VcfFile vcf_file(params->out_vcf_fn);
        PileupFile pileup_file(params->out_pileup_fn);
        LogFile log_file(params->out_log_fn);
        BamFiles bam(params->in_sam_fn, params->out_sam_fn);

        stats = new (std::nothrow) stats_t<T>(params, *bam.header);
        if (stats == nullptr) {
            fatal_error(
                "Allocation of the table for nucleotide statistics failed.\n");
            correctly_initialized = false;
            return;
        }

        info("Starting the main loop.\n");

        int32_t return_value;
        int64_t n_upd0 = 0;
        int64_t i_read = 0;

        while ((return_value = bam.read_alignment()) >= 0) {
            /* filtration on the alignment level */
            bool read_ok = check_read(bam.seqid, bam.flags, bam.mapq);
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

            // std::cerr << __PRETTY_FUNCTION__ << *rname << std::endl;

            /* iteration over individual bases */
            int32_t ref_pos = bam.mapping_pos;
            for (int32_t cigar_grp = 0, read_pos = 0; cigar_grp < bam.n_cigar;
                 cigar_grp++) {
                const int32_t op = bam_cigar_op(bam.cigar[cigar_grp]);
                const int32_t ol = bam_cigar_oplen(bam.cigar[cigar_grp]);

                const int32_t next_read_pos = read_pos + ol;
                switch (op) {
                    case BAM_CMATCH:
                    case BAM_CDIFF:
                    case BAM_CEQUAL:

                        for (; read_pos < next_read_pos;
                             ++read_pos, ++ref_pos) {
                            /* filtration on the level of base */
                            const uint8_t nt16 = bam_seqi(bam.seq, read_pos);
                            const uint8_t nt4  = nt16_nt4[nt16];
                            const int32_t bq   = bam.qual[read_pos];

                            if (bq != 0xff && bq < (stats->params->min_baseq)) {
                                continue;
                            }

                            if (nt4 == 0x4) {
                                continue;
                            }

                            /* updating counters */
                            pos_stats_uncompr_t psu;

                            // std::cerr << "\n" << read_pos << std::endl;
                            //_print_pos_stats(stats->seq_stats[seqid][ref_pos]);
                            psu.decompress(
                                stats->seq_stats[bam.seqid][ref_pos]);
                            //_print_pos_stats<T>(psu);
                            psu.increment(nt4);
                            // std::cerr << "       incr " << nt4_nt256[nt4]
                            //          << " at position " << ref_pos <<
                            //          std::endl;

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
                            if (stats->params->mode == mode_t::REALTIME) {
                                stats->call_consensus_position(
                                    vcf_file, pileup_file, bam.seqid, ref_pos,
                                    psu);
                            }

                            /* compressing the counters a putting them back to
                             * the statistics */
                            psu.compress(stats->seq_stats[bam.seqid][ref_pos]);
                            //_print_pos_stats(stats->seq_stats[seqid][ref_pos]);
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
            if (npos_low_cov > 0 && npos_high_cov < 0.5 * pseudo_rlen) {
                // at least 5% pos. low coverage => print read
                // if(npos_low_cov * 20 >= pseudo_rlen){
                int error_code = bam.print_alignment();
                if (error_code != 0) {
                    return_code = EXIT_FAILURE;
                    error("Writing SAM failed (error %d)", error_code);
                    break;
                }
            }

            /* Logging the number of updates from this alignment. */
            log_file.print(i_read, bam.rname, stats->params->n_upd - n_upd0);
            n_upd0 =
                stats->params->n_upd;  // todo: count automatically in the log

            i_read += 1;
        }  // while ((r = sam_read1

        if (return_value < -1) {
            error("Truncated BAM stream (error %" PRId32
                  "). Ococo will still try to print results.\n",
                  return_value);
            return_code = EXIT_FAILURE;
        }

        /*
         * Call final consensus and export stats.
         */

        if (params->out_stats_fn.size() > 0) {
            info("Saving the obtained statistics ('%s').\n",
                 params->out_stats_fn.c_str());

            int error_code = stats->export_stats(params->out_stats_fn);
            if (error_code != 0) {
                error("THe statistics could not be saved ('%s').\n",
                      params->out_stats_fn.c_str());
                return_code = EXIT_FAILURE;
            }
        }
    }
};

}  // namespace ococo
