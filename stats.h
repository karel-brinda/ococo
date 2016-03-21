#pragma once

#include "ococo.h"

#include "htslib/faidx.h"
#include "htslib/khash.h"
#include "htslib/kseq.h"
#include "htslib/kstring.h"
#include "htslib/sam.h"

/***********************
 *** Main statistics ***
 ***********************/

namespace ococo {

template <typename T, int counter_size, int refbase_size> struct stats_t {
    static_assert(8 * sizeof(T) >= 4 * counter_size + refbase_size,
                  "Too large counter size (does not fit into the main type).");

    int32_t n_seqs;
    bool *seq_active;
    int64_t *seq_len;
    std::string *seq_name;
    std::string *seq_comment;
    T **seq_stats;

    consensus_params_t params;

    // stats_t();
    stats_t(consensus_params_t params, bam_hdr_t &h);
    ~stats_t();

    /*******
     * I/O *
     *******/

    int import_stats(const std::string &stats_fn);
    int export_stats(const std::string &stats_fn) const;

    // Call consensus probabilistically.
    int call_consensus(FILE *vcf_file, FILE *pileup_file);
    int call_consensus_position(FILE *vcf_file, FILE *pileup_file,
                                int32_t seqid, int64_t pos);

    // Loader header from a BAM.
    int load_headers_bam_hdr(const bam_hdr_t &h);
    // Load header and data from a FASTA file and initialize statistics.
    int load_fasta(const std::string &fasta_fn);
    int save_fasta(const std::string &fasta_fn) const;

    int print_vcf_header(FILE *vcf_file, std::string cmd,
                         std::string fasta) const;
    int print_vcf_substitution(FILE *vcf_file, int32_t seqid, int64_t pos,
                               char old_base, char new_base,
                               const pos_stats_uncompr_t &psu) const;

    int print_pileup_line(FILE *pileup_file, int32_t seqid, int64_t pos,
                          const pos_stats_uncompr_t &psu) const;

    /*************************
     * Statistics & counters *
     *************************/

    inline int set_nucl_nt256(int32_t seqid, int64_t pos, const char &nt256);
    inline int get_nucl_nt256(int32_t seqid, int64_t pos, char &nt256) const;

    static T compress_position_stats(const pos_stats_uncompr_t &psu);
    static void decompress_position_stats(T psc, pos_stats_uncompr_t &psu);
    static T increment(T psc, nt4_t nt4);

    /***********************
     * Debuging & checking *
     ***********************/

    // Check if everything was initialized.
    bool check_allocation() const;

    // Check if a BAM header corresponds to the stats.
    bool check_headers_bam_hdr(const bam_hdr_t &h) const;

    void debug_print_counters() const;

    std::string debug_str_counters(int32_t seqid, int64_t pos) const;
};

}

