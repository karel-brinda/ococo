#ifndef _OCOCO_H_
#define _OCOCO_H_

#pragma once

#include "ococo_misc.h"

#include <iostream>
#include <sstream>
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <string>
#include <exception>
#include <unistd.h>
#include <cmath>
#include <cassert>
#include <zlib.h>
#include <cinttypes>
#include <ctime>

#include "htslib/sam.h"
#include "htslib/faidx.h"
#include "htslib/kstring.h"
#include "htslib/khash.h"
#include "htslib/kseq.h"

namespace ococo {
    
    KSEQ_INIT(gzFile, gzread);
    
    /*********************
     *** Configuration ***
     *********************/
    
    const int  fasta_line_l = 50;
    const int stats_delim_l = 10;

    enum mode_t {
        BATCH,
        REALTIME
    };
    enum strategy_t {
        STOCHASTIC,
        MAJORITY
    };
    
    typedef uint8_t nt4_t;
    typedef uint8_t nt16_t;
    typedef uint8_t nt256_t;
    
    /**************************
     *** Translation tables ***
     **************************/
    
    static const uint8_t nt256_nt4[] = {
        0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
    };
    
    static const uint8_t nt16_nt4[] = {
        4, 0, 1, 4,
        2, 4, 4, 4,
        3, 4, 4, 4,
        4, 4, 4, 4
    };
    
    static const uint8_t nt256_nt16[] = {
        15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
        15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
        15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
        1, 2, 4, 8, 15,15,15,15, 15,15,15,15, 15, 0 /*=*/,15,15,
        15, 1,14, 2, 13,15,15, 4, 11,15,15,12, 15, 3,15,15,
        15,15, 5, 6,  8,15, 7, 9, 15,10,15,15, 15,15,15,15,
        15, 1,14, 2, 13,15,15, 4, 11,15,15,12, 15, 3,15,15,
        15,15, 5, 6,  8,15, 7, 9, 15,10,15,15, 15,15,15,15,
        
        15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
        15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
        15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
        15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
        15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
        15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
        15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
        15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15
    };
    
    static const uint8_t nt16_nt256[] = "=ACMGRSVTWYHKDBN";
    
    static const uint8_t nt4_nt256[] = "ACGTN";
    
    static const uint8_t nt4_nt16[] = {1,2,4,8,15};

    
    
    
    /******************
     *                *
     *   Structures   *
     *                *
     ******************/
    
    
    /*****************
     *** Auxiliary ***
     *****************/

    struct single_seq_serial_t {
        bool    seq_active;
        int64_t seq_len;
        char    seq_name[1000];
        char    seq_comment[1000];
    };
    
    
    /************************************
     *** Single position uncompressed ***
     ************************************/
    
    struct pos_stats_uncompr_t {
        nt16_t nt16;
        
        int32_t counters[4];
        int32_t sum;
    };
    
    
    /****************************
     *** Consensus parameters ***
     ****************************/
    
    struct consensus_params_t {
        mode_t mode;
        strategy_t strategy;
        
        int32_t min_mapq;
        int32_t min_baseq;
        
        int32_t init_ref_weight;
        
        FILE *vcf_fo;
        FILE *fasta_cons_fo;
        
        consensus_params_t();
    };
    
    
    /***********************
     *** Main statistics ***
     ***********************/
    
    template<typename T, int counter_size, int refbase_size>
    struct stats_t {
        static_assert( 8*sizeof(T) >= 4*counter_size + refbase_size, "Too large counter size (does not fit into the main type)." );
        
        int32_t      n_seqs;
        bool         *seq_active;
        int64_t      *seq_len;
        std::string  *seq_name;
        std::string  *seq_comment;
        T            **seq_stats;
        
        consensus_params_t params;
        
        //stats_t();
        stats_t(consensus_params_t params,bam_hdr_t &h);
        ~stats_t();
        
        
        /*******
         * I/O *
         *******/
        
        int import_stats(const std::string &stats_fn);
        int export_stats(const std::string &stats_fn) const;
        
        // Call consensus probabilistically.
        int call_consensus();
        int call_consensus_position(int32_t seqid, int64_t pos);
        
        // Loader header from a BAM.
        int load_headers_bam_hdr(const bam_hdr_t &h);
        // Load header and data from a FASTA file and initialize statistics.
        int load_fasta(const std::string &fasta_fn);
        int save_fasta() const;
        
        int print_vcf_header(std::string cmd, std::string fasta) const;
        int print_vcf_substitution(int32_t seqid, int64_t pos, char old_base, char new_base, const pos_stats_uncompr_t &psu) const;
        
        
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
    
    
    
    /***********************************
     *                                 *
     *   Consensus calling functions   *
     *                                 *
     ***********************************/
    
    char cons_call_stoch(const pos_stats_uncompr_t &psu);
    char cons_call_stoch_amb(const pos_stats_uncompr_t &psu);
    char cons_call_maj(const pos_stats_uncompr_t &psu);
    char cons_call_maj_amb(const pos_stats_uncompr_t &psu);
};
    
#include "ococo_impl.h"
    
#endif
