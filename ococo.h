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
    
    struct consensus_params_t {
        mode_t mode;
        strategy_t strategy;

        int32_t min_mapq;
        int32_t min_baseq;

        int32_t init_ref_weight;

        FILE *vcf_fo;
        FILE *fasta_cons_fo;

        consensus_params_t():
            mode(BATCH),
            strategy(MAJORITY),
            min_mapq(1),
            min_baseq(0),
            init_ref_weight(2),
            vcf_fo(nullptr),
            fasta_cons_fo(nullptr)
        {}
    };
    
    typedef uint8_t nt4_t;
    typedef uint8_t nt16_t;
    typedef uint8_t nt256_t;
    
    
    /****************
     *** Counters ***
     ****************/
    
    /*
     data type for a single position in the genome (=4 counters)
     - counter structure (merged cells): 0...0|cell_T|cell_G|cell_C|cell_A
     */

    //typedef uint16_t counter_t;
    
    //typedef uint16_t pos_stats_compr_t;
    
    
    struct pos_stats_uncompr_t {
        nt16_t nt16;
        
        int32_t counters[4];
        int32_t sum;
    };
    
    //const int counter_size=sizeof(counter_t);
    //const int cell_bits=4;
    
    //static_assert(cell_bits*4 <= counter_size*8,"Counter type is too small");
    
    //const counter_t cell_maxval=(1<<cell_bits)-1;
    //const counter_t cell_maxval_shifted=cell_maxval>>1;
    /*
    const counter_t counter_norm_mask=cell_maxval_shifted \
    | (cell_maxval_shifted<<cell_bits)        \
    | (cell_maxval_shifted<<2*cell_bits)      \
    | (cell_maxval_shifted<<3*cell_bits);
    */
    
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
    
    
    /***********************************
     *** Consensus calling functions ***
     ***********************************/
    
    // Generate randomly nucleotide with respect to given frequencies.
    inline uint8_t rand_nucl(const pos_stats_uncompr_t &psu){
        if(psu.sum==0){
            return 'N';
        }
        
        const int32_t prefsum[]={
            psu.counters[0],
            psu.counters[0]+psu.counters[1],
            psu.counters[0]+psu.counters[1]+psu.counters[2],
            psu.counters[0]+psu.counters[1]+psu.counters[2]+psu.counters[3]};
        assert(prefsum[3]==psu.sum);

        const int32_t rn=rand() % psu.sum;
        for(int32_t i=0;i<4;i++){
            if (rn < prefsum[i]) {
                return nt4_nt256[i];
            }
        }
        
        return 'n';
    }
    
    
    /********************************
     *** Structure for statistics ***
     ********************************/
    
    template<typename T, int counter_size, int refbase_size>
    struct stats_t {
        static_assert( 8*sizeof(T) >= counter_size + refbase_size, "Too large counter size (does not fit into the main type)." );
        
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
        
        
        /***********
         *** I/O ***
         ***********/
        
        int import_stats(const std::string &stats_fn);
        int export_stats(const std::string &stats_fn) const;
        
        // Call consensus probabilistically.
        int call_consensus();
        int call_consensus_position(int32_t seqid, int64_t pos);
        
        // Loader header from a BAM.
        int load_headers_bam_hdr(const bam_hdr_t &h);
        // Load header from a FASTA file and initialize statistics (to level).
        int load_fasta(const std::string &fasta_fn);
        int save_fasta() const;
        //int load_headers_fai(const std::string &fai_fn);
        
        void print_vcf_header() const;
        void print_vcf_substitution(int32_t seqid, int64_t pos, char old_base, char new_base, const pos_stats_uncompr_t &psu) const;
        
        
        
        /************************
         *** Checking headers ***
         ************************/
        
        
        // Check if everything was initialized.
        bool check_state() const;
        
        // Check if a FAI header corresponds to the stats.
        bool check_headers_fai(const std::string &fai_fn) const;
        
        // Check if a BAM header corresponds to the stats.
        bool check_headers_bam_hdr(const bam_hdr_t &h) const;
        
        /*****************************
         *** Statistics & counters ***
         *****************************/
        
        inline char get_nucl_nt256(int32_t seqid, int64_t pos) const;
        static T compress_position_stats(const pos_stats_uncompr_t &psu);
        static void decompress_position_stats(T psc, pos_stats_uncompr_t &psu);
        static T increment(T psc, nt4_t nt4);
        
        /****************
         *** Debuging ***
         ****************/
        
        void debug_print_counters() const;
        std::string debug_str_counters(int32_t seqid, int64_t pos) const;
    };
    
    
    struct single_seq_serial_t {
        bool    seq_active;
        int64_t seq_len;
        char    seq_name[1000];
        char    seq_comment[1000];
    };
    
    
    
    /**********************************
     *** Manipulating with counters ***
     **********************************/
    
    /*
    inline int _CELL_SHIFT(nt16_t nt16) {
        return cell_bits * seq_nt16_int[nt16];
    }
    
    inline counter_t _COUNTER_CELL_VAL(counter_t counter,nt16_t nt16) {
        return (counter>>_CELL_SHIFT(nt16)) & cell_maxval;
    }
    
    inline counter_t _COUNTER_CELL_SET(counter_t counter,nt16_t nt16,int value) {
        assert(value>>cell_bits==0);
        return
        (
         value << _CELL_SHIFT(nt16) \
         ) | (
              counter
              ^
              (counter & (cell_maxval << _CELL_SHIFT(nt16)))
              );
    }
    
    inline counter_t _COUNTER_CELL_INC_NODIV(counter_t counter,nt16_t nt16) {
        assert(_COUNTER_CELL_VAL(counter,nt16)!=cell_maxval);
        return _COUNTER_CELL_SET (counter,nt16,_COUNTER_CELL_VAL(counter,nt16)+1);
    }
    
    inline counter_t _COUNTER_NORMALIZE(counter_t counter,bool divide) {
        //if(divide) fprintf(stderr,"Shift counter: %04x\n", counter);
        return divide ? (counter >> 1) & counter_norm_mask : counter;
    }
    
    inline counter_t _COUNTER_CELL_INC(counter_t counter,nt16_t nt16) {
        return _COUNTER_CELL_INC_NODIV (
                                        _COUNTER_NORMALIZE (counter,_COUNTER_CELL_VAL(counter,nt16)==cell_maxval),
                                        nt16
                                        );
    }*/
};
    
#include "ococo_impl.h"
    
#endif
