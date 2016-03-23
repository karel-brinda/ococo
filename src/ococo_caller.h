#pragma once

#include "ococo_stats.h"
#include "ococo_params.h"


#ifdef DEBUGGING_MODE

#ifndef DEBUGGING_SEVERITY
#define DEBUGGING_SEVERITY trace
#endif

#include <boost/log/core.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/trivial.hpp>

namespace logging = boost::log;

//void boost_logging_init() {
//    logging::core::get()->set_filter(logging::trivial::severity >=
//                                     logging::trivial::DEBUGGING_SEVERITY);
//}

#endif


namespace ococo{
    
    template <typename T, int counter_size, int refbase_size>
    struct caller_t {
        bool correctly_initialized;
        int return_code;
        
        hts_itr_t *iter;
        
        bam1_t *b;
        bam_hdr_t *header;
        
        stats_t<T, counter_size, refbase_size> *stats;
        
        params_t &params;
        
        caller_t(params_t &params);
        ~caller_t();
        
        bool check_read(int32_t seqid, int32_t flags, int32_t mapq);
        void run();
    };
    
    template <typename T, int counter_size, int refbase_size>
    caller_t<T, counter_size, refbase_size>::caller_t(params_t &params):
    params(params)
    {
        
        
        /*
         * Read SAM headers.
         */
        
        ococo::info("Initialing SAM/BAM reader.\n");
        
        correctly_initialized=true;
        
        iter = nullptr;
        b = nullptr;
        header = nullptr;
        stats = nullptr;
        
        if(DEBUGGING_MODE){
            BOOST_LOG_TRIVIAL(info) << "SAM/BAM reader initialization: reading '"
            << params.sam_fn.c_str() << "'.";
        }
        
        params.sam_file = sam_open(params.sam_fn.c_str(), "r");
        if (params.sam_file == nullptr) {
            ococo::fatal_error("Problem with opening SAM/BAM file ('%s').\n",
                               params.sam_fn.c_str());
            correctly_initialized=false;
            return;
        }
        
        if ((header = sam_hdr_read(params.sam_file)) == 0) {
            ococo::fatal_error("SAM/BAM headers are missing or corrupted.\n");
            correctly_initialized=false;
            return;
        }
        
        stats = new (std::nothrow) stats_t<T, counter_size, refbase_size>(params, *header);
        if (stats == nullptr || !stats->check_allocation()) {
            ococo::fatal_error("Allocation of the main structure failed.\n");
            correctly_initialized=false;
            return;
        }
        
        /*
         * Load FASTA and stats.
         */
        
        if (!params.stats_in_fn.empty() && !params.fasta_in_fn.empty()) {
            ococo::fatal_error("Initial FASTA reference and input statistics "
                               "cannot be used at the same time.\n");
            correctly_initialized=false;
            return;
        }
        
        if (!params.stats_in_fn.empty()) {
            ococo::info("Loading statistics ('%s').\n", params.stats_in_fn.c_str());
            if(DEBUGGING_MODE){
                BOOST_LOG_TRIVIAL(info) << "Importing statistics: '" << params.stats_in_fn
                << "'.";
            }
            
            int error_code = stats->import_stats(params.stats_in_fn);
            if (error_code != 0) {
                ococo::fatal_error("Import of statistics failed (file '%s').\n",
                                   params.stats_in_fn.c_str());
                correctly_initialized=false;
                return;
            }
        } else {
            if(DEBUGGING_MODE){
                
                BOOST_LOG_TRIVIAL(info) << "No file with statistics provided.";
            }
            
            if (!params.fasta_in_fn.empty()) {
                ococo::info("Loading reference ('%s').\n", params.fasta_in_fn.c_str());
                
                if(DEBUGGING_MODE){
                    BOOST_LOG_TRIVIAL(info) << "Loading FASTA: '" << params.fasta_in_fn
                    << "'.";
                }
                
                int error_code = stats->load_fasta(params.fasta_in_fn);
                if (error_code != 0) {
                    ococo::fatal_error("Loading of FASTA failed (file '%s').\n",
                                       params.fasta_in_fn.c_str());
                    correctly_initialized=false;
                    return;
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
        
        if (params.vcf_fn.size() > 0) {
            ococo::info("Opening VCF stream ('%s').\n", params.vcf_fn.c_str());
            
            if(DEBUGGING_MODE){
                BOOST_LOG_TRIVIAL(info) << "Open VCF: '" << params.vcf_fn << "'.";
            }
            
            if (params.vcf_fn == std::string("-")) {
                params.vcf_file = stdout;
            } else {
                params.vcf_file = fopen(params.vcf_fn.c_str(), "w+");
                if (params.vcf_file == nullptr) {
                    ococo::fatal_error("Problem with opening VCF file '%s'.\n",
                                       params.vcf_fn.c_str());
                    correctly_initialized=false;
                    return;
                }
            }
            
            char buf[PATH_MAX + 1];
            char *res = realpath(params.fasta_in_fn.c_str(), buf);
            std::string fasta_full_path;
            if (res) {
                fasta_full_path = std::string(buf);
            } else {
                fasta_full_path = params.fasta_in_fn;
            }
            
            stats->print_vcf_header(params.vcf_file, params.command, fasta_full_path);
        } else {
            if(DEBUGGING_MODE){
                BOOST_LOG_TRIVIAL(info) << "No VCF file required.";
            }
        }
        
        /*
         * Open pileup file.
         */
        
        if (params.pileup_fn.size() > 0) {
            ococo::info("Opening pileup stream ('%s').\n", params.pileup_fn.c_str());
            
            if(DEBUGGING_MODE){
                BOOST_LOG_TRIVIAL(info) << "Open pileup: '" << params.pileup_fn << "'.";
            }
            
            if (params.pileup_fn == std::string("-")) {
                params.pileup_file = stdout;
            } else {
                params.pileup_file = fopen(params.pileup_fn.c_str(), "w+");
                if (params.pileup_file == nullptr) {
                    ococo::fatal_error("Problem with opening pileup file '%s'.\n",
                                       params.pileup_fn.c_str());
                    correctly_initialized=false;
                    return;
                }
            }
            
        } else {
            if(DEBUGGING_MODE){
                BOOST_LOG_TRIVIAL(info) << "No pileup file required.";
            }
        }
        
        /*
         * Open consensus FASTA file.
         */
        
        if (params.fasta_out_fn.size() > 0) {
            params.fasta_out_file = fopen(params.fasta_out_fn.c_str(), "w+");
            
            ococo::info("Opening consensus file ('%s').\n", params.fasta_out_fn.c_str());
            
            if(DEBUGGING_MODE){
                BOOST_LOG_TRIVIAL(info) << "Open FASTA for consensus: '" << params.fasta_out_fn
                << "'.";
            }
            
            params.fasta_out_file = fopen(params.fasta_out_fn.c_str(), "w+");
            
            if (params.fasta_out_file == nullptr) {
                ococo::fatal_error(
                                   "Problem with opening FASTA for consensus: '%s'.\n",
                                   params.fasta_out_fn.c_str());
                correctly_initialized=false;
                return;
            }
        } else {
            if(DEBUGGING_MODE)
                BOOST_LOG_TRIVIAL(info) << "No FASTA file for consensus required.";
        }
    }
    
    /*
     //////////////////////////////////////////////////////
     //////////////////////////////////////////////////////
     //////////////////////////////////////////////////////
     */
    
    template <typename T, int counter_size, int refbase_size>
    bool caller_t<T, counter_size, refbase_size>::check_read(int32_t seqid, int32_t flags, int32_t mapq) {
        /* TODO: return back
         if(DEBUGGING_MODE){
         
         BOOST_LOG_TRIVIAL(debug)
         << "Reading alignment: rname='" << rname << ", chrom=" << seqid
         << ", pos=" << mappping_pos << ", mapq=" << mapq
         << ", flags=" << flags;
         }*/
        
        if ((flags & BAM_FUNMAP) != 0) {
            if(DEBUGGING_MODE){
                BOOST_LOG_TRIVIAL(debug) << "Discarded: read is not aligned.";
            }
            return false;
        }
        
        if (!stats->seq_active[seqid]) {
            if(DEBUGGING_MODE){
                BOOST_LOG_TRIVIAL(debug)
                << "Discarded: consensus calling is off for this chromosome.";
            }
            return false;
        }
        
        if (mapq < stats->params.min_mapq) {
            if(DEBUGGING_MODE){
                BOOST_LOG_TRIVIAL(debug)
                << "Discarded: mapping quality is too low.";
            }
            return false;
        }
        
        return true;
    }
    
    
    
    template <typename T, int counter_size, int refbase_size>
    void caller_t<T, counter_size, refbase_size>::run() {
        
        /*
         * Process alignments.
         */
        ococo::info("Starting the main loop.\n");
        
        if(DEBUGGING_MODE){
            BOOST_LOG_TRIVIAL(info) << "Starting the main loop.";
        }
        
        int32_t r;
        b = bam_init1();
        while ((r = sam_read1(params.sam_file, header, b)) >= 0) {
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
            
            bool read_ok = check_read(seqid, flags, mapq);
            if(!read_ok){
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
                                if(DEBUGGING_MODE){
                                    
                                    BOOST_LOG_TRIVIAL(trace)
                                    << "Omitting base (too low base quality): chrom="
                                    << seqid << ", pos=" << ref_pos
                                    << ", nucl=" << nt256 << ", quality=" << bq << ".";
                                }
                                continue;
                            }
                            
                            if (nt4 == 0x4) {
                                if(DEBUGGING_MODE){
                                    BOOST_LOG_TRIVIAL(trace)
                                    << "Omitting base (ambiguous nucleotide): chrom="
                                    << seqid << ", pos=" << ref_pos
                                    << ", nucl=" << nt256 << ", quality=" << bq << ".";
                                }
                                continue;
                            }
                            
                            if(DEBUGGING_MODE){
                                
                                BOOST_LOG_TRIVIAL(trace)
                                << "Incrementing counter: chrom=" << seqid
                                << ", pos=" << ref_pos << ", nucl=" << nt256
                                << ", quality=" << bq << ". Old state: "
                                << stats->debug_str_counters(seqid, ref_pos) << ",";
                            }
                            
                            stats->seq_stats[seqid][ref_pos] =
                            stats->increment(stats->seq_stats[seqid][ref_pos], nt4);
                            
                            if(DEBUGGING_MODE){
                                BOOST_LOG_TRIVIAL(trace)
                                << "           ...new state: "
                                << stats->debug_str_counters(seqid, ref_pos) << ".";
                            }
                            
                            if (stats->params.mode == ococo::mode_t::REALTIME) {
                                stats->call_consensus_position(params.vcf_file, params.pileup_file,
                                                               seqid, ref_pos);
                                if(DEBUGGING_MODE){
                                    
                                    BOOST_LOG_TRIVIAL(trace)
                                    << "Consensus called. New state: "
                                    << stats->debug_str_counters(seqid, ref_pos) << ".";
                                }
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
            
            if(DEBUGGING_MODE) {
                BOOST_LOG_TRIVIAL(debug) << "Alignment of '" << rname
                << "' incorporated into statistics.";
            }
        }
        
        /*
         * Calling final consensus and export stats.
         */
        
        if (stats->params.mode == ococo::mode_t::BATCH) {
            if(DEBUGGING_MODE){
                BOOST_LOG_TRIVIAL(info) << "Calling consensus for the entire reference "
                "sequence (batch mode).";
            }
            
            stats->call_consensus(params.vcf_file, params.pileup_file);
            
            if (params.fasta_out_fn.size() > 0) {
                if(DEBUGGING_MODE){
                    
                    BOOST_LOG_TRIVIAL(info) << "Saving FASTA: '" << params.fasta_out_fn
                    << "'.";
                }
                
                int error_code = stats->save_fasta(params.fasta_out_fn);
                if (error_code != 0) {
                    ococo::error("FASTA '%s' could not be saved.\n",
                                 params.fasta_out_fn.c_str());
                    return_code = -1;
                }
            } else {
                if(DEBUGGING_MODE){
                    BOOST_LOG_TRIVIAL(info) << "FASTA not saved.";
                }
            }
        }
        
        if (params.stats_out_fn.size() > 0) {
            if(DEBUGGING_MODE){
                BOOST_LOG_TRIVIAL(info) << "Saving statistics: '" << params.stats_out_fn
                << "'.";
            }
            
            ococo::info("Saving statistics ('%s').\n", params.stats_out_fn.c_str());
            
            int error_code = stats->export_stats(params.stats_out_fn);
            if (error_code != 0) {
                ococo::error("Statistics could not be saved ('%s').\n",
                             params.stats_out_fn.c_str());
                return_code = -1;
            }
        } else {
            if(DEBUGGING_MODE){
                BOOST_LOG_TRIVIAL(info) << "Statistics not saved.";
            }
        }
    }
    
    
    
    
    /*
     //////////////////////////////////////////////////////
     //////////////////////////////////////////////////////
     //////////////////////////////////////////////////////
     */
    
    
    template <typename T, int counter_size, int refbase_size>
    caller_t<T, counter_size, refbase_size>::~caller_t(){
        
        if(DEBUGGING_MODE){
            BOOST_LOG_TRIVIAL(info) << "Freeing memory.";
        }
        
        hts_itr_destroy(iter);
        bam_destroy1(b);
        bam_hdr_destroy(header);
        
        if (stats != nullptr) {
            delete stats;
        }
        
        if (return_code == 0) {
            ococo::info("Ococo successfully finished. Bye.\n");
        }
        
        if(DEBUGGING_MODE){
            BOOST_LOG_TRIVIAL(info) << "Ococo finished.";
        }
        
    }
    
}
