#include "ococo.h"

#define BOOST_LOG_DYN_LINK

#include <boost/format.hpp>

#include <boost/program_options.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>

namespace logging = boost::log;

void boost_logging_init()
{
    logging::core::get()->set_filter
    (
     logging::trivial::severity >= logging::trivial::warning
     //logging::trivial::severity >= logging::trivial::trace
     );
}


int main(int argc, const char* argv[])
{
    
    boost_logging_init();
    
    /*
     BOOST_LOG_TRIVIAL(trace) << "A trace severity message";
     BOOST_LOG_TRIVIAL(debug) << "A debug severity message";
     BOOST_LOG_TRIVIAL(info) << "An informational severity message";
     BOOST_LOG_TRIVIAL(warning) << "A warning severity message";
     BOOST_LOG_TRIVIAL(error) << "An error severity message";
     BOOST_LOG_TRIVIAL(fatal) << "A fatal severity message";
     */
    
    BOOST_LOG_TRIVIAL(info) << "Ococo starting.";
    
    /*
     * Default configuration.
     */
    
    ococo::consensus_params_t tmp_params = ococo::consensus_params_t();
    
    std::string fasta0_fn;
    std::string stats_fn;
    std::string sam_fn;

    std::string vcf_fn;
    std::string fasta_cons_fn;
    
    
    /*
     * Parse command-line parameters.
     */
    
    BOOST_LOG_TRIVIAL(info) << "Parsing command-line parameters.";
    try
    {
        namespace po = boost::program_options;
        
        po::positional_options_description pos;
        pos.add("input-file", -1);
        
        po::options_description vol("Ococo: On-line consensus caller.");
        
        std::string strategy;
        std::string mode;
        
        vol.add_options()
        //("help", "Print help message")
        ("input,i", po::value<std::string>(&sam_fn)->required(), "Input SAM/BAM file (- for standard input).")
        ("fasta-ref,f", po::value<std::string>(&fasta0_fn), "Initial FASTA reference (if not provided, sequence of N's is considered as the reference).")
        ("stats,s", po::value<std::string>(&stats_fn), "File with up-to-date statistics.")
        ("strategy,S", po::value<std::string>(&strategy), "Strategy for updates: majority / randomized. [majority]")
        ("mode,m", po::value<std::string>(&mode), "Mode: real-time / batch. [batch]")
        //("counter-size,s", po::value<int>(&counterSize), "Size of counter per nucleotide in bits [3]")
        ("min-MQ,q", po::value<int32_t>(&tmp_params.min_mapq), "Minimal mapping quality to increment a counter. [1]")
        ("min-BQ,Q", po::value<int32_t>(&tmp_params.min_baseq), "Minimal base quality to increment a counter. [0]")
        //("min-coverage,c", po::value<int32_t>(&params.min_coverage), "Minimal coverage to update the reference. [3]")
        //("accept-level,l", po::value<float>(&acceptanceLevel), "Acceptance level [0.60]")
        ("vcf-cons,v", po::value<std::string>(&vcf_fn), "VCF file with updates of consensus.")
        ("fasta-cons,c", po::value<std::string>(&fasta_cons_fn), "FASTA file with consensus, which is continuously updated (WARNING: will be rewritten).")
        ;
        
        po::variables_map vm;
        try
        {
            po::store(po::command_line_parser(argc, argv).options(vol).positional(pos).run(),vm); // can throw
            po::notify(vm); // throws on error, so do after help in case there are any problems
            
            if (vm.count("strategy")) {
                if (strategy.compare("majority")==0){
                    tmp_params.strategy=ococo::strategy_t::MAJORITY;
                }
                else if (strategy.compare("stochastic")==0){
                    tmp_params.strategy=ococo::strategy_t::STOCHASTIC;
                }
                else {
                    fprintf(stderr,"Unknown strategy '%s'. Possible strategies are 'majority' and 'stochastic'\n",strategy.c_str());
                    return EXIT_FAILURE;
                }
            }
            
            if (vm.count("mode")) {
                if (mode.compare("batch")==0){
                    tmp_params.mode=ococo::mode_t::BATCH;
                }
                else if (mode.compare("real-time")==0){
                    tmp_params.mode=ococo::mode_t::REALTIME;
                }
                else {
                    fprintf(stderr,"Unknown mode '%s'. Possible modes are 'batch' and 'real-time'\n",mode.c_str());
                    return EXIT_FAILURE;
                }
            }
            
            
            /*if (vm.count("help")) {
             cout << vol << "\n";
             return EXIT_FAILURE;
             }*/
        }
        catch(po::error& e)
        {
            std::cout << vol << "\n";
            fprintf(stderr,"Error: %s.\n",e.what());
            return EXIT_FAILURE;
        }
        
    }
    catch(std::exception& e)
    {
        fprintf(stderr,"Unhandled Exception: %s.\n",e.what());
        return EXIT_FAILURE;
    }
    
    
    /*
     * Read SAM headers.
     */
    hts_itr_t *iter = nullptr;
    
    samFile *in = nullptr;
    bam1_t *b = nullptr;
    bam_hdr_t *header = nullptr;
    
    BOOST_LOG_TRIVIAL(info) << "SAM/BAM reader initialization: reading '" << sam_fn.c_str() << "'.";
    in = sam_open(sam_fn.c_str(), "r");
    if(in==nullptr) {
        ococo::error_exit("Problem with opening input ('%s').\n", sam_fn.c_str());
    }
    if ((header = sam_hdr_read(in)) == 0){
        ococo::error_exit("SAM headers are missing or corrupted.\n");
    }
    
    ococo::stats_t<uint16_t,3,4> stats(tmp_params,*header);
    assert(stats.check_state());
    
    /*
     * Load FASTA and stats.
     */
    if (stats_fn.size()>0 && ococo::file_exists(stats_fn)){
        BOOST_LOG_TRIVIAL(info) << "Importing statistics: '" << stats_fn << "'.";
        stats.import_stats(stats_fn);
    }
    else{
        BOOST_LOG_TRIVIAL(info) << "No file with statistics provided.";
        if (fasta0_fn.size()>0){
            BOOST_LOG_TRIVIAL(info) << "Loading FASTA: '" << fasta0_fn << "'.";
            stats.load_fasta(fasta0_fn);
        }
    }
    
    /*
     * Open VCF file.
     */
    if(vcf_fn.size()>0){
        BOOST_LOG_TRIVIAL(info) << "Open VCF: '" << vcf_fn << "'.";
        if (vcf_fn==std::string("-")){
            stats.params.vcf_fo=stdout;
        }
        else{
            stats.params.vcf_fo=fopen(vcf_fn.c_str(),"w+");
        }
        
        if(stats.params.vcf_fo==nullptr){
            ococo::error_exit("Problem with opening VCF file '%s'", vcf_fn.c_str());
        }
        
        stats.print_vcf_header();
    }
    else {
        BOOST_LOG_TRIVIAL(info) << "No VCF file required.";
    }

    /*
     * Open consensus FASTA file.
     */
    if(fasta_cons_fn.size()>0){
        BOOST_LOG_TRIVIAL(info) << "Open FASTA for consensus: '" << fasta_cons_fn << "'.";
        stats.params.fasta_cons_fo=fopen(fasta_cons_fn.c_str(),"w+");
        
        if(stats.params.fasta_cons_fo==nullptr){
            ococo::error_exit("Problem with opening FASTA for consensus: '%s'", fasta_cons_fn.c_str());
        }
    }
    else {
        BOOST_LOG_TRIVIAL(info) << "No FASTA file for consensus required.";
    }
    
    
    /*
     * Process alignments.
     */
    BOOST_LOG_TRIVIAL(info) << "Starting the main loop.";
    
    int32_t r;
    b = bam_init1();
    while ((r = sam_read1(in, header, b)) >= 0) { // read one alignment from `in'
        const char* rname=bam_get_qname(b);
        const uint8_t *seq=bam_get_seq(b);
        const uint8_t *qual=bam_get_qual(b);
        const uint32_t *cigar=bam_get_cigar(b);
        const int32_t n_cigar=b->core.n_cigar;
        //+b->core.l_qname
        const int32_t seqid=b->core.tid;
        const int64_t read_pos=b->core.pos;
        const int32_t mapq=b->core.qual;
        const int32_t flags=b->core.flag;
        
        //fprintf(stderr,"pos %d, chrom %d, map q %d, flag %d, name %s \n",pos,chrom,mapq, flags, rname);
        
        BOOST_LOG_TRIVIAL(debug) << "Reading alignment: rname='" << rname << ", chrom=" << seqid << ", pos=" << read_pos <<", mapq="<< mapq << ", flags=" << flags;
        
        
        if ((flags & BAM_FUNMAP)!=0){
            BOOST_LOG_TRIVIAL(debug) << "Discarded: read is not aligned.";
            continue;
        }
        
        if (!stats.seq_active[seqid]){
            BOOST_LOG_TRIVIAL(debug) << "Discarded: consensus calling is off for this chromosome.";
            continue;
        }
        
        if (mapq<stats.params.min_mapq){
            BOOST_LOG_TRIVIAL(debug) << "Discarded: mapping quality is too low.";
            continue;
        }
        
        int32_t ref_pos;
        for (int32_t k=0,i=0; k < n_cigar; k++)
        {
            const int32_t op = bam_cigar_op(cigar[k]);
            const int32_t ol = bam_cigar_oplen(cigar[k]);
            int ni=0;
            
            switch (op) {
                case BAM_CMATCH:
                case BAM_CDIFF:
                case BAM_CEQUAL:
                    for (ni=i+ol;i<ni;i++){
                        ref_pos=read_pos+i;
                        const uint8_t &nt16=bam_seqi(seq, i);
                        const uint8_t &bq=qual[i];
                        
                        if (bq<stats.params.min_baseq){
                            BOOST_LOG_TRIVIAL(trace) << "Omitting base (too low base quality): chrom=" << seqid << ", pos=" << ref_pos << ", nucl=" << ococo::nt16_nt256[nt16] << ", quality=" << (int32_t)bq << ".";
                            
                        }
                        else{
                            BOOST_LOG_TRIVIAL(trace) << "Incrementing counter: chrom=" << seqid << ", pos=" << ref_pos << ", nucl=" << ococo::nt16_nt256[nt16] << ", quality=" << (int32_t)bq << ". New state: refbase='" << stats.get_nucl_nt256(seqid, ref_pos) << "', counters: " << stats.debug_str_counters(seqid,ref_pos);
                            stats.seq_stats[seqid][ref_pos] = stats.increment(stats.seq_stats[seqid][ref_pos],nt16);
                            BOOST_LOG_TRIVIAL(trace) << "           new state: counters: " << stats.debug_str_counters(seqid,ref_pos);

                            if(stats.params.mode==ococo::mode_t::REALTIME){
                                stats.call_consensus_position(seqid, ref_pos);
                            }
                        }
                    }
                    break;
                    
                case BAM_CDEL:
                case BAM_CSOFT_CLIP:
                case BAM_CREF_SKIP:
                    i+=ol;
                    break;
                    
                case BAM_CBACK:
                    BOOST_LOG_TRIVIAL(warning) << "Backward operation in CIGAR strings is not supported.";
                    continue;
                case BAM_CINS:
                    break;
                    
                case BAM_CPAD:
                case BAM_CHARD_CLIP:
                    break;
                    
                    //update stats
                    //skip ref
                    
            }
            //printf("%d%c",ol,BAM_CIGAR_STR[op]);
            
            /*if(params.mode==ococo::mode_t::REALTIME){
                for(int pos=read_pos;pos<ref_pos;pos++){
                    stats.call_consensus_position(chrom, pos, params.print_vcf);
                }
            }*/
        }
        
        BOOST_LOG_TRIVIAL(debug) << "Alignment incorporated into statistics.";
    }
    
    /*
     * Calling final consensus and export stats.
     */
    if(stats.params.mode==ococo::mode_t::BATCH){
        BOOST_LOG_TRIVIAL(info) << "Calling consensus for the entire reference sequence (batch mode).";
        stats.call_consensus();
        if (stats.params.fasta_cons_fo){
            BOOST_LOG_TRIVIAL(info) << "Saving FASTA: '" << fasta_cons_fn << "'.";
            stats.save_fasta();
        }
        else {
            BOOST_LOG_TRIVIAL(info) << "FASTA not saved.";
        }
    }
    
    if (stats_fn.size()>0){
        BOOST_LOG_TRIVIAL(info) << "Saving statistics: '" << stats_fn << "'.";
        stats.export_stats(stats_fn);
    }
    else {
        BOOST_LOG_TRIVIAL(info) << "Statistics not saved.";
    }
    
    
    /*
     * Free memory.
     */
    BOOST_LOG_TRIVIAL(info) << "Freeing memory.";
    hts_itr_destroy(iter);
    bam_destroy1(b);
    bam_hdr_destroy(header);
    
    /*
     * Close files.
     */
    sam_close(in);
    if(stats.params.vcf_fo!=nullptr){
        fclose(stats.params.vcf_fo);
    }
    if(stats.params.fasta_cons_fo!=nullptr){
        fclose(stats.params.fasta_cons_fo);
    }
    
    BOOST_LOG_TRIVIAL(info) << "Ococo finished. Bye.";
    
    return 0;
}