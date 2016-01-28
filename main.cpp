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
    
    ococo::consensus_params_t params;
    
    string fasta0_fn;
    string fasta1_fn;
    string stats_fn;
    string sam_fn;
    
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
        
        string strategy, mode;
        
        vol.add_options()
        //("help", "Print help message")
        ("input,i", po::value<string>(&sam_fn)->required(), "Input SAM/BAM file (- for standard input).")
        ("fasta-ref,f", po::value<string>(&fasta0_fn), "Initial FASTA reference (if not provided, sequence of N's is considered as the reference).")
        ("fasta-cons,c", po::value<string>(&fasta1_fn), "Consensus (up-to-date FASTA reference).")
        ("stats,s", po::value<string>(&stats_fn), "File with up-to-date statistics.")
        ("strategy,S", po::value<string>(&strategy), "Strategy for updates: majority / randomized. [majority]")
        ("mode,M", po::value<string>(&strategy), "Mode: real-time / batch. [batch]")
        //("counter-size,s", po::value<int>(&counterSize), "Size of counter per nucleotide in bits [3]")
        ("min-MQ,q", po::value<int>(&params.min_mapq), "Minimal mapping quality to increment a counter. [1]")
        ("min-BQ,Q", po::value<int>(&params.min_baseq), "Minimal base quality to increment a counter. [0]")
        //("min-coverage,c", po::value<int>(&params.min_coverage), "Minimal coverage to update the reference. [3]")
        //("accept-level,l", po::value<float>(&acceptanceLevel), "Acceptance level [0.60]")
        ("no-vcf,v", "Do not print VCF.")
        ;
        
        po::variables_map vm;
        try
        {
            po::store(po::command_line_parser(argc, argv).options(vol).positional(pos).run(),vm); // can throw
            po::notify(vm); // throws on error, so do after help in case there are any problems
            
            params.print_vcf=!vm.count("no-vcf");
            
            if (vm.count("strategy")) {
                if (strategy=="majority"){
                    params.strategy=ococo::strategy_t::MAJORITY;
                }
                else if (strategy=="stochastic"){
                    params.strategy=ococo::strategy_t::STOCHASTIC;
                }
                else {
                    fprintf(stderr,"Unknown strategy '%s'. Possible strategies are 'majority' and 'stochastic'\n",strategy.c_str());
                    return EXIT_FAILURE;
                }
            }
            
            if (vm.count("mode")) {
                if (mode=="batch"){
                    params.mode=ococo::mode_t::BATCH;
                }
                else if (strategy=="real-time"){
                    params.mode=ococo::mode_t::REALTIME;
                }
                else {
                    fprintf(stderr,"Unknown mode '%s'. Possible strategies are 'batch' and 'real-time'\n",strategy.c_str());
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
            cout << vol << "\n";
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
    hts_itr_t *iter=NULL;
    
    samFile *in = NULL;
    bam1_t *b= NULL;
    bam_hdr_t *header = NULL;
    
    BOOST_LOG_TRIVIAL(info) << "SAM/BAM reader initialization: reading '" << sam_fn.c_str() << "'.";
    in = sam_open(sam_fn.c_str(), "r");
    if(in==NULL) {
        ococo::error_exit("Problem with opening input ('%s').\n", sam_fn.c_str());
    }
    if ((header = sam_hdr_read(in)) == 0){
        ococo::error_exit("SAM headers are missing or corrupted.\n");
    }
    
    ococo::stats_t stats(params,*header);
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
            stats.load_fasta(fasta0_fn,2);
        }
    }
    
    if(params.print_vcf){
        stats.print_vcf_header(true);
    }
    
    /*
     * Process alignments.
     */
    BOOST_LOG_TRIVIAL(info) << "Starting the main loop.";
    
    int r;
    b = bam_init1();
    while ((r = sam_read1(in, header, b)) >= 0) { // read one alignment from `in'
        const char* rname=bam_get_qname(b);
        const uint8_t *seq=bam_get_seq(b);
        const uint8_t *qual=bam_get_qual(b);
        const uint32_t *cigar=bam_get_cigar(b);
        const int n_cigar=b->core.n_cigar;
        //+b->core.l_qname
        const int chrom=b->core.tid;
        const int read_pos=b->core.pos;
        const int mapq=b->core.qual;
        const int flags=b->core.flag;
        
        //fprintf(stderr,"pos %d, chrom %d, map q %d, flag %d, name %s \n",pos,chrom,mapq, flags, rname);
        
        BOOST_LOG_TRIVIAL(debug) << "Reading alignment: rname='" << rname << ", chrom=" << chrom << ", pos=" << read_pos <<", mapq="<< mapq << ", flags=" << flags;
        
        
        if ((flags & BAM_FUNMAP)!=0){
            BOOST_LOG_TRIVIAL(debug) << "Discarded: read is not aligned.";
            continue;
        }
        
        if (!stats.seq_used[chrom]){
            BOOST_LOG_TRIVIAL(debug) << "Discarded: consensus calling is off for this chromosome.";
            continue;
        }
        
        if (mapq<params.min_mapq){
            BOOST_LOG_TRIVIAL(debug) << "Discarded: mapping quality is too low.";
            continue;
        }
        
        int32_t ref_pos;
        for (int k=0,i=0; k < n_cigar; k++)
        {
            const int op = bam_cigar_op(cigar[k]);
            const int ol = bam_cigar_oplen(cigar[k]);
            int ni=0;
            
            switch (op) {
                case BAM_CMATCH:
                case BAM_CDIFF:
                case BAM_CEQUAL:
                    for (ni=i+ol;i<ni;i++){
                        ref_pos=read_pos+i;
                        const uint8_t &nt16=bam_seqi(seq, i);
                        const uint8_t &bq=qual[i];
                        
                        if (bq<params.min_baseq){
                            BOOST_LOG_TRIVIAL(trace) << "Omitting base (too low base quality): chrom=" << chrom << ", pos=" << ref_pos << ", nucl=" << ococo::nt16_nt256[nt16] << ", quality=" << (int)bq << ".";
                            
                        }
                        else{
                            BOOST_LOG_TRIVIAL(trace) << "Incrementing counter: chrom=" << chrom << ", pos=" << ref_pos << ", nucl=" << ococo::nt16_nt256[nt16] << ", quality=" << (int)bq << ". New state: refbase='" << stats.get_nucl(chrom, ref_pos) << "', counters: " << stats.debug_str_counters(chrom,ref_pos);
                            STATS_INCREMENT(stats,chrom,ref_pos,nt16);
                            BOOST_LOG_TRIVIAL(trace) << "           new state: counters: " << stats.debug_str_counters(chrom,ref_pos);

                            if(params.mode==ococo::mode_t::REALTIME){
                                stats.call_consensus_position(chrom, ref_pos, params.print_vcf);
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
    if(params.mode==ococo::mode_t::BATCH){
        BOOST_LOG_TRIVIAL(info) << "Calling consensus for the entire reference sequence (batch mode).";
        stats.call_consensus(params.print_vcf);
    }
    
    if (fasta1_fn.size()>0){
        BOOST_LOG_TRIVIAL(info) << "Saving FASTA: '" << fasta1_fn << "'.";
        stats.save_fasta(fasta1_fn);
    }
    else {
        BOOST_LOG_TRIVIAL(info) << "FASTA not saved.";
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
    sam_close(in);
    
    BOOST_LOG_TRIVIAL(info) << "Ococo finished. Bye.";
    
    return 0;
}