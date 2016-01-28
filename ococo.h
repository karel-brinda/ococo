#ifndef _OCOCO_H_
#define _OCOCO_H_

#pragma once


#include <iostream>
#include <sstream>
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <string>
#include <exception>
#include <unistd.h>
#include <math.h>
#include <assert.h>
#include <zlib.h>

#include "htslib/sam.h"
#include "htslib/faidx.h"
#include "htslib/kstring.h"
#include "htslib/khash.h"
#include "htslib/kseq.h"

using namespace std;

namespace ococo {
    
    KSEQ_INIT(gzFile, gzread);
    
    void error_exit(const char * format, ...){
        va_list args;
        va_start (args, format);
        vfprintf (stderr, format, args);
        va_end (args);
        exit(-1);
    }
    
    
    bool file_exists(const string &fn){
        return access( fn.c_str(), F_OK ) != -1;
    }
    
    
    /*********************
     *** Configuration ***
     *********************/
    
    const int  fasta_line_l = 50;
    const int stats_delim_l = 10;

    enum mode_t {BATCH, REALTIME};
    enum alg_t {STOCHASTIC, MAJORITY};
    
    struct consensus_params_t {
        mode_t mode;
        alg_t alg;

        int min_mapq;
        int min_baseq;

        int min_vote;

        bool print_vcf;

        consensus_params_t():
            mode(BATCH),
            alg(MAJORITY),
            min_mapq(1),
            min_baseq(0),
            print_vcf(true)
        {}
    };
    
    
    /****************
     *** Counters ***
     ****************/
    
    /*
     data type for a single position in the genome (=4 counters)
     - counter structure (merged cells): 0...0|cell_T|cell_G|cell_C|cell_A
     */
    typedef uint16_t counter_t ;
    
    const int counter_size=sizeof(counter_t);
    const int cell_bits=4;
    
    static_assert(cell_bits*4 <= counter_size*8,"Counter type is too small");
    
    const counter_t cell_maxval=(1<<cell_bits)-1;
    const counter_t cell_maxval_shifted=cell_maxval>>1;
    const counter_t counter_norm_mask=cell_maxval_shifted \
    | (cell_maxval_shifted<<cell_bits)        \
    | (cell_maxval_shifted<<2*cell_bits)      \
    | (cell_maxval_shifted<<3*cell_bits);
    
    
    /**************************
     *** Translation tables ***
     **************************/
    
    typedef uint8_t nt4_t;
    typedef uint8_t nt16_t;
    
    const uint8_t nt256_nt4[] = {
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
    
    const uint8_t nt16_nt4[] = {
        4, 0, 1, 4,
        2, 4, 4, 4,
        3, 4, 4, 4,
        4, 4, 4, 4
    };
    
    const uint8_t nt256_nt16[] = {
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
    
    const char nt16_nt256[] = "=ACMGRSVTWYHKDBN";
    
    const char nt4_nt256[] = "ACGTN";
    
    
    /***************************
     *** Auxiliary functions ***
     ***************************/
    
    
    // Generate randomly nucleotide with respect to given frequencies.
    inline char rand_nucl(int a, int c, int g, int t){
        const char nucls[] = {'A','C','G','T'};
        
        const int sum=a+c+g+t;
        const int max_val=max({a,c,g,t});
        const int vec[]={a,c,g,t};
        const int prefsum[]={a,a+c,a+c+g,a+c+g+t};
        
        //int rn=randint(0,sum);
        const int rn=rand() % sum;
        char nucl;
        int count=0;
        
        //printf("rn %d\n",rn);
        for(int i=0;i<4;i++){
            if (prefsum[i] > rn) {
                nucl=nucls[i];
                count=vec[i];
                //printf("selected i %d\n",i);
                break;
            }
        }
        
        if(count!=max_val){
            nucl=tolower(nucl);
        }
        
        //printf("%c ",nucl);
        
        return nucl;
    }
    
    // Print error message and exit with code -1.
    void error_exit(const char * format, ...);
    
    // Test if a file exists.
    bool file_exists(const string &fn);
    
    
    /********************************
     *** Structure for statistics ***
     ********************************/
    
    struct stats_t {
        int16_t   n_seqs;
        bool      *seq_used;
        int32_t   *seq_len;
        int32_t   *seq_comprseqlen;
        char      **seq_name;
        char      **seq_comment;
        uint8_t   **seq_comprseq;
        counter_t **counters;
        
        consensus_params_t params;
        
        //stats_t();
        stats_t(consensus_params_t params,bam_hdr_t &h);
        ~stats_t();
        
        
        /***********************
         *** Loading headers ***
         ***********************/
        
        // Load headers from a FAI index.
        int load_headers_fai(const string &fai_fn);
        
        // Load header from a FASTA file and initialize statistics (to level).
        int load_fasta(const string &fasta_fn, uint16_t initial_weight);
        
        // Loader header from a BAM.
        int load_headers_bam_hdr(const bam_hdr_t &h);
        
        
        /************************
         *** Checking headers ***
         ************************/
        
        // Check if everything was initialized.
        bool check_state() const;
        
        // Check if a FAI header corresponds to the stats.
        bool check_headers_fai(const string &fai_fn) const;
        
        // Check if a BAM header corresponds to the stats.
        bool check_headers_bam_hdr(const bam_hdr_t &h) const;
        
        
        /***********
         *** I/O ***
         ***********/
        
        int import_stats(const string &stats_fn);
        int export_stats(const string &stats_fn) const;
        
        // Call consensus probabilistically.
        int call_consensus(bool print_vcf);
        int call_consensus_position(int ref, int pos, bool print_vcf);
        
        int save_fasta(const string &fasta_fn) const;
        
        void print_vcf_header() const;
        void print_vcf_substitution(int ref, int pos, unsigned char old_base, unsigned char new_base) const;
        
        
        /*****************
         *** Auxiliary ***
         *****************/
        
        inline char get_nucl(int ref, int pos) const;
        inline void set_nucl(int ref, int pos, unsigned char nucl);
        //inline counter_t get_counter_value(int ref, int pos);
        void get_counters_values(int ref, int pos, counter_t &a, counter_t &c, counter_t &g, counter_t &t) const;
        
        /****************
         *** Debuging ***
         ****************/
        
        void debug_print_counters() const;
        string debug_str_counters(int ref, int pos) const;
    };
    
    /**********************************
     *** Manipulating with counters ***
     **********************************/
    
    /*
     const int nt16_A = 0x1;
     const int nt16_C = 0x2;
     const int nt16_G = 0x4;
     const int nt16_T = 0x8;
     const int nt16_N = 0xf;
     */
    
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
    }
    
    inline void STATS_INCREMENT(stats_t &stats,int seqid,int pos,nt16_t nt16) {
        assert(stats.seq_len[seqid]>pos);
        //fprintf(stderr,"Going to update stats: chrom %d, pos %d, nucl %d\n", seqid, pos, nt16);
        //fprintf(stderr,"Counter at pos %d: %04x\n", pos, stats.counters[seqid][pos]);
        stats.counters[seqid][pos] = _COUNTER_CELL_INC (stats.counters[seqid][pos],nt16);
        //fprintf(stderr,"Counter at pos %d: %04x\n", pos, stats.counters[seqid][pos]);
    }
    
    
    /**********************
     *** Implementation ***
     **********************/
    
    stats_t::stats_t(consensus_params_t params,bam_hdr_t &h):
    n_seqs(h.n_targets),
    seq_used(new bool[n_seqs]()),
    seq_len(new int32_t[n_seqs]()),
    seq_comprseqlen(new int32_t[n_seqs]()),
    seq_name(new char*[n_seqs]()),
    seq_comment(new char*[n_seqs]()),
    seq_comprseq(new uint8_t*[n_seqs]()),
    counters(new counter_t*[n_seqs]()),
    params(params)
    {
        for (int s=0;s<n_seqs;s++){
            seq_len[s]=h.target_len[s];
            seq_comprseqlen[s]=(int32_t)ceil(seq_len[s]/4.0);
            printf("haha %d\n",seq_comprseqlen[s]);
            printf("len %d\n",h.target_len[s]);
            seq_used[s]=true;
            //fprintf(stderr,"allocating %d chars\n",seq_len[i]);
            const int seq_len_name=strlen(h.target_name[s]);
            seq_name[s]=new char[seq_len_name+1];
            //printf("name: %s\n",h.target_name[i]);
            memcpy(seq_name[s], h.target_name[s],seq_len_name+1);
            
            seq_comprseq[s]=new uint8_t[seq_comprseqlen[s]]();
            counters[s]=new counter_t[seq_len[s]]();
            //printf("seq %s, len %d\n",stats->seqstats[i].name,stats->seqstats[i].length);
            //fprintf(stderr,"ok\n");
        }
    }
    
    
    stats_t::~stats_t(){
        for (int i=0;i<n_seqs;i++){
            delete[] seq_name[i];
            delete[] seq_comment[i];
            delete[] seq_comprseq[i];
            delete[] counters[i];
        }
        delete[] seq_used;
        delete[] seq_len;
        delete[] seq_comprseqlen;
        delete[] seq_name;
        delete[] seq_comment;
        delete[] seq_comprseq;
        delete[] counters;
    }
    
    
    int stats_t::load_fasta(const string &fasta_fn, uint16_t initial_weight) {
        gzFile fp;
        kseq_t *seq;
        int l;
        fp = gzopen(fasta_fn.c_str(), "r");
        seq = kseq_init(fp);
        for(int s=0;(l = kseq_read(seq)) >= 0;s++) {
            assert(strcmp(seq->name.s,seq_name[s])==0);
            assert((int)seq->seq.l == seq_len[s]);
            if (seq->comment.l && seq_comment[s]==NULL){
                seq_comment[s]=new char[seq->comment.l+1];
                memcpy(seq_comment[s], seq->comment.s,seq->comment.l+1);
            }
            //fprintf(stderr,"seq: %s\n", seq->seq.s);
            
            for(uint32_t i=0;i<seq->seq.l;i++){
                assert(counters[s][i]==0);
                
                char &nucl = seq->seq.s[i];
                set_nucl(s,i,nucl);
                const nt16_t nt16=nt256_nt16[(int)nucl];
                
                assert(counters[s][i]==0);
                counters[s][i]=_COUNTER_CELL_SET(0,nt16,initial_weight);
            }
            //if (seq->qual.l) printf("qual: %s\n", seq->qual.s);
        }
        kseq_destroy(seq); // STEP 5: destroy seq
        gzclose(fp); // STEP 6: close the file handler
        return 0;
    }
    
    
    int stats_t::save_fasta(const string &fasta_fn) const {
        assert(check_state());
        
        FILE *fp;
        fp = fopen(fasta_fn.c_str(), "w+");
        
        char fasta_buffer[fasta_line_l];
        for(int s=0;s<n_seqs;s++){
            //printf("%s\n",seq_name[s]);
            if(seq_comment[s]){
                fprintf(fp,">%s %s\n",seq_name[s],seq_comment[s]);
            }
            else{
                fprintf(fp,">%s\n",seq_name[s]);
            }
            
            for (int i=0,j=0;i<seq_len[s];i++,j++){
                fasta_buffer[j]=get_nucl(s,i);
                
                if(j==fasta_line_l-1 || i==seq_len[s]-1){
                    fwrite(fasta_buffer,1,j+1,fp);
                    fwrite("\n",1,1,fp);
                    j=-1;
                }
            }
        }
        
        fclose(fp);
        return 0;
    }
    
    
    bool stats_t::check_state() const {
        if(n_seqs==0) return false;
        if(seq_used==NULL || seq_len==NULL || seq_name==NULL || seq_comment==NULL || counters==NULL)
            return false;
        
        for(int i=0;i<n_seqs;i++){
            if(seq_name[i]==NULL || counters[i]==NULL){
                return false;
            }
        }
        
        return true;
    }
    
    
    bool stats_t::check_headers_fai(const string &fai_fn) const {
        if (!check_state()) return false;
        
        //todo
        return true;
    }
    
    
    bool stats_t::check_headers_bam_hdr(const bam_hdr_t &h) const {
        if (!check_state()) return false;
        
        //if (strcmp(seq->name.s,seq_name[s])!=0) return false;
        //if ((int)seq->seq.l != seq_len[s]) return false;
        
        for(int s=0;s<n_seqs;s++){
            if(seq_len[s]!=(int)h.target_len[s]) return false;
            if (strcmp(h.target_name[s],seq_name[s])!=0) return false;
        }
        
        
        /*for(int i=0;i<n_seqs;i++){
         if (strcmp(seq->name.s,seq_name[s])!=0) return false;
         if ((int)seq->seq.l != seq_len[s]) return false;
         }*/
        
        return true;
    }
    
    
    int stats_t::import_stats(const string &stats_fn){
        assert(check_state());
        
        FILE *fo=fopen(stats_fn.c_str(),"r");
        
        char delim[stats_delim_l]={};
        int16_t seq_name_l;
        int16_t seq_comment_l;
        
        int16_t n_seqs_loaded;
        bool    seq_used_loaded;
        int32_t seq_len_loaded;
        int16_t seq_name_l_loaded;
        int16_t seq_comment_l_loaded;
        
        /* number of seqs */
        fread(&n_seqs_loaded,sizeof(int16_t),1,fo);
        assert(n_seqs_loaded = n_seqs);
        
        for(int i=0;i<n_seqs;i++){
            
            /* delimiter */
            fread(delim,sizeof(char),stats_delim_l,fo);
            assert(delim[0]=='\255');
            assert(delim[stats_delim_l-1]=='\255');
            
            /* lengts */
            seq_name_l=strlen(seq_name[i]);
            seq_comment_l=seq_comment[i]==NULL ? 0 : strlen(seq_comment[i]);
            fread(&seq_used_loaded,sizeof(bool),1,fo);
            assert(seq_used_loaded==seq_used[i]);
            fread(&seq_len_loaded,sizeof(int32_t),1,fo);
            assert(seq_len_loaded==seq_len[i]);
            fread(&seq_name_l_loaded,sizeof(int16_t),1,fo);
            assert(seq_name_l_loaded==seq_name_l);
            fread(&seq_comment_l_loaded,sizeof(int16_t),1,fo);
            assert(seq_comment_l_loaded==seq_comment_l);
            
            /* strings */
            // values of strings are not checked
            char seq_name_loaded[seq_name_l+1];
            fread(seq_name_loaded,sizeof(char),seq_name_l+1,fo);
            char seq_comment_loaded[seq_comment_l+1];
            fread(seq_comment_loaded,sizeof(char),seq_comment_l+1,fo);
            
            /* counters */
            fread(counters[i],sizeof(counter_t),seq_len[i],fo);
        }
        
        fclose(fo);
        return 0;
    }
    
    int stats_t::export_stats(const string &stats_fn) const {
        assert(check_state());
        
        FILE *fo=fopen(stats_fn.c_str(),"w+");
        
        char delim[stats_delim_l]={};
        int16_t seq_name_l;
        int16_t seq_comment_l;
        
        delim[0]='\255';
        delim[stats_delim_l-1]='\255';
        
        /* number of seqs */
        fwrite(&n_seqs,sizeof(int16_t),1,fo);
        
        for(int i=0;i<n_seqs;i++){
            
            /* delimiter */
            fwrite(delim,sizeof(char),stats_delim_l,fo);
            
            /* lengts */
            seq_name_l=strlen(seq_name[i]);
            seq_comment_l=seq_comment[i]==NULL ? 0 : strlen(seq_comment[i]);
            
            fwrite(&seq_used[i],sizeof(bool),1,fo);
            fwrite(&seq_len[i],sizeof(int32_t),1,fo);
            fwrite(&seq_name_l,sizeof(int16_t),1,fo);
            fwrite(&seq_comment_l,sizeof(int16_t),1,fo);
            
            /* strings */
            assert(seq_name[i] != NULL);
            fwrite(seq_name[i],sizeof(char),seq_name_l+1,fo);
            if(seq_comment[i] != NULL){
                fwrite(seq_comment[i],sizeof(char),seq_comment_l+1,fo);
            }
            
            /* counters */
            assert(counters[i] != NULL);
            fwrite(counters[i],sizeof(counter_t),seq_len[i],fo);
        }
        
        fclose(fo);
        return 0;
    }
    
    
    int stats_t::call_consensus(bool print_vcf) {
        assert(check_state());
        
        for(int s=0;s<n_seqs;s++){
            for (int i=0;i<seq_len[s];i++){
                call_consensus_position(s, i, print_vcf);
            }
        }
        
        return 0;
    }
    
    int stats_t::call_consensus_position(int ref, int pos, bool print_vcf) {
        //BOOST_LOG_TRIVIAL(debug) << "Calling consensus for position " << pos;
        
        counter_t A,C,G,T;
        get_counters_values(ref, pos, A, C, G, T);
        
        const char new_base=(A+C+G+T<params.min_vote) ? 'N' : rand_nucl(A,C,G,T);
        const char old_base=get_nucl(ref,pos);
        
        if(old_base!=new_base){
            if(print_vcf){
                print_vcf_substitution(ref,pos,old_base,new_base);
            }
            set_nucl(ref,pos,new_base);
        }
        
        return 0;
    }
    
    void stats_t::get_counters_values(int ref, int pos, counter_t &a, counter_t &c, counter_t &g, counter_t &t) const {
        a=_COUNTER_CELL_VAL(counters[ref][pos],seq_nt16_table[(int)'A']);
        c=_COUNTER_CELL_VAL(counters[ref][pos],seq_nt16_table[(int)'C']);
        g=_COUNTER_CELL_VAL(counters[ref][pos],seq_nt16_table[(int)'G']);
        t=_COUNTER_CELL_VAL(counters[ref][pos],seq_nt16_table[(int)'T']);
    }
    
    void stats_t::print_vcf_header() const {
        assert(check_state());
        
        //todo: date
        printf(
               "##fileformat=VCFv4.3\n"
               "##fileDate=20150000\n"
               //"##reference=%s\n"
               "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
               );
        
    }
    
    void stats_t::print_vcf_substitution(int ref, int pos, unsigned char old_base, unsigned char new_base) const {
        assert(check_state());
        
        printf(
               "%s\t%d\t.\t%c\t%c\t100\tPASS\t.\n",
               seq_name[ref],
               pos+1,
               old_base,
               new_base
               );
        
    }
    
    string stats_t::debug_str_counters(int ref, int pos) const {
        counter_t a,c,g,t;
        get_counters_values(ref, pos, a, c, g, t);
        std::stringstream ss;
        ss << "(" << a << "," << c << "," << g << "," << t << ")";
        return ss.str();
    }
    
    void stats_t::debug_print_counters() const {
        for(int s=0;s<n_seqs;s++){
            fprintf(stderr,"%s\n",seq_name[s]);
            for (int i=0;i<seq_len[s];i++){
                fprintf(stderr,"%8d %04x \n",i,counters[s][i]);
            }
        }
    }
    
    inline char stats_t::get_nucl(int ref, int pos) const {
        if (counters[ref][pos]==0){
            // empty statistics => unknown nucleotide
            return 'N';
        }
        else {
            // non-empty statistics => look into comprseq
            
            const uint32_t comprseq_coor_1 = pos >> 2;
            const uint32_t comprseq_coor_2 = 3- (pos & 0x3);
            
            const nt4_t &nt4 = (seq_comprseq[ref][comprseq_coor_1] >> comprseq_coor_2) & 0x3;
            return nt4_nt256[nt4];
        }
    }
    
    inline void stats_t::set_nucl(int ref, int pos, unsigned char nucl){
        const nt4_t nt4=nt256_nt4[(int)nucl] & 0x3;
        
        const uint32_t comprseq_coor_1 = pos >> 2;
        const uint32_t comprseq_coor_2 = 3 - (pos & 0x3);
        
        uint8_t &cell=seq_comprseq[ref][comprseq_coor_1];
        cell ^=	cell & (0x3 << comprseq_coor_2);
        cell |= nt4 << comprseq_coor_2;
    }
    
};


#endif
