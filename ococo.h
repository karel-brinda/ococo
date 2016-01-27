#ifndef _OCOCO_H_
#define _OCOCO_H_

#pragma once 


#include <iostream>
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <string>
#include <exception>
#include <unistd.h>
#include <math.h>
#include <assert.h>
#include <zlib.h>  

#define BOOST_LOG_DYN_LINK

#include <boost/format.hpp>

#include <boost/program_options.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>

#include "htslib/sam.h"
#include "htslib/faidx.h"
#include "htslib/kstring.h"
#include "htslib/khash.h"
#include "htslib/kseq.h"


using namespace std;
namespace logging = boost::log;


/*********************
 *** Configuration ***
 *********************/

const int  fasta_line_l = 50;
const int      min_vote =  2;
const int stats_delim_l = 10;


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

typedef uint8_t nt4_t ;
typedef uint8_t nt16_t ;

typedef struct {
	int min_mapq;
	int min_baseq;
} cons_params_t;

//const uint8_t nt4_A = 0x0;
//const uint8_t nt4_C = 0x1;
//const uint8_t nt4_G = 0x2;
//const uint8_t nt4_T = 0x3;
//const uint8_t nt4_N = 0x4;

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

	//printf(" a %d, c %d, g %d, t %d \n",a,c,g,t);
	if (sum<min_vote){
		return 'N';
	}

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

// Init Boost logging,
void boost_logging_init();


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

	//stats_t();
	stats_t(bam_hdr_t &h);
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

inline void STATS_UPDATE(stats_t &stats,int seqid,int pos,nt16_t nt16) {
	assert(stats.seq_len[seqid]>pos);
	//fprintf(stderr,"Going to update stats: chrom %d, pos %d, nucl %d\n", seqid, pos, nt16);
	//fprintf(stderr,"Counter at pos %d: %04x\n", pos, stats.counters[seqid][pos]);
	stats.counters[seqid][pos] = _COUNTER_CELL_INC (stats.counters[seqid][pos],nt16);
	//fprintf(stderr,"Counter at pos %d: %04x\n", pos, stats.counters[seqid][pos]);
}

#endif
