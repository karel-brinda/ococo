#ifndef _OCOCO_H_
#define _OCOCO_H_


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


/****************************
 *** Auxiliary data types ***
 ****************************/

typedef uint16_t counter_t ;
typedef unsigned char nt16_t ;

typedef struct {
	int min_mapq;
	int min_baseq;
} cons_params_t;



/***************************
 *** Auxiliary functions ***
 ***************************/

// Generate randomly nucleotide with respect to given frequencies.
inline char rand_nucl(int a, int c, int g, int t);

// Print error message and exit with code -1.
void error_exit(const char * format, ...);

// Test if a file exists.
bool file_exists(const string &fn);

// Init Boost logging,
void boost_logging_init()
{
	logging::core::get()->set_filter
	(
		logging::trivial::severity >= logging::trivial::warning
	);
}


/********************************
 *** Structure for statistics ***
 ********************************/

struct stats_t {
	int16_t   n_seqs;
	bool      *seq_used;
	int32_t   *seq_len;
	char      **seq_name;
	char      **seq_comment;
	counter_t **counters;

	stats_t();
	stats_t(bam_hdr_t &h);
	~stats_t();


	/***********************
	 *** Loading headers ***
	 ***********************/

	// Load headers from a FAI index.
	int load_headers_fai(const string &fai_fn);

	// Load header from a FASTA file and initialize statistics (to level).
	int load_headers_fa(const string &fasta_fn, int level=0);

	// Loader header from a BAM.
	int load_headers_bam_hdr(const bam_hdr_t &h);


	/************************
	 *** Checking headers ***
	 ************************/

	// Check if everything was initialized.
	bool check_state();

	// Check if a FAI header corresponds to the stats.
	bool check_headers_fai(const string &fai_fn);

	// Check if a BAM header corresponds to the stats.
	bool check_headers_bam_hdr(const bam_hdr_t &h);


	/***********
	 *** I/O ***
	 ***********/

	int import_stats(const string &stats_fn);
	int export_stats(const string &stats_fn);

	// Generate consensus probabilistically.
	int generate_consensus(const string &fasta_fn);


	/****************
	 *** Debuging ***
	 ****************/

	void debug_print_counters();
};



/**********************************
 *** Manipulating with counters ***
 **********************************/


inline int _CELL_SHIFT(nt16_t nt16) {
		return 4 * seq_nt16_int[nt16];
}

inline counter_t _COUNTER_CELL_VAL(counter_t counter,nt16_t nt16) {
		return (counter>>_CELL_SHIFT(nt16)) & 0x0f;
}

inline counter_t _COUNTER_CELL_SET(counter_t counter,nt16_t nt16,int value) {
	return
		(
			value << _CELL_SHIFT(nt16) \
		) | (
			counter
			^
			(counter & (0x0f << _CELL_SHIFT(nt16)))
		);
}

inline counter_t _COUNTER_CELL_INC_NODIV(counter_t counter,nt16_t nt16) {
	return _COUNTER_CELL_SET (counter,nt16,_COUNTER_CELL_VAL(counter,nt16)+1);
}

inline counter_t _COUNTER_NORMALIZE(counter_t counter,bool divide) {
	if(divide) fprintf(stderr,"Shift counter: %04x\n", counter);
		return (counter >> divide) & 0x7777;
}

inline counter_t _COUNTER_CELL_INC(counter_t counter,nt16_t nt16) { \
	return _COUNTER_CELL_INC_NODIV (
		_COUNTER_NORMALIZE (counter,_COUNTER_CELL_VAL(counter,nt16)==0x0f),
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
