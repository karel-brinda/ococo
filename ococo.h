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



#include <boost/format.hpp>

#include <boost/program_options.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>


#include "htslib/sam.h"
#include "htslib/faidx.h"
#include "htslib/kstring.h"
#include "htslib/khash.h"
#include "htslib/kseq.h"


using namespace std;

const int fasta_line_l=50;

const int min_vote=2;

const char CMD_FLUSH[]="flush";

//extern const unsigned char seq_nt16_table[256];
//extern const char seq_nt16_str[16];
//extern const unsigned char seq_nt16_int[256];
//static char bam_nt16_nt4_table[] = { 4, 0, 1, 4, 2, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4 };

typedef uint16_t counter_t ;
typedef unsigned char nt16_t ;

struct stats_t {
	int n_seqs;
	int *seq_len;
	char **seq_name;
	char **seq_comment;
	counter_t **counters;

	stats_t();
	stats_t(bam_hdr_t &h);
	~stats_t();
	int load_fasta(string fasta_fn, int weight);
	int generate_fasta(string fasta_fn);
	void debug_print_counters();
};


inline char rand_nucl(int a, int c, int g, int t);



/*
	nt16: A=1, C=2, G=4, T=8
*/

/*
#define _NUCL_SHIFT(nt16) { \
		(4 * bam_nt16_nt4_table[nt16]) \
	}

#define _COUNTER_VAL(counters,nt16) { \
		(counters>>_NUCL_SHIFT(nt16)) & 0x0f \
	}

#define _COUNTER_SET(counters,nt16,value) { \
		( \
			value << _NUCL_SHIFT(nt16) \
		) | ( \
			counters \
			^ \
			(counters & (0x0f << _NUCL_SHIFT(nt16))) \
		) \
	}


#define _COUNTER_INC_NODIV(counters,nt16) { \
		 _COUNTER_SET (counters,nt16,_COUNTER_VAL(counters,nt16)+1) \
	}

#define _COUNTERS_NORMALIZE(counters,divide) {\
		(counters >> divide) & 0x7777 \
	}

#define _COUNTER_INC(counters,nt16) { \
		_COUNTER_INC_NODIV ( \
			_COUNTERS_NORMALIZE (counters,_COUNTER_VAL(counters,nt16)==0x0f), \
			nt16 \
		) \
	}

#define STATS_UPDATE(stats,seqid,pos,nts16) { \
		((stats->seqstats)[seqid].counters)[pos] = \
		_COUNTER_INC (((stats->seqstats)[seqid].counters)[pos],nts16), \
		(void)0 \
	}
*/


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
		//fprintf(stderr,"Going to update stats: chrom %d, pos %d, nucl %d\n", seqid, pos, nt16);
		//fprintf(stderr,"Counter at pos %d: %04x\n", pos, stats.counters[seqid][pos]);
		stats.counters[seqid][pos] = _COUNTER_CELL_INC (stats.counters[seqid][pos],nt16);
		//fprintf(stderr,"Counter at pos %d: %04x\n", pos, stats.counters[seqid][pos]);
	}

#endif
