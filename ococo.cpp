#include "ococo.h"
#include <zlib.h>  
#include <boost/utility/binary.hpp>

#include <stdint.h>
#include <stdlib.h>
#include <cstring>
#include <cstdio>
#include <unistd.h>
#include <math.h>
#include <inttypes.h>
#include <stdbool.h>
#include <assert.h>
#include "htslib/sam.h"
#include "htslib/faidx.h"
#include "htslib/kstring.h"
#include "htslib/khash.h"
#include "htslib/kseq.h"

KSEQ_INIT(gzFile, gzread);

const char CMD_FLUSH[]="flush";

static char bam_nt16_nt4_table[] = { 4, 0, 1, 4, 2, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4 };

typedef uint16_t counter_t ;
typedef unsigned char nt16_t ;


//enum ConsAlg { CA_MAJORITY , CA_STOCHASTIC };

/*! @typedef
 @abstract Structure for parameters for consensus calling..
 @field min_mapq    Minimum mapping quality to increase counter.
 @field min_baseq   Minimum base quality to increase counter.
 */
typedef struct {
	int min_mapq;
	int min_baseq;
}  cons_params_t;

/*
typedef struct {
	char *name;
	int length;
	counters_t *counters;
} seqstat_t;

typedef struct {
	int nseqs;
	seqstat_t *seqstats;
} stats_t;
*/

struct stats_t {
	int n_seqs;
	int *seq_len;
	char **seq_name;
	counter_t **counters;

	stats_t():
		n_seqs(0),seq_len(NULL), seq_name(NULL), counters(NULL)
	{		
	}

	stats_t(bam_hdr_t &h):
		n_seqs(h.n_targets), seq_len(new int[n_seqs]()), seq_name(new char*[n_seqs]()), counters(new counter_t*[n_seqs]())
	{
		for (int i=0;i<n_seqs;i++){
			seq_len[i]=h.target_len[i];
			fprintf(stderr,"allocating %d chars\n",seq_len[i]);
			const int seq_len_name=strlen(h.target_name[i]);
			seq_name[i]=new char[seq_len_name+1];
			//printf("name: %s\n",h.target_name[i]);
     		memcpy(seq_name[i], h.target_name[i],seq_len_name+1);

			counters[i]=new counter_t[seq_len[i]]();
			//printf("seq %s, len %d\n",stats->seqstats[i].name,stats->seqstats[i].length);
			fprintf(stderr,"ok\n");
		}
	}

	~stats_t(){
		for (int i=0;i<n_seqs;i++){
			delete[] seq_name[i];
			delete[] counters[i];
		}
		delete[] seq_len;
		delete[] seq_name; 
		delete[] counters; 
	}


	//int generate_fasta(char *fasta_fn, stats_t *stats, cons_params_t cps){
	//	return 0;
	//}

	int load_fasta(string fasta_fn)
	{
		gzFile fp;  
	    kseq_t *seq;  
	    int l;  
	    fp = gzopen(fasta_fn.c_str(), "r");
	    seq = kseq_init(fp);
	    while ((l = kseq_read(seq)) >= 0) {
	        fprintf(stderr,"name: %s\n", seq->name.s);  
	        if (seq->comment.l) printf("comment: %s\n", seq->comment.s);  
	        fprintf(stderr,"seq: %s\n", seq->seq.s);  
	        if (seq->qual.l) printf("qual: %s\n", seq->qual.s);  
	    }  
	    kseq_destroy(seq); // STEP 5: destroy seq  
	    gzclose(fp); // STEP 6: close the file handler
	    return 0;
	}
};


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
		return 4 * bam_nt16_nt4_table[nt16];
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
		return (counter >> divide) & 0x7777;
	}

inline counter_t _COUNTER_CELL_INC(counter_t counter,nt16_t nt16) { \
		return _COUNTER_CELL_INC_NODIV (
			_COUNTER_NORMALIZE (counter,_COUNTER_CELL_VAL(counter,nt16)==0x0f),
			nt16
		);
	}

inline void STATS_UPDATE(stats_t &stats,int seqid,int pos,nt16_t nt16) {
		fprintf(stderr,"Going to update stats: chrom %d, pos %d, nucl %d\n", seqid, pos, nt16);
		fprintf(stderr,"Counter at pos %d: %x\n", pos, stats.counters[seqid][pos]);
		stats.counters[seqid][pos] = _COUNTER_CELL_INC (stats.counters[seqid][pos],nt16);
		fprintf(stderr,"Counter at pos %d: %x\n", pos, stats.counters[seqid][pos]);
	}


/*int init_stats(bam_hdr_t *header, stats_t *stats){
	stats = new stats_t;
	stats->nseqs=header->n_targets;
	stats->seqstats=new seqstat_t[stats->nseqs];

	for (int i=0;i<stats->nseqs;i++){
		stats->seqstats[i].length=header->target_len[i];
		stats->seqstats[i].name=strdup(header->target_name[i]);
		stats->seqstats[i].counters=new counters_t[stats->seqstats[i].length];
		printf("seq %s, len %d\n",stats->seqstats[i].name,stats->seqstats[i].length);
	}
	return 0;
}*/

/*int free_seq_stat(seqstat_t *seqstat) {
	free(seqstat->name);
	free(seqstat->counters);
	free(seqstat);
	return 0;
}

int free_stats(stats_t *stats) {
	for (int i=0;i<stats->nseqs;i++){
		free_seq_stat(&stats->seqstats[i]);
	}
	free(stats);
	return 0;
}
*/

int main(int argc, const char* argv[])
{

	/*
		Default configuration.
	*/

	cons_params_t cps = {
		1, // min_mapq
		0, // min_baseq
	};

	bool debug=false;
	string fasta_fn;
	//string sam_fn("-");
	string sam_fn("BWA-MEM.bam"); //debuging

	/*
		Parse command-line parameters.
	*/

	try
	{
		namespace po = boost::program_options;

		po::positional_options_description pos;
		pos.add("input-file", -1);

		po::options_description vol("OPTIONS description");

		vol.add_options()
				("reference,r", po::value<string>(&fasta_fn), "FASTA reference file (with existing FAI index)")
				//("algorithm,a", po::value<string>(&alg), "Algorithm for updates: majority / randomized [majority]")
				//("counter-size,s", po::value<int>(&counterSize), "Size of counter per nucleotide in bits [3]")
				//("min-coverage,c", po::value<int>(&minCoverage), "Minimal coverage [3]")
				("min-map-qual,m", po::value<int>(&cps.min_mapq), "Minimal mapping quality [1]")
				("min-base-qual,b", po::value<int>(&cps.min_baseq), "Minimal base quality [0]")
				//("accept-level,l", po::value<float>(&acceptanceLevel), "Acceptance level [0.60]")
				("debug,d", "Debugging")
		;

		po::variables_map vm;
		try
		{
			po::store(po::command_line_parser(argc, argv).options(vol).positional(pos).run(),vm); // can throw

			if(vm.count("debug")){
				debug=true;
			}

			po::notify(vm); // throws on error, so do after help in case there are any problems

		}
		catch(po::error& e)
		{
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
		Read SAM headers.
	*/

	hts_itr_t *iter=NULL;

	samFile *in = NULL;
	bam1_t *b= NULL;
	bam_hdr_t *header = NULL;

	in = sam_open(sam_fn.c_str(), "r");
	if(in==NULL) {
		fprintf(stderr,"Problem with opening input ('%s').\n", sam_fn.c_str());
		return -1;
	}
	if ((header = sam_hdr_read(in)) == 0){
		fprintf(stderr,"SAM headers are missing or corrupted.\n");
		return -1;
	}

	stats_t stats(*header);

	/*
		Load FASTA.
	*/
	if (debug){
		if(fasta_fn.size()>0){
			fprintf(stderr, "Loading FASTA: %s.\n",fasta_fn.c_str());
		}
		else{
			fprintf(stderr, "No FASTA provided.\n");
		}
	}
	if (fasta_fn.size()>0){
		stats.load_fasta(fasta_fn);
	}


	/*
		Read alignments (main loop).
	*/
	if (debug){
		fprintf(stderr, "Reading alignments started.\n");
	}
	int r;
	b = bam_init1();
    while ((r = sam_read1(in, header, b)) >= 0) { // read one alignment from `in'
		const char* rname=bam_get_qname(b);
		const uint8_t *seq=bam_get_seq(b);
		const uint8_t *qual=bam_get_qual(b);
		const uint32_t *cigar = bam_get_cigar(b);
		const int n_cigar=b->core.n_cigar;
		//+b->core.l_qname
		const int chrom=b->core.tid;
		const int pos=b->core.pos;
		const int mapq=b->core.qual;
		const int flags=b->core.flag;

		//fprintf(stderr,"pos %d, chrom %d, map q %d, flag %d, name %s \n",pos,chrom,mapq, flags, rname);

		if (strcmp(rname,CMD_FLUSH)==0){
			// todo: flush fasta
			break;
		}

    	if (((b->core.flag & BAM_FUNMAP)==0) && (b->core.qual)>=cps.min_mapq){

	        for (int k=0,i=0; k < n_cigar; k++)
			{
				const int op = bam_cigar_op(cigar[k]);
				const int ol = bam_cigar_oplen(cigar[k]);
				int ni=0;
				uint8_t nt16=0;

				switch (op) {
					case BAM_CMATCH:
					case BAM_CDIFF:
					case BAM_CEQUAL:
						for (ni=i+ol;i<ni;i++){
							nt16=bam_seqi(seq, i);
							// opravdu &stats?
							//fprintf(stderr,"Going to update stats: chrom %d, pos %d, nucl %d\n", chrom, pos+i, nt16);
							STATS_UPDATE(stats,chrom,pos+i,nt16);
						}
						break;

					case BAM_CDEL:
					case BAM_CSOFT_CLIP:
					case BAM_CREF_SKIP:
						i+=ol;
						break;

					case BAM_CBACK:
						fprintf(stderr,"Backward operation in CIGAR string is not supported.\n");
						exit(1);
					case BAM_CINS:
						break;

					case BAM_CPAD:
					case BAM_CHARD_CLIP:
						break;

					//update stats
						//skip ref

				}
				//printf("%d%c",ol,BAM_CIGAR_STR[op]);

			}
	        printf("\n");

    	}
    	else {
    		if (debug){
				fprintf(stderr, "Read '%s' is not used for updating counters\n", rname);
    		}
    	}
    }
	if (debug){
		fprintf(stderr, "Reading alignments finished.\n");
	}


	/*
		Free memory.
	*/

	hts_itr_destroy(iter);
	bam_destroy1(b);
	bam_hdr_destroy(header);
	sam_close(in);

	return 0;
}
