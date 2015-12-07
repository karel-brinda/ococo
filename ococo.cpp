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

const char CMD_FLUSH[]="flush";

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

typedef struct {
	char *name;
	int length;
	uint16_t *counters;
} seqstat_t;

typedef struct {
	int nseqs;
	seqstat_t *seqstats;
} stats_t;

#define FULL_L 0xf0
#define FULL_R 0x0f

/*
	nt16: A=1, C=2, G=4, T=8
*/
#define COUNTER_NUCL_SHIFT(nt16)	(bam_nt16_nt4_table[nt16]*4)

#define COUNTER_VAL(counter,nt16) ((counter>>COUNTER_NUCL_POS[nt16]) & BOOST_BINARY( 1111 ))

#define COUNTER_INC_NO_SHIFT(counter,nt16) { \
		( \
			( \
				(COUNTER_VAL(counter,nt16)+1)<< COUNTER_NUCL_SHIFT[nt16] \
			) \
			| \
			( \
				counter \
				^ \
				(counter & (BOOST_BINARY(1111) << COUNTER_NUCL_SHIFT[nt16])) \
			) \
	}

#define COUNTER_SHIFT(counter,shift_bool) {\
		(counter >> shift_bool) & BOOST_BINARY( 0111011101110111 ) \
	}

#define COUNTER_INC(counter,nt16) { \
		COUNTER_INC_NO_SHIFT( \
			COUNTER_SHIFT(counter,COUNTER_VAL(counter,nt16)==BOOST_BINARY(1111)) \
		) \
	}


int flush_fasta(char *fasta_fn, stats_t *stats, cons_params_t cps){
	return 0;
}

KSEQ_INIT(gzFile, gzread)

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
	string sam_fn("-");

	/*
		Parse command-line parameters
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
		Read SAM headers
	*/

	hts_itr_t *iter=NULL;

	samFile *in = NULL;
	bam1_t *b= NULL;
	bam_hdr_t *header = NULL;

	in = sam_open(sam_fn.c_str(), "r");
	if(in==NULL) {
		fprintf(stderr,"Problem with opening input ('%s').", sam_fn.c_str());
		return -1;
	}
	if ((header = sam_hdr_read(in)) == 0){
		fprintf(stderr,"SAM headers are missing or corrupted.");
		return -1;
	}

	b = bam_init1();

	/*
		Load FASTA
	*/
	gzFile fp;  
    kseq_t *seq;  
    int l;  
    fp = gzopen(fasta_fn.c_str(), "r");
    seq = kseq_init(fp);
    while ((l = kseq_read(seq)) >= 0) {
        printf("name: %s\n", seq->name.s);  
        if (seq->comment.l) printf("comment: %s\n", seq->comment.s);  
        printf("seq: %s\n", seq->seq.s);  
        if (seq->qual.l) printf("qual: %s\n", seq->qual.s);  
    }  
    kseq_destroy(seq); // STEP 5: destroy seq  
    gzclose(fp); // STEP 6: close the file handler  

	/*
		Read alignments (main loop)
	*/
	if (debug){
		fprintf(stderr, "Reading alignments started.");
	}
	int r;
    while ((r = sam_read1(in, header, b)) >= 0) { // read one alignment from `in'

		printf("pos %d, chrom %d, map q %d, flag %d, name %s \n",
			b->core.pos+1,b->core.tid, b->core.qual, b->core.flag, b->data+b->core.l_qname-2);

		const char* rname=(char*)b->data+b->core.l_qname;
		if (strcmp(rname,CMD_FLUSH)==0){
			// todo: flush fasta
			break;
		}

    	if (((b->core.flag & BAM_FUNMAP)==0) && (b->core.qual)>=cps.min_mapq){

    	}
    	else {
    		if (debug){
				fprintf(stderr, "Read '%s' is not used for updating counters\n", rname);
    		}
    	}

		const auto cigar = bam_get_cigar(b);

        for (int k = 0; k < b->core.n_cigar; k++)
		{
			const int op = bam_cigar_op(cigar[k]);
			const int ol = bam_cigar_oplen(cigar[k]);

			printf("%d%c",ol,BAM_CIGAR_STR[op]);

		}
        printf("\n");
    }
	if (debug){
		fprintf(stderr, "Reading alignments finished.");
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
