#include "ococo.h"
#include <zlib.h>  


KSEQ_INIT(gzFile, gzread);


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
	int count;


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

stats_t::stats_t():
		n_seqs(0),seq_len(NULL), seq_name(NULL), seq_comment(NULL), counters(NULL)
{
}

stats_t::stats_t(bam_hdr_t &h):
		n_seqs(h.n_targets), seq_len(new int[n_seqs]()), seq_name(new char*[n_seqs]()), seq_comment(new char*[n_seqs]()), counters(new counter_t*[n_seqs]())
{
	for (int i=0;i<n_seqs;i++){
		seq_len[i]=h.target_len[i];
		//fprintf(stderr,"allocating %d chars\n",seq_len[i]);
		const int seq_len_name=strlen(h.target_name[i]);
		seq_name[i]=new char[seq_len_name+1];
		//printf("name: %s\n",h.target_name[i]);
		memcpy(seq_name[i], h.target_name[i],seq_len_name+1);

		counters[i]=new counter_t[seq_len[i]]();
		//printf("seq %s, len %d\n",stats->seqstats[i].name,stats->seqstats[i].length);
		//fprintf(stderr,"ok\n");
	}
}

stats_t::~stats_t(){
	for (int i=0;i<n_seqs;i++){
		delete[] seq_name[i];
		delete[] seq_comment[i];
		delete[] counters[i];
	}
	delete[] seq_len;
	delete[] seq_name; 
	delete[] seq_comment; 
	delete[] counters; 
}


int stats_t::load_fasta(string fasta_fn, int weight) {
	gzFile fp;  
	kseq_t *seq;  
	int l;  
	fp = gzopen(fasta_fn.c_str(), "r");
	seq = kseq_init(fp);
	for(int s=0;(l = kseq_read(seq)) >= 0;s++) {
		/*printf("%s\n",seq_name[s]);
		printf("%s\n",seq->name.s);
		printf("\n");
		printf("%d\n",seq_len[s]);
		printf("%d\n",seq->seq.l);
		printf("\n");*/
		assert(strcmp(seq->name.s,seq_name[s])==0);
		assert((int)seq->seq.l == seq_len[s]);
		if (seq->comment.l && seq_comment[s]==NULL){
			seq_comment[s]=new char[seq->comment.l+1];
			memcpy(seq_comment[s], seq->comment.s,seq->comment.l+1);
		}
		//fprintf(stderr,"seq: %s\n", seq->seq.s);  

		for(unsigned int i=0;i<seq->seq.l;i++){
			assert(counters[s][i]==0);
			const nt16_t nt16=seq_nt16_table[(int)seq->seq.s[i]];
			counters[s][i]=_COUNTER_CELL_SET(0,nt16,weight);
		}
		//if (seq->qual.l) printf("qual: %s\n", seq->qual.s);  
	}  
	kseq_destroy(seq); // STEP 5: destroy seq  
	gzclose(fp); // STEP 6: close the file handler
	return 0;
}

int stats_t::generate_fasta(string fasta_fn) {
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
			//fasta_buffer[j]='A';
			const int counter_a=_COUNTER_CELL_VAL(counters[s][i],seq_nt16_table['A']);
			const int counter_c=_COUNTER_CELL_VAL(counters[s][i],seq_nt16_table['C']);
			const int counter_g=_COUNTER_CELL_VAL(counters[s][i],seq_nt16_table['G']);
			const int counter_t=_COUNTER_CELL_VAL(counters[s][i],seq_nt16_table['T']);
			fasta_buffer[j]=rand_nucl(counter_a,counter_c,counter_g,counter_t);
			if(j==fasta_line_l || i==seq_len[s]-1){
				fwrite(fasta_buffer,1,j,fp);
				fwrite("\n",1,1,fp);
				j=0;
			}
		}

	}

	fclose(fp); // STEP 6: close the file handler
	return 0;
}

void stats_t::debug_print_counters(){
	for(int s=0;s<n_seqs;s++){
		fprintf(stderr,"%s\n",seq_name[s]);
		for (int i=0;i<seq_len[s];i++){
			fprintf(stderr,"%8d %04x \n",i,counters[s][i]);
		}
	}
}



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
		Parse command-line parameters.
	*/

	try
	{
		namespace po = boost::program_options;

		po::positional_options_description pos;
		pos.add("input-file", -1);

		po::options_description vol("OPTIONS description");

		vol.add_options()
				("sam,s", po::value<string>(&sam_fn), "Input SAM file (- for standard input) [-].")
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
		stats.load_fasta(fasta_fn,2);
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

		if (((flags & BAM_FUNMAP)==0) && mapq>=cps.min_mapq){

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


	//printf("writing_file\n");
	stats.generate_fasta("a.fa");
	//stats.debug_print_counters();
	//printf("ok\n");

	hts_itr_destroy(iter);
	bam_destroy1(b);
	bam_hdr_destroy(header);
	sam_close(in);

	return 0;
}
