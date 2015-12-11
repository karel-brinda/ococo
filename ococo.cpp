#include "ococo.h"
#include <zlib.h>  


KSEQ_INIT(gzFile, gzread);

int main(int argc, const char* argv[])
{

	/*
	 * Default configuration.
	 */

	cons_params_t cps = {
		1, // min_mapq
		0, // min_baseq
	};

	bool debug=false;
	string fasta0_fn;
	string fasta1_fn;
	string stats_fn;
	string sam_fn;

	/*
	 * Parse command-line parameters.
	 */

	try
	{
		namespace po = boost::program_options;

		po::positional_options_description pos;
		pos.add("input-file", -1);

		po::options_description vol("OPTIONS description");

		vol.add_options()
				("input,i", po::value<string>(&sam_fn)->required(), "Input SAM/BAM file (- for standard input).")
				("fa0,f", po::value<string>(&fasta0_fn), "Initial FASTA reference file.")
				("fa1,g", po::value<string>(&fasta1_fn), "Up-to-date FASTA reference.")
				("stats,s", po::value<string>(&stats_fn), "File with up-to-date statistics.")
				//("algorithm,a", po::value<string>(&alg), "Algorithm for updates: majority / randomized [majority]")
				//("counter-size,s", po::value<int>(&counterSize), "Size of counter per nucleotide in bits [3]")
				//("min-coverage,c", po::value<int>(&minCoverage), "Minimal coverage [3]")
				("min-map-qual,m", po::value<int>(&cps.min_mapq), "Minimal mapping quality. [1]")
				("min-base-qual,b", po::value<int>(&cps.min_baseq), "Minimal base quality. [0]")
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
	 * Read SAM headers.
	 */

	hts_itr_t *iter=NULL;

	samFile *in = NULL;
	bam1_t *b= NULL;
	bam_hdr_t *header = NULL;

	in = sam_open(sam_fn.c_str(), "r");
	if(in==NULL) {
		error_exit("Problem with opening input ('%s').\n", sam_fn.c_str());
	}
	if ((header = sam_hdr_read(in)) == 0){
		error_exit("SAM headers are missing or corrupted.\n");
	}

	stats_t stats(*header);

	/*
	 * Load FASTA and stats.
	 */
	if (debug){
		if(fasta0_fn.size()>0){
			fprintf(stderr, "Loading FASTA: %s.\n",fasta0_fn.c_str());
		}
		else{
			fprintf(stderr, "No FASTA provided.\n");
		}
	}
	if (fasta0_fn.size()>0){
		stats.load_headers_fa(fasta0_fn,2);
	}


	/*
	 * Read alignments (main loop).
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
	 * Save FASTA and stats.
	 */

	/*
	 * Free memory.
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
		n_seqs(0),
		seq_used(NULL),
		seq_len(NULL),
		seq_name(NULL),
		seq_comment(NULL),
		counters(NULL)
{
}

stats_t::stats_t(bam_hdr_t &h):
		n_seqs(h.n_targets),
		seq_used(new bool[n_seqs]()),
		seq_len(new uint16_t[n_seqs]()),
		seq_name(new char*[n_seqs]()),
		seq_comment(new char*[n_seqs]()),
		counters(new counter_t*[n_seqs]())
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
	delete[] seq_used;
	delete[] seq_len;
	delete[] seq_name; 
	delete[] seq_comment; 
	delete[] counters; 
}


int stats_t::load_headers_fa(const string &fasta_fn, int weight) {
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

bool stats_t::check_headers_fai(const string &fai_fn){
	//todo
	return true;
}
bool stats_t::check_headers_bam_hdr(const bam_hdr_t &h){
	//todo
	return true;
}

int stats_t::import_stats(const string &stats_fn){
	FILE *fo=fopen(stats_fn.c_str(),"r");

	char delim[stats_delim_l]={};
	uint16_t seq_name_l;
	uint16_t seq_comment_l;

	uint16_t n_seqs_loaded;
	bool     seq_used_loaded;
	uint16_t seq_len_loaded;
	uint16_t seq_name_l_loaded;
	uint16_t seq_comment_l_loaded;

	assert(seq_used!=NULL);
	assert(seq_len!=NULL);
	assert(seq_name!=NULL);
	assert(seq_comment!=NULL);

	/* number of seqs */
	fread(&n_seqs_loaded,sizeof(uint16_t),1,fo);
	assert(n_seqs_loaded = n_seqs);

	for(int i=0;i<n_seqs;i++){

		/* delimiter */
		fread(delim,sizeof(char),stats_delim_l,fo);
		assert(delim[0]=='\255');
		assert(delim[stats_delim_l-1]=='\255');

		/* lengts */
		seq_name_l=strlen(seq_name[i]);
		seq_comment_l=strlen(seq_comment[i]);
		fread(&seq_used_loaded,sizeof(bool),1,fo);
		assert(seq_used_loaded==seq_used[i]);
		fread(&seq_len_loaded,sizeof(uint16_t),1,fo);
		assert(seq_len_loaded==seq_len[i]);
		fread(&seq_name_l_loaded,sizeof(uint16_t),1,fo);
		assert(seq_name_l_loaded==seq_name_l);
		fread(&seq_comment_l_loaded,sizeof(uint16_t),1,fo);
		assert(seq_comment_l_loaded==seq_comment_l);

		/* strings */
		assert(seq_name[i] != NULL);
		char seq_name_loaded[seq_name_l+1];
		fread(seq_name_loaded,sizeof(char),seq_name_l+1,fo);
		assert(seq_comment[i] != NULL);
		char seq_comment_loaded[seq_comment_l+1];
		fread(seq_comment_loaded,sizeof(char),seq_comment_l+1,fo);

		/* counters */
		assert(counters[i]!=NULL);
		fread(counters[i],sizeof(counter_t),seq_len[i],fo);
	}


	fclose(fo);
	return 0;
}

int stats_t::export_stats(const string &stats_fn){
	FILE *fo=fopen(stats_fn.c_str(),"w+");

	char delim[stats_delim_l]={};
	uint16_t seq_name_l;
	uint16_t seq_comment_l;

	delim[0]='\255';
	delim[stats_delim_l-1]='\255';

	/* number of seqs */
	fwrite(&n_seqs,sizeof(uint16_t),1,fo);

	for(int i=0;i<n_seqs;i++){

		/* delimiter */
		fwrite(delim,sizeof(char),stats_delim_l,fo);

		/* lengts */
		seq_name_l=strlen(seq_name[i]);
		seq_comment_l=strlen(seq_comment[i]);
		fwrite(&seq_used[i],sizeof(bool),1,fo);
		fwrite(&seq_len[i],sizeof(uint16_t),1,fo);
		fwrite(&seq_name_l,sizeof(uint16_t),1,fo);
		fwrite(&seq_comment_l,sizeof(uint16_t),1,fo);

		/* strings */
		assert(seq_name[i] != NULL);
		fwrite(seq_name[i],sizeof(char),seq_name_l+1,fo);
		assert(seq_comment[i] != NULL);
		fwrite(seq_comment[i],sizeof(char),seq_comment_l+1,fo);

		/* counters */
		assert(counters[i] != NULL);
		fwrite(counters[i],sizeof(counter_t),seq_len[i],fo);
	}

	fclose(fo);
	return 0;
}


int stats_t::generate_fasta(const string &fasta_fn) {
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
			if(j==fasta_line_l-1 || i==seq_len[s]-1){
				fwrite(fasta_buffer,1,j+1,fp);
				fwrite("\n",1,1,fp);
				j=-1;
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
