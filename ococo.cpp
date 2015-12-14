#include "ococo.h"

KSEQ_INIT(gzFile, gzread);

void boost_logging_init()
{
	logging::core::get()->set_filter
	(
		logging::trivial::severity >= logging::trivial::warning
	);
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


/*
stats_t::stats_t():
		n_seqs(0),
		seq_used(NULL),
		seq_len(NULL),
		seq_name(NULL),
		seq_comment(NULL),
		counters(NULL)
{}
*/


stats_t::stats_t(bam_hdr_t &h):
		n_seqs(h.n_targets),
		seq_used(new bool[n_seqs]()),
		seq_len(new int32_t[n_seqs]()),
		seq_name(new char*[n_seqs]()),
		seq_comment(new char*[n_seqs]()),
		counters(new counter_t*[n_seqs]())
{
	for (int s=0;s<n_seqs;s++){
		seq_len[s]=h.target_len[s];
		printf("len %d\n",h.target_len[s]);
		seq_used[s]=true;
		//fprintf(stderr,"allocating %d chars\n",seq_len[i]);
		const int seq_len_name=strlen(h.target_name[s]);
		seq_name[s]=new char[seq_len_name+1];
		//printf("name: %s\n",h.target_name[i]);
		memcpy(seq_name[s], h.target_name[s],seq_len_name+1);

		counters[s]=new counter_t[seq_len[s]]();
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


bool stats_t::check_state(){
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


bool stats_t::check_headers_fai(const string &fai_fn){
	if (!check_state()) return false;

	//todo
	return true;
}


bool stats_t::check_headers_bam_hdr(const bam_hdr_t &h){
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


int stats_t::export_stats(const string &stats_fn){
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


int stats_t::generate_consensus(const string &fasta_fn) {
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
			//fasta_buffer[j]='A';
			const int counter_a=_COUNTER_CELL_VAL(counters[s][i],seq_nt16_table[(int)'A']);
			const int counter_c=_COUNTER_CELL_VAL(counters[s][i],seq_nt16_table[(int)'C']);
			const int counter_g=_COUNTER_CELL_VAL(counters[s][i],seq_nt16_table[(int)'G']);
			const int counter_t=_COUNTER_CELL_VAL(counters[s][i],seq_nt16_table[(int)'T']);
			fasta_buffer[j]=rand_nucl(counter_a,counter_c,counter_g,counter_t);
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


void stats_t::debug_print_counters(){
	for(int s=0;s<n_seqs;s++){
		fprintf(stderr,"%s\n",seq_name[s]);
		for (int i=0;i<seq_len[s];i++){
			fprintf(stderr,"%8d %04x \n",i,counters[s][i]);
		}
	}
}
