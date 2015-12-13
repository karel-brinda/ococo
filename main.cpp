#include "ococo.h"


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

		po::options_description vol("Ococo: On-line consensus calling.");

		vol.add_options()
			    //("help", "Print help message")
				("input,i", po::value<string>(&sam_fn)->required(), "Input SAM/BAM file (- for standard input).")
				("init-fasta,f", po::value<string>(&fasta0_fn), "Initial FASTA reference.")
				("consensus-fasta,c", po::value<string>(&fasta1_fn), "Consensus (up-to-date FASTA reference).")
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
			po::notify(vm); // throws on error, so do after help in case there are any problems

			if(vm.count("debug")){
				debug=true;
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

	in = sam_open(sam_fn.c_str(), "r");
	if(in==NULL) {
		error_exit("Problem with opening input ('%s').\n", sam_fn.c_str());
	}
	if ((header = sam_hdr_read(in)) == 0){
		error_exit("SAM headers are missing or corrupted.\n");
	}

	stats_t stats(*header);
	assert(stats.check_state());

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
	if (stats_fn.size()>0 && file_exists(stats_fn)){
		stats.import_stats(stats_fn);
	}
	else{
		if (fasta0_fn.size()>0){
			stats.load_headers_fa(fasta0_fn,2);
		}
	}

	/*
	 * Process alignments.
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

		if (((flags & BAM_FUNMAP)==0) && stats.seq_used[chrom] && mapq>=cps.min_mapq){

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
	 * Generate consensus and export stats.
	 */
	if (fasta1_fn.size()>0){
		stats.generate_consensus(fasta1_fn);
	}

	if (stats_fn.size()>0){
		stats.export_stats(stats_fn);
	}


	/*
	 * Free memory.
	 */
	hts_itr_destroy(iter);
	bam_destroy1(b);
	bam_hdr_destroy(header);
	sam_close(in);

	return 0;
}