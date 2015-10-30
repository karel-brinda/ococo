#include "ococo.h"

int main(int argc, const char* argv[])
{
	string alg="";
	string oldRefFn="";
	//string pileupFn="";
	int minBaseQ=0;
	int minCoverage=3;
	float acceptanceLevel=0.6;
	//string vcfFn="";

    bool debug=false;

    calling_method method=CM_UNKNOWN;

	try
	{
		namespace po = boost::program_options;

		po::options_description mandatory("Mandatory options");
		mandatory.add_options()
                ("calling-alg,a", po::value<string>(&alg)->required(), "parikh")
                ("reference,r", po::value<string>(&oldRefFn)->required(), "FASTA reference file")
		;
		//("pileup,p", po::value<string>(&pileupFn)->required(), "Pileup file")
		

		po::options_description voluntary("Voluntary options");
		voluntary.add_options()
                ("min-coverage,c", po::value<int>(&minCoverage), "Minimal coverage [3]")
                ("min-base-qual,b", po::value<int>(&minBaseQ), "Minimal base quality [0]")
                ("accept-level,l", po::value<float>(&acceptanceLevel), "Acceptance level [0.60]")
                //("vcf,v", po::value<string>(&vcfFn), "VCF file for called variants")
                ("debug,d", "Debug (print information about every position)")
		;

		po::options_description all("OPTIONS");
		all.add(mandatory).add(voluntary);

		po::variables_map vm;
		try
		{
			po::store(po::parse_command_line(argc, argv, all),vm); // can throw

            if(vm.count("debug")){
                debug=true;
                cerr << "DEBUG MODE" << endl;
            }

            po::notify(vm); // throws on error, so do after help in case there are any problems

		}
		catch(po::error& e)
		{
			std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
			std::cerr << all << std::endl;
			return EXIT_FAILURE;
		}

	}
	catch(std::exception& e)
    {
        error_message_exit(str(boost::format("Unhandled Exception reached the top of main: %1%,\napplication will now exit") % e.what()) );
	}


	ifstream &iFa=openIfStream(oldRefFn);
    //ofstream *oVcf;
    /*if(vcfFn != ""){
        oVcf = new ofstream(vcfFn.c_str());
        if (oVcf->fail()) {
            error_message_exit("Unable to open file '" + vcfFn + "'");
        }
        *oVcf << vcf_header();
    }*/
    cout << vcf_header();

	istream &iPileup=cin;

	if (alg.compare("parikh")==0)
	{
		method=CM_PARIKH;
	}
	else
	{
		error_message_exit("Unknown method '" + alg +"'");
	}

	/* loading the first pileup line */
	PileupLineType pileupLine;
	string seqName = "";
	int okIPileup = load_next_pileup_line(pileupLine, minBaseQ, iPileup);
    if(debug){
        cerr << "DEBUG: " << get_debug_info_pileup_line(pileupLine) << endl;
    }

	/* output buffer */
	int oFaPosition = 0;
	string oFaLine = string(FASTA_WIDTH, ' ');

	/* reading input fasta */
	string iFaLine;
	while (getline(iFa, iFaLine))
	{
		if (iFaLine.length() > 0)
		{
			if (iFaLine[0] == '>')
			{
				/*
				    fasta header
				*/

                auto spacePos=iFaLine.find_first_of(" \t");
                if (spacePos==string::npos){
                    spacePos=iFaLine.length();
                }

				seqName = iFaLine.substr(1, spacePos - 1);
				oFaPosition = 0;
				//cout << ">" << seqName << endl;
			}
			else
			{
				/*
				    fasta body
				*/
				for (int i = 0; i < (int)iFaLine.length(); i++)
				{
					/* is a nucleotide?  */
					if (iFaLine[i] != '\r' && iFaLine[i] != '\n' && iFaLine[i] != ' ')
					{
                        if(!is_genomic_char(iFaLine[i])){
                            error_message_exit(string() + "Unknown character in Fasta.\n"
                            + "Error in Function '" + __FUNCTION__
                                    + "' (" __FILE__ ", " + BOOST_PP_STRINGIZE(__LINE__)
                                        + ")");
                        }

						/*
							is there a corresponding position in pileup?
							(1+ ... because it is 0-based)
						*/
						if (okIPileup && pileupLine.position == 1 + oFaPosition && pileupLine.chr.compare(seqName) == 0)
						{
							/************
								yes => try to update, load next pileup line
							**************/
							char oldChar=iFaLine[i];
							char newChar=' ';
							switch (method)
							{
							case CM_PARIKH:
								newChar=call_variant_parikh(iFaLine[i], pileupLine, minCoverage, acceptanceLevel);
								break;
							default:
                                error_message_exit(string() + "An unsupported method.\n"
                                        + "Error in Function '" + __FUNCTION__
                                        + "' (" __FILE__ ", " + BOOST_PP_STRINGIZE(__LINE__)
                                        + ")");
								break;

							}

							if(tolower(newChar) != tolower(oldChar))
							{
                                /*
                                    There was an update.
                                 */
                                //cerr << "update at pos " << oFaPosition << endl;
                                /*if(oVcf!= 0){
							        *oVcf << vcf_line(seqName, oFaPosition, oldChar, newChar) << endl;
                                }*/
						        cout << vcf_line(seqName, oFaPosition, oldChar, newChar) << endl;
							}

							/* modify corresponding part of buffer, load next line from pileup */
							oFaLine[(oFaPosition++) % FASTA_WIDTH] = newChar;
							okIPileup = load_next_pileup_line(pileupLine, minBaseQ, iPileup);
                            if(debug){
                                cerr << "DEBUG: " << get_debug_info_pileup_line(pileupLine) << endl;
                            }
						}
						else
						{
							/* no, print fasta char */
							oFaLine[(oFaPosition++) % FASTA_WIDTH] = tolower(iFaLine[i]);
						}

						/* fasta columning :) */
						if (oFaPosition % FASTA_WIDTH == 0)
						{
							//cout << oFaLine << endl;
						}
					}
				} // for i
			} // if(iFaLine[0]=='>')
		}
	}
	/* the rest of the last line (flushing the buffer) */
	//cout << oFaLine.substr(0, oFaPosition % FASTA_WIDTH) << endl;

    if (!iPileup.eof()){
        cerr<<"Some lines of the pileup were not read, e.g.,:" << endl;
        for(int i=1;i<4;i++){
            string line;
            getline(iPileup, line);
            cerr << line << endl;
            if (iPileup.eof())
            {
                break;
            }
        }
        error_message_exit("Error appered.");
    }
	iFa.close();
    //oVcf->close();

	return 0;
}

inline char call_variant_parikh(const char &originalBase, const PileupLineType& pileupLine, const int &minCoverage, const float &majority)
{
	//TODO: assert
	if (pileupLine.strictCoverage >= minCoverage)
	{

		float majority = (parikhMinRate * pileupLine.strictCoverage);

		char candidat = tolower(originalBase);

		if (pileupLine.aSize + pileupLine.ASize > majority)
		{
			candidat = 'A';
		}
		else if (pileupLine.cSize + pileupLine.CSize > majority)
		{
			candidat = 'C';
		}
		else if (pileupLine.gSize + pileupLine.GSize > majority)
		{
			candidat = 'G';
		}
		else if (pileupLine.tSize + pileupLine.TSize > majority)
		{
			candidat = 'T';
		}

		if (tolower(candidat) == tolower(originalBase))
		{
			/* no update */
			return tolower(originalBase);
		}
		else
		{
			/* update */
			return candidat;
		}
	}
	else
	{
		/* no update */
		return tolower(originalBase);
	}
}

string vcf_header ()
{
	/*
		##contig=<ID=chr1,length=1075152>\n\
		##INFO=<ID=mt,Number=1,Type=String,Description=\"Variant Type: SUBSTITUTE/INSERT/DELETE\">\n\
	*/

	return "\
##fileformat=VCFv4.1\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
";
}

string vcf_line (const string &chr, int pos, char oldNucl, char newNucl)
{
	stringstream ss;
	ss << chr << "\t" << pos+1 << "\t" << ".\t" << oldNucl << "\t"  << newNucl
	   <<  "\t100\tPASS\tmt=SUBSTITUTE";
	return ss.str();
}


string get_debug_info_pileup_line(const PileupLineType &pileupLine)
{

	string AValues;
	string aValues;
	string CValues;
	string cValues;
	string GValues;
	string gValues;
	string TValues;
	string tValues;

	{
		stringstream ss;
		for (int i = 0; i < pileupLine.ASize; i++)
			ss << pileupLine.A[i] << " ";
		AValues = ss.str();
	}
	{
		stringstream ss;
		for (int i = 0; i < pileupLine.aSize; i++)
			ss << pileupLine.a[i] << " ";
		aValues = ss.str();
	}
	{
		stringstream ss;
		for (int i = 0; i < pileupLine.CSize; i++)
			ss << pileupLine.C[i] << " ";
		CValues = ss.str();
	}
	{
		stringstream ss;
		for (int i = 0; i < pileupLine.cSize; i++)
			ss << pileupLine.c[i] << " ";
		cValues = ss.str();
	}
	{
		stringstream ss;
		for (int i = 0; i < pileupLine.GSize; i++)
			ss << pileupLine.G[i] << " ";
		GValues = ss.str();
	}
	{
		stringstream ss;
		for (int i = 0; i < pileupLine.gSize; i++)
			ss << pileupLine.g[i] << " ";
		gValues = ss.str();
	}
	{
		stringstream ss;
		for (int i = 0; i < pileupLine.TSize; i++)
			ss << pileupLine.T[i] << " ";
		TValues = ss.str();
	}
	{
		stringstream ss;
		for (int i = 0; i < pileupLine.tSize; i++)
			ss << pileupLine.t[i] << " ";
		tValues = ss.str();
	}

	string output = (boost::format(
				     "\n---------------------------------------\n\
chr:             %1%\n\
bases:           %2%\n\
qualities:       %3%\n\
position:        %4%\n\
coverage:        %5%\n\
strict coverage: %6%\n\
\n\
A (%7%):         %8%\n\
a (%9%):         %10%\n\
C (%11%):         %12%\n\
c (%13%):         %14%\n\
G (%15%):         %16%\n\
g (%17%):         %18%\n\
T (%19%):         %20%\n\
t (%21%):         %22%\n\
---------------------------------------\n")
			     % pileupLine.chr
			     % pileupLine.bases
			     % pileupLine.qualities
			     % pileupLine.position
			     % pileupLine.coverage
			     % pileupLine.strictCoverage
			     % pileupLine.ASize % AValues
			     % pileupLine.aSize % aValues
			     % pileupLine.CSize % CValues
			     % pileupLine.cSize % cValues
			     % pileupLine.GSize % GValues
			     % pileupLine.gSize % gValues
			     % pileupLine.TSize % TValues
			     % pileupLine.tSize % tValues).str();

	return output;
}

inline int load_next_pileup_line(PileupLineType &pileupLine, const int &minBaseQuality, istream &input = cin)
{
	/* read line, if end of file => stop */
	string line;
	getline(input, line);
	if (input.eof())
	{
		return 0;
	}
	char *cstr;
	cstr = new char[line.size() + 1];
	strcpy(cstr, line.c_str());

	/* save new data into pileupLine*/
	pileupLine.chr = strtok(cstr, "\t");
	pileupLine.position = atoi(strtok(NULL, "\t"));
	char refBase = strtok(NULL, "\t")[0];
	pileupLine.coverage = atoi(strtok(NULL, "\t"));

	/* delete old data from pileupLine */
	pileupLine.ASize = 0;
	pileupLine.aSize = 0;
	pileupLine.CSize = 0;
	pileupLine.cSize = 0;
	pileupLine.GSize = 0;
	pileupLine.gSize = 0;
	pileupLine.TSize = 0;
	pileupLine.tSize = 0;
	pileupLine.astSize = 0;
	pileupLine.bases = "";
	pileupLine.qualities = "";
	pileupLine.strictCoverage = 0;

	/* continue only when coverage is enough */
	if (pileupLine.coverage > 0)
	{

		pileupLine.bases = strtok(NULL, "\t");
		pileupLine.qualities = strtok(NULL, "\t");
		/* load bases */
		int posQ = 0;
		/*
			state:
				0 ... normal
				1 ... reading number after + or -
		 */
		int state = 0;
		int numberStartsAt=0;
		for (int posB = 0; posB < (int)pileupLine.bases.length(); posB++)
		{
			/* reading number? */
			if (state==1)
			{
				/* not reading number? */
				if (!isdigit(pileupLine.bases[posB]))
				{
					/* set normal state, skip enough positions */
					state = 0;
					posB += atoi(pileupLine.bases.substr(numberStartsAt, posB - numberStartsAt).c_str());
				}
			}
			/* already normal reading? */
			if (state==0)
			{
				int currentQuality = QUALITY_TO_INT(pileupLine.qualities, posQ);
				char currentChar = pileupLine.bases[posB];
				if (currentChar=='.')
				{
					currentChar=toupper(refBase);
				}
				else if (currentChar==',')
				{
					currentChar=tolower(refBase);
				}
				switch (currentChar)
				{
				case '\n':
					break;
				case '\r':
					break;
				case '$':
					break;
				case '^':
					//cout << "huraa" << endl;
					/* skip the quality information */
					++posB;
					break;
				case ' ':
					break;
				case 'A':
					if (minBaseQuality <= currentQuality)
					{
						pileupLine.A[pileupLine.ASize++] = currentQuality;
						pileupLine.strictCoverage++;
					}
					posQ++;
					break;
				case 'a':
					if (minBaseQuality <= currentQuality)
					{
						pileupLine.a[pileupLine.aSize++] = currentQuality;
						pileupLine.strictCoverage++;
					}
					posQ++;
					break;
				case 'C':
					if (minBaseQuality <= currentQuality)
					{
						pileupLine.C[pileupLine.CSize++] = currentQuality;
						pileupLine.strictCoverage++;
					}
					posQ++;
					break;
				case 'c':
					if (minBaseQuality <= currentQuality)
					{
						pileupLine.c[pileupLine.cSize++] = currentQuality;
						pileupLine.strictCoverage++;
					}
					posQ++;
					break;
				case 'G':
					if (minBaseQuality <= currentQuality)
					{
						pileupLine.G[pileupLine.GSize++] = currentQuality;
						pileupLine.strictCoverage++;
					}
					posQ++;
					break;
				case 'g':
					if (minBaseQuality <= currentQuality)
					{
						pileupLine.g[pileupLine.gSize++] = currentQuality;
						pileupLine.strictCoverage++;
					}
					posQ++;
					break;
				case 'T':
					if (minBaseQuality <= currentQuality)
					{
						pileupLine.T[pileupLine.TSize++] = currentQuality;
						pileupLine.strictCoverage++;
					}
					posQ++;
					break;
				case 't':
					if (minBaseQuality <= currentQuality)
					{
						pileupLine.t[pileupLine.tSize++] = currentQuality;
						pileupLine.strictCoverage++;
					}
					posQ++;
					break;
				case '*':
					if (minBaseQuality <= currentQuality)
					{
						pileupLine.ast[pileupLine.astSize++] = currentQuality;
						pileupLine.strictCoverage++;
					}
					posQ++;
					break;
				case '-':
				case '+':
					/* must be followed by a number */
					state = 1;
					numberStartsAt = posB + 1;
					break;

				}
			} // switch (state)
		} // for posB
		assert(posQ == (int)pileupLine.qualities.length() && "Length of qualities and bases are not same. Probably a bug of the program.");
	} // if (pileupLine.coverage>0)

	//printPileupLine(pileupLine,cerr);

	return 1;
}
