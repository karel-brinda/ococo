#ifndef _CALL_VARIANTS
#define _CALL_VARIANTS

#define FASTA_WIDTH 60
#define MAX_COVERAGE 1000
#define QUALITY_TO_INT(__stringOfQualities__,__position__)	(int(__stringOfQualities__[__position__]-33))

float parikhMinRate = 0.6;

#include <fstream>
#include <sstream>
#include <list>

#include <iostream>
#include <cstdlib>
#include <cstring>
#include <string>
#include <exception>

#include <boost/format.hpp>

#include <boost/program_options.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

using namespace std;


enum calling_method {
    CM_UNKNOWN = 0, CM_PARIKH
};

inline void error_message_exit(const string error) {
    cerr << "Error: " << error << endl;
    exit(EXIT_FAILURE);
}


inline bool is_genomic_char(const char &character){
    switch (character){
        case 'a':
        case 'A':
        case 'c':
        case 'C':
        case 'g':
        case 'G':
        case 't':
        case 'T':
        case 'n':
        case 'N':
            return true;
        default:
            return false;
    }
}

ifstream& openIfStream(string filename) {
    ifstream* is = new ifstream(filename.c_str());
    if (is->fail()) {
        error_message_exit("Unable to open file '" + filename + "'");
    }
    return *is;
}

typedef struct PileupLineType
{
	string chr;
	string bases;
	string qualities;
	int position;
	int coverage;

	int aSize;
	int ASize;
	//int AaSize;
	int cSize;
	int CSize;
	//int CcSize;
	int gSize;
	int GSize;
	//int GgSize;
	int tSize;
	int TSize;
	//int TtSize;
	int astSize;

	/* = aSize + ASize + cSize + ...*/
	int strictCoverage;

	int A[MAX_COVERAGE];
	int a[MAX_COVERAGE];
	int C[MAX_COVERAGE];
	int c[MAX_COVERAGE];
	int G[MAX_COVERAGE];
	int g[MAX_COVERAGE];
	int T[MAX_COVERAGE];
	int t[MAX_COVERAGE];
	int ast[MAX_COVERAGE];

} PileupLineType;

inline char call_variant_parikh(const char &originalBase, const PileupLineType& pileupLine, const int &minCoverage, const float &majority);

/*
    Read a new line from a input and save info about the pileup into pileupLine
*/
inline int load_next_pileup_line(PileupLineType &pileupLine, const int &minBaseQuality, istream &input);

/*
    Get debug information about a pileup line
*/
string get_debug_info_pileup_line(const PileupLineType &pileupLine);

/*
    Get a VCF header
*/
string vcf_header ();

/*
    Get a VCF line
*/
string vcf_line (const string &chr, int pos, char oldNucl, char newNucl);

#endif

