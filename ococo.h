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

#endif

