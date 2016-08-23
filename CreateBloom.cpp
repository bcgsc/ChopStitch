#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <cstdlib>
#include <getopt.h>
#include <cassert>
#include <cmath>

#include "nthash.hpp"
#include "Uncompress.h"

#ifdef _OPENMP
# include <omp.h>
#endif


using namespace std;

namespace opt {
unsigned nThrd=1;
unsigned kmLen=64;
unsigned rBuck=4194304;
unsigned rBits=22	;
unsigned sBits=16;
size_t totKmer=0;
bool canon=false;
bool samH=true;
}

static const char shortopts[] = "t:s:r:k:c";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
    { "threads",	required_argument, NULL, 't' },
    { "kmer",	required_argument, NULL, 'k' },
    { "rbit",	required_argument, NULL, 'r' },
    { "sbit",	required_argument, NULL, 's' },
    { "canonical",	no_argument, NULL, 'c' },
    { "help",	no_argument, NULL, OPT_HELP },
    { "version",	no_argument, NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};

void getEfq(size_t &parkCount, uint8_t *m_Counter, std::ifstream &in) {
    bool good = true;
    for(string seq, hseq; good;) {
        good = getline(in, seq);
        good = getline(in, hseq);
        good = getline(in, hseq);
        if(good && seq.length()>=opt::kmLen) 
			ntRead(seq, parkCount, m_Counter);      
        good = getline(in, hseq);
    }
}

void getEfa(size_t &parkCount, uint8_t *m_Counter, std::ifstream &in) {		
	bool good = true;
    for(string seq, hseq; good;) {		
		string line;
		good = getline(in, seq);
        while(good&&seq[0]!='>') {
			line+=seq;
			good = getline(in, seq);
		}
        if(line.length()>=opt::kmLen) 
			ntRead(line, parkCount, m_Counter);
    }
}

void loadSeqr(BloomFilter & myFilter, const string & seq) {
    if (seq.size() < opt::kmerLen) return;

    string kmer = seq.substr(0,opt::kmerLen);
    uint64_t fhVal, rhVal;
	if(!dbFilter.contains(kmer.c_str()))
		dbFilter.insert(kmer.c_str());
	else
		sbFilter.insert(kmer.c_str());
	
    //myFilter.insert(kmer.c_str(), fhVal, rhVal);
    for (size_t i = 1; i < seq.size() - opt::kmerLen + 1; i++) {
        myFilter.insertF(fhVal, seq[i-1], seq[i+opt::kmerLen-1]);
    }
}

void querySeqr(BloomFilter & myFilter, const string & seq) {
    if (seq.size() < opt::kmerLen) return;
    string kmer = seq.substr(0,opt::kmerLen);
    uint64_t fhVal;
    if(myFilter.containsF(kmer.c_str(), fhVal)) {
        #pragma omp atomic
        ++fHit;
    }
    for (size_t i = 1; i < seq.size() - opt::kmerLen + 1; i++) {
        if(myFilter.containsF(fhVal, seq[i-1], seq[i+opt::kmerLen-1])) {
            #pragma omp atomic
            ++fHit;
        }
    }
}

void loadBfs(BloomFilter &dbFilter, BloomFilter &sbFilter, const char* faqFile) {
    getFtype(faqFile);
    ifstream uFile(faqFile);
    bool good = true;
    #pragma omp parallel
    for(string line; good;) {
        #pragma omp critical(uFile)
        good = getSeq(uFile, line);
        if(good) {
		loadSeqr(myFilter, line);
        }
    }
    uFile.close();
}


int main(int argc, char** argv) {
	double sTime = omp_get_wtime();
	
    bool die = false;
    for (int c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) {
        case '?':
            die = true;
            break;
        case 't':
            arg >> opt::nThrd;
            break;
        case 's':
            arg >> opt::sBits;
            break;
        case 'r':
            arg >> opt::rBits;
            break;
        case 'k':
            arg >> opt::kmLen;
            break;
        case 'c':
            opt::canon=true;
            break;
        case OPT_HELP:
            std::cerr << USAGE_MESSAGE;
            exit(EXIT_SUCCESS);
        case OPT_VERSION:
            std::cerr << VERSION_MESSAGE;
            exit(EXIT_SUCCESS);
        }
        if (optarg != NULL && !arg.eof()) {
            std::cerr << PROGRAM ": invalid option: `-"
                      << (char)c << optarg << "'\n";
            exit(EXIT_FAILURE);
        }
    }
    if (argc - optind < 1) {
        std::cerr << PROGRAM ": missing arguments\n";
        die = true;
    }
    if (die) {
        std::cerr << "Try `" << PROGRAM << " --help' for more information.\n";
        exit(EXIT_FAILURE);
    }
    vector<string> inFiles;
    for (int i = optind; i < argc; ++i) {
        string file(argv[i]);
        if(file[0]=='@') {
            string inName;
            ifstream inList(file.substr(1,file.length()).c_str());
            while(getline(inList,inName))
                inFiles.push_back(inName);
        }
        else
            inFiles.push_back(file);
    }

	size_t dbfSize,sbfSize;
	getCardinality(dbfSize,sbfSize);
	BloomFilter dbFilter(dbfSize, opt::nhash, opt::kmerLen);
	BloomFilter sbFilter(sbfSize, opt::nhash, opt::kmerLen);
	
	popBf(dbFilter,sbFilter,inFiles);
	
	



       
#ifdef _OPENMP
    omp_set_num_threads(opt::nThrd);
#endif

cerr << "time(sec): " <<setprecision(4) << fixed << omp_get_wtime() - sTime << "\n";
    return 0;
}
