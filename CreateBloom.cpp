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

#include "BloomFilter.hpp"
//#include "RollingHashIterator.h"

#include "Uncompress.h"

#ifdef _OPENMP
# include <omp.h>
#endif

#define PROGRAM "CreateBloom"

static const char VERSION_MESSAGE[] =
    PROGRAM " Version 1.0.0 \n"
    "Written by Hamid Mohamadi.\n"
    "Copyright 2016 Canada's Michael Smith Genome Science Centre\n";

static const char USAGE_MESSAGE[] =
    "Usage: " PROGRAM " [OPTION]... FILES...\n"
    "Estimates the number of k-mers in FILES(>=1).\n"
    "Accepatble file formats: fastq, fasta, sam, bam, gz, bz, zip.\n"
    "\n"
    " Options:\n"
    "\n"
    "  -t, --threads=N	use N parallel threads [1]\n"
    "  -k, --kmer=N	the length of kmer [64]\n"
    "  -c, --canonical	get the estimate for cannonical form\n"
    "      --help	display this help and exit\n"
    "      --version	output version information and exit\n"
    "\n"
    "Report bugs to hmohamadi@bcgsc.ca.\n";

using namespace std;

namespace opt {
unsigned nThrd=1;
unsigned ibits=8;
unsigned nhash=3;
unsigned kmLen=50;
}

static const char shortopts[] = "t:k:b:h:c";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
    { "threads",	required_argument, NULL, 't' },
    { "kmer",	required_argument, NULL, 'k' },
    { "ibit",	required_argument, NULL, 'b' },
    { "hash",	required_argument, NULL, 'h' },
    { "help",	no_argument, NULL, OPT_HELP },
    { "version",	no_argument, NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};

unsigned getftype(std::ifstream &in) {
    std::string hseq;
    getline(in,hseq);
    if(hseq[0]=='>') {
        return 1;
    }
    if(hseq[0]=='@') {
        return 0;
    }
    return 2;
}

inline void seqLoad(BloomFilter &dbFilter, BloomFilter &sbFilter, const string &seq) {
    if (seq.size() < opt::kmLen) return;
    uint64_t hVal, fhVal=0, rhVal=0;
    for(unsigned seqIndex=0; seqIndex<seq.length()-opt::kmLen+1;) {
        bool rollFlag=(seqIndex==0)?false:true;
        if(seedTab[seq[seqIndex+opt::kmLen-1]]==seedN) {
            seqIndex+=opt::kmLen;
            rollFlag=false;
        }
        if(!rollFlag) {
            bool hGood=false;
            while(!hGood && seqIndex<seq.length()-opt::kmLen+1) {
                string kmer = seq.substr(seqIndex, opt::kmLen);
                unsigned locN=0;
                hGood = NTPC64(kmer.c_str(), opt::kmLen, fhVal, rhVal, hVal, locN);
                seqIndex+=locN+1;
            }
            if(hGood) {
				if(!dbFilter.insert_make_change(hVal))sbFilter.insert(hVal);
            }
        }
        else {
            hVal=NTPC64(fhVal, rhVal, seq[seqIndex-1], seq[seqIndex-1+opt::kmLen], opt::kmLen);
            seqIndex++;            		
			if(!dbFilter.insert_make_change(hVal))sbFilter.insert(hVal);
        }
    }
}

void loadBFfq(BloomFilter &dbFilter, BloomFilter &sbFilter, std::ifstream &in) {
    bool good = true;

    //good = getline(in, seq);
    //good = getline(in, hseq);
    //good = getline(in, hseq);
    //if(good) seqLoad(dbFilter, sbFilter, seq);
    
    
        
    #pragma omp parallel
    for(string seq, hseq; good;) {
        #pragma omp critical(in)
        {
            good = getline(in, hseq);
            good = getline(in, seq);
            good = getline(in, hseq);
            good = getline(in, hseq);
        }
        if(good) {
			 seqLoad(dbFilter, sbFilter, seq);
			/*RollingHashIterator itr(seq, opt::nhash, opt::kmLen);			
			while (itr != itr.end()) {
			//	if(!dbFilter.contains(*itr))
					dbFilter.insert(*itr);
				//else
					//sbFilter.insert(*itr);
				itr++;
			}*/
		 }
    }
}


bool getSeq(std::ifstream &uFile, std::string &line) {
    bool good=false;
    std::string hline;
    line.clear();
    do {
            good=getline(uFile, hline);
            if(hline[0]=='>'&&!line.empty()) break;// !line.empty() for the first rec
            if(hline[0]!='>')line+=hline;
        } while(good);
    
    if(!good&&!line.empty())
            good=true;

    return good;
}

/*void loadBFfa(BloomFilter &dbFilter, BloomFilter &sbFilter, std::ifstream &in) {
    bool good = true;
    #pragma omp parallel
    for(string seq, hseq; good;) {
        #pragma omp critical(in)
        good = getSeq(in,seq);
        //if(good) seqLoad(dbFilter, sbFilter, seq);
        if(good) {
			 //seqLoad(dbFilter, sbFilter, seq);
			RollingHashIterator itr(seq, opt::nhash, opt::kmLen);			
			while (itr != itr.end()) {
				if(!dbFilter.contains(*itr))
					dbFilter.insert(*itr);
				else
					sbFilter.insert(*itr);
				itr++;
			}
		 }
    }
}*/

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
        case 'k':
            arg >> opt::kmLen;
            break;
        case 'b':
            arg >> opt::ibits;
            break;
        case 'h':
            arg >> opt::nhash;
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

#ifdef _OPENMP
    omp_set_num_threads(opt::nThrd);
#endif



    size_t dbfSize=11000000000,sbfSize=4000000000; // Human_smallerBF
    //size_t dbfSize=22000000000,sbfSize=4000000000; // Human
    //size_t dbfSize=300000000,sbfSize=200000000; // C. elegans
    
    BloomFilter dbFilter(dbfSize*opt::ibits, opt::nhash, opt::kmLen);
    BloomFilter sbFilter(sbfSize*opt::ibits, opt::nhash, opt::kmLen);
    
    for (unsigned file_i = 0; file_i < inFiles.size(); ++file_i) {
        std::ifstream in(inFiles[file_i].c_str());
		loadBFfq(dbFilter, sbFilter, in);		
		in.close();
	}
	
	
    cerr << "Load time(sec): " <<setprecision(4) << fixed << omp_get_wtime() - sTime << "\n";
    
	sbFilter.storeFilter("sfilter.bf");

	cout << "h= " << opt::nhash << "\n";
	cout << "k= " << opt::kmLen << "\n";
	cout << "b= " << opt::ibits << "\n";
	cout << "FPR= " << setprecision(4) << fixed << pow(1.0*sbFilter.getPop()/sbfSize/opt::ibits,opt::nhash) << "\n";
	cout << "popcnt of dbFilter and sbFilter: " << dbFilter.getPop() << "\t" << sbFilter.getPop() << "\n";
	
    cerr << "time(sec): " <<setprecision(4) << fixed << omp_get_wtime() - sTime << "\n";
    return 0;
}
