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

#include "nthist.hpp"
#include "BloomFilter.hpp"
#include "ntHashIterator.hpp"

#include "Uncompress.h"

#ifdef _OPENMP
# include <omp.h>
#endif

#define PROGRAM "CreateBloom"

static const char VERSION_MESSAGE[] =
    PROGRAM " Version 1.0.0 \n"
    "Written by Hamid Mohamadi and Hamza Khan.\n"
    "Copyright 2017 Canada's Michael Smith Genome Science Centre\n";

static const char USAGE_MESSAGE[] =
    "Usage: " PROGRAM " [OPTION]... FILES...\n"
    "Creates a Bloom filter (BF) to find exon-exon junctions.\n"
    "Accepatble file formats: fastq, fasta, sam, bam, gz, bz, zip.\n"
    "\n"
    " Options:\n"
    "\n"
    "  -t, --threads=N	use N parallel threads [1]\n"
    "  -k, --kmer=N	the length of kmer [50]\n"
    "  -d, --fpr1=N	primary BF fpr [0.01]\n"
    "  -s, --fpr2=N	secondary BF fpr [0.01]\n"
    "  -r, --ref	using FASTA reference as input instead of FASTQ reads. Don't use fpr2 in this case\n"
    "      --help	display this help and exit\n"
    "      --version	output version information and exit\n"
    "\n"
    "Report bugs to hmohamadi@bcgsc.ca or hkhan@bcgsc.ca \n";

using namespace std;

namespace opt {
unsigned nThrd=1;
unsigned nhash1;
unsigned nhash2;
unsigned kmLen=50;
size_t m1;
size_t m2;
double fpr1=0.01;
double fpr2=0.01;
bool ref=false;
size_t dbfSize=0;
size_t sbfSize=0;
}

static const char shortopts[] = "t:k:d:s:H:h:r:c";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
    { "threads",	required_argument, NULL, 't' },
    { "kmer",	required_argument, NULL, 'k' },
    { "hash1",	required_argument, NULL, 'H' },
    { "hash2",	required_argument, NULL, 'h' },
    { "fpr1",	required_argument, NULL, 'd' },
    { "fpr2",	required_argument, NULL, 's' },
    { "ref",	no_argument, NULL, 'r' },
    { "help",	no_argument, NULL, OPT_HELP },
    { "version",	no_argument, NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};

bool getFAseq(std::ifstream &uFile, std::string &line) {
    bool good=false;
    std::string hline;
    line.clear();
    do {
        good=static_cast<bool>(getline(uFile, hline));
        if(hline[0]=='>'&&!line.empty())
            break;
        if(hline[0]!='>')
            line+=hline;
    } while(good);
    if(!good&&!line.empty())
        good=true;
    return good;
}

void loadBFfa(std::ifstream &in, BloomFilter &dbFilter) {
    bool good = true;
    #pragma omp parallel
    for(string seq, hseq; good;) {
        #pragma omp critical(in)
        good = getFAseq(in,seq);
        if(good) {
            ntHashIterator itr(seq, opt::nhash1, opt::kmLen);
            while (itr != itr.end()) {
                dbFilter.insert(*itr);
                ++itr;
            }
        }
    }
}

void genBFref(const vector<string> &inFiles) {
    cerr << "k-mer length: " << opt::kmLen << "\n";
    cerr << "Number of distinct k-mers: " << opt::dbfSize << "\n";
    cerr << "Number of k-mers with freq>1: " << opt::sbfSize << "\n";
    cerr << "BF fpr: " << opt::fpr1 << "\n";
    cerr << "BF bits: " << opt::m1 << "\t bits/kmer= " << opt::m1/opt::dbfSize << "\n";
    cerr << "BF hashes: " << opt::nhash1 << "\n";

    BloomFilter dbFilter(opt::m1, opt::nhash1, opt::kmLen);
    for (unsigned file_i = 0; file_i < inFiles.size(); ++file_i) {
        std::ifstream in(inFiles[file_i].c_str());
        loadBFfa(in, dbFilter);
        in.close();
    }
    dbFilter.storeFilter("Bfilter.bf");

    ofstream bfinfo("Bfilter.inf");
    bfinfo << opt::m1 << "\n" << opt::nhash1 << "\n" <<opt::kmLen << "\n" << pow(1.0*dbFilter.getPop()/opt::m1,opt::nhash1);
    bfinfo.close();

    cerr << "BF actual fpr: " << setprecision(4) << fixed << pow(1.0*dbFilter.getPop()/opt::m1,opt::nhash1) << "\n";
    cerr << "Popcnt of bf: " << dbFilter.getPop() << "\n";
}

void loadBFfq(std::ifstream &in, BloomFilter &dbFilter, BloomFilter &sbFilter) {
    bool good = true;
    unsigned maxHash= std::max(opt::nhash1,opt::nhash2);

    #pragma omp parallel
    for(string seq, hseq; good;) {
        #pragma omp critical(in)
        {
            good =static_cast<bool>(getline(in, hseq));
            good = static_cast<bool>(getline(in, seq));
            good = static_cast<bool>(getline(in, hseq));
            good = static_cast<bool>(getline(in, hseq));
        }
        if(good) {
            ntHashIterator itr(seq, maxHash, opt::kmLen);
            while (itr != itr.end()) {
                if(!dbFilter.insert_make_change(*itr))
                    sbFilter.insert(*itr);
                ++itr;
            }
        }

    }

}

void genBFseq(const vector<string> &inFiles) {
    cerr << "k-mer length: " << opt::kmLen << "\n";
    cerr << "Number of distinct k-mers: " << opt::dbfSize << "\n";
    cerr << "Number of k-mers with freq>1: " << opt::sbfSize << "\n";
    cerr << "Primary BF fpr: " << opt::fpr1 << "\n";
    cerr << "Secondary BF fpr: " << opt::fpr2 << "\n";
    cerr << "Primary BF bits: " << opt::m1 << "\t bits/kmer= " << opt::m1/opt::dbfSize << "\n";
    cerr << "Secondary BF bits: " << opt::m2 << "\t bits/kmer= " << opt::m2/opt::sbfSize << "\n";
    cerr << "Primary BF hashes: " << opt::nhash1 << "\n";
    cerr << "Secondary BF hashes: " << opt::nhash2 << "\n";

    BloomFilter dbFilter(opt::m1, opt::nhash1, opt::kmLen);
    BloomFilter sbFilter(opt::m2, opt::nhash2, opt::kmLen);
    for (unsigned file_i = 0; file_i < inFiles.size(); ++file_i) {
        std::ifstream in(inFiles[file_i].c_str());
        loadBFfq(in, dbFilter, sbFilter);
        in.close();
    }
    sbFilter.storeFilter("Bfilter.bf");

    ofstream bfinfo("Bfilter.inf");
    bfinfo << opt::m1 << "\n" << opt::nhash1 << "\n" <<opt::kmLen << "\n" << pow(1.0*dbFilter.getPop()/opt::m1,opt::nhash1) << "\n";
    bfinfo << opt::m2 << "\n" << opt::nhash2 << "\n" <<opt::kmLen << "\n" << pow(1.0*sbFilter.getPop()/opt::m2,opt::nhash2);
    bfinfo.close();

    cerr << "Primary BF actual fpr: " << setprecision(4) << fixed << pow(1.0*dbFilter.getPop()/opt::m1,opt::nhash1) << "\n";
    cerr << "Secondary BF actual fpr: " << setprecision(4) << fixed << pow(1.0*sbFilter.getPop()/opt::m2,opt::nhash2) << "\n";
    cerr << "Popcnt of pbf and sbf: " << dbFilter.getPop() << "\t" << sbFilter.getPop() << "\n";
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
        case 'k':
            arg >> opt::kmLen;
            break;
        case 'd':
            arg >> opt::fpr1;
            break;
        case 's':
            arg >> opt::fpr2;
            break;
        case 'r':
            opt::ref=true;
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

    getHist(opt::dbfSize, opt::sbfSize, opt::kmLen, opt::nThrd, inFiles);

    opt::m1 = -1*log(opt::fpr1)/log(2)/log(2)*opt::dbfSize;
    opt::m2 = -1*log(opt::fpr2)/log(2)/log(2)*opt::sbfSize;

    opt::nhash1 = opt::m1/opt::dbfSize*log(2);
    opt::nhash2 = opt::m2/opt::sbfSize*log(2);
    if(opt::nhash1==0)
        opt::nhash1 = 1;
    if(opt::nhash2==0)
        opt::nhash2 = 1;

    if(opt::ref)
        genBFref(inFiles);
    else
        genBFseq(inFiles);

    cerr << "time(sec): " <<setprecision(4) << fixed << omp_get_wtime() - sTime << "\n";
    return 0;
}
