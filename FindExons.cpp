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
#include <utility>
#include <algorithm>
#include <cassert>
#include <getopt.h>
#include <iostream>
#include <cstring>
#include "BloomFilter.hpp"
#include "pstream.h"
#include "Options.h"
#include "FastaConcat.h"


#define PROGRAM "FindExons"

static const char VERSION_MESSAGE[] =
    PROGRAM " Version 1.0.0 \n"
    "Written by Hamza Khan, Hamid Mohamadi.\n"
    "Copyright 2016 Canada's Michael Smith Genome Science Centre\n";

static const char USAGE_MESSAGE[] =
    "Usage: " PROGRAM " [OPTION]... FILES...\n"
    "Find putative exons in TransAbySS Transcriptome assembly file\n"
    "Accepatble file formats: fasta\n"
    "\n"
    " Options:\n"
    "\n"
    "  -i, --input-bloom=FILE     load bloom filter from FILE\n"
    "  -l, --leniency=N           leniency for exon-exon juction detection [10]\n"
    "  -f, --lfactor=N            leniency calculated as ceil(FPR*lfactor*k) \n"
    "  -s, --lsplicesignals=csv   Comma separated 5' splicesignals \n"
    "  -r, --rsplicesignals=csv   Comma separated 3' splicesignals \n"
    "      --allexons             Also output exons on either ends of contigs\n"
    "      --help	              display this help and exit\n"
    "      --version	          output version information and exit\n"
    "\n"
    "Report bugs to hmohamadi@bcgsc.ca; hkhan@bcgsc.ca\n";

using namespace std;

namespace opt {
string lsplicesignals="";
string rsplicesignals="";
unsigned leniency=10;
unsigned lfactor=0;
unsigned k=50;
unsigned nhash;
static string inputBloomPath;
size_t m;
bool internalexons = true;
int allExons=0;
double FPR=0.0;
}

static const char shortopts[] = "k:l:i:h:f:s:r:a:v:";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
    { "leniency",	required_argument, NULL, 'l' },
    { "input-bloom",	required_argument, NULL, 'i' },
    { "lfactor",	required_argument, NULL, 'f' },
    { "allexons",	no_argument, &opt::allExons, 1},
    { "lsplicesignals",	required_argument, NULL, 's' },
    { "rsplicesignals",	required_argument, NULL, 'r' },
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

void findJunctions(BloomFilter& bloom, int optind, char** argv) {
    FastaReader reader(argv[optind], FastaReader::FOLD_CASE);
    BloomFilter& bf = bloom;
    std::ofstream bfile("boundaries.csv");
    std::ofstream efile("exons.fa");
    std::ofstream cefile("confident_exons.fa");
    string read1, read2, final_Seq;

    for (FastaRecord rec; reader >> rec;) {
        if(int((rec.seq).length())<(int(opt::k)+5))
            continue;
        unsigned pos = (opt::k)-1,
                 myflag=0,
                 start = 0,
                 end = 0,
                 min_exon=5,
                 unmatch_count=0,
                 snp_chance_A=0,
                 snp_chance_T=0,
                 snp_chance_G=0,
                 snp_chance_C=0,
                 snp_nochance=0;

        int first_match=0, last_match_pos=0,
            big1=0, small1=0, big2=0;
        std::string current_sec = rec.seq;
        vector<int> arr;
        vector<int> carr;

        if(opt::leniency==10 && opt::lfactor!=0)
            opt::leniency = ceil(opt::lfactor*opt::FPR*opt::k);

        int str_length = current_sec.length();
        for (unsigned int x=0; x < (str_length-((opt::k)-1)-2); x++) {
            string it = current_sec.substr(x,(opt::k));
            pos=pos+1;
            if(bf.contains(it.c_str()))
                last_match_pos=pos;
            if(!(bf.contains(it.c_str())) && (myflag==0) && (first_match==1)) {
                myflag=1;
                start = pos;
                unmatch_count=1;
                snp_chance_A=0, snp_chance_C=0, snp_chance_G=0, snp_chance_T=0, snp_nochance=0;
            }

            if(!(bf.contains(it.c_str())) && (myflag==1)) {
                if(unmatch_count<=opt::leniency) {
                    string it_a = (it).replace(((it).length())-unmatch_count, 1, "A");
                    string it_t = (it).replace(((it).length())-unmatch_count, 1, "T");
                    string it_g = (it).replace(((it).length())-unmatch_count, 1, "G");
                    string it_c = (it).replace(((it).length())-unmatch_count, 1, "C");
                    if(bf.contains(it_a.c_str())) {
                        snp_chance_A++;
                        unmatch_count+=1;
                        goto SNP_found;
                    }
                    else if(bf.contains(it_t.c_str())) {
                        snp_chance_T++;
                        unmatch_count+=1;
                        goto SNP_found;
                    }
                    else if(bf.contains(it_g.c_str())) {
                        snp_chance_G++;
                        unmatch_count+=1;
                        goto SNP_found;
                    }

                    else if(bf.contains(it_c.c_str())) {
                        snp_chance_C++;
                        unmatch_count+=1;
                        goto SNP_found;
                    }
                    else {
                        snp_nochance++;
                    }
                    unmatch_count+=1;
                }
                if(unmatch_count == opt::leniency+1) {
                    if((snp_chance_A  >= 2*snp_nochance) && (snp_chance_A > (snp_chance_T+snp_chance_G+snp_chance_C)))
                        current_sec.replace(pos-(1+opt::leniency),1,"A");
                    if((snp_chance_T  >= 2*snp_nochance) && (snp_chance_T > (snp_chance_A+snp_chance_G+snp_chance_C)))
                        current_sec.replace(pos-(1+opt::leniency),1,"T");
                    if((snp_chance_G  >= 2*snp_nochance) && (snp_chance_G > (snp_chance_A+snp_chance_T+snp_chance_C)))
                        current_sec.replace(pos-(1+opt::leniency),1,"G");
                    if((snp_chance_C  >= 2*snp_nochance) && (snp_chance_C > (snp_chance_A+snp_chance_G+snp_chance_T)))
                        current_sec.replace(pos-(1+opt::leniency),1,"C");
                    unmatch_count+=1;
                }
            }
SNP_found:
            if((bf.contains(it.c_str())) && (myflag==1) && (first_match==1)) {
                unsigned int next_kmer = x;
                if (!(bf.contains((current_sec.substr (++next_kmer,(opt::k))).c_str())))
                    unmatch_count+=1;
            }
            unsigned int a = x;
            if((bf.contains(it.c_str())) && (myflag==1) &&  (bf.contains((current_sec.substr (++a,(opt::k))).c_str())) && (first_match==1)) {
                if(bf.contains((current_sec.substr (++a,(opt::k))).c_str())) {
                    myflag = 0;
                    end = pos;
                    if(end-start > ((opt::k)+min_exon)) {
                        int minus_counter = 0, snp_chance_A=0, snp_chance_C=0, snp_chance_G=0, snp_chance_T=0, snp_nochance=0;
                        for(int a_minus=(x-1); a_minus>=((int)(x-opt::leniency)); a_minus-- ) {
                            string it_rev = current_sec.substr (a_minus,(opt::k));
                            string it_a = (it_rev).replace((minus_counter), 1, "A");
                            string it_t = (it_rev).replace((minus_counter), 1, "T");
                            string it_g = (it_rev).replace((minus_counter), 1, "G");
                            string it_c = (it_rev).replace((minus_counter), 1, "C");
                            minus_counter = minus_counter + 1;

                            if(bf.contains(it_a.c_str()))
                                snp_chance_A++;
                            else if(bf.contains(it_t.c_str()))
                                snp_chance_T++;
                            else if(bf.contains(it_g.c_str()))
                                snp_chance_G++;
                            else if(bf.contains(it_c.c_str()))
                                snp_chance_C++;
                            else
                                snp_nochance++;
                        }

                        if(snp_chance_A  >= 2*snp_nochance || snp_chance_T  >= 2*snp_nochance || snp_chance_G  >= 2*snp_nochance || snp_chance_C  >= 2*snp_nochance ) {
                            bfile << "," << start-1 << "\n" << (rec.id) << ","<< start ;
                            arr.push_back(start-1);
                            arr.push_back(start);
                        }
                        else {
                            bfile << "," << start-1 << "\n" << (rec.id) <<","<< start << "," << (end-opt::k)-1<< "\n" << (rec.id) <<","<< end-opt::k ;
                            arr.push_back(start-1);
                            arr.push_back(start);
                            arr.push_back((end-opt::k)-1);
                            arr.push_back(end-opt::k);
                        }
                    }

                    if((end - start >= (((opt::k)-opt::leniency)-2)) && ((end-start)<= ((opt::k)+min_exon))) {
                        bfile << "," << start-1 << "\n" << (rec.id) <<","<< end-opt::k ;
                        arr.push_back(start-1);
                        arr.push_back(end-opt::k);
                    }
                }
            }

            //Testing the first matching kmer in the bloom filter
            if(bf.contains(it.c_str()) && first_match==0) {
                unsigned int next_it = x;
                //Added extra checks to avoid False positives
                if(bf.contains((current_sec.substr (++next_it,(opt::k))).c_str()) && bf.contains((current_sec.substr (++next_it,(opt::k))).c_str())) {
                    bfile << (rec.id)<< ",";
                    bfile << pos-((opt::k)-1);
                    arr.push_back(pos-((opt::k)-1));
                    first_match=1;
                }
            }
       }
        if(first_match==1) {
            bfile << "," << last_match_pos << "\n";
            arr.push_back(last_match_pos);
        }

        if(opt::internalexons) {
            int check_flag=0;
            carr=arr;
            if(!carr.empty() && carr.size()==2) {
                if (*(carr.begin())!=1 && *(carr.end()-1) < str_length-2) {
                    string single_read = rec.seq.substr(((*(carr.begin()))-1),((*(carr.end()-1)-(*(carr.begin()))+1)));
                    cefile << ">" << rec.id << "_" << *(carr.begin())<<"_"<< *(carr.end()-1)<< std::endl;
                    cefile << single_read<<endl;
                    std::flush(cefile);
                }
            }
            if(!carr.empty() && (carr.size())>2) {
                int beGIN = *(carr.begin());
                int enD = *(carr.end()-1);
                check_flag=1;
                if(beGIN == 1)
                    carr.erase(carr.begin(),carr.begin()+2);
                if(enD >= str_length-2)
                    carr.erase (carr.end()-2,carr.end());
            }
            if(!carr.empty() && carr.size()==2 && check_flag==1) {
                string single_read = rec.seq.substr(((*(carr.begin()))-1),((*(carr.end()-1)-(*(carr.begin()))+1)));
                cefile << ">" << rec.id << "_" << *(carr.begin())<<"_"<< *(carr.end()-1)<< std::endl;
                cefile << single_read<<endl;
                std::flush(cefile);
                check_flag=0;
            }
            if(!carr.empty() && (carr.size())>2) {
                int beGIN = *(carr.begin());
                int enD = *(carr.end()-1);
                int Big1=0, Small1=0, Big2=0;
                string Read1, Read2;
                for (std::vector<int>::const_iterator w = (carr.begin())+1; w != (carr.end())-1; w+=2) {
                    if(*w < *(w+1)) {
                        Big1 = *(w+1);
                        Small1 = *w;
                    }
                    else {
                        Small1 = *(w+1);
                        Big1 = *w;
                    }
                    if(w == (carr.begin())+1) {
                        Read1 = rec.seq.substr(beGIN-1,(Small1-beGIN)+1);
                        cefile << ">" << rec.id << "_" << beGIN<<"_"<< Small1<< std::endl;
                        cefile << Read1<<endl;
                        std::flush(cefile);
                        Big2 = Big1;
                    }
                    else {
                        Read2 = rec.seq.substr(Big2-1,(Small1-Big2)+1);
                        cefile << ">" << rec.id << "_" << Big2<<"_"<< Small1<< std::endl;
                        cefile << Read2<<endl;
                        std::flush(cefile);
                        Big2 = Big1;
                    }
                }
                if(!carr.empty() && (carr.size())>2) {
                    string End_Read = rec.seq.substr(Big2-1,(enD-Big2)+1);
                    cefile << ">" << rec.id << "_" << Big2<<"_"<< enD<< std::endl;
                    cefile << End_Read<<endl;
                    std::flush(cefile);
                }
            }
        }
        if(!arr.empty() && arr.size()==2) {
            string single_read = rec.seq.substr(((*(arr.begin()))-1),((*(arr.end()-1)-(*(arr.begin()))+1)));
            efile << ">" << rec.id << "_" << *(arr.begin())<<"_"<< *(arr.end()-1)<< std::endl;
            efile << single_read<<endl;
            std::flush(efile);
        }
        if(!arr.empty() && (arr.size())>2) {
            int bEgin = *(arr.begin());
            int eNd = *(arr.end()-1);
            for (std::vector<int>::const_iterator q = (arr.begin())+1; q != (arr.end())-1; q+=2) {
                if(*q < *(q+1)) {
                    big1 = *(q+1);
                    small1 = *q;
                }
                else {
                    small1 = *(q+1);
                    big1 = *q;
                }
                if(q == (arr.begin())+1) {
                    read1 = rec.seq.substr(bEgin-1,(small1-bEgin)+1);
                    efile << ">" << rec.id << "_" << bEgin<<"_"<< small1<< std::endl;
                    efile << read1<<endl;
                    std::flush(efile);
                    big2 = big1;
                }
                else {
                    read2 = rec.seq.substr(big2-1,(small1-big2)+1);
                    efile << ">" << rec.id << "_" << big2<<"_"<< small1<< std::endl;
                    efile << read2<<endl;
                    std::flush(efile);
                    big2 = big1;
                }
            }
            if(!arr.empty() && (arr.size())>2) {
                string end_read = rec.seq.substr(big2-1,(eNd-big2)+1);
                efile << ">" << rec.id << "_" << big2<<"_"<< eNd<< std::endl;
                efile << end_read<<endl;
                std::flush(efile);
            }
        }
    }
    efile.close();
}

string check_modified_kmer(string it, BloomFilter& bf, unsigned int &increments){
  string it_a = "A"+it;
  string it_t = "T"+it;
  string it_g = "G"+it;
  string it_c = "C"+it;
  if(bf.contains(it_a.c_str())){
    increments+=1;
    return it_a;
  }
  if(bf.contains(it_g.c_str())){
     increments+=1;
     return it_g;
  }
  if(bf.contains(it_t.c_str())){
    increments+=1;
    return it_t;
  }

  if(bf.contains(it_c.c_str())){
    increments+=1;
    return it_c;
  }
  else
    return it;
}


string check_modified_kmer_last(string it, BloomFilter& bf, unsigned int &increments){
  string it_a = it+"A";
  string it_t = it+"T";
  string it_g = it+"G";
  string it_c = it+"C";
  if(bf.contains(it_a.c_str())){
    increments+=1;
    return it_a;
  }
  if(bf.contains(it_g.c_str())){
     increments+=1;
     return it_g;
  }
  if(bf.contains(it_t.c_str())){
    increments+=1;
    return it_t;
  }

  if(bf.contains(it_c.c_str())){
    increments+=1;
    return it_c;
  }
  else
    return it;
}


string PostProcess_end(BloomFilter& bloom, string& current_sec, vector<string> rss) {
     BloomFilter& bf = bloom;
     string last_kmer = current_sec.substr(current_sec.length()-(opt::k),(opt::k));
     //cout << "last_kmer = " << last_kmer << endl;
     unsigned int increments = 0;
     for (unsigned int x=0; x<6; x++) {
         last_kmer = check_modified_kmer_last(last_kmer, bf, increments);
         last_kmer = last_kmer.substr((last_kmer.length()-(opt::k)),(opt::k));
         //cout << "last_kmer-after = " << last_kmer << endl;
         //std::cout << "increments = " << increments << std::endl;
   }
     if(increments>1) {
       for (unsigned int y=0; y<increments-1; y++){
         //std::cout << "dimer=" << last_kmer.substr(((opt::k)+y-increments),2) << std::endl;
         bool testcond = false;
         string tail_plus_rec, tail;
         for(unsigned int i = 0; i < rss.size(); ++i){
           if(last_kmer.substr(((opt::k)+y-increments),2) == rss[i]){
             //std::cout << "End Dimer present" << std::endl;
             tail = last_kmer.substr(opt::k-increments, y);
             tail_plus_rec =current_sec.substr(current_sec.length()-(opt::k),(opt::k))+tail;
             //std::cout << "tail=" << tail<< "\ntail_plus_rec=" << tail_plus_rec << std::endl;
             testcond=true;
             break;
           }
         }
         if(testcond==true){
           current_sec = current_sec+tail;
           //std::cout << "current_sec_after_testcond_true = " << current_sec << std::endl;
           break;
         }
       }
      }
      return current_sec;
    }



void PostProcess_begin(BloomFilter& bloom, vector<string> lss, vector<string> rss ) {
    FastaReader reader("confident_exons.fa", FastaReader::FOLD_CASE);
    std::ofstream pfile("processed_exons.fa");
    BloomFilter& bf = bloom;
    int count = 0;
    for (FastaRecord rec; reader >> rec;) {
        //cout << rec.seq << "\n";
        std::string current_sec = rec.seq;
        string processed_seq;
        if(current_sec.length()<(opt::k)){
          pfile << ">" << rec.id << "\n" << rec.seq << endl;
          continue;
        }

        string first_kmer = current_sec.substr(0,(opt::k));
        //cout << "First Kmer = " << first_kmer << endl;
        unsigned int increments = 0;
        for (unsigned int x=0; x<6; x++) {
            first_kmer = check_modified_kmer(first_kmer, bf, increments);
            first_kmer = first_kmer.substr(0,(opt::k));
            //cout << "Kmer-after = " << first_kmer << endl;
            //std::cout << "increments" << increments << std::endl;
      }
        if(increments>1) {
          for (unsigned int y=increments; y>1; y--){
            //std::cout << "dimer=" << first_kmer.substr(y-2,2) << std::endl;
            bool testcond = false;
            string head_plus_rec;
            for(unsigned int i = 0; i < lss.size(); ++i){
              if(first_kmer.substr(y-2,2) == lss[i]){
                //std::cout << "Dimer present" << std::endl;
                count+=1;
                string head = first_kmer.substr(y,(increments-y));
                head_plus_rec = head + current_sec.substr(0,(opt::k));
                //std::cout << "head=" << head << "\nhead_plus_rec=" << head_plus_rec << std::endl;
                testcond=true;
                break;
            }
          }
            if(testcond==true){
                current_sec = head_plus_rec + current_sec.substr((opt::k), (current_sec.length()-(opt::k)));
                //std::cout << "current_sec_after_testcond_true = " << current_sec << std::endl;
                break;
              }
          }
      }
        processed_seq = PostProcess_end(bloom, current_sec, rss);
        pfile << ">" << rec.id << "\n" << processed_seq << endl;
      }
      //std::cout << "Count = " << count << std::endl;
    }

template<typename Out>
void split(const std::string &s, char delim, Out result){
    std::stringstream ss;
    ss.str(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        *(result++) = item;
    }
}
std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, std::back_inserter(elems));
    return elems;
}

int main(int argc, char** argv) {

    bool die = false;
    for (int c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) {
        case '?':
            die = true;
            break;
        case 'l':
            arg >> opt::leniency;
            break;
        case 'f':
            arg >> opt::lfactor;
            break;
        case 'i':
            arg >> opt::inputBloomPath;
            break;
        case 's':
            arg >> opt::lsplicesignals;
            break;
        case 'r':
            arg >> opt::rsplicesignals;
            break;
        case 'a':
            opt::allExons=1;
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
    if (argc - optind < 0) {
        std::cerr << PROGRAM ": missing arguments\n";
        die = true;
    }
    if (die) {
        std::cerr << "Try `" << PROGRAM << " --help' for more information.\n";
        exit(EXIT_FAILURE);
    }

    std::vector<std::string> lss = split(opt::lsplicesignals, ',');
    std::vector<std::string> rss = split(opt::rsplicesignals, ',');
    //std::cout << "lsplicesignals = " << std::endl;
    //for(unsigned int i = 0; i < lss.size(); ++i)
    //cout << lss[i] << endl;

    //std::cout << "rsplicesignals = " << std::endl;
    //for(unsigned int i = 0; i < rss.size(); ++i)
    //cout << rss[i] << endl;

    if (!opt::inputBloomPath.empty()) {

        if (opt::verbose)
            std::cerr << "Loading Bloom filter from `"
                      << opt::inputBloomPath << "'...\n";

        string infoPath;
        size_t pos = opt::inputBloomPath.rfind("/");
        if (pos == std::string::npos)
            infoPath = "Bfilter.inf";
        else {
            infoPath=opt::inputBloomPath.substr(0, pos) + "/Bfilter.inf";
        }

        cerr << infoPath << "\n";

        ifstream bfinfo(infoPath.c_str());
        string bfLine;
        while(bfinfo >> opt::m >> opt::nhash >> opt::k >> opt::FPR) {}
        cerr << opt::m << "\t" << opt::nhash<< "\t" <<opt::k <<"\t" << opt::FPR << "\n";
        bfinfo.close();

        BloomFilter bloom(opt::m, opt::nhash, opt::k, opt::inputBloomPath.c_str());
        cerr << bloom.getPop() << "\n";

        findJunctions(bloom, optind, argv);
        if(!opt::allExons)
            remove("exons.fa");
        
        remove("boundaries.csv");
   
        if(opt::lsplicesignals!=""||opt::rsplicesignals!=""){
             PostProcess_begin(bloom, lss, rss);
             remove("confident_exons.fa");
        }
    }
    else
        cerr << "No input file\n";


}
