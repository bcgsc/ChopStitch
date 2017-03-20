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

#include "ntHashIterator.hpp"
#include "Uncompress.h"

#ifdef _OPENMP
# include <omp.h>
#endif

using namespace std;

namespace opt {
unsigned t=1;
unsigned k=64;
unsigned rBuck=4194304;
unsigned rBits=22	;
unsigned sBits=16;
size_t totKmer=0;
bool canon=true;
bool samH=true;
}

unsigned getftype(std::ifstream &in, std::string &samSeq) {
    std::string hseq;
    getline(in,hseq);
    if(hseq[0]=='>') {
        return 1;
    }
    if(hseq[0]=='@') {
        if( (hseq[1]=='H'&& hseq[2]=='D') ||
                (hseq[1]=='S'&& hseq[2]=='Q') ||
                (hseq[1]=='R'&& hseq[2]=='G') ||
                (hseq[1]=='P'&& hseq[2]=='G') ||
                (hseq[1]=='C'&& hseq[2]=='O') ) {
            return 2;
        }
        else
            return 0;
    }
    opt::samH=false;
    samSeq=hseq;
    return 2;
}

inline void ntComp(const uint64_t hVal, uint8_t *m_Counter) {
    if(hVal>>(63-opt::sBits) == 1) {
        size_t shVal=hVal&(opt::rBuck-1);
        if(m_Counter[shVal]<255)
            //#pragma omp atomic
            ++m_Counter[shVal];
    }
}

inline void ntRead(const string &seq, size_t &parkCount, uint8_t *m_Counter) {
    ntHashIterator itr(seq, 1, opt::k);
    while (itr != itr.end()) {
        ntComp((*itr)[0],m_Counter);
        ++itr;
        ++parkCount;
    }
}

void getEfq(size_t &parkCount, uint8_t *m_Counter, std::ifstream &in) {
    bool good = true;
    for(string seq, hseq; good;) {
        good = static_cast<bool>(getline(in, seq));
        good = static_cast<bool>(getline(in, hseq));
        good = static_cast<bool>(getline(in, hseq));
        if(good && seq.length()>=opt::k)
            ntRead(seq, parkCount, m_Counter);
        good = static_cast<bool>(getline(in, hseq));
    }
}

void getEfa(size_t &parkCount, uint8_t *m_Counter, std::ifstream &in) {
    bool good = true;
    for(string seq, hseq; good;) {
        string line;
        good = static_cast<bool>(getline(in, seq));
        while(good&&seq[0]!='>') {
            line+=seq;
            good = static_cast<bool>(getline(in, seq));
        }
        if(line.length()>=opt::k)
            ntRead(line, parkCount, m_Counter);
    }
}

void getEsm(size_t &parkCount, uint8_t *m_Counter, std::ifstream &in, const std::string &samSeq) {
    std::string samLine,seq;
    std::string s1,s2,s3,s4,s5,s6,s7,s8,s9,s11;
    if(opt::samH) {
        while(getline(in,samLine))
            if (samLine[0]!='@') break;
    }
    else
        samLine=samSeq;
    do {
        std::istringstream iss(samLine);
        iss>>s1>>s2>>s3>>s4>>s5>>s6>>s7>>s8>>s9>>seq>>s11;
        if(seq.length()>=opt::k)
            ntRead(seq, parkCount, m_Counter);
    } while(getline(in,samLine));
}

bool getHist(size_t &dbsize, size_t &sbsize, const unsigned klen, const unsigned tnum, const vector<string> &inFiles) {
    double sTime = omp_get_wtime();

    opt::rBuck = ((unsigned)1) << opt::rBits;

    uint8_t *t_Counter = new uint8_t [opt::rBuck];
    for(size_t i=0; i<opt::rBuck; i++) t_Counter[i]=0;

    opt::t = tnum;
    opt::k = klen;

#ifdef _OPENMP
    omp_set_num_threads(opt::t);
#endif

    #pragma omp parallel
    {
        size_t parkCount=0;

        uint8_t *m_Counter = new uint8_t [opt::rBuck];
        for(size_t i=0; i<opt::rBuck; i++) m_Counter[i]=0;

        #pragma omp for schedule(dynamic) nowait
        for (unsigned file_i = 0; file_i < inFiles.size(); ++file_i) {
            std::ifstream in(inFiles[inFiles.size()-file_i-1].c_str());
            std::string samSeq;
            unsigned ftype = getftype(in,samSeq);
            if(ftype==0)
                getEfq(parkCount, m_Counter, in);
            else if (ftype==1)
                getEfa(parkCount, m_Counter, in);
            else if (ftype==2)
                getEsm(parkCount, m_Counter, in, samSeq);
            in.close();
        }
        #pragma omp critical(cmrg)
        {
            for(size_t i=0; i<opt::rBuck; i++) {
                if(t_Counter[i]+m_Counter[i]<255)
                    t_Counter[i] += m_Counter[i];
            }
        }
        delete [] m_Counter;

        #pragma omp atomic
        opt::totKmer+=parkCount;
    }


    unsigned singlton=0,totton=0;
    unsigned x[256];
    for(size_t i=0; i<256; i++) x[i]=0;

    for(size_t i=0; i<opt::rBuck; i++) {
        ++x[t_Counter[i]];
        if(t_Counter[i]==1) ++singlton;
        if(t_Counter[i]) ++totton;
    }

    delete [] t_Counter;

    double f[256];
    for(size_t i=0; i<256; i++) f[i]=0;
    f[1]= -1.0*x[1]/(x[0]*(log(x[0])-opt::rBits*log(2)));
    for(size_t i=2; i<256; i++) {
        double sum=0.0;
        for(size_t j=1; j<i; j++)
            sum+=j*x[i-j]*f[j];
        f[i]=-1.0*x[i]/(x[0]*(log(x[0])-opt::rBits*log(2)))-sum/(i*x[0]);
    }
    double F0= (opt::rBits*log(2)-log(x[0])) * 1.0* ((size_t)1<<(opt::sBits+opt::rBits));

    dbsize = (long long)F0;
    sbsize = dbsize - (long long)(f[1]*F0);

    //cerr << "Total #k-mer:" << opt::totKmer << "\n";

    cerr << "time(sec): " <<setprecision(4) << fixed << omp_get_wtime() - sTime << "\n";
    return true;
}
