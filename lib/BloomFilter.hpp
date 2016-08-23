/*
 *
 * BloomFilter.hpp
 * Author: Hamid Mohamadi
 * Genome Sciences Centre,
 * British Columbia Cancer Agency
 */


#ifndef BLOOMFILTER_H_
#define BLOOMFILTER_H_
#include <string>
#include <stdint.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <sys/stat.h>
#include <cstring>
#include <cassert>
#include <cstdlib>
#include <stdio.h>
#include <cstring>
#include "nthash.hpp"

using namespace std;

inline unsigned popCnt(unsigned char x) {
    return ((0x876543210 >>
             (((0x4332322132212110 >> ((x & 0xF) << 2)) & 0xF) << 2)) >>
            ((0x4332322132212110 >> (((x & 0xF0) >> 2)) & 0xF) << 2))
           & 0xf;
}

class BloomFilter {
public:

    BloomFilter(size_t filterSize, unsigned hashNum, unsigned kmerSize, const char * fPath):
        m_size(filterSize), m_hashNum(hashNum), m_kmerSize(kmerSize) {
        m_filter = new unsigned char [(m_size + 7)/8];
        ifstream myFile(fPath, ios::in | ios::binary);
        myFile.seekg (0, ios::beg);
        myFile.read ((char *)m_filter, (m_size + 7)/8);
        myFile.close();
    }

    BloomFilter(size_t filterSize, unsigned hashNum, unsigned kmerSize):
        m_size(filterSize), m_hashNum(hashNum), m_kmerSize(kmerSize) {
        m_filter = new unsigned char [(m_size + 7)/8];
        for(size_t i = 0; i < (m_size + 7)/8; i++)
            m_filter[i]=0;
    }

    void insert(const uint64_t hVal) {
        for (unsigned i = 0; i < m_hashNum; i++) {
            uint64_t mhVal = hVal * (i ^ m_kmerSize * multiSeed);
            mhVal ^= mhVal >> multiShift;
            size_t hLoc = mhVal % m_size;
            __sync_or_and_fetch(&m_filter[hLoc / 8], (1 << (7 - hLoc % 8)));
        }
    }

    void insert(const char* kmer) {
        uint64_t hVal = NTC64(kmer, m_kmerSize);
        for (unsigned i = 0; i < m_hashNum; i++) {
            uint64_t mhVal = hVal * (i ^ m_kmerSize * multiSeed);
            mhVal ^= mhVal >> multiShift;
            size_t hLoc = mhVal % m_size;
            __sync_or_and_fetch(&m_filter[hLoc / 8], (1 << (7 - hLoc % 8)));
        }
    }

    void insert(const char * kmer, uint64_t& fhVal, uint64_t& rhVal) {
        fhVal = getFhval(kmer, m_kmerSize);
        rhVal = getRhval(kmer, m_kmerSize);
        uint64_t hVal = fhVal ^ rhVal;
        for (unsigned i = 0; i < m_hashNum; i++) {
            uint64_t mhVal = hVal * (i ^ m_kmerSize * multiSeed);
            mhVal ^= mhVal >> multiShift;
            size_t hLoc = mhVal % m_size;
            __sync_or_and_fetch(&m_filter[hLoc / 8], (1 << (7 - hLoc % 8)));
        }
    }

    void insert(uint64_t& fhVal, uint64_t& rhVal, const char charOut, const char charIn) {
        fhVal = rol(fhVal, 1)^rol(seedTab[charOut], m_kmerSize)^seedTab[charIn];
        rhVal = ror(rhVal, 1)^ ror(seedTab[(unsigned char)charOut + cpOff], 1)^ rol(seedTab[(unsigned char)charIn+cpOff],m_kmerSize-1);
        uint64_t hVal = fhVal ^ rhVal;
        for (unsigned i = 0; i < m_hashNum; i++) {
            uint64_t mhVal = hVal * (i ^ m_kmerSize * multiSeed);
            mhVal ^= mhVal >> multiShift;
            size_t hLoc = mhVal % m_size;
            __sync_or_and_fetch(&m_filter[hLoc / 8], (1 << (7 - hLoc % 8)));
        }
    }

	bool contains(const uint64_t hVal) const {
        for (unsigned i = 0; i < m_hashNum; i++) {
            uint64_t mhVal = hVal * (i ^ m_kmerSize * multiSeed);
            mhVal ^= mhVal >> multiShift;
            size_t hLoc = mhVal % m_size;
            if ((m_filter[hLoc / 8] & (1 << (7 - hLoc % 8))) == 0)
                return false;
        }
        return true;
    }

    bool contains(const char* kmer) const {
        uint64_t hVal = getFhval(kmer, m_kmerSize) ^ getRhval(kmer, m_kmerSize);
        for (unsigned i = 0; i < m_hashNum; i++) {
            uint64_t mhVal = hVal * (i ^ m_kmerSize * multiSeed);
            mhVal ^= mhVal >> multiShift;
            size_t hLoc = mhVal % m_size;
            if ((m_filter[hLoc / 8] & (1 << (7 - hLoc % 8))) == 0)
                return false;
        }
        return true;
    }

    bool contains(const char * kmer, uint64_t& fhVal, uint64_t& rhVal) {
        fhVal = getFhval(kmer, m_kmerSize);
        rhVal = getRhval(kmer, m_kmerSize);
        uint64_t hVal = fhVal ^ rhVal;
        for (unsigned i = 0; i < m_hashNum; i++) {
            uint64_t mhVal = hVal * (i ^ m_kmerSize * multiSeed);
            mhVal ^= mhVal >> multiShift;
            size_t hLoc = mhVal % m_size;
            if ((m_filter[hLoc / 8] & (1 << (7 - hLoc % 8))) == 0)
                return false;
        }
        return true;
    }

    bool contains(uint64_t& fhVal, uint64_t& rhVal, const char charOut, const char charIn) {
        fhVal = rol(fhVal, 1)^ rol(seedTab[charOut], m_kmerSize)^ seedTab[charIn];
        rhVal = ror(rhVal, 1)^ ror(seedTab[(charOut + cpOff)], 1)^ rol(seedTab[(charIn + cpOff)],m_kmerSize - 1);
        uint64_t hVal = (rhVal < fhVal) ? rhVal : fhVal;
        for (unsigned i = 0; i < m_hashNum; i++) {
            uint64_t mhVal = hVal * (i ^ m_kmerSize * multiSeed);
            mhVal ^= mhVal >> multiShift;
            size_t hLoc = mhVal % m_size;
            if ((m_filter[hLoc / 8] & (1 << (7 - hLoc % 8))) == 0)
                return false;
        }
        return true;
    }    
    
    void storeFilter(const char * fPath) const {
        ofstream myFile(fPath, ios::out | ios::binary);
        myFile.write(reinterpret_cast<char*>(m_filter), (m_size + 7)/8);
        myFile.close();
    }

    size_t getPop() const {
        size_t i, popBF=0;
        #pragma omp parallel for reduction(+:popBF)
        for(i=0; i<(m_size + 7)/8; i++)
            popBF = popBF + popCnt(m_filter[i]);
        return popBF;
    }

    unsigned getHashNum() const {
        return m_hashNum;
    }

    unsigned getKmerSize() const {
        return m_kmerSize;
    }

    ~BloomFilter() {
        delete[] m_filter;
    }

private:
    BloomFilter(const BloomFilter& that); //to prevent copy construction
    unsigned char * m_filter;
    size_t m_size;
    unsigned m_hashNum;
    unsigned m_kmerSize;
};

#endif /* BLOOMFILTER_H_ */