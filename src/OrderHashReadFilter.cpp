#include "OrderHashReadFilter.h"
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <random>
#include <set>

orderhash::ReadFilter::~ReadFilter() {}

void orderhash::MinHashReadFilter::initialize(ReadData &rD) {
    this->rD = &rD;
    numReads = rD.getNumReads();
    readPos = &rD.getReadPos();
    readPosSorted.clear();
    for (read_t r = 0; r < numReads; ++r) {
        readPosSorted.push_back(std::make_pair((*readPos)[r], r));
    }
    std::sort(readPosSorted.begin(), readPosSorted.end());

    std::vector<kMer_t> sketches(n * numReads);

    generateRandomNumbers(n);
//    for(int i=0;i<n;i++)
//    {
//        std::cout<<randNumbers[i]<<" ";
//    }std::cout<<std::endl;

    size_t maxNumkMers;
    if (rD.maxReadLen < k-1)
        maxNumkMers = 0;
    else
        maxNumkMers = rD.maxReadLen - k + 1;
//    std::cout<<"l="<<l<<std::endl;
//    std::string readStr1,readStr2;
//    rD.getRead(0, readStr1);
//    rD.getRead(1, readStr2);
//    std::vector<mer_info> kMersVec1(maxNumkMers),kMersVec2(maxNumkMers);
//    std::vector<kMer_t> hashesVec1(n),hashesVec2(n);
//    std::cout<<readStr1<<std::endl;
//    std::cout<<readStr2<<std::endl;
//    string2Sketch(readStr1, sketches.data() , kMersVec1, hashesVec1);
//    std::cout<<"second"<<std::endl;
//    string2Sketch(readStr2, sketches.data()+n , kMersVec2, hashesVec2);
//    for(int i=0;i<2*n;i++)
//    {
//        std::cout<<sketches[i]<<"*";
//        if(i==n-1)std::cout<<std::endl;
//    }
//    std::cout<<std::endl;

#pragma omp parallel
    {
        // We define these vectors here rather than allocate inside string2Sketch
        // to avoid thread contention during repeated allocation and deallocation.
        // Note that memory allocation typically leads to waits when multiple threads
        // do it at the same time.
        std::vector<mer_info> kMersVec(maxNumkMers);
        std::vector<kMer_t> hashesVec(n);
        std::string readStr;
#pragma omp for
        for (read_t i = 0; i < numReads; ++i){
            rD.getRead(i, readStr);
            string2Sketch(readStr, sketches.data() + i * n, kMersVec, hashesVec);
        }
    } // pragma omp parallel

    populateHashTables(sketches);
}

void orderhash::MinHashReadFilter::generateRandomNumbers(size_t n) {
    if (randNumbers)
        delete[] randNumbers;
    randNumbers = new kMer_t[n];

    std::random_device rd;
    std::mt19937_64 gen(rd());

    /* This is where you define the number generator for unsigned long long: */
    std::uniform_int_distribution<unsigned long long> dis;

    for (size_t i = 0; i < n; ++i) {
        randNumbers[i] = dis(gen);
    }
}

void orderhash::MinHashReadFilter::getFilteredReads(kMer_t sketch[],
                                         std::vector<read_t> &results) {
    std::vector<read_t> matches;
    results.clear();
    for (size_t sketchIndex = 0; sketchIndex < n; ++sketchIndex) {
        kMer_t curHash = sketch[sketchIndex];
        hashTables[sketchIndex].pushMatchesInVector(curHash,matches);
    }
    std::sort(matches.begin(), matches.end());
    auto end = matches.end();
    auto next = matches.begin();
    for (auto it = matches.begin(); it != end; it = next) {
        next = std::upper_bound(it, end, *it);

        if (next - it >= (ssize_t)overlapSketchThreshold) {
            results.push_back(*it);
        }
    }
}

void orderhash::MinHashReadFilter::getFilteredReads(const std::string &s,
                                         std::vector<read_t> &results) {
    results.clear();
    std::vector<kMer_t> sketch(n);
    size_t numKmers;
    if (s.size() <= k - 1) {
        numKmers = 0;
        return;
    }
    else
        numKmers = s.size() - k + 1;
    std::vector<mer_info> kMersVec(numKmers); // preallocation
    std::vector<kMer_t>hashesVec(n);
    string2Sketch(s, sketch.data(), kMersVec, hashesVec);
    getFilteredReads(sketch.data(), results);
}

orderhash::MinHashReadFilter::MinHashReadFilter() {}

kMer_t orderhash::MinHashReadFilter::kMerToInt(const std::string &s) {
    size_t l = s.length();
    kMer_t result = 0;
    for (size_t i = 0; i < l; ++i) {
        result <<= 2;
        result |= baseToInt(s[i]);
    }
    return result;
}

// Using the bit operations version of this function provides a 13X improvement
// in speed
char orderhash::MinHashReadFilter::baseToInt(const char base) {
    return (base & 0b10) | ((base & 0b100) >> 2);
}

void orderhash::MinHashReadFilter::string2Sketch(const std::string &s, kMer_t *sketch, std::vector<mer_info> &kMers, std::vector<kMer_t> &hashes) {

    ssize_t numKMers = s.length() - k + 1;
    if (numKMers < 0)
        return;
//    std::cout<<s.substr(0,5)<<" "<<"yes"<<std::endl;
    string2KMers(s, k, kMers);
    //为每个hash函数挑选前l个最小的mer_info
//    std::cout<<"get kmer_info"<<std::endl;
    xxhash hash;
    for(size_t i=0;i<n;i++) {
        //通过xxhash计算hash值
        for (size_t j = 0; j < numKMers; j++) {
            //初始化随机种子
            hash.reset(randNumbers[i]);
            //拆入kmer和occ
            hash.update(&kMers[j].kmer, sizeof(kMers[j].kmer));
            hash.update(&kMers[j].occ, sizeof(kMers[j].occ));
            //计算hash值
            kMers[j].hash = hash.digest();
        }

        //按hash值排序
        std::partial_sort(kMers.begin(), kMers.begin() + l, kMers.begin() + numKMers,
                          [&](const mer_info &x, const mer_info &y) { return x.hash < y.hash; });

        //前l个按pos排序
        std::sort(kMers.begin(), kMers.begin() + l,
                  [&](const mer_info &x, const mer_info &y) { return x.pos < y.pos; });
        //把前l个mer_info哈希成一个值，等效于minhash中的某一个最小hash值
        hash.reset(randNumbers[i]);
        for (size_t j = 0; j < l; j++) {
            hash.update(&kMers[j].kmer, sizeof(kMers[j].kmer));
            hash.update(&kMers[j].occ, sizeof(kMers[j].occ));
//            std::cout<<kMers[j].kmer<<" "<<kMers[j].occ<<std::endl;
        }
        sketch[i] = hash.digest();
    }
}

void orderhash::MinHashReadFilter::hashKMer(kMer_t kMer, std::vector<kMer_t> &hashes) {
    for (size_t l = 0; l < n; l++)
        hashes[l] = hasher(kMer^randNumbers[l]);
}

//void MinHashReadFilter::string2KMers(const std::string &s, const size_t k,
//                                     std::vector<kMer_t> &kMers) {
//    ssize_t maxI = s.length() - k + 1;
//    if (maxI <= 0)
//        return;
//    kMer_t currentKMer = kMerToInt(s.substr(0, k));
//    kMers[0] = currentKMer;
//    const unsigned long long mask = (1ull << (2 * k)) - 1;
//    for (size_t i = 1; i < (size_t)maxI; ++i) {
//        currentKMer =
//            ((currentKMer << 2) | MinHashReadFilter::baseToInt(s[i + k - 1])) &
//            mask;
//        kMers[i] = currentKMer;
//    }
//}
void orderhash::MinHashReadFilter::string2KMers(const std::string &s, const size_t k, std::vector<mer_info> &KMers_info) {
    ssize_t maxI = s.length() - k + 1;
    if (maxI <= 0)
        return;
    std::unordered_map<kMer_t, unsigned> occurrences;
    kMer_t currentKMer = kMerToInt(s.substr(0, k));
    KMers_info[0].init(0,0,currentKMer);
    occurrences[currentKMer]=1;
    const unsigned long long mask = (1ull << (2 * k)) - 1;

    for (size_t i = 1; i < (size_t)maxI; ++i) {
        currentKMer =
                ((currentKMer << 2) | MinHashReadFilter::baseToInt(s[i + k - 1])) &
                mask;
        auto occ = occurrences[currentKMer]++;
        KMers_info[i].init(i,occ,currentKMer);
    }
}

orderhash::MinHashReadFilter::~MinHashReadFilter() {
    if (randNumbers)
        delete[] randNumbers;
}

void orderhash::MinHashReadFilter::populateHashTables(const std::vector<kMer_t> &sketches) {
    std::cout << "Starting to populate hash tables" << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    hashTables.resize(n);
#pragma omp parallel for
    for (size_t i = 0; i < n; ++i)
        hashTables[i].initialize(n,numReads,sketches.data(),i,tempDir);

    auto end = std::chrono::high_resolution_clock::now();
    auto duration =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "finished populating hash tables" << std::endl;
    std::cout << duration.count() << " milliseconds passed" << std::endl;
}
