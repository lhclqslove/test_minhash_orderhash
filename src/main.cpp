#include <iostream>
#include <map>
#include <set>

#include "ReadData.h"
#include "ReadFilter.h"
#include "Test_minhash.h"
#include "omp.h"
#include "Test_Orderminhash.h"
#include "OrderHashReadFilter.h"

int main(int argc, char **argv) {
    /*
     * 参数说明
     * 1：data.fasta 2: dataid 3:线程数 4:n 5：k 6:t 7:L
     */
    //读入比对paf文件，和fastq文件
    omp_set_num_threads(atoi(argv[3]));
    {
    std::string filename=argv[1];
//    std::string paf_file=argv[2];
    ReadData rD;
    rD.tempDir = "/home/lihc/projects/test_minhash_orderhash/cmake-build-release/tmp";
    rD.loadFromFile(filename.c_str(),ReadData::Filetype::FASTQ);

    //读入datasetid文件形成标准集合
    int *belong=new int[rD.getNumReads()+5];
    std::ifstream datasetid;
    datasetid.open(argv[2]);
    int setid=0;
    int readid=0;
    std::string id;
    std::vector<std::string>idtmp;
    std::vector<int>setsize;
    int set_size=0;
    int su=0;
    while(getline(datasetid,id))
    {
        if(id==">")
        {
            setid++;
            std::cout<<set_size<<" ";
            su+=set_size;
            setsize.push_back(set_size);
            set_size=0;
        }
        else
        {
            belong[readid++]=setid;
            set_size++;
        }
    }
    std::cout<<su<<"\n";
    //测试minhash
//    minhash::MinHashReadFilter filter;
//    filter.n=atoi(argv[4]);
//    filter.k=atoi(argv[5]);
//    filter.overlapSketchThreshold=atoi(argv[6]);
//    filter.tempDir="/home/lihc/projects/test_minhash_orderhash/cmake-build-release/tmp";
//    filter.initialize(rD);
//    Testminhash testminhash;
//    testminhash.init(filter,belong,setsize);
//    testminhash.test(atoi(argv[4]),atoi(argv[5]),atoi(argv[6]));


    //测试orderminhash
    orderhash::MinHashReadFilter filter1;
    filter1.n=atoi(argv[4]);
    filter1.k=atoi(argv[5]);
    filter1.overlapSketchThreshold=atoi(argv[6]);
    filter1.l=atoi(argv[7]);
    filter1.tempDir="/home/lihc/projects/test_minhash_orderhash/cmake-build-release/tmp";
    filter1.initialize(rD);
    Test_Orderhash testOrderhash;

    testOrderhash.init(filter1,belong,setsize);

    testOrderhash.test(atoi(argv[4]),atoi(argv[5]),atoi(argv[6]),atoi(argv[7]));
    return 0;
    }
}
