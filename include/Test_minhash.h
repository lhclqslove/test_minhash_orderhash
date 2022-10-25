//
// Created by USER on 2022/4/6.
//
#ifndef TEST_MINHASH_ORDERHASH_TEST_MINHASH_H
#define TEST_MINHASH_ORDERHASH_TEST_MINHASH_H
#include "ReadFilter.h"
#include "ReadData.h"
#include "string"
class Testminhash
{
    public:
//        ReadData *rD;
        minhash::MinHashReadFilter *minHashReadFilter;
        int *belong;
        std::vector<int> *setsize;
        void init(minhash::MinHashReadFilter &minHashReadFilter,int *belong,std::vector<int> &setsize)
        {
            this->minHashReadFilter=&minHashReadFilter;
            this->setsize=&setsize;
            this->belong=belong;
        }
        void test(int  n,int k,int t)
        {
            std::vector<int>sec[setsize->size()];
            std::vector<int>fail[setsize->size()];
            std::vector<double>sec_bi_sum;
            std::vector<double>fail_bi_sum;
            std::vector<double>sec_bi_setsize;

            int id_pos=0;
             int sec_count=0;
             int fail_count=0;
            std::string read;
            for(int i=0;i<setsize->size();i++)
            {
                 sec_count=0;
                 fail_count=0;
//                std::cout<<setsize[i]<< std::endl;
                for(int j=0;j<(*setsize)[i];j++)
                {
                    int id=id_pos++;
//                    std::cout<<id<<"id done here"<<std::endl;
                    minHashReadFilter->rD->getRead(id,read);
                    std::vector<read_t> result_ids;
                    int offset=minHashReadFilter->rD->avgReadLen/4;
                    std::set<read_t> res;

                    for(int st=0;st<read.size();st+=offset)
                    {
                        std::string sub_read=read.substr(st,offset);
                        minHashReadFilter->getFilteredReads(sub_read,result_ids);
                        for(auto  &x:result_ids)
                        {
                            res.insert(x);
                        }
                    }
                    for(auto x:res){
                        if(belong[x]==belong[id])sec_count++;
                        else fail_count++;
                    }
                }
//                std::cout<<i<<" done here"<<std::endl;
                sec[i].push_back(sec_count);
                fail[i].push_back(fail_count);
//                if(sec_count+fail_count<=0)
//                {
//                    printf("chu 0");
//                }
                sec_bi_sum.push_back(1.0*sec_count/(sec_count+fail_count));
                fail_bi_sum.push_back(1.0*fail_count/(sec_count+fail_count));
                sec_bi_setsize.push_back(1.0*sec_count/((*setsize)[i]*(*setsize)[i]));
            }
//            std::cout<<"yes done here"<<std::endl;
            //打印日志文件
            std::ofstream logfile;
            logfile.open("n"+std::to_string(n)+"k"+std::to_string(k)+"t"+std::to_string(t)+"logfile");

            for(int i=0;i<setsize->size();i++){
                logfile<<"sec count set:"<<i<<"\n";
                std::cout<<"sec count set:"<<i<<"\n";
                for(auto  &v:sec[i]) {
                    logfile<< v << " ";
                    std::cout<< v << " ";
                }
                logfile<<"\nfail count set:"<<i<<"\n";
                std::cout<<"\nfail count set:"<<i<<"\n";
                for(auto  &v:fail[i]) {
                    logfile << v << " ";
                    std::cout << v << " ";
                } logfile << "\n ";
                std::cout << "\n";
            }

            auto printlog=[&](std::vector<double> &info,std::string str){
                logfile<<str;
                std::cout<<str;
                double su=0;
                for(auto &v:info)
                {
                    su+=v;
                    std::cout<<v<<" ";
                    logfile<<v<<" ";
                }
                std::cout<<"\n";
                logfile<<"\n";
                logfile<<"average  rate:"<<su/info.size()<<"\n";
                std::cout<<"average  rate:"<<su/info.size()<<"\n";
            };
            printlog(sec_bi_sum,"sec_bi_sum\n");
            printlog(fail_bi_sum,"fail_bi_sum\n");
            printlog(sec_bi_setsize,"sec_bi_setsize\n");
        }

};
#endif //TEST_MINHASH_ORDERHASH_TEST_MINHASH_H
