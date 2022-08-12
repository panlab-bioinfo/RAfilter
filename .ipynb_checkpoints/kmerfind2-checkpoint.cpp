# include <cstdio>
# include <cstdlib>
# include <iostream>
# include <algorithm>
# include "parallel_hashmap/phmap.h"
# include <unordered_set>
# include <string>
# include <vector>
# include <fstream>
#include <mutex>
# include "thread_pool.hpp"

using namespace  std;

vector<phmap::flat_hash_set<uint32_t>> kmer; //存储kmer的hash>
mutex m;
// vector<unordered_set<uint32_t>> kmer;

void print_bin(int n){
    int l = sizeof(n)*8;//总位数。
    int i;
    if(i == 0){
         printf("0");
         return;
     }
    for(i = l-1; i >= 0; i --){
        if(n&(1<<i)) break;
    }
    for(;i>=0; i --)
        printf("%d", (n&(1<<i)) != 0);
    printf("\n");
}
// 读取k-mer为二进制整数
char ctoi(char c){
    switch (c)
    {
        case 'A': return 0;
        case 'T': return 3;
        case 'C': return 1;
        case 'G': return 2;
        default:
            return -1;
    }
}

// 
uint64_t tobin(string *s){
    uint64_t k_value=0;
    int i=0;
    for(auto c=s->begin();c<s->end();c++){
        // cout<<*c;
        int v=ctoi(*c);
        if (v!=-1){
            k_value = k_value | (v<<i);
            i+=2;
        }
    }
    // cout<<endl;
    return k_value;
}
//reverse the bin of kmer
uint64_t reversebin(uint64_t x){
    uint64_t mid = x & 0x300000;
    x = (x & 0x3ffffc00000) >> 22 | (x & 0xfffff) << 22;
    x = (x & 0x3ff000ffc00) >> 10 | (x & 0xffc003ff) << 10;
    unit64_t quarter = x & 0x300c00c030;
    x = (x & 0x3C0F00F03C0) >> 6 | (x & 0xF03C03C0F)<<6;
    x = (x & 0x30CC30C330C) >> 2 | (x & 0xC330C30CC3)<<2;
    x = x | mid | quarter;
    return x;
}


int creat_dict(const char *k_path,uint64_t step_len, int idx){
    ifstream k_mer(k_path);
    phmap::flat_hash_set<uint32_t> dict;
    int i=0;
    int v;
    uint64_t k_value=0;
    // printf("%s",s);
    char *tmp=new char[step_len+2];
    *tmp={'\0'};
    k_mer.seekg(idx*(step_len),k_mer.beg);
    k_mer.read(tmp,step_len);
    for(char *p=tmp;p<tmp+step_len;p++){
        v=ctoi(*p);
        if (v!=-1){
            k_value = k_value | (v<<i);
            i+=2;
            // printf(" [%c,%d] ",*p,v);
        }
        else if(v == -1 && k_value != 0){
            dict.insert(k_value);
            // printf("%u\n",k_value);
            // cout<<tmp<<","<<endl;
            k_value=0;
            i=0;
        }
    }
    m.lock();
    kmer.push_back(dict);
    m.unlock();
    delete tmp;
    k_mer.close();
    return 1;
}

int read_kmer(const char *k_path,int t){
    /*
    将文件分割为与线程数一样多的字符串块；
    每块文件获得一个线程，对其进行二进制转换并存入unorderset集合
    */
    ThreadPool pool(t);
    pool.init();
    ifstream k_mer(k_path);

    // 分块
    k_mer.seekg(0,k_mer.end);
    uint64_t total_len=k_mer.tellg();
    k_mer.seekg(0,k_mer.beg);
    string s;
    getline(k_mer, s);
    uint64_t line_len = s.size();
    k_mer.close();
    uint64_t row_cnt=total_len/line_len/t;  //每块的长度
    if(row_cnt*t*line_len<total_len){
        row_cnt++;
    }
    // printf("%lu\n%lu",total_len,step_len);
    // string tmp;
    
    for(int i=0;i<t;i++){
        // unordered_set<uint32_t> d;
        // kmer.push_back(d);
        // printf(tmp);
        //creat_dict(tmp,step_len);
        pool.submit(creat_dict,k_path, row_cnt*line_len, i);
    }
    pool.shutdown();
    return 1;
}


/* int findindex(string seq,string name_path){
}
*/ 

int main(int argc,char* argv[]){
    // const char *kmerfile="../kmer.sample.dump";

    const char *kmerfile="/data/panweihua/human_chm13/hifi/fastqs/chm13.asm.p_ctg.compressed.kmers.1.dump";
    string fafile="chm13.asm.p_ctg.compressed.kmers.1.dump";
    read_kmer(kmerfile,1);
    ifstream fa(fafile);
    string line,name,qkmer;
    uint64_t len,ikmer,rkmer;
    while (getline(fa,line)){
        if(line[0] == '>'){
            name = line.substr(1,20);
        }
        else{
            len = line.length();
            for(int i = 0;i < len - 20;i++){
                qkmer = line.substr(i,21);
                ikmer=tobin(&(line.substr(i,21)));
                rkmer=
                for(auto d:kmer){
                    if(d.count()){
                        cout<<name<<"\t"<<len<<"\t"<<i<<"\t"<<qkmer<<"\t+"<<endl;
                        break;
                    }
                
                    else if(hashmap.find(kmer_rc) != hashmap.end()){
                        cout<<name<<"\t"<<len<<"\t"<<i<<"\t"<<kmer_rc<<"\t-"<<endl;
                    }
                }
            }
        }
    }
    fa.close();
    }
/*
    char s[25];
    while(printf("Input seq: ") && scanf("%s",s)){
        uint32_t ibin=tobin(s);
        for (auto d:kmer){
            if (d.count(ibin))
                cout<<"yes"<<endl;
        }    
    }
    
    printf("Reading finished!");
    */
    return 0;
}
