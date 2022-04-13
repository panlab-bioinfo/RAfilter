# include <cstdio>
# include <ctime>
# include <cstdlib>
# include <iostream>
# include <string>
# include <vector>
# include <fstream>
# include <mutex>
# include "thread_pool.hpp"
# include "uthash.h"

using namespace  std;

//vector<phmap::flat_hash_set<uint64_t>> kmer; //存储kmer的hash>
mutex m;
//vector<unordered_set<uint64_t>> kmer;
//

typedef struct my_struct {
    uint64_t key;                    /* key */
    //char name;
    UT_hash_handle hh;         /* makes this structure hashable */
} hash_node;

class hash_set{
public:
    hash_node *dict=NULL;
    
    void insert(uint64_t *key){
        hash_node *s; 
        HASH_FIND_INT(dict, key, s);
        if (s == NULL)
        {
            s = (struct my_struct*)malloc(sizeof *s);
            s->key = *key;
            HASH_ADD_INT(dict, key, s);
            // cout<<"insert finished!"<<endl;
        }
    }

    hash_node *find(uint64_t *key)
    {
        hash_node *tmp=NULL;
        HASH_FIND_INT(dict, key, tmp); 
        return tmp;
    }
};

vector<hash_set> kmer;

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
            return 4;
    }
}

// 
uint64_t tobin(string *s){
    uint64_t k_value=0;
    int i=0;
    for(auto c=s->begin();c<s->end();c++){
        // cout<<*c;
        uint64_t v=ctoi(*c);
        if (v!=4){
            k_value = k_value | (v<<i);
            i+=2;
        }
    }
    // cout<<endl;
    return k_value;
}
//reverse the bin of kmer
uint64_t reversebin(uint64_t x){
    uint64_t mid = x & 0x300000ull;
    x = (x & 0x3ffffc00000ull) >> 22|((x & 0xfffffull) << 22);
    x = (x & 0x3ff000ffc00ull) >> 10 | (x & 0xffc003ffull) << 10;
    uint64_t quarter = x & 0x300c00c030ull;
    x = (x & 0x3C0F00F03C0ull) >> 6 | (x & 0xF03C03C0Full) << 6;
    x = (x & 0x30CC30C330Cull) >> 2 | (x & 0xC330C30CC3ull) << 2;
    x = x | mid | quarter;
    return (~x) & 0x3ffffffffffull;
}



int creat_dict(const char *k_path,uint64_t step_len, int idx){
    ifstream k_mer(k_path);
    //phmap::flat_hash_set<uint64_t> dict;
    //unordered_set<uint64_t> dict;
    hash_set dict;
    int i=0;
    uint64_t v;
    uint64_t k_value=0;
    // printf("%s",s);
    char *tmp=new char[step_len+2];
    *tmp={'\0'};
    k_mer.seekg(idx*(step_len),k_mer.beg);
    k_mer.read(tmp,step_len);
    for(char *p=tmp;p<tmp+step_len;p++){
        v=ctoi(*p);
        if (v!=4){
            k_value = k_value | (v<<i);
            i+=2;
            // printf(" [%c,%d] ",*p,v);
        }
        else if(v == 4 && k_value != 0){
            dict.insert(&k_value);
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

int main(int argc,char* argv[])
{
    // const char *kmerfile="../kmer.sample.dump";

    const char *kmerfile="/data/panweihua/human_chm13/hifi/fastqs/chm13.asm.p_ctg.compressed.kmers.1.0_100M.dump";
    string fafile="/data/panweihua/human_chm13/hifi/fastqs/chm13.asm.p_ctg.compressed.fa";
	int t=1;
    read_kmer(kmerfile,t);
    ifstream fa(fafile);
    string line,name,qkmer;
    uint64_t kbin,rbin;
    while (getline(fa,line)){
        
        
        if(line[0] == '>'){
            name = line.substr(1,20);
        }
        else{
            
			//cout<<"swich_s: "<<clock()<<endl;
			uint64_t base;
			uint32_t len=line.length();
			string qkmer=line.substr(0,21);
			uint64_t sbin=tobin(&qkmer);
			uint64_t rbin=reversebin(sbin);
			//cout<<"sbin: "<<sbin<<"	rbin:	"<<rbin<<endl;
            for(uint32_t i = 21;i < len;i++){

                for(int j=0;j<t;j++){
                    auto start=clock();
                    if(kmer[j].find(&sbin)!=NULL){
                        printf("%s\t%ld\t%d\t%s\t+\n",name,len,i-20,sbin);
                        // cout<<name<<"\t"<<len<<"\t"<<i<<"\t"<<sbin<<"\t+"<<endl;
                        break;
                    }
                    else if(kmer[j].find(&rbin)!=NULL){
                        printf("%s\t%ld\t%d\t%s\t-\n",name,len,i-20,rbin);
                        // cout<<name<<"\t"<<len<<"\t"<<i<<"\t"<<rbin<<"\t-"<<endl;
						break;
                    }
                }
				base=ctoi(line[i]);
        		sbin = sbin >> 2 | base << 40;
		        rbin = (rbin << 2 | (3ull-base)) & 0x3ffffffffffull;
            }
            
        }
        
    }
    fa.close();
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
