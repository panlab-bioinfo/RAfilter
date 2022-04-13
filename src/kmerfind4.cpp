# include <cstdio>
# include <ctime>
# include <cstdlib>
# include <iostream>
# include <string>
# include <vector>
# include <fstream>
# include <mutex>
# include <bitset>
# include "thread_pool.hpp"
# include "uthash.h"
# include "kmerfind4.hpp"

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
        }
    }
    hash_node *find(uint64_t *key)
    {
        hash_node *tmp=NULL;
        HASH_FIND_INT(dict, key, tmp); 
        return tmp;
    }
    int size() {
        return HASH_COUNT(dict);
    }
};


vector<hash_set> kmer;
// 读取k-mer为二进制整数
uint64_t ctoi(char c){
    switch (c)
    {
        case 'A': return 0ull;
        case 'T': return 3ull;
        case 'C': return 1ull;
        case 'G': return 2ull;
        default:
            return 4ull;
    }
}

// 
uint64_t tobin(string *s){
    uint64_t k_value=0;
    int i=0;
    for(auto c=s->begin();c<s->end();c++){
        uint64_t v=ctoi(*c);
        if (v!=4){
            k_value = k_value | (v<<i);
            i+=2;
        }
    }
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
    uint64_t i=0;
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
    cout<<dict.size()<<"\n";
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
    uint64_t line_len = s.size()+1;
    k_mer.close();
    uint64_t row_cnt=total_len/line_len/t;  //每块的长度
    if(row_cnt*t*line_len<total_len){
        row_cnt++;
    }    
    for(int i=0;i<t;i++){
        pool.submit(creat_dict,k_path, row_cnt*line_len, i);
    }
    pool.shutdown();
    return 1;
}


/* store the sbin(bin of kmer),rbin(reverse bin of kmer), pos at contig and contig name*/  
class KMER{           
public:
    uint64_t sbin,rbin;
    uint32_t pos;
    int nameid;
    int is_ukmer(int i,FILE *fp){
            if(kmer[i].find(&sbin)!=NULL){
                        //cout<<nameid<<"\t"<<pos<<"\t"<<sbin<<"\t+"<<endl;
                        fprintf(fp,"%-8x\t%-11lx\n",pos,sbin);
                        // return '+';
                    }
            else if(kmer[i].find(&rbin)!=NULL){
                        // cout<<nameid<<"\t"<<"\t"<<i<<"\t"<<rbin<<"\t-\n";
                        fprintf(fp,"%-8x\t%-11lx\n",pos,rbin);
						// return '-';
                    }
            return 0;
    }
};

int search_kmer(string line, int t1, string name, FILE *file[], bool mask[], int n){
    //Find a free file
    FILE *fp;
    int index;
    m.lock();
    for (int i=0;i<n;i++){
        if (mask[i]==0){
            index=i;
            fp = file[i];
            mask[i]=1;
            break;
        }
    }
    m.unlock();
    KMER k;
    fprintf(fp,"@%s\n",name.c_str());
    uint64_t base;
    uint32_t len=line.length();
    string qkmer=line.substr(0,21);
    k.sbin=tobin(&qkmer);
    k.rbin=reversebin(k.sbin);
    //cout<<k.rbin<<" "<<k.sbin<<"\n";
    for(uint32_t i = 21;i < len;i++){
        k.pos=i-20;
        for(int j=0;j<t1;j++){
            k.is_ukmer(j,fp);
        }
        base=ctoi(line[i]);
        k.sbin = k.sbin >> 2 | base << 40;
        k.rbin = ((k.rbin << 2 | (3ull-base)) & 0x3ffffffffffull);
    }
    for(int j=0;j<t1;j++){
            k.is_ukmer(j,fp);
        }
    m.lock();
    mask[index]=0;
    m.unlock();
    return 1;
}

int build_pos(const char  *fasta_file, string out_path, bool type, int t1, int t2)
{
    auto start = clock();
    cout<<"Number of threads: "<<t2<<"\n"<<"Reading the kmer to dicts\n";
    ThreadPool fpool(t2);
    fpool.init();
    ifstream fa(fasta_file);
    string line,name;

    FILE *file[t2];

    /* 
        Creat mutiple files for mutiple threads;
        File pointers are store in a FILE * array.
    */ 
    char outname[15];
	if (type)
		strcpy(outname,"/kmerAA_r.pos");
	else
		strcpy(outname,"/kmerAA_q.pos");
	string mk="mkdir -p "+out_path;
	system(mk.c_str());

    for (int i=0;i<t2;i++){
		char *outfile = (char*)malloc((out_path.size()+25)*sizeof(char));
		strcpy(outfile,out_path.c_str());strcat(outfile,outname);	
        file[i]=fopen(outfile,"w");
        if(outname[5]=='Z')
            outname[5]='a';
        else
            outname[5]=outname[5]+1;
        if (outname[5] == 'z'){
            if(outname[6]=='Z')
                outname[6]='a';
            else
                outname[6]++;
        } 
    }
    bool mask[t2]={0};
    // auto size=sizeof(KMER);
	cout<<"Searching kmer in fasta file!\n";
    while (getline(fa,line)){
        if(line[0] == '>'){
            name=line.substr(1,line.size());
        }
        else
            fpool.submit(search_kmer, line, t1, name, &file[0], &mask[0], t2);
    }
    fpool.shutdown();
    for(int i=0;i<t2;i++){
        fclose(file[i]);
    }
    fa.close();
	string cmd,cmd2;
	if(type){
		cmd="cat " + out_path +"/* > "+ out_path + "/ref.pos";
		cmd2="rm "+out_path + "/*_r.pos";
	}
	else{
		cmd="cat " + out_path +"/* > "+ out_path + "/query.pos";
		cmd2="rm "+out_path + "/*_q.pos";
	}
	// system(cmd.c_str());
	// system(cmd2.c_str());
    return 1;
}
