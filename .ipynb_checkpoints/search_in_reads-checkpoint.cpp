#include <iostream>
#include <fstream> 
#include <hash_map>
#include <string>

using namespace std;
using namespace __gnu_cxx;

struct str_hash
{
    size_t operator()(const string &s) const
    {
        return __stl_hash_string(s.c_str());
    }
};

struct str_compare
{
    int operator()(const string &a, const string &b) const
    {
        return (a==b);
    }
};

typedef hash_map<string, string, str_hash, str_compare> StrMap;

string getrc(string kmer){
    char kmer_rc[21] = {0};
    for(int i = 20;i >= 0;i--){
        char c = kmer[i];
        if(c == 'A' || c == 'a')
            kmer_rc[21 - i - 1] = 'T';
        else if(c == 'T' || c == 't')
            kmer_rc[21 - i - 1] = 'A';
        else if(c == 'C' || c == 'c')
            kmer_rc[21 - i - 1] = 'G';
        else if(c == 'G' || c == 'g')
            kmer_rc[21 - i - 1] = 'C';
        else{
            cout<<"ERROR!!!!!letters out of ATCG"<<endl;
            exit(0);
        }
    }
    return kmer_rc;
}

int main(int argc,char* argv[]){
    StrMap hashmap;
    ifstream kmerfile(argv[1]); 
    string temp,kmer,name,kmer_rc;
    int len;
    while(getline(kmerfile,temp)){
        kmer = temp.substr(0,21);
        hashmap[kmer] = " ";
    }
    kmerfile.close();
    ifstream fafile("/data/panweihua/human_chm13/hifi/fastqs/chm13.asm.p_ctg.compressed.fa"); 
    while(getline(fafile,temp)){
        if(temp[0] == '>'){
            name = temp.substr(1,20);
        }
        else{
            len = temp.length();
            for(int i = 0;i < len - 20;i++){
                kmer = temp.substr(i,21);
                kmer_rc = getrc(kmer);
                if(hashmap.find(kmer) != hashmap.end()){
                    cout<<name<<"\t"<<len<<"\t"<<i<<"\t"<<kmer<<"\t+"<<endl;
                }
                else if(hashmap.find(kmer_rc) != hashmap.end()){
                    cout<<name<<"\t"<<len<<"\t"<<i<<"\t"<<kmer_rc<<"\t-"<<endl;
                }
            }
        }
    }
    fafile.close();
}
