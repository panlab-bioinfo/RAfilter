//# include "src/thread_pool.hpp"
# include <mutex>
# include <vector>
# include <unordered_set>
# include <iostream>
//# include "src/parallel_hashmap/phmap.h"
# include <algorithm>
# include <cstdio>
using namespace std;
/*
mutex mt;
vector<phmap::flat_hash_set<uint32_t>> arrays;
int put_array(vector<uint32_t> a){
    phmap::flat_hash_set<uint32_t> d;
    d.insert(a.cbegin(),a.cend());
    mt.lock();
    arrays.push_back(d);
    mt.unlock();
    return 1;
}
*/
void print_bin(uint64_t n){
    int l = sizeof(n)*8;//总位数。
    cout<<"total_len: "<<l<<endl;
    int i;
    if(n == 0){
         printf("0");
         return;
     }
    // for(i = l-1; i >= 0; i --){
    //     if(n&(1<<i)) break;
    // }
    for(i=41;i>=0; i --)
        printf("%d", (n&(1<<i)) != 0);
    printf("\n");
}


int ctoi(char c){
    switch (c)
    {
        case 'A': return 0;
        case 'T': return 1;
        case 'C': return 2;
        case 'G': return 3;
        default:
            return -1;
    }
}

uint64_t tobin(string *s){
    uint64_t k_value=0;
    int i=0;
    for(auto c=s->begin();c!=s->end();c++){
        // cout<<*c;
        int v=ctoi(*c);
        if (v!=-1){
            k_value = k_value | (v<<i);
            i+=2;
        }
        else
            break;
    }
    // cout<<endl;
    return k_value;
}
uint64_t reversebin(uint64_t x){
    uint64_t mid = x & 0x300000;
    x = ((x & 0x3ffffc00000) >> 22) | ((x & 0xfffff) << 22);
    x = (x & 0x3ff000ffc00) >> 10 | (x & 0xffc003ff) << 10;
    uint64_t quarter = x & 0x300c00c030u;
    x = (x & 0x3C0F00F03C0) >> 6 | (x & 0xF03C03C0F)<<6;
    x = (x & 0x30CC30C330C) >> 2 | (x & 0xC330C30CC3)<<2;
    x = x | mid | quarter;
    return x;
}
int main(){
    string s="AATTCGGTCGTTGTCGTCGGC";
    cout<<s<<endl;
    uint64_t sbin=tobin(&s);
    print_bin(sbin);
    uint64_t rbin=reversebin(sbin);
    print_bin(rbin);
    cout<<endl;
    cout<<sbin<<"   "<<rbin<<endl;
    reverse(s.begin(),s.end());
    sbin=tobin(&s);
    cout<<sbin<<endl;
    print_bin(sbin);
    return 0;
}

//1110011101011110011111100101000011100111010111100111111001010000
//AATTCGGTCGTTGTCG