#include<bits/stdc++.h>
using namespace std;

typedef list<uint64_t> KLIST;
typedef pair<uint32_t,uint64_t> KMER;
typedef vector<KMER> KMER_LIST;
typedef unordered_map<string,pair<uint32_t,uint32_t>> lib;
lib ref_lib;
lib reads_lib;
KMER_LIST ref_data;

void split(string q,string item[])
{
    stringstream ss;
    int i=0;
    ss.clear();
    ss.str(q); 
    while(!ss.fail())
    {
        ss>>item[i];
        i++;
    }
}

//initize hash
void lib_init(string path,lib &library)
{
    string buf,name;
    uint32_t fposstart,fposend;
    int flag=0;
    ifstream  posfile(path);
    while(getline(posfile,buf))
    {
        if(buf[0]=='@'){
        if(flag==0)
        {
            name=buf.substr(1,buf.size());
            fposstart=posfile.tellg();
            flag=1;
         } 
        else
        {
            library[name]={fposstart,fposend};
            name=buf.substr(1,buf.size());
            fposstart=posfile.tellg();
        }
        }
        fposend=posfile.tellg();
    }
    library[name]={fposstart,fposend};
    cout<<"library success"<<endl;
}

void ref_lib_init(string name ,string path)
{
    uint32_t p,fposstart,fposend;
    fposstart=ref_lib[name].first;
    fposend=ref_lib[name].second;
    uint64_t k;
    string buf;
    FILE *reffile;
    reffile=fopen(path.c_str(),"r");
    fseek(reffile,fposstart,SEEK_SET);
    while(ftell(reffile)<fposend)
    {
        fscanf(reffile,"%8x %11lx",&p,&k);
        ref_data.push_back(make_pair(p,k));
    }
}



inline bool compare(KMER a,uint32_t b)
{
    return a.first<b; 
}

unordered_map<uint64_t,uint32_t> searchinref(int startpos,int endpos)
{
    unordered_map<uint64_t,uint32_t> klist;
    auto iters= lower_bound(ref_data.begin(),ref_data.end(),startpos,compare);
    //cout<<iters->first<<endl;
    auto itere = lower_bound(ref_data.begin(),ref_data.end(),endpos,compare);
    while(iters!=itere)
    {
        klist[iters->second]=iters->first;
        iters++;
    }
    return klist;
}


KLIST binary_searchinreads(string name,int startpos,int endpos,FILE *fp)
{
    
    KLIST klist;
    uint32_t p,fposstart,fposend;
    uint64_t k;
    int num;
    fposstart=reads_lib[name].first;
    fposend=reads_lib[name].second;
    num=(fposend - fposstart)/21;
    int lo = 0, hi = num;
    while (lo < hi) {
        int mid = (lo + hi) / 2;
        fseek(fp,fposstart+mid*21,SEEK_SET);
        fscanf(fp,"%8x%11lx",&p,&k);
        if (p < startpos) {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    fseek(fp,fposstart+lo*21,SEEK_SET);
    while(ftell(fp)<fposend)
    {
        fscanf(fp,"%8x%11lx",&p,&k);
        if(p>endpos)    break;
        klist.push_back(k);
    }
    return klist;
}



int lengthOfLIS(vector<int>& nums) {
    int len = 1, n = (int)nums.size();
    if (n == 0) {
        return 0;
    }
    vector<int> d(n + 1, 0);
    d[len] = nums[0];
    for (int i = 1; i < n; ++i) {
        if (nums[i] > d[len]) {
            d[++len] = nums[i];
        } else {
            int l = 1, r = len, pos = 0; 
            while (l <= r) {
                int mid = (l + r) >> 1;
                if (d[mid] < nums[i]) {
                    pos = mid;
                    l = mid + 1;
                } else {
                    r = mid - 1;
                }
            }
            d[pos + 1] = nums[i];
        }
    }
    return len;
}
int longestCommonSubsequence(unordered_map<uint64_t,uint32_t> &list1,KLIST &list2)
{
    vector<int> dp;
    auto iter=list2.begin();
    while(iter!=list2.end())
    {
        auto a=list1.find(*iter);
        if(a!=list1.end())
        {
            dp.push_back(a->second);
        }
        iter++;
    }
    return lengthOfLIS(dp);
}

int evaluation(string b[],FILE *fp)
{
    int v;
    float p;
    unordered_map<uint64_t,uint32_t> ref_list;
    KLIST reads_list;
    ref_list=searchinref(atoi(b[7].c_str()),atoi(b[8].c_str()));
    reads_list=binary_searchinreads(b[0],atoi(b[2].c_str()),atoi(b[3].c_str()),fp);
    if(b[4]=="-")   reads_list.reverse();
    int common=0;
    common=longestCommonSubsequence(ref_list,reads_list);
    //cout<<common<<endl;
    p=1-(float)(common)/min(ref_list.size(),reads_list.size());
    v=-10 * log10(p);
    return v;
}
int read_file(string paffiles,string reffile,string readsfile)
{
    string temp; 
    string name=" ";
    string b[100];
    int value;
    ifstream paffile(paffiles);
    FILE *fw=NULL,*fp=NULL;
    string path=paffiles+".fliter";
    fw=fopen(path.c_str(),"w");
    fp=fopen(readsfile.c_str(),"r");
    if (fp == NULL)
    {
        perror("file fopen error!");
        exit(0);
    }
    lib_init(reffile,ref_lib);
    lib_init(readsfile,reads_lib);
    while(getline(paffile,temp))
    {
        split(temp,b);
         if(name!=b[5])
        {
            ref_data.clear();
            ref_lib_init(b[5],reffile);
            name=b[5];
        }
        value=evaluation(b,fp);
        fprintf(fw,"%s\t%d\n",temp.c_str(),value);

    }
    fclose(fp);
    fclose(fw);
    return 1;
}
