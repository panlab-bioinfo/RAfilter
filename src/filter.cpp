#include<bits/stdc++.h>
#include <sam.h>
using namespace std;

typedef list<uint64_t> KLIST;
typedef pair<uint32_t,uint64_t> KMER;
typedef vector<KMER> KMER_LIST;
typedef unordered_map<string,pair<uint64_t,uint64_t>> lib;
lib ref_lib;
lib reads_lib;
KMER_LIST ref_data;
int static_threshold;

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
    uint64_t fposstart,fposend;
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
    uint32_t p;
	uint64_t fposstart,fposend;
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
    uint32_t p;
    uint64_t k,fposstart,fposend;
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

uint32_t getQueryStart(bam1_t *src)
{
    uint32_t * cigar_p;
    uint32_t start_offset = 0;
    uint32_t k, op;
    cigar_p =bam_get_cigar(src);
    uint32_t len=src->core.n_cigar;
   for(int k=0;k<len;k++)
    {
        op = cigar_p[k] & BAM_CIGAR_MASK;
        if (op == BAM_CHARD_CLIP){
            if (start_offset != 0 && start_offset != src->core.l_qseq)
            {
                cerr << "ERROR: Wrong base with value "<<start_offset<< endl ;
            }
        }   
        else if (op == BAM_CSOFT_CLIP) {
            start_offset += cigar_p[k] >> BAM_CIGAR_SHIFT;
        }
        else    break;
    }
    return start_offset;
}

uint32_t getQueryEnd(bam1_t *src)
{    
    uint32_t * cigar_p = bam_get_cigar(src);
    uint32_t end_offset = src->core.l_qseq;
    uint32_t k, op;
    uint32_t len=src->core.n_cigar;
    // if there is no sequence, compute length from cigar string
    if( end_offset == 0){
        for(int k=0;k<len;k++)
        {
            op = cigar_p[k] & BAM_CIGAR_MASK;
            if (op == BAM_CMATCH ||
               op == BAM_CINS ||
               op == BAM_CEQUAL ||
               op == BAM_CDIFF ||
              (op == BAM_CSOFT_CLIP && end_offset == 0))
                end_offset += cigar_p[k] >> BAM_CIGAR_SHIFT;
        }
    }
    else{
        //walk backwards in cigar string
        for(int k=len-1;k>=1;k--)
        {   
            op = cigar_p[k] & BAM_CIGAR_MASK;
            if (op == BAM_CHARD_CLIP)
            {   
                if (end_offset != src->core.l_qseq)
                    cerr << "ERROR: Wrong base with value "<< end_offset << endl ;
            }
            else if(op == BAM_CSOFT_CLIP)
            {
                end_offset -= cigar_p[k] >> BAM_CIGAR_SHIFT;
            }
            else    break;
        }
    }
    return end_offset;
}


int paf_evaluation(string item[],FILE *fp)
{
    int v;
    float p;
    unordered_map<uint64_t,uint32_t> ref_list;
    KLIST reads_list;
    ref_list=searchinref(atoi(item[7].c_str()),atoi(item[8].c_str()));
    reads_list=binary_searchinreads(item[0],atoi(item[2].c_str()),atoi(item[3].c_str()),fp);
    if(item[4]=="-")   reads_list.reverse();
    int common=0;
    common=longestCommonSubsequence(ref_list,reads_list);
	int fenmu=min(ref_list.size(),reads_list.size());
	// cout<<ref_list.size()<<" "<<reads_list.size()<<endl;
    p=(fenmu-common)*1.0/fenmu;
	if(fenmu==0){
		v = -1;
	}else if(p<=1e-6){
		v = 60;
	}
	else{
		v=-10 * log10(p);
	}
    return v;
}

int sam_evaluation(int r_posstart,int r_posend,string queryname,int q_posstart,int q_posend ,FILE *fp)
{
    int v;
    float p;
    unordered_map<uint64_t,uint32_t> ref_list;
    KLIST reads_list;
    ref_list=searchinref(r_posstart,r_posend);
    reads_list=binary_searchinreads(queryname,q_posstart,q_posend,fp);
    int common=0;
    common=longestCommonSubsequence(ref_list,reads_list);
	int fenmu=min(ref_list.size(),reads_list.size());
	// cout<<ref_list.size()<<" "<<reads_list.size()<<endl;
    p=(fenmu-common)*1.0/fenmu;
	if(fenmu==0){
		v = -1;
	}else if(p<=1e-6){
		v = 60;
	}
	else{
		v=-10 * log10(p);
	}
    return v;
}


int paf_filter(string paf,string r_pos, FILE *q_pos, FILE *fo){
	string temp;
	string name = " ";
	string b[100];
	int value;
	int filter = 0, total = 0;
	ifstream fp(paf);
	while(getline(fp,temp)){
		split(temp,b);
		if(name!=b[5])
		{
			ref_data.clear();
			ref_lib_init(b[5],r_pos);
			name=b[5];
		}
		value=paf_evaluation(b,q_pos);
		total++;
		if (value > static_threshold || value < 0){
			fprintf(fo,"%s\n",temp.c_str());
			filter++;
		}
	}
	cout.width(30);
	cout<<"Total alignments count:\t";
	cout<<total<<"\n";
	cout<<"Filter alignments count:\t";
	cout<<filter<<endl;
	fp.close();
	return 1;
}


int sam_filter(string inpath,string r_pos, FILE *q_pos, FILE *fo,string outpath) {
		htsFile *in,*out;
		sam_hdr_t *hdr;
		bam1_t *b;
		hts_idx_t *idx = NULL;
		hts_itr_t *iter = NULL;
		int ret;
        int value;
		int filter = 0, total = 0;
		string name="/";
        string outflie=outpath+"filter.sam";
		if ((in = hts_open(inpath.c_str(), "rb")) == NULL) {
		fprintf(stderr, "Error opening '%s'\n", inpath.c_str());
		return -3;
		}
		if ((out = hts_open(outpath.c_str(), "w")) == NULL) {
			fprintf(stderr, "Error opening '%s'\n", outpath.c_str());
			return -3;
		}
		if ((hdr = sam_hdr_read(in)) == NULL) {
			fprintf(stderr, "[E::%s] couldn't read header for inputs'\n", __func__);
			return -1;
		}
		if ((b = bam_init1()) == NULL) {
			fprintf(stderr, "[E::%s] Out of memory allocating BAM struct.\n", __func__);
			goto fail;
		}
		if (sam_hdr_write(out,hdr) < 0) {
			fprintf(stderr, "[E::%s] Error writing alignments.\n", __func__);
			goto fail;
		}
		while ((ret = sam_read1(in, hdr, b)) >= 0) {
			int r_startpos = b->core.pos ;
			int r_endpos=bam_endpos(b);
			string refname = "*";
			if (b->core.tid != -1)
				refname = hdr->target_name[b->core.tid];          
			string queryname = bam_get_qname(b);
			int q_startpos=getQueryStart(b);
			int q_endpos=getQueryEnd(b);
			if(name!=refname)
			{
                ref_data.clear();
				ref_lib_init(refname,r_pos);
                name=refname;
			}
            value=sam_evaluation(r_startpos,r_endpos,queryname,q_startpos,q_endpos,q_pos);
			total++;
			if (value > 12 || value < 0){
				filter++;
				if (sam_write1(out, hdr, b) < 0) {
					fprintf(stderr, "[E::%s] Error writing alignments.\n", __func__);
					goto fail;
				}
			}
			
		}

		cout.width(30);
		cout<<"Total alignments count:\t";
		cout<<total<<"\n";
		cout<<"Filter alignments count:\t";
		cout<<filter<<endl;
		if (ret < -1) {
			fprintf(stderr, "[E::%s] Error parsing input.\n", __func__);
			goto fail;
		}
		bam_destroy1(b);
		sam_hdr_destroy(hdr);
		if ((ret = hts_close(out)) < 0) {
			fprintf(stderr, "Error closing output.\n");
			return  -3;
		}
		if ((ret = hts_close(in)) < 0) {
			fprintf(stderr, "Error closing input.\n");
			return  -3;
		}
		return 0;
	fail:
		if (iter) sam_itr_destroy(iter);
		if (b) bam_destroy1(b);
		if (idx) hts_idx_destroy(idx);
		if (hdr) sam_hdr_destroy(hdr);
		if ((ret = hts_close(out)) < 0) {
			fprintf(stderr, "Error closing output.\n");
			return  -3;
		}
		if ((ret = hts_close(in)) < 0) {
			fprintf(stderr, "Error closing input.\n");
			return  -3;
		}
		return 1;
}


int read_file(string align_file,string r_pos,string q_pos,bool fmt,string out_path,int thre)
{
	static_threshold = thre;
    lib_init(r_pos,ref_lib);
    lib_init(q_pos,reads_lib);
	string outfile = out_path+"/filter.paf";
	cout<<outfile<<endl;
	FILE *fo = fopen(outfile.c_str(),"w");
	FILE *fq = fopen(q_pos.c_str(),"r"); 
	cout<<"readfile: "<< fq <<"\n"<<"outputfile:"<< fo <<endl;
	if (fmt){								//The format of alignment file is BAM
		;
        sam_filter(align_file,r_pos,fq,fo,out_path);
		// bamfile handle
	}	
	else{	//The format of alignment file is PAF
		string sortcmd = "sort -k6,6 " + align_file + " > "+out_path+"/aln.sort.paf";
		system (sortcmd.c_str());
		string sort_aln_file = out_path + "/aln.sort.paf";
		paf_filter(sort_aln_file, r_pos, fq, fo);
	}
	fclose(fq);
	fclose(fo);
    return 1;
}

// int main(){
// 	string path = "/data/yangjinbao/data/repeat/sim/human_cut/ont/";
// 	string paf = path + "sort.aln.paf";
// 	string rpos = path + "result/ref.pos";
// 	string qpos = path + "result/query.pos";
// 	string output = path+"result/";
// 	read_file(paf,rpos,qpos,0,output);
// 	return 0;
// }
