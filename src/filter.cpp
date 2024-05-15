/* RAfilter filter */
// false positive alignments are filter out in this step of program.

#include <bits/stdc++.h>
#include "htslib/htslib/sam.h"

#ifndef FILE_SPLIT_
#define FILE_SPLIT_
#include "file_processing.hpp"
#endif // file_processing.hpp

#ifndef MULTIPLE_PROCESS_
#define MULTIPLE_PROCESS_
#include "thread_pool.hpp"
#endif // thread_pool.hpp

using namespace std;

typedef list<uint64_t> KLIST;
typedef pair<uint32_t,uint64_t> KMER;
typedef vector<KMER> KMER_LIST;
typedef unordered_map<string,pair<uint64_t,uint64_t>> lib;
lib ref_lib;
lib reads_lib;
KMER_LIST ref_data;
int static_threshold;

class Ref{
    public:
      string refname_;
      KMER_LIST ref_data_;
      void ref_lib_init(string path);
      unordered_map<uint64_t, uint32_t> searchinref(int startpos, int endpos);
};

vector<string> split(const string& str, const char& ch) 
{
    string next;
    vector<string> result;

    // For each character in the string
    for (string::const_iterator it = str.begin(); it != str.end(); it++) {
        // If we've hit the terminal character
        if (*it == ch) {
            // If we have some characters accumulated
            if (!next.empty()) {
                // Add them to the result vector
                result.push_back(next);
                next.clear();
            }
        } else {
            // Accumulate the next character into the sequence
            next += *it;
        }
    }
    if (!next.empty())
         result.push_back(next);
    return result;
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

void Ref::ref_lib_init(string path)
{
    uint32_t p;
    uint64_t fposstart, fposend;
    fposstart = ref_lib[this->refname_].first;
    fposend = ref_lib[this->refname_].second;
    uint64_t k;
    string buf;
    FILE *reffile;
    reffile = fopen(path.c_str(), "r");
    assert(reffile);

    fseek(reffile, fposstart, SEEK_SET);
    while (ftell(reffile) < fposend)
    {
        fscanf(reffile, "%8x %11lx", &p, &k);
        this->ref_data_.push_back(make_pair(p, k));
    }
    fclose(reffile);
}

inline bool compare(KMER a,uint32_t b)
{
    return a.first<b; 
}

unordered_map<uint64_t, uint32_t> Ref::searchinref(int startpos, int endpos)
{
    unordered_map<uint64_t, uint32_t> klist;
    auto iters = lower_bound(this->ref_data_.begin(), this->ref_data_.end(), startpos, compare);
    // cout<<iters->first<<endl;
    auto itere = lower_bound(this->ref_data_.begin(), this->ref_data_.end(), endpos, compare);
    while (iters != itere)
    {
        klist[iters->second] = iters->first;
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

int paf_evaluation(vector<string> item, FILE *fp, Ref target)
{
    int v;
    float p;
    unordered_map<uint64_t, uint32_t> ref_list;
    KLIST reads_list;
    ref_list = target.searchinref(atoi(item[7].c_str()), atoi(item[8].c_str()));
    reads_list = binary_searchinreads(item[0], atoi(item[2].c_str()), atoi(item[3].c_str()), fp);
    if (item[4] == "-")
        reads_list.reverse();
    int common = 0;
    common = longestCommonSubsequence(ref_list, reads_list);
    int fenmu = min(ref_list.size(), reads_list.size());
	// cout<<ref_list.size()<<" "<<reads_list.size()<<endl;
    p = (fenmu - common) * 1.0 / fenmu;
    if (fenmu == 0)
    {
		v = -1;
    }
    else if (p <= 1e-6)
    {
		v = 60;
	}
    else
    {
        v = -10 * log10(p);
	}
    return v;
}

int sam_evaluation(int r_posstart, int r_posend, string queryname,
    int q_posstart, int q_posend, FILE *fp, Ref target,int flag)
{
    int v;
    float p;
    unordered_map<uint64_t, uint32_t> ref_list;
    KLIST reads_list;
    ref_list = target.searchinref(r_posstart, r_posend);
    reads_list = binary_searchinreads(queryname, q_posstart, q_posend, fp);
    if (flag & 16 ==1 )
        reads_list.reverse();
    int common = 0;
    common = longestCommonSubsequence(ref_list, reads_list);
    int fenmu = min(ref_list.size(), reads_list.size());
	// cout<<ref_list.size()<<" "<<reads_list.size()<<endl;
    p = (fenmu - common) * 1.0 / fenmu;
    if (fenmu == 0)
    {
		v = -1;
    }
    else if (p <= 1e-6)
    {
		v = 60;
	}
    else
    {
        v = -10 * log10(p);
	}
    return v;
}

int paf_filter(string paf, string r_pos, string q_pos,string outfile)
{
    FILE *fo = fopen(outfile.c_str(), "w");
    assert(fo);
    FILE *fq = fopen(q_pos.c_str(), "r");
    assert(fq);

    vector<string> b;
	string temp;
    Ref target;
    target.refname_ = " ";
	int value;
    int correct_alignment = 0, total = 0;
	ifstream fp(paf);
    while (getline(fp, temp))
    {
        b=split(temp,'\t');
        if (target.refname_ != b[5])
		{
            target.ref_data_.clear();

            target.refname_ = b[5];
            target.ref_lib_init(r_pos);
            // cout<<"ref read\n";

        }
        value = paf_evaluation(b, fq, target);
		total++;
		if (value > static_threshold || value < 0){
            //printf("%s\n", temp.c_str());
            fprintf(fo, "%s\t%d\n", temp.c_str(),value);
            correct_alignment++;
		}
	}
	cout.width(30);
    cout << "Total alignments count:\t";
    cout << total << "\n";
    cout << "correct  alignments count:\t";
    cout << correct_alignment << endl;
	fp.close();
    fclose(fo);
    fclose(fq);
	return 1;
}

int sam_filter(string inpath, string r_pos, FILE *q_pos, string outpath)
{
    htsFile *in, *out;
		sam_hdr_t *hdr;
		bam1_t *b;
		hts_idx_t *idx = NULL;
		hts_itr_t *iter = NULL;
		int ret;
        int value;
    int correct_alignment = 0, total = 0;
    Ref target;
    target.refname_ = "/";
    string outflie = outpath + "filter.sam";
    if ((in = hts_open(inpath.c_str(), "rb")) == NULL)
    {
		fprintf(stderr, "Error opening '%s'\n", inpath.c_str());
		return -3;
		}
    if ((out = hts_open(outpath.c_str(), "w")) == NULL)
    {
			fprintf(stderr, "Error opening '%s'\n", outpath.c_str());
			return -3;
		}
    if ((hdr = sam_hdr_read(in)) == NULL)
    {
			fprintf(stderr, "[E::%s] couldn't read header for inputs'\n", __func__);
			return -1;
		}
    if ((b = bam_init1()) == NULL)
    {
			fprintf(stderr, "[E::%s] Out of memory allocating BAM struct.\n", __func__);
			goto fail;
		}
    if (sam_hdr_write(out, hdr) < 0)
    {
			fprintf(stderr, "[E::%s] Error writing alignments.\n", __func__);
			goto fail;
		}
    while ((ret = sam_read1(in, hdr, b)) >= 0)
    {
        int r_startpos = b->core.pos;
        int flag=b->core.flag;
        int r_endpos = bam_endpos(b);
			string refname = "*";
			if (b->core.tid != -1)
				refname = hdr->target_name[b->core.tid];          
			string queryname = bam_get_qname(b);
        int q_startpos = getQueryStart(b);
        int q_endpos = getQueryEnd(b);
        if (target.refname_ != refname)
			{
            target.ref_data_.clear();
            target.refname_ = refname;
            target.ref_lib_init(r_pos);
        }
        value = sam_evaluation(r_startpos, r_endpos, queryname, q_startpos,
            q_endpos, q_pos, target,flag);
			total++;
        if (value > static_threshold || value < 0)
        {
            correct_alignment++;
            if (sam_write1(out, hdr, b) < 0)
            {
					fprintf(stderr, "[E::%s] Error writing alignments.\n", __func__);
					goto fail;
				}
			}
			
		}

		cout.width(30);
    cout << "Total alignments count:\t";
    cout << total << "\n";
    cout << "Correct alignments count:\t";
    cout << correct_alignment << endl;
    if (ret < -1)
    {
			fprintf(stderr, "[E::%s] Error parsing input.\n", __func__);
			goto fail;
		}
		bam_destroy1(b);
		sam_hdr_destroy(hdr);
    if ((ret = hts_close(out)) < 0)
    {
			fprintf(stderr, "Error closing output.\n");
        return -3;
		}
    if ((ret = hts_close(in)) < 0)
    {
			fprintf(stderr, "Error closing input.\n");
        return -3;
		}
		return 0;
	fail:
    if (iter)
        sam_itr_destroy(iter);
    if (b)
        bam_destroy1(b);
    if (idx)
        hts_idx_destroy(idx);
    if (hdr)
        sam_hdr_destroy(hdr);
    if ((ret = hts_close(out)) < 0)
    {
			fprintf(stderr, "Error closing output.\n");
        return -3;
		}
    if ((ret = hts_close(in)) < 0)
    {
			fprintf(stderr, "Error closing input.\n");
        return -3;
		}
		return 1;
}




int read_file(string align_file, string r_pos, string q_pos, bool fmt,
             string out_path, int threshold, int thread_num){
    static_threshold = threshold;
    lib_init(r_pos, ref_lib);
    lib_init(q_pos, reads_lib);
    string outfile = out_path + "/rafiltered.paf";
    if (fmt) { // The format of alignment file is BAM
        FILE *fq = fopen(q_pos.c_str(), "r");

        assert(fq);

        sam_filter(align_file, r_pos, fq, out_path);
        fclose(fq);
    } // bamfile processing
    else { // The format of alignment file is PAF
		string sortcmd = "sort -k6,6 -T"+ out_path +" "+ align_file + " > "+out_path+"/aln.sort.paf";
		system (sortcmd.c_str());
		string sort_aln_file = out_path + "/aln.sort.paf";
        splitpaf (sort_aln_file, thread_num, out_path); 
        ThreadPool pool(thread_num);
        pool.init();
        vector<string> files;
        get_dir_file(out_path, "subpaf*", files);
        cout << "Split the alignments file.\n";
        int i = 0; // Tail marker of mutiple output file
        for (auto subpaf : files){
            string sub_outfile =out_path + "/subfilter" + to_string(++i);
            cout << "filter processing : " << subpaf <<"\n"
            <<"sub outputfile: " << sub_outfile;
            pool.submit(paf_filter, subpaf, r_pos, q_pos, sub_outfile);
	}
        pool.shutdown();
        string cmd = "cat " + out_path + "/subfilter* >" + out_path
                +"/rafiltered.paf; rm " +out_path + "/sub*";
        system(cmd.c_str());
    }
    cout << "Alignments filter finished" << endl;
    return 1;
}

int splitfile(string file){

;

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
