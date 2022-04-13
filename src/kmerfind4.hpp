#include <string>


//  t1 is number of dicts; t2 is number of threads for finding kmers.
int build_pos(const char  *fasta_file, std::string out_path, bool type, int t1, int t2);

// Read the kmerfile and store in hashtable
int read_kmer(const char *k_path,int t);
