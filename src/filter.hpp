#include<string>
using namespace std;

//read file and builing dict,paffile,reffile,query file;
int read_file(std::string align_file, std::string r_pos, std::string q_pos, bool fmt,
             std::string out_path, int threshold, int thread_num);
//calculate the kmap-Q item is line of paffile with spilt ,fp is query file's pointer 
int evaluation(std::string item[], FILE *fp);

//Traveral the aligment file(paf) and calculate the kmap
// int main();
