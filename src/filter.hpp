#include<string>
using namespace std;

//read file and builing dict,paffile,reffile,query file;
int read_file(string align_file,string r_pos,string q_pos,bool fmt,string out_path);

//calculate the kmap-Q item is line of paffile with spilt ,fp is query file's pointer 
int evaluation(string item[],FILE *fp);

//Traveral the aligment file(paf) and calculate the kmap
// int main();
