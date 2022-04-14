#include<string>
using namespace std;

//read file and builing dict,paffile,reffile,query file;
int read_file(string paffile,string reference,string queryfile);

//calculate the kmap-Q item is line of paffile with spilt ,fp is query file's pointer 
int evaluation(string item[],FILE *fp);

//Traveral the aligment file(paf) and calculate the kmap
