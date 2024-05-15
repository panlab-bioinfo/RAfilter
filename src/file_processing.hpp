// A function of spliting paf file and sam/bam file by a constant
// the function can make an application on mutiple processing

#include <stdio.h>
#include <string>
#include <vector>
// #ifndef HTS_LIB_
// #define HTS_LIB_
// #include "htslib/htslib/sam.h"
// #endif

#define M 0

using namespace std;


// Get matching files with a certain pattern from a given path
void get_dir_file(const string path, const string pattern, vector<string> &files){
	string cmd = "ls " + path + "/" + pattern;
	FILE *ls_result = popen(cmd.c_str(), "r");
	char tmp[100] = {"\0"};
	string file;
	while (!feof(ls_result)) {
		fscanf(ls_result, "%s\n", tmp);
		printf("%s\n", tmp);
		file = tmp;
		files.push_back(file);
	}
	pclose(ls_result);
}

// split the paf file to mutiple files
int splitpaf(const string paf_file, const int thread, const string out_path)
{
	int linenum; // Line number
	string cmd = "wc -l " + paf_file;
	FILE *cmdResult = popen(cmd.c_str(), "r");
	char result[512];
	if (cmdResult != NULL){
		fscanf(cmdResult, "%d\t%s\n", &linenum, result);
		pclose(cmdResult);
	}
	int block = linenum / thread + 1;
	// Use command line "split" to split the paf file to a Given folder
	cmd = "split -l " + to_string(block) + " " + paf_file.c_str() + " " +
			out_path.c_str() + "/subpaf";
	system(cmd.c_str());
	return 0;
}

#if M // Debug switch
int main(int argc, char *argv[])
{
	vector<string> files;
	get_dir_file(argv[1], "*cpp", files);
	for (auto i : files){
		printf("%s\n", i.c_str());
	}
	// int t = atoi(argv[2]);
	// splitpaf(file, t, argv[3]);
	return 0;
}
#endif
