/*
Name: RAfilter
description: a programer for correct alignments.
	RAfilter can filter the alignments file with the PAF/BAM file by calculation
of k-mer lists of reference and query. The programer has two sections and at
first, a k-mer with its pos database is set up. Secondly, alignments with the
KMAPQ achieved the threshold are leave. The port to set the threshold is exposed
in options.
*/


#include <iostream>
#include <thread>
#include "kmerfind4.hpp"
#include "filter.hpp"
#ifndef CLIPP_
#define CLIPP_
#include "clipp.h"
#endif //clipp.h

using namespace clipp;
using std::cout;
using std::endl;
using std::string;
using std::thread;

string program_name = "rafilter";
string version = "V1.0.0";

int main(int argc, char *argv[])
{
	bool help = false;
	char mode = 'N';
    string kmer_file = "";
    string ref_file = "";
	string query_file = "";
	string output_path = "./";
	string align = "";
	string rpos = "", qpos = "";
	int threshold = -1;
	bool align_fmt = 0;
	// t1:[1,2,3,4,5]
	int t1 = 3, t2 = 8;
	auto build = "build:"
			% (option("-t", "--threads")/*.doc("Number of threads") */& value("threads", t2),
		option("-o", "--out-put").doc("Diractory of output") & value("out_path", output_path),
			option("-q", "--query").doc("Query/Reads fasta file") & value("query", query_file),
		option("-r", "--reference").doc("Reference fasta file") & value("reffile", ref_file),
			value("kmerfile", kmer_file).doc("Kmer file with creating by jellyfish"));
	auto filter = "filter:"
		% (option("-t", "--threads")/*.doc("Number of threads")*/ & value("threads", t2),
		option("-o", "--out-put").doc("Diractory of output") & value("out_path", output_path),
		(option("-p") | option("-b").set(align_fmt, true)).doc("Format of alignment file [p:paf/b:bam]"),
		option("--threshold").doc("Threshold of KMAPQ to filter") & value("threshold", threshold),
		value("r_pos", rpos).doc("Ref pos file built with build"),
		value("q_pos", qpos).doc("Query pos file built with build"),
		value("paf/bam", align).doc("Alignment file"));

	auto cli = ((command("build").set(mode, 'b'), build)
			| (command("filter").set(mode, 'f'), filter),
			option("-t", "--threads").doc("Number of threads"),
			option("-h", "--help").set(help, true).doc("Show this page"),
			option("-V", "--version")([](){printf("%s\n", version.c_str());}));

	auto fmt = doc_formatting{}
		.first_column(8)                           //left border column for text body
		.doc_column(30)                            //column where parameter docstring starts
		.last_column(100)                          //right border column for text body
		.indent_size(4)                            //indent of documentation lines for children of a documented group
		.line_spacing(0)                           //number of empty lines after single documentation lines
		.paragraph_spacing(1)                      //number of empty lines before and after paragraphs
		.flag_separator(", ")                      //between flags of the same parameter
		.param_separator(" ")                      //between parameters 
		.group_separator(" ")                      //between groups (in usage)
		.alternative_param_separator("/")          //between alternative flags 
		.alternative_group_separator(" | ")        //between alternative groups 
		.surround_group("(", ")")                  //surround groups with these 
		.surround_alternatives("(", ")")           //surround group of alternatives with these
		.surround_alternative_flags("", "")        //surround alternative flags with these
		.surround_joinable("(", ")")               //surround group of joinable flags with these
		.surround_optional("[", "]")               //surround optional parameters with these
		.surround_repeat("", "...")                //surround repeatable parameters with these
		.empty_label("")                           //used if parameter has no flags and no label
		.max_flags_per_param_in_usage(1)           //max. # of flags per parameter in usage
		.max_flags_per_param_in_doc(32)            //max. # of flags per parameter in detailed documentation
		.split_alternatives(true)                  //split usage into several lines for large alternatives
		.alternatives_min_split_size(3)            //min. # of parameters for separate usage line
		.merge_alternative_flags_with_common_prefix(false)  //-ab(cdxy|xy) instead of -abcdxy|-abxy
		.ignore_newline_chars(false)               //ignore '\n' in docstrings
		;


    if(!parse(argc, const_cast<char **>(argv), cli)) {
		cout << "Usage:\n" << usage_lines(cli, "rafilter", fmt)
     << "\nOptions:\n" << documentation(cli, fmt) << "\nERROR: Required parameter missing\n";
		// throw "Division by zero condition!";
		exit(0);
	}
	if (help){
		cout << "Usage:\n" << usage_lines(cli, "rafilter", fmt)
     << "\nOptions:\n" << documentation(cli, fmt) << '\n';
		// throw "Division by zero condition!";
		return 0;
	}




	if (output_path!="./"){
		string cmd = "mkdir -p "+output_path;
		system(cmd.c_str());
	}
	
	/* 
	 *	Build the kmer pos for certain fasta file
	*/
	if (mode=='b'){
		cout<<"Building kmer pos for input fa\nReading the kmer file\n.........\n";
		if (fopen(kmer_file.c_str(),"r")){
			read_kmer(kmer_file.c_str(),t1);
			cout<<"Saved the kmer in hash table\n";
		}
		else{
			perror("Kmer file not exist!");
			exit(0);
		}
		if (fopen(ref_file.c_str(),"r"))
			build_pos(ref_file.c_str(), output_path, 1, t1, t2);
		if (fopen(query_file.c_str(),"r"))
			build_pos(query_file.c_str(),output_path, 0, t1, t2);
		cout<<"Building pos finished!\n";
	}


	if (mode=='f'){
		if (align != ""){
			// align_fmt will be transmited later.
			cout<<"Thread number : "<< t2
			<< "\nOutput dictionary : "<<output_path
			<<"\nAlignments file : "<< align
			<<"\nAlignment format : "<<align_fmt<<endl;
			read_file (align, rpos, qpos, align_fmt, output_path, threshold, t2);
	}
	}
	return 0;
}
