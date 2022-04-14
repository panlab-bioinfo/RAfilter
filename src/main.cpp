#include "clipp.h"
#include <iostream>
#include "kmerfind4.hpp"
#include "filter.hpp"
#include <thread>
using namespace clipp;using std::cout;using std::string;using std::thread;

int main(int argc, char* argv[]) { 
    string kmer_file = "";
    string ref_file = "";
	string query_file = "";
	string output_path = "./";
	string align = "";
	bool align_fmt = 0;
    int t1=4,t2=8;
    auto cli = (
        option("-t", "--threads").doc("Number of threads") & value("threads", t2),
		option("-o", "--out-put").doc("Diractory of output") & value("out_path", output_path),
		option("-q", "--query").doc("Query/Reads fasta file") & value("query",query_file),
		option("-r", "--reference").doc("Reference fasta file") & value("reffile", ref_file),
		(option("-p").set(align_fmt) | option("-b").set(align_fmt=1)).doc("fmt of alignment file p:paf/b:bam") & value("paf/bam", align),
        value("kmerfile", kmer_file).doc("Kmer file with creating by jellyfish")
    );

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
		.alternative_param_separator("|")          //between alternative flags 
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
		cout << "Usage:\n" << usage_lines(cli, "KAfilter", fmt)
     << "\nOptions:\n" << documentation(cli, fmt) << '\n';
		cout << "\nerror: Too few arguments!\n";
		// throw "Division by zero condition!";
		exit(0);
	};
	// system("mkdir ");
	read_kmer(kmer_file.c_str(),t1);
    cout << "Create dict finished!\n";
	if (ref_file!="")
		build_pos(ref_file.c_str(), output_path, 1, t1, t2);
	if (query_file!="")
		build_pos(query_file.c_str(), output_path, 0, t1, t2);
	cout<<"Building pos finished!\n";
	if (align != ""){
		// align_fmt will be transmited later.
		string refpos=output_path+"/ref.pos";
		string querypos=output_path+"/query.pos";
		read_file(align,refpos,querypos);
	}
	return 0;
}