# <div align=center>RAFILTER</div>

## Introduction

**rafilter** is a cli software of filter to obtain the correct alignments and it can only use on linux(HPC) because of the cost of resource consumption. The program uses ordered kemr list to evaluate alignments and decides if them should be filtered.

## Install

**rafilter** requests [**jellyfish**](https://github.com/gmarcais/Jellyfish) to get the kmers of sequence and [htslib](https://github.com/samtools/htslib) to handle .BAM file. We consider integrating  feature of jellyfish in **rafilter** in the future. The htslib have be put in our workflow, so you don't need additional installation. But **jellyfish** is needed to install by yourself.

1. The installation method of jullyfish is following <https://github.com/gmarcais/Jellyfish>.
   Please see <https://github.com/samtools/htslib> about htslib.
2. Get and install rafilter  
   `git clone https://github.com/panlab-bioinfo/RAfilter.git`  
   `cd RAfilter/src`  
   `make`
  
Note: If the following error: "rafilter: error while loading shared libraries: libhts.so.3: cannot open shared object file: No such file or directory" occurs on runing rafilter, please add the htslib path to the `LD_LIBRARY_PATH`  

## Usage

The workflow of rafilter is a two-step process. The first is build kmer pos library and the second step is filter alignments.  

```shell
 $ rafilter
Usage:
        rafilter build [-t <threads>] [-o <out_path>] [-q <query>] [-r <reffile>] <kmerfile>
        rafilter filter [-o <out_path>] [-p/-b] [--threshold <threshold>] <r_pos> <q_pos> <paf/bam>
        rafilter -h
        rafilter -V
Options:
        build:
            -t, --threads     Number of threads
            -o, --out-put     Diractory of output
            -q, --query       Query/Reads fasta file
            -r, --reference   Reference fasta file
            <kmerfile>        Kmer file with creating by jellyfish

        filter:
            -o, --out-put     Diractory of output
            -p/-b             Format of alignment file [p:paf/b:bam]
            --threshold       Threshold of KMAPQ to filter
            <r_pos>           Ref pos file built with build
            <q_pos>           Query pos file built with build
            <paf/bam>         Alignment file

        -h, --help            Show this page
        -V, --version         Version
```

### **Input file**

**rafilter** requests the result of **jellyfish** and the format as follow. On the workflow, kmer length must be 21. The detail see https://github.com/gmarcais/Jellyfish/blob/master/doc/Readme.md

```shell{}
jellyfish count -m 21 -s 1G -t 16 -C reference.fasta

# for HiFi reads alignments, we recommand use unique kmer.
jellyfish dump -c -U 1 mer_counts.jf > unique.kmer   # HiFi alignment

# For ONT reads alignments, we recommand use the rare kmer with the frequency Lower than 4.
jellyfish dump -c -U 3 mer_counts.jf > rare.kmer   # ONT alignment
```

NOTE: The parameter -c is necessary for rafilter because of input format.  

Fasta format with one line ACGT sequence for one only be support for input sequence by **rafilter bulid -q/-r**.  
For example:

```shell{.line_numbers}
>S1_1
CTAGCTCCAGTCCCACCCCGGCCTGCAGAGTGGCTGGGCTGCAGGCATACCCCA
>S1_2
TGCCACAGCGGAGCTTGGATGAGCAAAAGAGGAAGTGGAGCATCTGAACTCTT
>S1_3
GATCTCGAACCTACTCATCTTGTGTAACAAAACTTTATACCCTTTGAACAGTCACC
>S1_4
AATGTTTTTAAAGTGGCCATACTGCCCAAAGCAGTTTATAGATTCAATGCTATTCCT
...
```

And you can get the format file using the follow command.

```shell{}
sed -n '1~4s/^@/>/p;2~4p' reads.fastq > query.fa  # For fastq


# fasta to fasta with single line sequence
awk '/^>/{print n $1; n = "\n"} !/^>/{printf "%s",$0}' test.fasta > reference.fa  # For fasta

```

### **Build library**  

As follow command illustration, the -t is number of threads and its value should be lower than 54. At least, it is a nessary parameter to one or both of -p and -r. If the \<out_path\> provided is not exist, it will be created.  Otherwise, it better be empty. if it doesn't be provided, it will be ./. The kmerfile is a necessary parameter obtained in the previous step by jellyfish.  

```shell{}
# Example for build
rafilter build -t 32 -q query.fa -r reference.fa -o test/ unique.kmer
```

The kmer poses library file will get in this step named `ref.pos and query.pos` depend on the kmers of reference and as input in filter step. PLease make sure your storage space is enough to store the all k-mer pos.
> ref.pos and quer.pos both are generated by the rare k-mer of reference sequnce. 

### **Filter**

In this step, alignment file and pos file creating in [**build**](#build-library) and format of alignment can be .paf or .bam.  

```shell{}
# For example
rafilter filter -o test/ --threshold 12 -p ref.pos query.pos alignment.paf
```
This step will generate the last filtered result using orgin format and a simple filter report. The sample data and scripts are avaliable from [ruoyu1123/rafilter_test](https://github.com/ruoyu1123/rafilter_test).
> We recommend apply PAF format of minimap2 as the input file.

## Declaration and cite

This project uses MIT open source protocol. Please click [here]() to see the detail. If you will use the program to complete your study, please cite the paper with:[Yang, J. et al. RAfilter: an algorithm for detecting and filtering false-positive alignments in repetitive genomic regions. Horticulture Research 10, uhac288 (2023)](https://academic.oup.com/hr/article/doi/10.1093/hr/uhac288/6965242).  
