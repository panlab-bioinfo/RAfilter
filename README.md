# <center>kmfilter

## Introduction

**Kmfilter** is a filter cli software to obtain the correct alignments and it can only use on linux because of the cost of resource consumption. The program uses  ordered kemr list to evaluate alignments and decides if them should be filtered.

## Install

**kmfilter** requests [**jellyfish**](https://github.com/gmarcais/Jellyfish) to get the kmers of sequence and [htslib](https://github.com/samtools/htslib) to handle .BAM file. We consider integrating  feature of jellyfish in **kmfilter** in the future. For now, one jellyfish available is placed in our codebase.

1. The installation method of jullyfish is following <https://github.com/gmarcais/Jellyfish>.
   Please see <https://github.com/samtools/htslib> about htslib.
2. Get and install kmfilter  
   `git clone https://github.com/ruoyu1123/kmfilter.git`  
   `cd kmfilter/src`  
   `make`
3. Add the path to environment variables.  
   `echo export PATH = $PATH:$install_PATH/kmfilter/src/`  

Note: If the following error: "kmfilter: error while loading shared libraries: libhts.so.3: cannot open shared object file: No such file or directory" occurs on runing kmfilter, please add the htslib path to the `LD_LIBRARY_PATH`  

## Run

The workflow of kmfilter is a two-step process. The first is build kmer pos library and the second step is filter alignments.  

### **Input file**

**kmfilter** requests the result of **jellyfish** and the format as follow. On the workflow, kmer length must be 21. The detail see https://github.com/gmarcais/Jellyfish/blob/master/doc/Readme.md

```shell{}
jellyfish count -m 21 -s 1G -t 16 -C reads.fasta

# for HiFi reads alignments, we recommand use unique kmer.
jellyfish dump -c -U 1 mer_counts.jf > unique.kmer   # HiFi reads

# For ONT reads alignments, we recommand use the rare kmer with the frequency Lower than 4.
jellyfish dump -c -U 3 mer_counts.jf > rare.kmer   # ONT reads
```

NOTE: The parameter -c is necessary for kmfilter because of input format.  

Fasta format with one line ACGT sequence for one only be support for input sequence by **kmfilter bulid -q/-r**.  
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
sed -n '1~4s/^@/>/p;2~4p' test.fastq > test.fa  # For fastq
# fasta 转单行fasta



```


### **Build library**  

1. As follow command illustration, the -t is number of threads and its value should be lower than 54. At least, it is a nessary parameter to -p/-r. If the \<out_path\> provided is not exist, it will be created.  Otherwise, it better be empty. if it doesn't be provided, it will be ./. 

```shell{.line-numbers}
 $ kmfilter
Usage:
        Kmfilter build [-t <threads>] [-o <out_path>] [-q <query>] [-r <reffile>] <kmerfile>
        Kmfilter filter [-o <out_path>] [-p/-b] <r_pos> <q_pos> <paf/bam>
        Kmfilter -h
        Kmfilter -V
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
            <r_pos>           Ref pos file built with build
            <q_pos>           Query pos file built with build
            <paf/bam>         Alignment file

        -h, --help            Show this page
        -V, --version         Version
```
