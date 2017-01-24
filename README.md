# ChopStitch 1.0.0.
####Exon annotation and splice graph reconstruction using transcriptome assembly and whole genome sequencing data
                     
ChopStitch is a new method for finding putative exons and constructing splice graphs using an assembled transcriptome and whole genome shotgun sequencing (WGSS) data. ChopStitch identifies exon-exon boundaries in *de novo* assembled RNA-seq data with the help of a Bloom filter that represents the *k*-mer spectrum of WGSS reads. The algorithm also detects base substitutions in transcript sequences corresponding to sequencing or assembly errors, haplotype variations, or putative RNA editing events. The primary output of our tool is a FASTA file containing putative exons. Further, exon edges are interrogated for alternative exon-exon boundaries to detect transcript isoforms, which are reported as splice graphs in dot output format.

###Requirements:
Install [pip](https://pip.pypa.io/en/latest/installing/) 
         
Install requirements by running:
```
pip install -r requirements.txt
```
Install [Graphviz version 2.4.0](http://www.graphviz.org/Download..php) command line tools
              
              
###Install ChopStitch:
When installing ChopStitch from GitHub source the following tools are required:

* [Autoconf](http://www.gnu.org/software/autoconf)
* [Automake](http://www.gnu.org/software/automake)

To generate the configure script and make files:

```
./autogen.sh
```

To compile and install ChopStitch in /usr/local:

```
$ ./configure
$ make
$ sudo make install
```
To install ChopStitch in a specified directory:

```
$ ./configure --prefix=/opt/ChopStitch
$ make 
$ make install 
```

ChopStitch uses OpenMP for parallelization, which requires a modern compiler such as GCC 4.2 or greater. If you have an older compiler, it is best to upgrade your compiler if possible. If you have multiple versions of GCC installed, you can specify a different compiler:

```
$ ./configure CC=gcc-xx CXX=g++-xx 
```

For the best performance of ChopStitch, pass `-O3` flag:  

```
$ ./configure CFLAGS='-g -O3' CXXFLAGS='-g -O3' 
```

To run ChopStitch, its executables, `CreateBloom` and `FindExons`, should be found in your PATH. If you installed ChopStitch in /opt/ChopStitch, add /opt/ChopStitch/bin to your PATH:

```
$ PATH=/opt/ChopStitch/bin:$PATH
```


###Run CreateBloom

```
Usage: CreateBloom [OPTION]... FILES...
Creates a Bloom filter (BF) to be used for FindExons.
Acceptable file formats: fastq, fasta, sam, bam, gz, bz, zip.

 Options:

  -t, --threads=N  use N parallel threads [1]
  -k, --kmer=N	the length of kmer [50]
  -d, --fpr1=N	primary BF fpr [0.01]
  -s, --fpr2=N	secondary BF fpr [0.01]
  -r, --ref	using FASTA reference as input instead of FASTQ reads. Don't use fpr2 in this case
      --help	display help and exit
      --version	output version information and exit

```
Example:
```
./CreateBloom -t 32 -k 50 --fpr1 0.01 --fpr2 0.01 <FASTQ1> <FASTQ2>
```
OR
```
./CreateBloom --ref -t 32 -k 50 --fpr1 0.01  <REFERENCE FASTA> 
```
               
Output:
            
Bfilter.bf : Bloom filter file
             
Bfilter.inf : Info file required for FindExons 
        
             
###Run FindExons
Find putative exons in TransAbySS Transcriptome assembly file
Acceptable file formats: FASTA
```
  Options:
   -i, --input-bloom=FILE     load bloom filter from FILE
   -l, --leniency=N           leniency for exon-exon juction detection [10]
   -f, --lfactor=N            leniency calculated as ceil(FPR*lfactor*k) 
       --help	                display this help and exit
       --version	            output version information and exit

````
           
Example:
```
./FindExons -i Bfilter.bf <Transcriptome assembly file (TransABySS FASTA file)>
```
Output:
A FASTA file of exons with headers in this format - 
```
>TranscriptName_startCoordinate_Endcoordinate
```         
              
###Run MakeSplicegraph.py with the putative exons FASTA file outputted by FindExons(confident-exons.fa)
    
Example:
```
python MakeSplicegraph.py -i <Putative exon in FASTA format> -o <Splicegraph-outputfile>
```
       
###Run Graphviz ccomps to obtain Splice sub-graphs
    
Example:   
```
ccomps <Splicegraph DOT file from MakeSplicegraph.py> -o <splicegraph_subgraph>
```
         
