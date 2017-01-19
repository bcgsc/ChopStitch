# ChopStitch 1.0.0.
Exon annotation and splice graph reconstruction using transcriptome assembly and whole genome sequencing data
ChopStitch, a new method for finding putative exons and constructing splice graphs using an assembled transcriptome and whole genome shotgun sequencing (WGSS) data. ChopStitch identifies exon-exon boundaries in *de novo* assembled RNA-seq data with the help of a Bloom filter that represents the *k*-mer spectrum of WGSS reads. The algorithm also detects base substitutions in transcript sequences corresponding to sequencing or assembly errors, haplotype variations, or putative RNA editing events. The primary output of our tool is a FASTA file containing putative exons. Further, exon edges are interrogated for alternative exon-exon boundaries to detect transcript isoforms, which are reported as splice graphs in dot output format.

###Requirements:
Install [pip](https://pip.pypa.io/en/latest/installing/) 
Install requirements by running:
```
pip install -r requirements.txt
```
Install [Graphviz version 2.4.0](http://www.graphviz.org/Download..php) command line tools
              
              
###Install ChopStitch:
``
./configure
make
make install
```
            
###Run CreateBloom

```
Usage: CreateBloom [OPTION]... FILES...
Creates a Bloom filter (BF) to find exon-exon junctions.
Accepatble file formats: fastq, fasta, sam, bam, gz, bz, zip.

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
Accepatble file formats: FASTA

  Options:
    -i, --input-bloom=FILE  load bloom filter from FILE
    -l, --leniency=N        Calculated as ceil(FPR*k*leniency factor) [10]
        --help	            display help and exit
        --version           output version information and exit

Example:
```
./FindExons -i Bfilter.bf <Transcriptome assembly file (TransABySS FASTA file)>
```
Output:
A FASTA file of exons with headers in this format - 
>TranscriptName_startCoordinate_Endcoordinate
         
            
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
         
