#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <cstdlib>
#include <getopt.h>
#include <cassert>
#include <cmath>
//#include "config.h"
#include <utility>
#include <algorithm>
#include <vector>
#include "lib/pstream.h"
#include <cassert>
#include <getopt.h>
#include <iostream>
#include <cstring>
#include <algorithm> 

#include "lib/BloomFilter.hpp"


#define PROGRAM "Chopstitch"

static const char VERSION_MESSAGE[] =
    PROGRAM " Version 1.0.0 \n"
    "Written by Hamza Khan, Hamid Mohamadi.\n"
    "Copyright 2016 Canada's Michael Smith Genome Science Centre\n";

static const char USAGE_MESSAGE[] =
    "Usage: " PROGRAM " [OPTION]... FILES...\n"
    "Find putative exons in TransAbySS Transcriptome assembly file\n"
    "Accepatble file formats: fasta\n"
    "\n"
    " Options:\n"
    "\n"
    "  -k, --kmer=N	              the length of kmer [50]\n"
    "  -i, --input-bloom=FILE     load bloom filter from FILE\n"
    "  -l, --leniency=N           Calculated as ceil(FPR*) [10]\n"   
    "  -X, --internalexons        output internal confident exons in a fasta file \n"
    "  -v, --verbose              display verbose output\n"
    "      --help	              display this help and exit\n"
    "      --version	          output version information and exit\n"
    "\n"
    "Report bugs to hmohamadi@bcgsc.ca; hkhan@bcgsc.ca\n";

using namespace std;

namespace opt {
unsigned leniency=0;
unsigned k=50;
unsigned nhash=1;
static string inputBloomPath;
unsigned ibits=8;   
bool internalexons = false;
unsigned verbose=0;
}

static const char shortopts[] = "k:l:i:h:Xv";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
    { "kmer",	required_argument, NULL, 'k' },
    { "leniency",	required_argument, NULL, 'l' },
    { "input-bloom",	required_argument, NULL, 'i' },
    { "hash",	required_argument, NULL, 'h' },
    { "help",	no_argument, NULL, OPT_HELP },
    { "version",	no_argument, NULL, OPT_VERSION },
    { "internalexons",    no_argument, NULL, 'X' },
	{ "verbose",          no_argument, NULL, 'v' },
    { NULL, 0, NULL, 0 }
};


bool getSeq(std::ifstream &uFile, std::string &line) {
    bool good=false;
    std::string hline;
    line.clear();
    do {
            good=getline(uFile, hline);
            if(hline[0]=='>'&&!line.empty()) break;// !line.empty() for the first rec
            if(hline[0]!='>')line+=hline;
        } while(good);
    
    if(!good&&!line.empty())
            good=true;

    return good;
}



/*

void findJunctions(BloomFilter& bloom, int optind, char** argv)

{
   FastaReader reader(argv[optind], FastaReader::FOLD_CASE);

   BloomFilter& bf = bloom;
 

   std::ofstream bfile("boundaries.csv", std::ios::app);
   std::ofstream efile("exons.fa", std::ios::app);
   std::ofstream cefile("confident_exons.fa", std::ios::app);


   for (FastaRecord rec; reader >> rec;) 
     {
	   
     if(int((rec.seq).length())<(int(opt::k)+5)){continue;}
     
	 unsigned pos = (opt::k)-1,
              myflag=0, 
              start = 0, 
              end = 0, 
              n = 10, 
              min_exon=5, 
              unmatch_count=0, 
              snp_chance_A=0,
              snp_chance_T=0,
              snp_chance_G=0,
              snp_chance_C=0,
              snp_nochance=0;
	 
	 int first_match=0, last_match_pos=0,
	     big1=0, small1=0, big2=0;
	 std::string current_sec = rec.seq;
	 vector<int> arr;
	 vector<int> carr;  
	     
     double FDR = 0.02;
     
     if(opt::leniency==0)
     {
    	 opt::leniency = ceil(n*FDR*(opt::k));
     }
 
     int str_length = current_sec.length();

     for (unsigned int x=0; x < (str_length-((opt::k)-1)-2); x++) 
            {
             //Kmer it = Kmer(current_sec.substr (x,(opt::k)));
    	       string it = current_sec.substr (x,(opt::k));
    	       getCanon(it);
    	       
		       pos=pos+1;
		       if((bf[it]))
		           {
		               last_match_pos=pos;
		           }
		        
              	
		        if((!bf[it]) && (myflag==0) && (first_match==1))
		           {
		             myflag=1;
		             start = pos;
                     unmatch_count=1;
                     snp_chance_A=0, snp_chance_C=0, snp_chance_G=0, snp_chance_T=0, snp_nochance=0;                     
		           }
		        
		         if((!bf[it]) && (myflag==1))
		            {
                    if(unmatch_count<=opt::leniency)
                        {  
                         
                         string it_a = (it.str()).replace(((it.str()).length())-unmatch_count, 1, "A"); 
                         string it_t = (it.str()).replace(((it.str()).length())-unmatch_count, 1, "T");                        
                         string it_g = (it.str()).replace(((it.str()).length())-unmatch_count, 1, "G");                        
                         string it_c = (it.str()).replace(((it.str()).length())-unmatch_count, 1, "C");
                        
                         if(bf[Kmer(it_a)]){snp_chance_A++;unmatch_count+=1;
                                                 goto SNP_found;}

                         else if(bf[Kmer(it_t)]){snp_chance_T++;unmatch_count+=1;
                                                 goto SNP_found;}
                         
                         else if(bf[Kmer(it_g)]){snp_chance_G++;unmatch_count+=1;
                                                 goto SNP_found;}
                         
                         else if(bf[Kmer(it_c)]){snp_chance_C++;unmatch_count+=1;
                                                 goto SNP_found;}
                         
                         else{snp_nochance++;}
                         unmatch_count+=1;
                       }
                      
                      if(unmatch_count == opt::leniency+1)
                        {

                         if((snp_chance_A  >= 2*snp_nochance) && (snp_chance_A > (snp_chance_T+snp_chance_G+snp_chance_C)))
                             {
                              current_sec.replace(pos-(1+opt::leniency),1,"A");

                             }   
                         if((snp_chance_T  >= 2*snp_nochance) && (snp_chance_T > (snp_chance_A+snp_chance_G+snp_chance_C)))
                             {
                              current_sec.replace(pos-(1+opt::leniency),1,"T");

                             }   
                         if((snp_chance_G  >= 2*snp_nochance) && (snp_chance_G > (snp_chance_A+snp_chance_T+snp_chance_C)))
                             {
                              current_sec.replace(pos-(1+opt::leniency),1,"G");

                              }   
                         if((snp_chance_C  >= 2*snp_nochance) && (snp_chance_C > (snp_chance_A+snp_chance_G+snp_chance_T)))
                             {
                              current_sec.replace(pos-(1+opt::leniency),1,"C");
   
                             }   
                          unmatch_count+=1;
                       }
                                                 
		            }                
                SNP_found:

                if((bf[it]) && (myflag==1) && (first_match==1))
                    {
                       unsigned int next_kmer = x;
                       if (!(bf[Kmer(current_sec.substr (++next_kmer,(opt::k)))]))
	                        {  
                             unmatch_count+=1;
                           }
                    }


              unsigned int a = x;                           

		        if((bf[it]) && (myflag==1) &&  (bf[Kmer(current_sec.substr (++a,(opt::k)))]) && (first_match==1))
	               {  
                       
                        if((bf[Kmer(current_sec.substr (++a,(opt::k)))]))
                            {
                               
			                      myflag = 0;
			                      end = pos;
			                    
                               if(end-start > ((opt::k)+min_exon))
                                    {     
                                          int minus_counter = 0, snp_chance_A=0, snp_chance_C=0, snp_chance_G=0, snp_chance_T=0, snp_nochance=0;
                                          for(int a_minus=(x-1); a_minus>=(x-opt::leniency); a_minus-- )
                                          {
                                           
                                          string it_rev = Kmer(current_sec.substr (a_minus,(opt::k))).str();                                                                       
                                          string it_a = (it_rev).replace((minus_counter), 1, "A");
                                          string it_t = (it_rev).replace((minus_counter), 1, "T");                                      
                                          string it_g = (it_rev).replace((minus_counter), 1, "G");                                          
                                          string it_c = (it_rev).replace((minus_counter), 1, "C");
                                         
                                          minus_counter = minus_counter + 1;

                                         if(bf[Kmer(it_a)]){snp_chance_A++;}

                                         else if(bf[Kmer(it_t)]){snp_chance_T++;}
                         
                                         else if(bf[Kmer(it_g)]){snp_chance_G++;}
                         
                                         else if(bf[Kmer(it_c)]){snp_chance_C++;}
                         
                                         else{snp_nochance++;}

                                         }

                                         if(snp_chance_A  >= 2*snp_nochance || snp_chance_T  >= 2*snp_nochance || snp_chance_G  >= 2*snp_nochance || snp_chance_C  >= 2*snp_nochance )
                                             {
                                          
                                           bfile << "," << start-1 << "\n" << (rec.id) << ","<< start ;
                                           arr.push_back(start-1);
                                           arr.push_back(start);    
                                             }
                                         else 
                                             {
                                       
                                           bfile << "," << start-1 << "\n" << (rec.id) <<","<< start << "," << (end-opt::k)-1<< "\n" << (rec.id) <<","<< end-opt::k ;
                                           arr.push_back(start-1);
                                           arr.push_back(start);
                                           arr.push_back((end-opt::k)-1);
                                           arr.push_back(end-opt::k);
                                             }     
                                        }
                    
			                      if((end - start >= (((opt::k)-opt::leniency)-2)) && ((end-start)<= ((opt::k)+min_exon)))
                                     {
			                                 
			                                  bfile << "," << start-1 << "\n" << (rec.id) <<","<< end-opt::k ;
                                           arr.push_back(start-1);
                                           arr.push_back(end-opt::k);
			                            }
                             }
		             }

              //Testing the first matching kmer in the bloom filter
              if((bf[it]) && first_match==0)
		              {
                          unsigned int next_it = x;
                          //Added extra checks to avoid False positives
                          if(bf[Kmer(current_sec.substr (++next_it,(opt::k)))] && bf[Kmer(current_sec.substr (++next_it,(opt::k)))])
                            { 
                             bfile << (rec.id)<< ",";
	    		                 bfile << pos-((opt::k)-1);
                             arr.push_back(pos-((opt::k)-1));
			                    first_match=1;
                            }
		              }
        
		    }  
      if(first_match==1)
         {
		      	bfile << "," << last_match_pos << "\n";
               arr.push_back(last_match_pos);
         } 

     

  if(opt::internalexons)
         {
             int check_flag=0;
             carr=arr;

             if(!carr.empty() && carr.size()==2)
                 {
                   if (*(carr.begin())!=1 && *(carr.end()-1) < str_length-2)
                      {
                      string single_read = rec.seq.substr(((*(carr.begin()))-1),((*(carr.end()-1)-(*(carr.begin()))+1)));
                      cefile << ">" << rec.id << "_" << *(carr.begin())<<"_"<< *(carr.end()-1)<< std::endl;
                      cefile << single_read<<endl; 
                      std::flush(cefile);             
                      }              
                  }

            if(!carr.empty() && (carr.size())>2) 
                 {   
                  int beGIN = *(carr.begin());
                  int enD = *(carr.end()-1);
                  check_flag=1;

                  if(beGIN == 1)
                     {
                       carr.erase(carr.begin(),carr.begin()+2);
                     }
                 
                  if(enD >= str_length-2)
                     {
                       
                       carr.erase (carr.end()-2,carr.end());
                     }
                  }

           
             if(!carr.empty() && carr.size()==2 && check_flag==1)
                 {

                      string single_read = rec.seq.substr(((*(carr.begin()))-1),((*(carr.end()-1)-(*(carr.begin()))+1)));
                      cefile << ">" << rec.id << "_" << *(carr.begin())<<"_"<< *(carr.end()-1)<< std::endl;
                      cefile << single_read<<endl; 
                      std::flush(cefile);             
                      check_flag=0;             
                  }


            if(!carr.empty() && (carr.size())>2) 
                 { 
                     int beGIN = *(carr.begin());
                     int enD = *(carr.end()-1);
                     int Big1=0, Small1=0, Big2=0;
                     string Read1, Read2;

                     for (std::vector<int>::const_iterator w = (carr.begin())+1; w != (carr.end())-1; w+=2)
                         {
                         if(*w < *(w+1))
                             {
                                  Big1 = *(w+1);
                                  Small1 = *w;
                             }
                         else
                             {
                                  Small1 = *(w+1);
                                  Big1 = *w;   
                             }
              
                         if(w == (carr.begin())+1) 
                             {
                                 Read1 = rec.seq.substr(beGIN-1,(Small1-beGIN)+1);
                                 cefile << ">" << rec.id << "_" << beGIN<<"_"<< Small1<< std::endl;
                                 cefile << Read1<<endl;
                                 std::flush(cefile); 
                                 Big2 = Big1; 
                             } 
                         else 
                             {
                                 Read2 = rec.seq.substr(Big2-1,(Small1-Big2)+1);
                                 cefile << ">" << rec.id << "_" << Big2<<"_"<< Small1<< std::endl;
                                 cefile << Read2<<endl;
                                 std::flush(cefile); 
                                 Big2 = Big1;  
                             }       
                         }
        
                         if(!carr.empty() && (carr.size())>2) 
                             {                    
                                 string End_Read = rec.seq.substr(Big2-1,(enD-Big2)+1);
                                 cefile << ">" << rec.id << "_" << Big2<<"_"<< enD<< std::endl;
                                 cefile << End_Read<<endl;
                                 std::flush(cefile);                   
                             } 
                  }


         }


    if(!arr.empty() && arr.size()==2)
        {
             
               string single_read = rec.seq.substr(((*(arr.begin()))-1),((*(arr.end()-1)-(*(arr.begin()))+1)));
               efile << ">" << rec.id << "_" << *(arr.begin())<<"_"<< *(arr.end()-1)<< std::endl;
               efile << single_read<<endl; 
               std::flush(efile);             
              
        }


     if(!arr.empty() && (arr.size())>2) 
        {   
              int bEgin = *(arr.begin());
              int eNd = *(arr.end()-1);

              for (std::vector<int>::const_iterator q = (arr.begin())+1; q != (arr.end())-1; q+=2)
                 {
                   if(*q < *(q+1))
                      {
                        big1 = *(q+1);
                        small1 = *q;
 
                      }
                  else
                      {
                         small1 = *(q+1);
                         big1 = *q;
                        
                      }
              
                  if(q == (arr.begin())+1) 
                      {
                           read1 = rec.seq.substr(bEgin-1,(small1-bEgin)+1);
                           efile << ">" << rec.id << "_" << bEgin<<"_"<< small1<< std::endl;
                           efile << read1<<endl;
                           std::flush(efile); 
                           big2 = big1; 
                           
                      } 
                 else 
                      {
                          
                            read2 = rec.seq.substr(big2-1,(small1-big2)+1);
                            efile << ">" << rec.id << "_" << big2<<"_"<< small1<< std::endl;
                            efile << read2<<endl;
                            std::flush(efile); 
                            big2 = big1;  

                   } 

              
           }
        
         if(!arr.empty() && (arr.size())>2) 
               {                    
                     string end_read = rec.seq.substr(big2-1,(eNd-big2)+1);
                     efile << ">" << rec.id << "_" << big2<<"_"<< eNd<< std::endl;
                     efile << end_read<<endl;
                     std::flush(efile); 
                   
               } 


  }
}

 efile.close();     

}	

*/

int main(int argc, char** argv) {

    bool die = false;
    for (int c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) {
        case '?':
            die = true;
            break;
        case 'k':
            arg >> opt::k;
            break;
        case 'l':
            arg >> opt::leniency;
            break;
        case 'i':
         	arg >> opt::inputBloomPath; 
       	    break;
        case 'h':
            arg >> opt::nhash;
            break;
        case 'b':
            arg >> opt::ibits;
		case 'X':
		    	opt::internalexons = true; break;  
		case 'v':
			opt::verbose++; break;
        case OPT_HELP:
            std::cerr << USAGE_MESSAGE;
            exit(EXIT_SUCCESS);
        case OPT_VERSION:
            std::cerr << VERSION_MESSAGE;
            exit(EXIT_SUCCESS);
        }
        if (optarg != NULL && !arg.eof()) {
            std::cerr << PROGRAM ": invalid option: `-"
                      << (char)c << optarg << "'\n";
            exit(EXIT_FAILURE);
        }
    }
    if (argc - optind < 0) {
        std::cerr << PROGRAM ": missing arguments\n";
        die = true;
    }
    if (die) {
        std::cerr << "Try `" << PROGRAM << " --help' for more information.\n";
        exit(EXIT_FAILURE);
    }

    
	if (!opt::inputBloomPath.empty()) {

		if (opt::verbose)
			std::cerr << "Loading Bloom filter from `"
				<< opt::inputBloomPath << "'...\n";
		
		const char* inputPath = opt::inputBloomPath.c_str();
		cerr << inputPath <<"\n";
		
		size_t dbfSize=3000000000,sbfSize=2000000000;
		BloomFilter bloom(sbfSize*opt::ibits, opt::nhash, opt::k,inputPath);
		//BloomFilter bloom(inputPath);
		cerr << bloom.getPop() << "\n";
        //string it = "ATCGCTGATGATCGCTGATGATCGCTGATGATCGCTGATGATCGCTGATG";
		//string it = "CTCTTCTTGCTCAAAGTATTGTTATGCTCATCTGTATGATTTTGATGCTG";
		//string it = "AAATAAATGCTAAATTTTCTGGCCTGATTTAATGTAGAAAAATAAAATCT";
		  string it = "GTGATGTCTGCATTCAAGTCACAGAGTTGAACATTGCCTTTCATAGAGCA";
		  if(bloom.contains(it.c_str()))
			  cerr << "kmer is in filter\n";
		  else
			  cerr << "kmer NOT in filter\n";
		//bloom = new BloomFilter(sbfSize*opt::ibits, opt::nhash, opt::k,inputPath);
       // findJunctions(*bloom, int optind, char** argv);

    }
	else
		cerr << "No input file\n";	
	
}
