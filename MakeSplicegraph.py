#!/usr/bin/env python
"""
Create a Splice graph from confident_exons.fa 
"""

import argparse, datetime
import sys, getopt
from pprint import pprint
import networkx 
from networkx.algorithms.components.connected import connected_components


ts = datetime.datetime.now()

__title__ = 'MakeSplicegraph:Splice graph Construction'
__version__ = '0'
__description__ = "A tool to easily create a splice graph from confident-exons.fa"
__author__ = 'Hamza Khan'
__license__ = 'GPL license'
__author_email__ = "hkhan@bcgsc.ca"
epi = "Licence: %s by %s <%s>\n\n" % (__license__,
__author__,
__author_email__)
__doc__ = "***************************************************************\
          \n %s v%s - %s \n************************************************\
***************" % (__title__,
__version__,
__description__)


kmer_dict = {}
kmer_dict_rt = {}
directed_graph = {}
seq_dict = {}
new_nodes = {}
seqlen_dict={}
first_last_exon_dict={}


def to_graph(l):
    ''' Including networkx connected components alogorithm
    '''
    G = networkx.Graph()
    for part in l:
        # each sublist is a bunch of nodes
        G.add_nodes_from(part)
        # it also implies a number of edges:
        G.add_edges_from(to_edges(part))
    return list(connected_components(G))


def to_edges(l):
    """ 
        treat `l` as a Graph and returns it's edges 
        to_edges(['a','b','c','d']) -> [(a,b), (b,c),(c,d)]
    """
    it = iter(l)
    last = next(it)

    for current in it:
        yield last, current
        last = current    


def get_next_line(reader):
  """(file open for reading) -> string or None
  
  Read the next non-blank line from a file and return it as a string after
  stripping all leading and trailing whitespace.  A "non-blank line" is a line
  which contains some character(s) other than whitespace.  None is returned
  if there are no more lines to read.
  """

  # Read lines from the file until a line is found which is not empty.
  line = ''
  while not line:
    # Read the next line.
    line = reader.readline()
    # If an empty string is returned, there are no more lines.
    if not line:
      return None
    else:
      # Otherwise remove whitespace from the line.  If the line contains only 
      # whitespace (eg: it would appear blank in a text editor), this operation
      # will result in an empty string.
      line = line.strip()

  return line


def split_line(line, sep = ','):
  """(string, string) -> list:

  Break a string into a list at some separator character(s).  Resulting elements
  are stripped of leading and trailing whitespace and then returned in a list
  of strings.
  """
  
  # Break line into elements at the separator character.
  elements = line.split(sep)
  
  # Remove any leading or trailing whitespace.
  for i in range(len(elements)):
    elements[i] = elements[i].strip()
    
  return elements
  


def spew_kmers(name,st, k, n):
  """ (string, int, int) -> kmer_list"""  

  if(len(st)>=(k+n-1)):  

     Kmers=[]
     for i in xrange(0, n): Kmers.append(st[i:i+k])
  
     for kmer in Kmers:
        if kmer in kmer_dict:
           temp = []
           temp.append(name)
           kmer_dict[kmer].extend(temp)
        else:
            temp=[]
            temp.append(name)
            kmer_dict[kmer] = temp


def spew_kmers_rt(name,st, k, n):
  """ (string, int, int) -> kmer_list"""
 
  if(len(st)>=(k+n-1)): 

     Kmers=[]
     for i in range(len(st),len(st)-n, -1): Kmers.append(st[(i-k):(i)])
   
     for kmer in Kmers:
        if kmer in kmer_dict_rt:
           temp = []
           temp.append(name)
           kmer_dict_rt[kmer].extend(temp)
        else:
           temp=[]
           temp.append(name)
           kmer_dict_rt[kmer] = temp



def read_fasta_records(inputfile):
   """(string) -> data_list_of_lists"""

     
   with open(inputfile, 'r') as reader:

     next_line = get_next_line(reader)
    
     previous_name = ''         
     previous_start = ''
     previous_stop = ''
     similar_header = 0
 
     while next_line:

       if str(next_line).startswith(">"):
          
            #print "\nnext_line = ", next_line
            header = next_line.split('_')
            #print "\nheader = ", header
            name = str([(''.join(header[:-2]))])[3:-2]
            #print "\nname = ", name
            start = str(header[-2])
            #print "\nstart = ", start
            stop = str(header[-1])
            #print "\nstop = ", stop
            
            #Constructing the sequence dictionary(key->fasta header, value->sequence)

            identifier = (name+"_"+start+"_"+stop)
            seq = get_next_line(reader)
            seq_dict[identifier]= seq.rstrip()   
            #print seq_dict
            seqlen_dict[identifier]=len(seq.rstrip())
            #print seqlen_dict

            #print "name = ", name, "\nstart = ", start, "\n stop = ", stop, "\n"
            #print "\nSimilar Header = ", similar_header, ", previous name = ", previous_name    
                

            if(name==previous_name and similar_header == 0):

               first_last_exon_dict[(previous_name+"_"+previous_start+"_"+previous_stop)]=1
              
               #print (previous_name+"_"+previous_start+"_"+previous_stop)
               directed_graph[(previous_name+"_"+previous_start+"_"+previous_stop)] = (name+"_"+start+"_"+stop)

               # Ignoring the starting k-mers of the first exon in a transcript. Uncomment this line to include
               #spew_kmers((previous_name+"_"+previous_start+"_"+previous_stop) , seq_dict[(previous_name+"_"+previous_start+"_"+previous_stop)], 30, 10)
               spew_kmers_rt((previous_name+"_"+previous_start+"_"+previous_stop) , seq_dict[(previous_name+"_"+previous_start+"_"+previous_stop)], 30, 10)
               spew_kmers((name+"_"+start+"_"+stop) , seq_dict[(name+"_"+start+"_"+stop)], 30, 10)

               last_pos = reader.tell()
               next_header_line = reader.readline()
               #print "\nsim_next_header_line = ", next_header_line
               next_header = next_header_line.split('_')
               #print ", sim_next_header = ", next_header
               next_name = str([(''.join(next_header[:-2]))])[3:-2]
               #print ", sim_next_name = ", next_name
               
               if(next_name == name):
                      
                   spew_kmers_rt((name+"_"+start+"_"+stop) , seq_dict[(name+"_"+start+"_"+stop)], 30, 10)
                   first_last_exon_dict[(name+"_"+start+"_"+stop)]=0 

               else:
                  
                  first_last_exon_dict[(name+"_"+start+"_"+stop)]=1 
              
               reader.seek(last_pos)

               similar_header = 1


            #This statement deletes the le        
            #elif(name!=previous_name and similar_header == 1):


            elif(name==previous_name and similar_header == 1):
            
               #print (previous_name+"_"+previous_start+"_"+previous_stop)
               directed_graph[(previous_name+"_"+previous_start+"_"+previous_stop)] = (name+"_"+start+"_"+stop)
               spew_kmers((name+"_"+start+"_"+stop) , seq_dict[(name+"_"+start+"_"+stop)], 30, 10) 

               last_pos = reader.tell()
               next_header_line = reader.readline()
               #print "\n\nnext_header_line = ", next_header_line
               next_header = next_header_line.split('_')
               #print ",next_header = ", next_header
               next_name = str([(''.join(next_header[:-2]))])[3:-2]
               #print ", next_name = ", next_name
               
               if(next_name == name):
   
                   spew_kmers_rt((name+"_"+start+"_"+stop) , seq_dict[(name+"_"+start+"_"+stop)], 30, 10)  
                   first_last_exon_dict[(name+"_"+start+"_"+stop)]=0 


               else:
                  
                  first_last_exon_dict[(name+"_"+start+"_"+stop)]=1 
                          

               reader.seek(last_pos)

            else:
               similar_header = 0

            previous_name = str([(''.join(header[:-2]))])[3:-2]          
            previous_start = str(header[-2])
            previous_stop = str(header[-1])      


       next_line = get_next_line(reader)
   #print("\n\nDirected Graph = \n") 
   #pprint(directed_graph)
   #print("\n\n kmer_dict = \n")
   #pprint(kmer_dict)
   #print("\n\n kmer_dict_rt = \n")
   #pprint(kmer_dict_rt)
        
     

def data_lists(inputfile, outputfile):

     data_list =read_fasta_records(inputfile)
     #print data_list[0:1]



def kmer_dict_to_list(kmer_dict):
    
   nodes = []
   for key in kmer_dict:
       if(type(kmer_dict[key]) is list and len(set(kmer_dict[key]))>1):
           #list(set(kmer_dict[key]))
           #print (list(set(kmer_dict[key])), kmer_dict[key])
           #print  key,kmer_dict[key]
           nodes.append(kmer_dict[key])
   return(nodes)




def new_node_names(g, cc):

   #print "Directed Graph\n\n", g
   #print "\n\ncc = ", cc

   for item in cc:

       max_length=0
       biggest_exon=''
       for elements in item:
          if(int(seqlen_dict[elements]) > max_length):
             max_length = int(seqlen_dict[elements])
             biggest_exon = elements

       dissimilar_exons_flag = 0
       for elements in item:
           if int(seqlen_dict[elements]) not in range(max_length-5,max_length+5):
               dissimilar_exons_flag =1 
               break

       if(dissimilar_exons_flag==1):
          new_node=''
          if(first_last_exon_dict[biggest_exon]==0):
                 new_node = new_node+"_OR_"+str(biggest_exon)

          for x in item:
              if(first_last_exon_dict[x]==1):
                 new_node = new_node+"_OR_"+str(x)                
          
          for x in item:
              if(first_last_exon_dict[x]==1):
                 new_nodes[x] = new_node[4:]     
     
          new_nodes[biggest_exon] = new_node[4:]

       if(dissimilar_exons_flag==0):
          new_node=''
          for x in item: 
             new_node = new_node+"_OR_"+str(x)
          for x in item:
             new_nodes[x] = new_node[4:] 
   #print "New_node = ", new_nodes
   return new_nodes



def write_dot(directed_graph, node_names_dict, outputfile):
      
   with open(outputfile, 'w') as writer:
      writer.write("digraph graphname {\n")

      for key in directed_graph:
          if(key in node_names_dict):
             writer.write(str(node_names_dict[key]))
          else:
             writer.write(key)
       
          writer.write("->")      

          if(directed_graph[key] in node_names_dict):
             writer.write(str(node_names_dict[(directed_graph[key])]))
          else:
             writer.write(str(directed_graph[key]))
       
          writer.write(";\n")  
     
      writer.write("\n}")


def main(argv):
   if(len(argv)==0):
        print '\nERROR!:No input provided\n\nUsage: python splicegraph.py -i <inputfile> -o <outputfile>'
        sys.exit(2)
   inputfile = ''
   outputfile = 'splicegraph'
   try:
      opts, args = getopt.getopt(argv,"hi:o:",["ifile=","ofile="])
   except getopt.GetoptError:
      print __doc__,"\n",epi,'Usage: python splicegraph.py -i <inputfile> -o <outputfile>'
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print __doc__,"\n",epi,'Usage: python splicegraph.py -i <inputfile> -o <outputfile>'
         sys.exit()
      elif opt in ("-i", "--ifile"):
         inputfile = arg
      elif opt in ("-o", "--ofile"):
         outputfile = arg
   print 'Input file is "', inputfile
   print 'Output file is "', outputfile
   data_lists(inputfile, outputfile)
   nodes_list_left = kmer_dict_to_list(kmer_dict)
   nodes_list_right = kmer_dict_to_list(kmer_dict_rt)
   #print nodes_list_left,"\n\n",nodes_list_right
   total_nodes = []
   total_nodes.extend(nodes_list_left)
   total_nodes.extend(nodes_list_right)
   cc = to_graph(total_nodes)


   node_names_dict = new_node_names(directed_graph, cc)
   write_dot(directed_graph, node_names_dict, outputfile)
   tf = datetime.datetime.now()
   print "\n Time required - ",tf-ts

   #pprint(first_last_exon_dict)

if __name__ == "__main__":
   main(sys.argv[1:])

