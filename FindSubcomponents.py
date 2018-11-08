
'''
A script to: 
1) Find connected components in ChopStitch Splicegraph
2) Create a transcript to gene map tsv file 
3) Create a putative exon to gene map tsv file
'''

import argparse
import networkx as nx
from collections import defaultdict
import csv

class findsubcom:
 
    def __init__(self,dotfile):
        self.dotfile = dotfile
        self.nodes_dict={}

    def findcc(self):
        ''' 
        Read input DOT file and return a list of 
        connected component sub graphs
        '''
         
        self.nodes_dict = {}
        self.segment_dict={}
        node_id = 0
        G = nx.MultiDiGraph()
        with open (self.dotfile) as df:
            for line in df:
                if "->" in line:
                    line = line.replace(";","").replace("}", "").replace("{","").replace("\"","")
                    line = line.strip().split("->")

                    #Assuming the user only gives the DOT file output from MakeSplicegraph.py
                    #which has only two nodes per line
                    node1 = line[0]
                    node2 = line[1]
     
                    if node1 not in self.segment_dict:
                         node_id+=1
                         self.nodes_dict[node_id]=node1
                         self.segment_dict[node1] = node_id
                         G.add_node(node_id, segments=node1)

                    if node2 not in self.segment_dict:
                         node_id+=1
                         self.nodes_dict[node_id]=node2
                         self.segment_dict[node2] = node_id
                         G.add_node(node_id, segments=node1)

                    G.add_edge(self.segment_dict[node1], self.segment_dict[node2])
                    
        G_undirected = G.to_undirected()
        cc = nx.connected_components(G_undirected)
        cc_subgraphs = []

        for x in cc:
            H = G.subgraph(x)
            cc_subgraphs.append(H)
        return cc_subgraphs

                   
    def writeoutput(self, cc_subgraphs, geneMap):
        '''
        Given a list of subcomponent graphs, 
        write a CSV file with GeneID and ExonID
        '''
         
        genemap=defaultdict(set)
        with open ('splice_subgraphs.csv', 'w') as fh:
            csv_writer = csv.writer(fh, delimiter=',')
            csv_writer.writerow(['ExonID','GeneID'])
            gene_count = 0
            for G in cc_subgraphs:
                gene_count+=1    
                for node in G.nodes():
                    csv_writer.writerow([self.nodes_dict[node], gene_count])                    
                    for transcript in self.nodes_dict[node].split("_OR_"):
                        genemap[gene_count].add("_".join(transcript.split("_")[0:-2]))

        #Check if the user wants the geneMap output as well
        if(geneMap):
            with open('geneMap.tsv', 'w') as gm:
                tsv_writer = csv.writer(gm, delimiter='\t')
                tsv_writer.writerow(['TranscriptID','GeneID'])   
                for gene in genemap:
                    for transcript in genemap[gene]:
                        tsv_writer.writerow([transcript, gene])    

     
    def writedot(self, cc_subgraphs):
        '''
        Given a list of subcomponents graphs, 
        write an output DOT file
        '''
        with open ("splice_subgraphs.dot", 'w') as fh:
            fh.write("digraph graphname {\n")
            for G in xrange(0, len(cc_subgraphs)):
                fh.write("\tsubgraph graphname_cc_"+ str(G) +" { \n")
                for node in cc_subgraphs[G].nodes():
                    for neighbor in cc_subgraphs[G].neighbors(node):
                        fh.write("\t\t" + self.nodes_dict[node] + " -> " + self.nodes_dict[neighbor] + ";\n")
                fh.write("\t}\n")


def parse_args():

    #Parse command line arguments
    parser = argparse.ArgumentParser(
        description = 'Find graph subcomponents and write output')

    #Positional arguments
    parser.add_argument('-g', '--dotfile',
                       default = 'None', required=True,
                       help = 'Graph DOT file from MakeSplicegraph.py')

    parser.add_argument('-m', '--geneMap', action="store_true",
                       default=False, help = 'Write a file with mappings of  transcripts to genes')

    parser.add_argument('-w', '--writesplicesubgraphs', action="store_true", 
                       default=False, help = 'Write splice subgraphs to DOT file')

    args = parser.parse_args()
    return args


def main():
    
    #Parse arguments
    args = parse_args()
    
    #Create an object
    obj = findsubcom (args.dotfile)
    
    #Find ccomponents
    cc_subgraph = obj.findcc()
    
    #Write CSV output
    obj.writeoutput(cc_subgraph, args.geneMap)
 
    #Check flag and generate DOT output
    if(args.writesplicesubgraphs):
        obj.writedot(cc_subgraph)

    

if __name__=='__main__':
    main()
