
'''
Created by Hamza Khan on 2018-10-17
A script to find connected components in ChopStitch Splicegraph
'''


import argparse
import networkx as nx
import csv

class findsubcom:
 
    def __init__(self,dotfile):
        self.dotfile = dotfile

    def findcc(self):
        ''' 
        Read input DOT file and return a list of 
        connected component sub graphs
        '''
        G = nx.drawing.nx_agraph.read_dot(self.dotfile)
        G_undirected = G.to_undirected()
        cc = nx.connected_components(G_undirected)
        cc_subgraphs = []
        for x in cc:
            H = G.subgraph(x)
            cc_subgraphs.append(H)
        return cc_subgraphs

                   
    def writeoutput(self, cc_subgraphs):
        '''
        Given a list of subcomponent graphs, 
        write a CSV file with GeneID and ExonID
        '''
        with open ("splice_subgraphs.csv", 'w') as fh:
            csv_writer = csv.writer(fh, delimiter=',')
            csv_writer.writerow(['GeneID','ExonID'])
            gene_count = 0
            for G in cc_subgraphs:
                gene_count+=1    
                for node in G.nodes():
                    csv_writer.writerow([gene_count,node])        
                
     
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
                        fh.write("\t\t" + node + " -> " + neighbor + ";\n")
                fh.write("\t}\n")


def parse_args():

    #Parse command line arguments
    parser = argparse.ArgumentParser(
        description = 'Find graph subcomponents and write output')

    #Positional arguments
    parser.add_argument('-g', '--dotfile',
                       default = 'None', required=True,
                       help = 'Graph DOT file from MakeSplicegraph.py')

    parser.add_argument('--writesplicesubgraphs', action="store_true", 
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
    obj.writeoutput(cc_subgraph)
 
    #Check flag and generate DOT output
    if(args.writesplicesubgraphs):
        obj.writedot(cc_subgraph)


if __name__=='__main__':
    main()
