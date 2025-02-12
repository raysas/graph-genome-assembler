import networkx as nx
import pandas as pd
from IPython.display import display
import matplotlib.pyplot as plt

from utils import *

class HAGraph(nx.DiGraph):
    '''
    Hamiltionian Assembly Graph class
        inherits from networkx.DiGraph (directed graph)
        takes a list of kmers while instantiating and transforms each into a Kmer object
    '''

    def __establish_nodes_links(self):
        # 1. vertices
        self.add_nodes_from(self.kmers)
        # 2. edges (O(n2))
        for kmer_1 in self.nodes:
            for kmer_2 in self.nodes: #directed graph, have to go through all of them
                if suffix(kmer_1) == prefix(kmer_2) and kmer_1 != kmer_2:
                    self.add_edge(kmer_1, kmer_2)

    def __init__(self, kmers):
        super().__init__()
        self.__kmers = kmers # transformation into Kmer object happens in __setattr__
        self.__k = len(self.__kmers[0])


    @property
    def k(self):
        return self.__k
    @k.setter
    def k(self, value):
        raise AttributeError('k is read-only, deduced from kmers')

    @property
    def kmers(self):
        return self.__kmers
    @kmers.setter
    def kmers(self, value):
        # this will invoke __setattr__ method where validation and kmerization happen
        self.__kmers=value

    def __setattr__(self, name, value):

        if name == '_HAGraph__kmers': # when list changes reset the graph
            value= kmerize(value)
            self.__k = len(value[0])
            self.clear()
            if not all(validate_kmer(kmer,self.__k) for kmer in value):
                raise ValueError('Invalid kmers (should be all of dna sequences of the same length)')
            super().__setattr__(name, value)
            self.__establish_nodes_links()

        else:
            return super().__setattr__(name, value)



    def __str__(self):
        stats={'kmers':self.number_of_nodes(),'overlaps':self.number_of_edges()}
        stats=pd.Series(stats)
        display(stats)

        self.viz()

        return 'Hamiltionian Assembly Graph'


    def viz(self,title='None'):
        if title=='None':
            title='Hamiltionian Assembly Graph'
        pos=nx.spring_layout(self, seed=42)
        fig, ax = plt.subplots()
        # color title
        ax.set_title(title, fontsize=10, color='#660033')
        fig.patch.set_alpha(0)  
        ax.set_facecolor("none")  
        nx.draw(self, pos, with_labels=True, edge_color='#660033', node_color='lightgrey', 
                font_size=10, node_size=1000, ax=ax)
        plt.show()

        
if __name__ == '__main__':
    dna_string = "ATGCGATGACCTGACT"
    k = 3
    a=split_into_kmers(dna_string, k)
    g=HAGraph(a)
    print(g)