import networkx as nx
import pandas as pd
from IPython.display import display
import matplotlib.pyplot as plt

from utils import *

class DBGraph(nx.MultiDiGraph):
    '''
    -- De Bruijn Graph class --
    ----------------------------
        inherits from networkx.MultiDiGraph (directed graph with multiple edges)
        takes a list of kmers while instantiating (no transformation into Kmer object because nodes are unique)

    * edges: k-mers  
    * nodes: prefix and suffix of k-mers (k-1 mers)
    * k: length of the kmers

    Multiple edges are allowed because here edges are k-mers, to account for multiple k-mers of the same sequence -> allowed duplicate edges
    '''

    def __establish_nodes_links(self):
        for kmer in self.kmers:
            pref, suff = prefix(kmer), suffix(kmer)
            self.add_node(pref)
            self.add_node(suff)
            self.add_edge(pref, suff, label=kmer)

    def __init__(self, kmers):
        super().__init__()
        self.__kmers = kmers
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
        self.__kmers = value

    def __setattr__(self, name, value):
        if name == '_DBGraph__kmers':
            self.__k = len(value[0])
            self.clear()
            if not all(validate_kmer(kmer,self.__k) for kmer in value):
                raise ValueError('Invalid kmers (should be all of dna sequences of the same length)')
            super().__setattr__(name, value)
            self.__establish_nodes_links()
        else:
            return super().__setattr__(name, value)
        
    def viz(self, title='De Bruijn Graph'):
        '''
        -- visualize the graph --
        '''
        pos=nx.spring_layout(self, seed=42)
        fig, ax = plt.subplots()
        pos = nx.spring_layout(self)
        ax.set_title(title, fontsize=10, color='#4E9BB9')
        fig.patch.set_alpha(0)  
        ax.set_facecolor("none")  
        labels = nx.get_edge_attributes(self, 'label')
        nx.draw(self, pos, with_labels=False, node_size=800, node_color="#91CBD7", connectionstyle='arc3, rad = 0.1', edge_color='lightgrey')
        nx.draw_networkx_edge_labels(self, pos, edge_labels=labels, font_color='#4E9BB9', font_size=10, bbox=dict(facecolor='none', edgecolor='none'))
        plt.show()

        
    def __str__(self):
        stats = {'kmers': self.number_of_nodes(), 'overlaps': self.number_of_edges()}
        stats = pd.Series(stats)
        display(stats)
        self.viz()
        return 'De Bruijn Graph'
    
if __name__ == '__main__':
    dna_string = "ATGCGATGACCTGACT"
    k = 3
    a=split_into_kmers(dna_string, k)
    g=DBGraph(a)
    print(g)