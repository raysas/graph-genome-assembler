import networkx as nx

def suffix(kmer):
    return kmer[1:]
def prefix(kmer):
    return kmer[:-1]
def validate_dna(dna):
    return all(c in 'ACGT' for c in dna)

class HAGraph(nx.DiGraph):
    '''
    Hamiltionian Assembly Graph class
        inherits from networkx.DiGraph
        takes a list of kmers while instantiating
    '''
    def __init__(self, kmers):
        super().__init__()
        self.add_nodes_from(kmers)
        


