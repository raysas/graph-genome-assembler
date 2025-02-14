'''Utility functions for graph construction and model eval'''

from Bio import pairwise2
import pandas as pd
from IPython.display import display

def suffix(kmer):
    return kmer[1:]
def prefix(kmer):
    return kmer[:-1]
def validate_dna(dna):
    return all(c in 'ACGT' for c in dna)
def validate_kmer(kmer,k):
    return validate_dna(kmer) and len(kmer) == k

class Kmer:
    '''
    -- Kmer class --
    ----------------
    
    simple helper class Kmer defined as data type for nodes (to allow for duplicate nodes in the graph)
    
    implements the following dunder methods to allow for pythonic semantics like len(kmer), kmer[1:], etc.

    * __init__: constructor
    * __repr__: string representation show from NodeView
    * __str__: string representation show from print 
    * __len__: length of the kmer
    * __getitem__: indexing and slicing

    '''
    def __init__(self, kmer):
        self.kmer = kmer

    # -- some useful magic methods to allow string-like operations
    def __repr__(self): #--> string representation show from NodeView
        return f'{len(self.kmer)}-mer {self.kmer}'
    def __str__(self): #--> string representation show from print
        return self.kmer
    def __len__(self):
        return len(self.kmer)
    def __getitem__(self, item): #--> spliceable
        return self.kmer[item]
    
# k1=Kmer('AUG')
# k2=Kmer('AUG')
# k1==k2
# k1[1:]
# prefix(k1)
# validate_kmer(k1,3)
# # -- successful


def kmerize(kmers):
    if isinstance(kmers, str):
        return Kmer(kmers)
    if all(isinstance(kmer, str) for kmer in kmers):
        return [Kmer(kmer) for kmer in kmers]
    if all(isinstance(kmer, Kmer) for kmer in kmers):
        return kmers


def split_into_kmers(dna_string, k):
    kmers = [dna_string[i:i+k] for i in range(len(dna_string) - k + 1)]
    return kmers

def evaluate_modeled_sequence(real_seq, graph_modeled_seq):
    '''
    evaluates how good the modeled seq is compared to the real seq (testing), using needleman-wunsch for alignment
    display a df of parameters used in accuracy computation
    param: 
        real_seq: str, true sequence
        graph_modeled_seq: str, predicted sequence
    return:
        accuracy: float, accuracy of the prediction which is defined based on length of the real sequence, gaps and mismatches
    '''
    len_real_seq = len(real_seq)
    len_graph_seq = len(graph_modeled_seq)

    aln = pairwise2.align.globalxx(real_seq, graph_modeled_seq, one_alignment_only=True)[0]
    aligned_real_seq, aligned_graph_seq, score, begin, end = aln

    gaps = aligned_real_seq.count('-') + aligned_graph_seq.count('-')
    mismatches = sum(1 for i in range(len(aligned_real_seq)) if aligned_real_seq[i] != aligned_graph_seq[i])
    accuracy = (len_real_seq - gaps - mismatches) / len_real_seq

    df=pd.DataFrame({'real_seq':[len_real_seq],'graph_seq':[len_graph_seq],'gaps':[gaps],'mismatches':[mismatches],'accuracy':[accuracy]})
    display(df)


    return accuracy