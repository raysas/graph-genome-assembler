'''Utility functions for graph construction'''


def suffix(kmer):
    return kmer[1:]
def prefix(kmer):
    return kmer[:-1]
def validate_dna(dna):
    return all(c in 'ACGT' for c in dna)
def validate_kmer(kmer,k):
    return validate_dna(kmer) and len(kmer) == k

class Kmer:
    '''class Kmer defined as data type for nodes (to allow for duplicate nodes in the multi-digraph)'''
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
# -- successful


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