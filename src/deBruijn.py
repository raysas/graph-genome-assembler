# ----------------------------------- useful functions -----------------------------------
#########################################################################################
def suffix(kmer):
    return kmer[1:]
def prefix(kmer):
    return kmer[:-1]
def validate_dna(dna):
    return all(c in 'ACGT' for c in dna)

# ----------------------------------- useful classes ----------------------------------- 
#########################################################################################

class DBGNode:
    '''
    ---- DeBruijn Graph Node ----
    ------------------------------
    Represents a node in a DeBruijn Graph which is a prefix or suffix of a k-mer
    General aim: solve the assembly problem efficiently

    This class is made to acount for the conditions to be met in order to get a Eulerian path of a graph
    As in the problem of the bridge of Königsberg, the graph must have 0 or 2 nodes with odd degree to have a Eulerian path (2 for the start and end nodes)  
    In the case of directed graphs, the in-degree and out-degree of a node must be equal - except for the start and end nodes - otherwsie the graph is not Eulerian

    Purpose is to save the in-degree and out-degree of each node in the graph for them to be easily accessible and checked
    '''
    
    def __init__(self, value):
        self.value = value
        self.in_degree = 0
        self.out_degree = 0

    def __str__(self):
        return f'''{self.value}'''
    
    def __repr__(self):
        return f'DBGNode({self.value})'
    
    def __eq__(self, other):
        if isinstance(other, DBGNode):
            return self.value == other.value
        elif isinstance(other, str):
            return self.value == other
        return False
    
    def __hash__(self):
        return hash(self.value)

class DBGEdge:
    '''
    ---- DeBruijn Graph Edge ----
    ------------------------------
    Represents an edge in a DeBruijn Graph which is a k-mer (solving teh Eulerian path for assembly)

    Since we will allow for multiple edges between two nodes, we will count the multiplicity of each edge through a decorator class CountMultiplicity
    Validation of the DNA sequence is done through the validate_dna function from initialisation (raises ValueError if not valid - allowed alphabet is {A, C, G, T})
    '''
    
    nodes=[]

    def __init__(self, kmer):
        self.kmer=kmer  
        self.from_node=DBGNode(prefix(kmer))
        self.to_node=DBGNode(suffix(kmer))
        self.multiplicity=1


    def __setattr__(self, name, value):
        if name=='kmer':
            if not validate_dna(value):
                raise ValueError('Invalid sequence alphabet - DNA should be composed of A, C, G, T') 
            self.__dict__[name] = value 
        elif name == 'from_node':
            if value not in self.nodes:
                value.out_degree += 1 
                self.nodes.append(value)
                self.__dict__[name] = value
            else:
                
                temp=self.nodes[self.nodes.index(value)]
                temp.out_degree += 1
                self.__dict__[name] = temp
        elif name == 'to_node':
            if value not in self.nodes:
                value.in_degree += 1
                self.nodes.append(value)
                self.__dict__[name] = value
            else:
                temp=self.nodes[self.nodes.index(value)]
                temp.in_degree += 1
                self.__dict__[name] = temp
        else:
            super().__setattr__(name, value)
                
    def __str__(self):
        print(f'{self.from_node}--{self.kmer}-->{self.to_node}')
        return self.kmer
    
    def __repr__(self):
        return f'DBGEdge({self.kmer})'
    
    def __eq__(self, other):
        return self.kmer == other.kmer
    
    def __hash__(self):
        return hash(self.kmer)


# ----------------------------------- main class ----------------------------------- #
######################################################################################
class DBG:
    '''
    -- DeBruijn Graph --
    ---------------------
    Directed graph with k-mers as edges, prefixes and suffixes of kmers as vertices  

    This will be a directed graph with the possibility of having multiple edges between two nodes  
    This will be encapsulated through a decorator class CountMultiplicity that will count the number of times an edge is added between two nodes
    '''
    def __init__(self, k=4):
        self.k = k
        self.edges = []
        self.nodes = []

    
    # -- graph construction --
    def __setattr__(self, name, value):
        '''this is to validate the added edge obeys the length of k'''
        if name=='edges':
            if not all(len(e.kmer) == self.k for e in value):
                raise ValueError(f'Invalid edge length - should be a k-mer of size {self.k}')
            self.__dict__[name] = value
        super().__setattr__(name, value)

    def update_nodes(self):
        self.nodes=DBGEdge.nodes
            
    def add_edge(self, kmer):
        if len(kmer) != self.k:
            print(f'[!] this k-mer ({kmer}) of size {len(kmer)} is not valid for this graph (k={self.k}) -- passing this edge')
            return

        edge = DBGEdge(kmer)#this has a multiplicity of 1
        edges_list=[e for e in self.edges] #this step is to be able to set the edges list explicitly in order to call the __setattr__ method (rather thna append)
        if edge in self.edges:
            # print('found it in te list')
            edge= edges_list[edges_list.index(edge)] #getting the isntance of the edge in the list
            edge.multiplicity += 1 
            # edges_list[edges_list.index(edge)]=edge #should be by reference so same instance will have its multiplicity incremented in the list
        else:
            # print('not found in the list')
            edges_list.append(edge)
        self.edges=edges_list #now calling __seattr__ :)
        self.update_nodes()

    # -- graph representation --

    # def plot(self):

    def __str__(self):
        return f'{self.edges}'
    
    def __repr__(self):
        return f'DBG({self.k})'
    
    # -- Eulerian path --

    # step1 1: verification
    def is_eulerian(self):
        '''
        checks if the graph has a Eulerian circuit

        Since this is a directed graph, we need to check if the in-degree and out-degree of each node is equal
        We can have a margin of 1 for the start and end nodes (in-degree=out-degree+1 for start and out-degree=in-degree+1 for end)
        '''
        source=[]; sink=[]
        for node in self.nodes:
            if node.in_degree - node.out_degree ==1:
                sink.append(node)
            elif node.out_degree - node.in_degree ==1:
                source.append(node)
            elif node.in_degree != node.out_degree:
                return False
        return len(source) == len(sink) == 0 or len(source) == len(sink) == 1
    
    def find_source(self):
        for i in range(len(self.nodes)):
            node=self.nodes[i]
            if node.out_degree - node.in_degree == 1:
                return i
        # if all are equal we can start at any node HAVING AN OUT EDGE
        for i in range(len(self.nodes)):
            node=self.nodes[i]
            if node.out_degree > 0:
                return i #gets the 1st node with an out edge
        return None
            
    # def find_sink(self):
    #     for node in self.nodes:
    #         if node.in_degree - node.out_degree == 1:
    #             return node
    #     return None