'''
simulated annealing on graph to solve hamiltonian path problem
'''

import networkx as nx
import numpy as np
import random
import matplotlib.pyplot as plt

from utils import *
from HAG import HAGraph

class SA:
    '''
    -- Simulated Annealing class --
    -------------------------------
    
    Simulated Annealing algorithm to solve the Hamiltonian Path problem on a graph
    
    * graph: a networkx.DiGraph instance
    * T: initial temperature
    * alpha: cooling rate
    * n_iter: number of iterations
    * verbose: print the progress of the algorithm
    '''
    
    def __init__(self, graph, T=1, alpha=0.99, n_iter=1000):
        self.graph = graph
        self.T = T
        self.alpha = alpha
        self.n_iter = n_iter
        # self.verbose = verbose
        
    def __acceptance_probability(self, old_cost, new_cost, T):
        '''
        formula: e^((old_cost - new_cost) / T)
        probability of accepting a new path (this will be considered when old < new, the larger the difference the less likely to accept)
        its used to get out of local minima
        '''
        delta_E =new_cost - old_cost

        return np.exp(
            -(delta_E) / T
            )

    def __cost(self, path):
        '''
        cost function: paths have to be complete (go through all nodes) 
        cost is on how many times we re revisiting a node and if the edge exists (not a very imp step when im picking from a set of valid edges)
        '''
        if len(path) < len(self.graph.nodes):  
            return float('inf')  # maximum penalty wehn incomplete
        
        visit_counts = {}  
        revisits = 0 
        
        for node in path:
            if node in visit_counts:
                visit_counts[node] += 1
                revisits += 1  # Count each revisit
            else:
                visit_counts[node] = 1

        edge_penalty = 0
        for i in range(len(path) - 1):
            if not self.graph.has_edge(path[i], path[i + 1]):
                edge_penalty += 10  

        return revisits * 5 + edge_penalty  

    def __anneal(self):
        '''
        -- simulated annealing --
        --------------------------
        
        main simulated annealing algorithm
        '''
        plots=[]
        # -- random initialization
        path = list(self.graph.nodes)
        random.shuffle(path)
        old_cost = self.__cost(path)
        best_path = path.copy()
        best_cost = old_cost
        T = self.T
        
        for _ in range(self.n_iter):
            new_path = path.copy()
            i, j = random.sample(range(len(path)), 2)
            new_path[i], new_path[j] = new_path[j], new_path[i]
            new_cost = self.__cost(new_path)
            
            if new_cost < old_cost or random.random() < self.__acceptance_probability(old_cost, new_cost, T): 
                path = new_path
                old_cost = new_cost
                
            if new_cost < best_cost:
                best_path = path.copy()
                best_cost = new_cost
                
            T=T*self.alpha
            
            print(f'iteration {_} - cost: {best_cost} - temperature: {T}')
            plot= self.viz_path(best_path,title=f'iteration {i} - cost: {best_cost} - temperature: {T}')
            plots.append(plot)
        
        return best_path, best_cost, plots
    
    def __call__(self, *args, **kwds):

        return self.__anneal()
    
    def viz_path(self, path,title='Hamiltonian Path'):
        '''
        -- visualize the path --
        '''
        pos = nx.spring_layout(self.graph)
        fig, ax = plt.subplots()
        ax.set_title('Hamiltonian Path', fontsize=10, color='#4E9BB9')
        fig.patch.set_alpha(0)  
        nx.draw(self.graph, pos, with_labels=True, edge_color='#91CBD7', node_color='lightgrey', 
                node_size=500, font_size=8, font_color='black', font_weight='bold', ax=ax)
        nx.draw_networkx_edges(self.graph, pos, edgelist=[(path[i], path[i+1]) for i in range(len(path) - 1)], edge_color='#ED6571', ax=ax)
        # plt.show()
        return plt
    
from matplotlib.animation import FuncAnimation
def dynamic_plots(plots):
    '''makes a movie of the plots'''
    fig, ax = plt.subplots()
    def update(frame):
        ax.clear()
        ax.imshow(plots[frame].gca().images[0].make_image())
        ax.axis('off')
    ani = FuncAnimation(fig, update, frames=len(plots), repeat=False)
    return ani

if __name__ == '__main__':
    dna_string = "ATGTAGTGAGTAGTGAAAATGAGTGTCCAGTGCTCGTCAGGTAGCTGACTAGCGTAGTCGCTA"
    k = 3
    a=split_into_kmers(dna_string, k)
    g=HAGraph(a)
    sa=SA(g, T=1, alpha=0.99, n_iter=10)
    best_path, best_cost, plots=sa()
    print(f'best path: {best_path} - best cost: {best_cost}')

    for i,plot in enumerate(plots):
        plot.savefig(f'./images/HAG/sa_{i}.png')