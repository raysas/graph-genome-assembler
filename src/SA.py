'''
simulated annealing on graph to solve hamiltonian path problem
'''

import networkx as nx
import numpy as np
import random
import matplotlib.pyplot as plt

from utils import *

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
        probability of accepting a new path (should be better than the old one <=> lower cost)
        if cost difference is positive, 
        '''
        diff = old_cost - new_cost
        return np.exp((diff) / T)
        # -- note while debugging: 
        #      1. we will only consider cases where diff is positive in the algo
        #      2. the greater the difference (more positive), the higher 
    
    def __cost(self, path):
        '''
        cost function: sum of the weights of the edges in the path 
        will also penalize going thoug a node several times
        '''
        total_cost = 0
        for i in range(len(path) - 1):
            edge_weight = self.graph[path[i]][path[i+1]]['weight']
            total_cost += edge_weight
        
        visited = set()
        penalty = 0
        penalty_weight = max(total_cost, 1) * 10 

        for node in path:
            if node in visited:
                penalty += penalty_weight 
            visited.add(node)

        return total_cost + penalty 

    def __anneal(self):
        '''
        -- simulated annealing --
        --------------------------
        
        main simulated annealing algorithm
        '''
        # -- random initialization
        path = list(self.graph.nodes)
        random.shuffle(path)
        old_cost = self.__cost(path)
        best_path = path.copy()
        best_cost = old_cost
        T = self.T
        
        for i in range(self.n_iter):
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
            
            print(f'iteration {i} - cost: {best_cost} - temperature: {T}')
        
        return best_path, best_cost