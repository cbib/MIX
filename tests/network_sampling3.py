"""
A set of functions to sample nodes of a graph.
 
Supported sampling methods:
- sample_with_replacements
- uniform_independent_node_sample
- degree_weighted_independent_node_sample
- random_walk

Supported estimators:
- estimate_relative_size
- estimate_mean


More details in:
M. Gjoka, M. Kurant, C. T. Butts and A. Markopoulou
"Walking in Facebook: A Case Study of Unbiased Sampling of OSNs"
INFOCOM 2010

Example:
>>> import networkx as nx
>>> G = nx.generators.wheel_graph(10)
>>> list(uniform_independent_node_sample(G, size=20))
[8, 6, 4, 7, 9, 5, 5, 2, 4, 5, 7, 9, 3, 6, 1, 6, 9, 8, 6, 4]
>>> sample = uniform_independent_node_sample(G, size=20)   #generator
>>> estimate_mean(sample, values=lambda x: 2*x)
8.9000000000000004
>>> sample = list(nx.random_walk(G, size=20))
>>> sample
[6, 7, 6, 7, 0, 3, 2, 3, 2, 0, 5, 0, 7, 6, 7, 6, 5, 6, 7, 6]
>>> estimate_relative_size(S, labeled=lambda x: x in [0,1,2,3,4] , weights=graph.degree)
0.50679999999999814
For more examples try 'nx.test_sampling()'

"""


#Author: Maciej Kurant



import networkx as nx
import random
import bisect
import itertools



__author__ = """Maciej Kurant"""
__all__ = ['sample_with_replacements', 'uniform_independent_node_sample', 'degree_weighted_independent_node_sample', 'random_walk', 'estimate_relative_size', 'estimate_mean', 'test_sampling']


############################
def sample_with_replacements(population, size=-1, weights=None):
    """
    sample_with_replacements(population, size=None, weights=None)
    
    Generates random items (with replacement) from a population, uniformly (by default) or with explicit sampling weights
    
    Parameters
    ----------  
    population:   - indexable container (e.g., a list) with items to sample from
    size:         - desired sample length (int). If -1 (default), then the generator never stops
    weights:      - function: item -> sampling_weight of this item
    """    
    
    n = len(population)
    _random = random.random
    
    if weights==None:  # uniform
        for c in itertools.count():
            if c==size:  return
            yield population[int(_random()*n)]
    
    else:    # weighted              
        cum_weights = [0]*n
        tot_weight = 0

        for i,v in enumerate(population):
            tot_weight += weights(v)
            cum_weights[i] = tot_weight
                    
        for c in itertools.count():
            if c==size:  return
            i = bisect.bisect_right(cum_weights, _random() * tot_weight)
            yield population[i]

        
        
####################
def uniform_independent_node_sample(graph, size=-1):
    """
    uniform_independent_node_sample(graph, size=-1)
    
    Generates a uniform independent sample (UIS) of nodes, with replacements.
    If size==-1 (default), then the generator never stops.
    """

    return sample_with_replacements(graph.nodes(), size)
        
####################
def degree_weighted_independent_node_sample(graph, size=-1):
    """
    degree_weighted_independent_node_sample(graph, size=-1)
    
    Generates an independent sample of nodes, with replacements, with sampling probabilities proportional to node degrees. 
    If size==-1 (default), then the generator never stops.
    """
    
    return sample_with_replacements(graph.nodes(), size=size, weights=graph.degree)
    
    #Faster:
    #for v in sample_with_replacements(V, size=size):
    #    yield random.choice(graph.neighbors(v))

                 

####################
def random_walk(graph, start_node=None, size=-1, metropolized=False):    
    """
    random_walk(G, start_node=None, size=-1):
    
    Generates nodes sampled by a random walk (classic or metropolized)
    
    Parameters
    ----------  
    graph:        - networkx.Graph 
    start_node    - starting node (if None, then chosen uniformly at random)
    size          - desired sample length (int). If -1 (default), then the generator only stops when hitting a sink node, or never stop for undirected graphs
    metropolized  - False (default): classic Random Walk
                    True:  Metropolis Hastings Random Walk (with the uniform target node distribution) 
    """
    
    # if type(graph) != nx.Graph:
    #     raise nx.NetworkXException("Graph must be a simple undirected graph!") 
        
    if start_node==None:
        start_node = random.choice(graph.nodes())
    
    v = start_node
    for c in itertools.count():
        if c==size:  return
        next_candidates = graph.neighbors(v)
        if len(next_candidates)==0: return 

        if metropolized:   # Metropolis Hastings Random Walk (with the uniform target node distribution) 
            candidate = random.choice(next_candidates)
            v = candidate if (random.random() < float(graph.degree(v))/graph.degree(candidate)) else v
        else:              # classic Random Walk
            v = random.choice(next_candidates)
            
        yield v
    

####################
def estimate_mean(sample, values, weights=None):
    """
    estimate_mean(sample, values, weights=None)

    Based on a sample, estimate and return the average value over all existing items. 

    Parameters
    ----------  
    sample:       - a sample of items (iterable)
    values:       - function: item -> value
    weights:      - function: item -> sampling_weight of this item   
    """           

    if weights==None:   # uniform sample
        weights = lambda x: 1
    
    up = down = 0.
    for v in sample:
        up += 1.*values(v)/weights(v) 
        down += 1./weights(v)
    return up/down


####################
def estimate_relative_size(sample, labeled, weights=None):
    """
    
    Based on a sample, estimate the relative number of labeled nodes. 
    
    Parameters
    ----------  
    sample:      - a sample of items (iterable)
    labeled:     - function: item -> True/False  (True if item is labeled)
    weights:     - function: item -> sampling_weight of this item     
    """
           
    values = lambda s: 1. if labeled(s) else 0.                
    return estimate_mean(sample, values,  weights)


####################
def test_sampling():
    
    graph = nx.generators.wheel_graph(10)
    print 
    print nx.info(graph)

    print """
Collect samples of N=10000 nodes with three methods:
1) uniform_independent_node_sample(graph, size=N)
2) random_walk(graph, size=N, metropolized=True)
3) random_walk(graph, size=N)"""

    N = 10000
    UIS_sample = list(uniform_independent_node_sample(graph, size=N))  # 'uniform' is the correct estimator
    MHRW_sample =list(random_walk(graph, size=N, metropolized=True)) # 'uniform' is the correct estimator
    RW_sample  = list(random_walk(graph, size=N))                      # 'random_walk' is the correct estimator

    
    print '\nSet node values to node numbers (0-9). The mean value is AV=4.5'
    values = lambda v: v
    
    print "Estimate AV using correct sampling weights: ",
    print "%0.2f, %0.2f, %0.2f" % (estimate_mean(UIS_sample, values), estimate_mean(MHRW_sample, values), estimate_mean(RW_sample, values, weights=graph.degree))
    
    print "Estimate AV using wrong sampling weights:   ",
    print "%0.2f, %0.2f, %0.2f" % (estimate_mean(UIS_sample, values, weights=graph.degree), estimate_mean(MHRW_sample, values, weights=graph.degree), estimate_mean(RW_sample, values))
        
    
    
    print "\nLabel nodes 0 and 2. The fraction of labeled nodes is F=0.2 " 
    labeled = lambda v: v in [0,2] 
    
    print "Estimate F using correct sampling weights: ", 
    print "%0.2f, %0.2f, %0.2f " % (estimate_relative_size(UIS_sample, labeled), estimate_mean(MHRW_sample, labeled), estimate_mean(RW_sample, labeled, weights=graph.degree))
    
    print "Estimate F using wrong sampling weights:   ",
    print "%0.2f, %0.2f, %0.2f " % (estimate_relative_size(UIS_sample, labeled, weights=graph.degree), estimate_mean(MHRW_sample, labeled, weights=graph.degree), estimate_mean(RW_sample, labeled))
