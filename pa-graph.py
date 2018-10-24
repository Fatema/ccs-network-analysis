#na3pa_in-degrees.py
#matthew johnson 19 january 2017

#####################################################

"""the following code creates a PA graph and is based on code from
http://www.codeskulptor.org/#alg_dpa_trial.py; """

import random

#first we need
class PATrial:
    """
    Used when each new node is added in creation of a PA graph.
    Maintains a list of node numbers with multiple instances of each number.
    The number of instances of each node number are in proportion to the
    probability that it is linked to.
    Uses random.choice() to select a node number from this list for each trial.
    """

    def __init__(self, num_nodes):
        """
        Initialize a PATrial object corresponding to a 
        complete graph with num_nodes nodes
        
        Note the initial list of node numbers has num_nodes copies of
        each node number
        """
        self._num_nodes = num_nodes
        self._node_numbers = [node for node in range(num_nodes) for dummy_idx in range(num_nodes)]


    def run_trial(self, num_nodes):
        """
        Conduct num_node trials using by applying random.choice()
        to the list of node numbers
        
        Updates the list of node numbers so that the number of instances of
        each node number is in the same ratio as the desired probabilities
        
        Returns:
        Set of nodes
        """       
        #compute the neighbors for the newly-created node
        new_node_neighbors = set()
        for dummy_idx in range(num_nodes):
            new_node_neighbors.add(random.choice(self._node_numbers))
        # update the list of node numbers so that each node number 
        # appears in the correct ratio
        self._node_numbers.append(self._num_nodes)
        self._node_numbers.extend(list(new_node_neighbors))        
        #update the number of nodes
        self._num_nodes += 1
        return new_node_neighbors
    
def make_complete_graph(num_nodes):
    """Takes the number of nodes num_nodes and returns a dictionary
    corresponding to a complete directed graph with the specified number of
    nodes. A complete graph contains all possible edges subject to the
    restriction that self-loops are not allowed. The nodes of the graph should
    be numbered 0 to num_nodes - 1 when num_nodes is positive. Otherwise, the
    function returns a dictionary corresponding to the empty graph."""
    #initialize empty graph
    complete_graph = {}
    #consider each vertex
    for vertex in range(num_nodes):
        #add vertex with list of neighbours
        complete_graph[vertex] = set([j for j in range(num_nodes) if j != vertex])
    return complete_graph
    
def make_PA_Graph(total_nodes, out_degree):
    """creates a PA_Graph on total_nodes where each vertex is iteratively
    connected to a number of existing nodes equal to out_degree"""
    #initialize graph by creating complete graph and trial object
    PA_graph = make_complete_graph(out_degree)
    trial = PATrial(out_degree)
    for vertex in range(out_degree, total_nodes):
        PA_graph[vertex] = trial.run_trial(out_degree)
    return PA_graph


def compute_in_degrees(digraph):
    """Takes a directed graph and computes the in-degrees for the nodes in the
    graph. Returns a dictionary with the same set of keys (nodes) and the
    values are the in-degrees."""
    #initialize in-degrees dictionary with zero values for all vertices
    in_degree = {}
    for vertex in digraph:
        in_degree[vertex] = 0
    #consider each vertex
    for vertex in digraph:
        #amend in_degree[w] for each outgoing edge from v to w
        for neighbour in digraph[vertex]:
            in_degree[neighbour] += 1
    return in_degree

def in_degree_distribution(digraph):
    """Takes a directed graph and computes the unnormalized distribution of the
    in-degrees of the graph.  Returns a dictionary whose keys correspond to
    in-degrees of nodes in the graph and values are the number of nodes with
    that in-degree. In-degrees with no corresponding nodes in the graph are not
    included in the dictionary."""
    #find in_degrees
    in_degree = compute_in_degrees(digraph)
    #initialize dictionary for degree distribution
    degree_distribution = {}
    #consider each vertex
    for vertex in in_degree:
        #update degree_distribution
        if in_degree[vertex] in degree_distribution:
            degree_distribution[in_degree[vertex]] += 1
        else:
            degree_distribution[in_degree[vertex]] = 1
    return degree_distribution

EX_GRAPH_PA1 = make_PA_Graph(27770, 13)

#plotting the normalized in_degree distribution it seems that, perhaps, we have found
#a better model for citation graphs

EX_DIST = in_degree_distribution(EX_GRAPH_PA1)
normalized_PA_distribution = {}
for degree in EX_DIST:
    normalized_PA_distribution[degree] = EX_DIST[degree] / 27770.0
    
#create arrays for plotting
xdata = []
ydata = []
for degree in normalized_PA_distribution:
    xdata += [degree]
    ydata += [normalized_PA_distribution[degree]]
    
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

#clears plot
plt.clf()

#plot degree distribution 
plt.xlabel('In-Degree')
plt.ylabel('Normalized Rate')
plt.title('In-Degrees of PA Graph on 27700 vertices with out-degree 13')
plt.loglog(xdata, ydata, marker='.', linestyle='None', color='b')
plt.savefig('na3pa.png')
