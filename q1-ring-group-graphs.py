import math
import random
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt


"""
The functions degree_distribution, normalized_degree_distribution and average_normalized_degree_distribution are taken from lecture 3 solution code (na3random_in-degree.py)
- average_normalized_degree_distribution: it was modified to fit with the ring group graph
"""


def make_ring_group_graph(m, k, p, q):
    """
    Returns a dictionary to a ring group graph with the specified number of nodes (m * k nodes)
    and edge probabilities.  The nodes of the graph are numbered 0 to
    m * k - 1.  For every pair of nodes, v and u, the pair is considered
    once: an edge (v,u) and an edge (u,v) with probability p if v and u are in the same group
    or if v and u are in adjacent groups (the labels for the nodes are 1 mod m apart),
    for all other cases the edges will be added with probability q.
    """
    num_vertices = m * k

    ## initialize the graph
    ring_group_graph = {}
    for i in range(num_vertices): ring_group_graph[i] = set()

    num_pass = 0
    num_fail = 0

    for v in range(num_vertices):
        v_group = v // k

        for u in range(v + 1, num_vertices):
            u_group = u // k
            random_number = random.random()
            if v_group == u_group or (abs(v_group - u_group) % m) == 1:
                num_pass += 1
                if random_number < p:
                    ring_group_graph[v].add(u)
                    ring_group_graph[u].add(v)
            else:
                # it seems that it is more likely that this condition will be selected making this a random graph
                # with size m*k and probability q
                num_fail += 1
                if random_number < q:
                    ring_group_graph[v].add(u)
                    ring_group_graph[u].add(v)
    return ring_group_graph


def degree_distribution(graph):
    """Takes a graph and computes the unnormalized distribution of the degrees of the nodes.
    Returns a dictionary whose keys correspond to degrees of nodes in the graph and values are the number of nodes
    with that degree.
    Degrees with no corresponding nodes in the graph are not included in the dictionary."""
    # initialize dictionary for degree distribution
    degree_distribution = {}
    # consider each vertex
    for vertex in graph:
        # update degree_distribution
        if len(graph[vertex]) in degree_distribution:
            degree_distribution[len(graph[vertex])] += 1
        else:
            degree_distribution[len(graph[vertex])] = 1
    return degree_distribution


def normalized_degree_distribution(graph, num_nodes):
    """Takes a graph and computes the normalized distribution of the degrees of the graph.
    Returns a dictionary whose keys correspond to in-degrees of nodes in the graph and values are the
    fraction of nodes with that in-degree.
    Degrees with no corresponding nodes in the graph are not included in the dictionary."""
    unnormalized_dist = degree_distribution(graph)
    normalized_dist = {}
    for degree in unnormalized_dist:
        normalized_dist[degree] = 1.0 * unnormalized_dist[degree] / num_nodes
    return normalized_dist


def average_normalized_degree_distribution(t, m, k, p ,q):
    """finds the average degree distribution of a random ring group graph on 2000 nodes
    in group of 100 of size 20 with edge probabilities 0.35 and 0.15 by taking k trials for each data point"""
    cumulative_dist = {}
    num_nodes = m * k
    for i in range(t):
        graph = make_ring_group_graph(m, k, p ,q)
        #find the distribution
        dist = normalized_degree_distribution(graph, num_nodes)
        for deg in dist:
            if deg in cumulative_dist:
                cumulative_dist[deg] += dist[deg]
            else:
                cumulative_dist[deg] = dist[deg]
    average_dist = {}
    for deg in cumulative_dist:
        average_dist[deg] = cumulative_dist[deg] / t
    return average_dist


def create_degree_distribution_plot(degree_distribution, plot_file_name, plot_name):
    # create arrays for plotting
    xdata = []
    ydata = []
    for degree in degree_distribution:
        xdata += [degree]
        ydata += [degree_distribution[degree]]

    # clears plot
    plt.clf()

    # plot degree distribution
    plt.xlabel('Degree')
    plt.ylabel('Normalized Rate')
    plt.title('Degree Distribution of ' + plot_name)
    plt.plot(xdata, ydata, marker='.', linestyle='None', color='b')
    plt.savefig('distributions/' + plot_file_name + '.png')


for i in range(1,25,5):
    p = round(0.25 + i/100, 2)
    q = round(0.25 - i/100, 2)
    m = 100
    k = 10
    print(m,k,p,q)
    mGreater = average_normalized_degree_distribution(100, m, k, p, q)
    create_degree_distribution_plot(mGreater, 'q1/' + str(m) + '-' + str(k) + '-' + str(p) + '-' + str(q), 'Ring Group Graph')
    kGreater = average_normalized_degree_distribution(100, k, m, p, q)
    create_degree_distribution_plot(kGreater, 'q1/' + str(k) + '-' + str(m) + '-' + str(p) + '-' + str(q), 'Ring Group Graph')