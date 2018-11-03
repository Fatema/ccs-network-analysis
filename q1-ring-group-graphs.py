import math
import random
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import networkx as nx


"""
The functions degree_distribution, normalized_degree_distribution and average_normalized_degree_distribution are taken from lecture 3 solution code (na3random_in-degree.py)
- average_normalized_degree_distribution: it was modified to fit with the ring group graph
"""

# TODO redo the graphs for ring_group


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

    # initialize the graph
    ring_group_graph = {}

    for i in range(num_vertices): ring_group_graph[i] = set()

    for v in range(num_vertices):
        v_group = v // k

        for u in range(v + 1, num_vertices):
            u_group = u // k
            random_number = random.random()
            if v_group == u_group or (abs(v_group - u_group) % m ) == 1 or (abs(v_group - u_group) % m ) == (m - 1):
                if random_number < p:
                    ring_group_graph[v].add(u)
                    ring_group_graph[u].add(v)
            else:
                # it seems that it is more likely that this condition will be selected making this a random graph
                # with size m*k and probability q if m >> k
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
    num_instances = t
    for i in range(t):
        graph = make_nx_ring_group_graph(m, k, p, q)
        if not nx.is_connected(graph):
            num_instances -= 1
            print('not connected')
            continue
        #find the distribution
        dist = normalized_degree_distribution(graph, num_nodes)
        for deg in dist:
            if deg in cumulative_dist:
                cumulative_dist[deg] += dist[deg]
            else:
                cumulative_dist[deg] = dist[deg]

    if num_instances == 0: return -1

    average_dist = {}
    for deg in cumulative_dist:
        average_dist[deg] = cumulative_dist[deg] / num_instances
    return average_dist


def create_degree_distribution_plot(degree_distribution, plot_file_name, plot_name):
    colours = ['b', 'g', 'r', 'c', 'm', 'y', 'k']

    i = 0
    num_colours = len(colours)

    # clears plot
    plt.clf()

    # plot degree distribution
    plt.xlabel('Degree')
    plt.ylabel('Normalized Rate')
    plt.title('Degree Distribution of ' + plot_name)
    for dist in degree_distribution:
        plt.plot(list(dist.keys()), list(dist.values()), marker='.', linestyle='None', color=colours[i % num_colours])
        i += 1
    plt.savefig('distributions/' + plot_file_name + '.png')


def make_nx_ring_group_graph(m, k, p, q):
    num_vertices = m * k

    ## initialize the graph
    ring_group_graph = nx.Graph()

    for v in range(num_vertices):
        v_group = v // k

        for u in range(v + 1, num_vertices):
            u_group = u // k
            random_number = random.random()
            if v_group == u_group or (abs(v_group - u_group) % m ) == 1 or (abs(v_group - u_group) % m ) == (m - 1):
                if random_number < p:
                    ring_group_graph.add_edge(v,u)
            else:
                # it seems that it is more likely that this condition will be selected making this a random graph
                # with size m*k and probability q if m >> k
                if random_number < q:
                    ring_group_graph.add_edge(v,u)
    return ring_group_graph


def average_ring_group_graph_diameter(t,m,k,p,q):
    diameter = 0
    num_instances = t
    for i in range(t):
        graph = make_nx_ring_group_graph(m,k,p,q)
        if not nx.is_connected(graph):
            num_instances -= 1
            print('not connected')
            continue
        diameter += nx.diameter(graph)
    if num_instances == 0: return -1
    return round(diameter/num_instances,2)


def make_p_diameter(t,m,k,q):
    p_diameter = {}
    p = round(q + 5 / 100, 5)
    while p < 1:
        val = average_ring_group_graph_diameter(t,m,k,p,q)
        if val is not -1:
            p_diameter[p] = val
            print(p,p_diameter[p])
        p *= 1.2
    return p_diameter


def create_diameter_p_plot(p_diameter, plot_file_name, plot_name):
    # clears plot
    plt.clf()

    colours = ['b','g','r','c','m','y','k']

    i = 0
    num_colours = len(colours)

    # plot degree distribution
    plt.xlabel('probability p')
    plt.ylabel('diameter')
    plt.title('Diameter vs p probability for ' + plot_name)
    for dist in p_diameter:
        plt.plot(list(dist.keys()), list(dist.values()), marker='.', linestyle='None', color=colours[i % num_colours])
        i += 1
    plt.savefig('distributions/' + plot_file_name + '.png')


def rgg(m,k,p,q):
    # the group itself, left and right groups then minus the vertex itself
    p_exp = p * (3 * k - 1)
    # the totak m*k - 3*k from above
    q_exp = q * k * (m - 3)
    total_exp = p_exp + q_exp
    return total_exp, p_exp, q_exp


m = 50
k = 10
q = 0.001

mGreater_dia = []
kGreater_dia = []

for i in range(3):
    mGreater_dia += [make_p_diameter(100, m, k, q)]
    kGreater_dia += [make_p_diameter(100, m, k, q)]
    q *= 10

    # rgg_dict = {}
    mGreater_dd = []
    kGreater_dd = []

    for j in range(5, int(round(100 - q*100)) + 5, 5):
        p = round(q + j / 100, 5)
        val1 = average_normalized_degree_distribution(100, m, k, p, q)
        if val1 is not -1:
            mGreater_dd += [val1]
        val2 = average_normalized_degree_distribution(100, k, m, p, q)
        if val2 is not -1:
            kGreater_dd += [val2]

    create_degree_distribution_plot(kGreater_dd, 'q1/prob-p/' + str(k) + '-' + str(m) + '-' + str(q), 'Ring Group Graph')

create_diameter_p_plot(mGreater_dia, 'q1/diameter/' + str(m) + '-' + str(k), 'Ring Group Graph')
create_diameter_p_plot(kGreater_dia, 'q1/diameter/' + str(k) + '-' + str(m), 'Ring Group Graph')
