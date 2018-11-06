import math
import random
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import networkx as nx
from sortedcontainers import SortedDict


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


def create_n_degree_distribution_plot(degree_distributions, plot_file_name, plot_name):
    colours = ['b', 'g', 'r', 'c', 'm', 'y', 'k']

    i = 0
    num_colours = len(colours)

    # clears plot
    plt.clf()

    plt.ylim(0,0.25)
    # plt.xlim(0,280)

    # plot degree distribution
    plt.xlabel('Degree')
    plt.ylabel('Normalized Rate')
    plt.title('Degree Distribution of ' + plot_name)
    for dist in degree_distributions:
        m, k, p, q, dd = dist
        xdata = []
        ydata = []
        for deg in dd:
            xdata += [deg]
            ydata += [dd[deg]]
        plt.plot(xdata, ydata, marker='.', linestyle='None', color=colours[i % num_colours], label='p = ' + str(p)
                                                                                                   + ', q = ' + str(q))
        i += 1
    plt.legend(loc='upper right', ncol=1)
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
        print('diameter calculations', i, p, q)
        graph = make_nx_ring_group_graph(m,k,p,q)
        if not nx.is_connected(graph):
            num_instances -= 1
            print('not connected')
            continue
        temp = nx.diameter(graph)
        print(temp)
        diameter += temp
    if num_instances == 0: return -1
    return round(diameter/num_instances,2)


def make_p_diameter(t,m,k,q):
    p_diameter = {}
    p = round(q + 1 / 10, 7)
    while p < 1:
        val = average_ring_group_graph_diameter(t,m,k,p,q)
        if val is not -1:
            p_diameter[p] = val
            print(p,p_diameter[p])
        p = round(p * 1.2,7)
    return p_diameter


def create_n_diameter_p_plot(p_diameter_dists, plot_file_name, plot_name):
    # clears plot
    plt.clf()

    colours = ['b','g','r','c','m','y','k']

    i = 0
    num_colours = len(colours)

    plt.ylim(0,12)

    # plot degree distribution
    plt.xlabel('probability p')
    plt.ylabel('diameter')
    plt.title('Diameter vs p probability for ' + plot_name)
    for dist in p_diameter_dists:
        q, p_diameter = dist
        p_diameter = SortedDict(p_diameter)
        xdata = []
        ydata = []
        for p in p_diameter:
            xdata += [p]
            ydata += [p_diameter[p]]
        print(p_diameter)
        plt.plot(xdata, ydata, marker='.', linestyle='-', color=colours[i % num_colours], label='q = '+ str(q))
        i += 1
    plt.legend(loc='upper right')
    plt.savefig('distributions/' + plot_file_name + '.png')


def rgg(m,k,p,q):
    # the group itself, left and right groups then minus the vertex itself
    p_exp = p * (3 * k - 1)
    # the totak m*k - 3*k from above
    q_exp = q * k * (m - 3)
    total_exp = p_exp + q_exp
    return total_exp, p_exp, q_exp

def main_dia():
    m = 50
    k = 20
    q = 0.0001

    mGreater_dia = []
    kGreater_dia = []

    for i in range(4):
        mGreater_dia += [(q, make_p_diameter(20, m, k, q))]
        kGreater_dia += [(q, make_p_diameter(20, k, m, q))]
        q *= 10
        print(mGreater_dia)
        print(kGreater_dia)

    create_n_diameter_p_plot(mGreater_dia, 'q1/diameter/' + str(m) + '-' + str(k), 'Ring Group Graph')
    create_n_diameter_p_plot(kGreater_dia, 'q1/diameter/' + str(k) + '-' + str(m), 'Ring Group Graph')


def main_dd():
    m = 4
    k = 100
    # q = 0.0001

    # for i in range(4):
    #     # rgg_dict = {}
    #     mGreater_dd = []
    #     kGreater_dd = []
    #
    #     for j in range(20, int(round(100 - q*100)), 15):
    #         p = round(q + j / 100, 7)
    #         val1 = average_normalized_degree_distribution(25, m, k, p, q)
    #         if val1 is not -1:
    #             mGreater_dd += [(m,k,p,q,val1)]
    #         # val2 = average_normalized_degree_distribution(25, k, m, p, q)
    #         # if val2 is not -1:
    #         #     kGreater_dd += [(m,k,p,q,val2)]
    #         print(i,j,p)
    #
    #     create_n_degree_distribution_plot(mGreater_dd, 'q1/prob-p/' + str(m) + '-' + str(k) + '-' + str(q), 'Ring Group Graph')
    #     # create_n_degree_distribution_plot(kGreater_dd, 'q1/prob-p/' + str(k) + '-' + str(m) + '-' + str(q), 'Ring Group Graph')
    #
    #     q *= 10

    q = 0.001

    for i in range(3):
        # rgg_dict = {}
        mGreater_dd = []
        kGreater_dd = []

        p = round(0.5 - q, 3)
        temp_q = round(q, 3)

        while p > temp_q:
            val1 = average_normalized_degree_distribution(25, m, k, p, temp_q)
            if val1 is not -1:
                mGreater_dd += [(m, k, p, temp_q, val1)]
            # val2 = average_normalized_degree_distribution(25, k, m, p, temp_q)
            # if val2 is not -1:
            #     kGreater_dd += [(m, k, p, temp_q, val2)]

            temp_q = round(temp_q + 0.05, 3)
            p = round(0.5 - temp_q, 3)
            print(i,p,temp_q)

        # create_n_degree_distribution_plot(kGreater_dd, 'q1/' + str(k) + '-' + str(m) + '-' + str(q), 'Ring Group Graph')
        create_n_degree_distribution_plot(mGreater_dd, 'q1/' + str(m) + '-' + str(k) + '-' + str(q), 'Ring Group Graph')

        q *= 10


main_dd()
# m = 50
# k = 20
#
# mGreater_dia = [(0.0001, {0.1001: 11.8, 0.12012: 10.65, 0.144144: 9.85, 0.1729728: 9.6, 0.2075674: 9.05, 0.2490809: 8.8, 0.2988971: 8.45, 0.3586765: 8.1, 0.4304118: 7.9, 0.5164942: 7.8, 0.619793: 7.85, 0.7437516: 7.6, 0.8925019: 7.45}),
#                 (0.001, {0.101: 8.05, 0.1212: 7.2, 0.14544: 6.3, 0.174528: 6.1, 0.2094336: 6.0, 0.2513203: 5.45, 0.3015844: 5.1, 0.3619013: 5.0, 0.4342816: 5.0, 0.5211379: 5.0, 0.6253655: 5.0, 0.7504386: 4.85, 0.9005263: 4.2}),
#                 (0.01, {0.11: 4.1, 0.132: 4.0, 0.3284582: 4.0, 0.19008: 4.0, 0.8173092: 3.0, 0.4729798: 3.35, 0.681091: 3.0, 0.1584: 4.0, 0.2737152: 4.0, 0.5675758: 3.1, 0.980771: 3.0, 0.3941498: 4.0, 0.228096: 4.0}),
#                 (0.1, {0.2: 3.0, 0.24: 3.0, 0.288: 2.8, 0.8599634: 2.0, 0.41472: 2.35, 0.7166362: 2.0, 0.5971968: 2.0, 0.3456: 2.5, 0.497664: 2.15})]
#
# kGreater_dia = [(0.0001, {0.1001: 6.9, 0.12012: 6.25, 0.144144: 6.05, 0.1729728: 6.0, 0.2075674: 5.75, 0.2490809: 5.3, 0.2988971: 5.1, 0.3586765: 5.1, 0.4304118: 5.0, 0.5164942: 5.1, 0.619793: 5.1, 0.7437516: 5.05, 0.8925019: 4.95}),
#                 (0.001, {0.101: 5.0, 0.1212: 5.0, 0.14544: 5.0, 0.174528: 4.95, 0.2094336: 4.15, 0.2513203: 4.0, 0.3015844: 4.0, 0.3619013: 4.0, 0.4342816: 4.0, 0.5211379: 4.0, 0.6253655: 4.0, 0.7504386: 3.35, 0.9005263: 3.0}),
#                 (0.01, {0.11: 4.0, 0.132: 4.0, 0.3284582: 3.0, 0.19008: 3.1, 0.8173092: 3.0, 0.4729798: 3.0, 0.681091: 3.0, 0.1584: 3.8, 0.2737152: 3.0, 0.5675758: 3.0, 0.980771: 3.0, 0.3941498: 3.0, 0.228096: 3.0}),
#                 (0.1, {0.2: 2.65, 0.24: 2.05, 0.288: 2.1, 0.8599634: 2.0, 0.41472: 2.0, 0.7166362: 2.0, 0.5971968: 2.0, 0.3456: 2.0, 0.497664: 2.0})]
#
# create_n_diameter_p_plot(mGreater_dia, 'q1/diameter/' + str(m) + '-' + str(k), 'Ring Group Graph')
# create_n_diameter_p_plot(kGreater_dia, 'q1/diameter/' + str(k) + '-' + str(m), 'Ring Group Graph')