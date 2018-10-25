import random
import networkx as nx
import matplotlib.pyplot as plt


def load_graph(file):
    """
    Loads a graph from a text file.
    Then returns the graph as a dictionary.
    """
    graph = open(file)

    read_vertices = False
    read_edges = False

    vertices_dict = {}

    answer_graph = {}

    edges = 0

    for line in graph:
        if read_vertices and line != '' and 'Edges' not in line:
            vertex_key, vertex_val = line.strip().split(' ',1)
            vertices_dict[vertex_key] = vertex_val

        if read_edges and line != '':
            node_v, node_u, weight = line.strip().split(' ')
            if node_v == node_u:
                continue
            answer_graph[node_v].add(node_u)
            answer_graph[node_u].add(node_v)

            edges += 1

        if 'Vertices' in line:
            read_vertices = True
            print('starting to read vertices')

        if 'Edges' in line:
            read_vertices = False
            read_edges = True
            for v in vertices_dict:
                answer_graph[v] = set()
            print('starting to read edges')

    print("Loaded graph with", len(answer_graph.keys()), "nodes")
    print('number of edges', edges)

    return vertices_dict, answer_graph


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

    edges = 0

    pedges = 0
    qedges = 0

    for i in range(num_vertices): ring_group_graph[i] = set()

    for v in range(num_vertices):
        v_group = v // k

        for u in range(v + 1, num_vertices):
            u_group = u // k
            random_number = random.random()
            if v_group == u_group or (abs(v_group - u_group) % m ) == 1 or (abs(v_group - u_group) % m ) == (m - 1):
                if random_number < p:
                    edges += 1
                    pedges +=1
                    ring_group_graph[v].add(u)
                    ring_group_graph[u].add(v)
            else:
                # it seems that it is more likely that this condition will be selected making this a random graph
                # with size m*k and probability q if m >> k
                if random_number < q:
                    edges += 1
                    qedges += 1
                    ring_group_graph[v].add(u)
                    ring_group_graph[u].add(v)
    print(edges,'p',pedges,'q',qedges)
    return ring_group_graph


# This code is based on networkx implementation and lecture3 implementation of directed PA graph
def make_pa_graph(total_nodes, num_nodes):
    repeated_nodes = []

    PA_graph = {}

    edges = 0

    for vertex in range(total_nodes):
        PA_graph[vertex] = set()

    for v in range(num_nodes):
        for u in range(v + 1, num_nodes):
            PA_graph[v].add(u)
            PA_graph[u].add(v)
        repeated_nodes += [v] * num_nodes
    for v in range(num_nodes, total_nodes):
        for dummy_idx in range(num_nodes):
            node = random.choice(repeated_nodes)
            PA_graph[v].add(node)
            PA_graph[node].add(v)
        repeated_nodes.extend(list(PA_graph[v]))
        repeated_nodes.extend([v] * len(PA_graph[v]))
        edges += len(PA_graph[v])
    return PA_graph


def vertex_brilliance(graph, source):
    """finds the distance (the length of the shortest path) from the source to
    every other vertex in the same component using breadth-first search, and
    returns the value of the largest distance found"""
    kstar = {}
    for v in graph[source]:
        kstar[v] = 0
    if not kstar.values(): return 0
    while True:
        for v in kstar:
            for u in graph[v]:
              if u in kstar:
                  kstar[v] -= 1
        minkstar = min(kstar.values())
        if minkstar == 0: break
        else:
            remove = []
            for k in kstar:
                if kstar[k] == minkstar and len(remove) == 0:
                    remove += [k]
                else:
                    kstar[k] = 0
            kstar.pop(remove[0], None)
    return len(kstar)


def vertex_brilliance_distribution(graph):
    """Takes a graph and computes the unnormalized distribution of the degrees of the nodes.
    Returns a dictionary whose keys correspond to degrees of nodes in the graph and values are the number of nodes
    with that degree.
    Degrees with no corresponding nodes in the graph are not included in the dictionary."""
    brilliance = {}
    for v in graph:
        brilliance[v] = vertex_brilliance(graph, v)
    # initialize dictionary for degree distribution
    brilliance_distribution = {}
    # consider each vertex
    for vertex in brilliance:
        # update degree_distribution
        if brilliance[vertex] in brilliance_distribution:
            brilliance_distribution[brilliance[vertex]] += 1
        else:
            brilliance_distribution[brilliance[vertex]] = 1
    return brilliance_distribution


def normalized_vertex_brilliance_distribution(graph, num_nodes):
    """Takes a graph and computes the normalized distribution of the degrees of the graph.
    Returns a dictionary whose keys correspond to in-degrees of nodes in the graph and values are the
    fraction of nodes with that in-degree.
    Degrees with no corresponding nodes in the graph are not included in the dictionary."""
    unnormalized_dist = vertex_brilliance_distribution(graph)
    normalized_dist = {}
    for degree in unnormalized_dist:
        normalized_dist[degree] = 1.0 * unnormalized_dist[degree] / num_nodes
    return normalized_dist


def create_vertex_brilliance_distribution_plot(brilliance_distribution, plot_file_name, plot_name):
    # create arrays for plotting
    xdata = []
    ydata = []
    for vertex in brilliance_distribution:
        xdata += [vertex]
        ydata += [brilliance_distribution[vertex]]

    # clears plot
    plt.clf()

    # plot vertex distribution
    plt.xlabel('Vertex Brilliance')
    plt.ylabel('Normalized Rate')
    plt.title('Vertex Brilliance Distribution of ' + plot_name)
    plt.plot(xdata, ydata, marker='.', linestyle='None', color='b')
    plt.savefig('distributions/q2/' + plot_file_name + '.png')


def compute_degrees(graph):
    """Takes a directed graph and computes the in-degrees for the nodes in the
    graph. Returns a dictionary with the same set of keys (nodes) and the
    values are the in-degrees."""
    #initialize in-degrees dictionary with zero values for all vertices
    degree = {}
    for vertex in graph:
        degree[vertex] = 0
    #consider each vertex
    for vertex in graph:
        #amend degree[w] for each outgoing edge from v to w
        for neighbour in graph[vertex]:
            degree[neighbour] += 1
    return degree


def degree_distribution(graph):
    """Takes a directed graph and computes the unnormalized distribution of the
    in-degrees of the graph.  Returns a dictionary whose keys correspond to
    in-degrees of nodes in the graph and values are the number of nodes with
    that in-degree. In-degrees with no corresponding nodes in the graph are not
    included in the dictionary."""
    #find in_degrees
    degree = compute_degrees(graph)
    #initialize dictionary for degree distribution
    degree_distribution = {}
    #consider each vertex
    for vertex in degree:
        #update degree_distribution
        if degree[vertex] in degree_distribution:
            degree_distribution[degree[vertex]] += 1
        else:
            degree_distribution[degree[vertex]] = 1
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


def average_normalized_distribution(distributions):
    """finds the average distribution of a list of distributions"""
    cumulative_dist = {}

    for dist in distributions:
        for deg in dist:
            if deg in cumulative_dist:
                cumulative_dist[deg] += dist[deg]
            else:
                cumulative_dist[deg] = dist[deg]
    average_dist = {}
    for deg in cumulative_dist:
        average_dist[deg] = cumulative_dist[deg] / len(distributions)
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
    plt.savefig('distributions/q2/' + plot_file_name + '.png')


# graph = {1: {3,2}, 2: {3,1}, 3: {1,2,4,6}, 4: {3,5}, 5: {4,7,6}, 6: {3,7,5,8}, 7: {6,5,8,9}, 8: {6,7}, 9: {7}}
#
# print(vertex_brilliance_distribution(graph))


# vertices_dict, coauthorship_graph = load_graph("coauthorship.txt")
# vertex_brilliance_distribution_coauthorship = normalized_vertex_brilliance_distribution(coauthorship_graph,
#                                                                                             len(coauthorship_graph.keys()))
# create_vertex_brilliance_distribution_plot(vertex_brilliance_distribution_coauthorship,
#                                                 'coauthorship-vertex-brilliance','coauthorship')


# n = 1559
# m = 30
#
# t = 100
#
# degree_distributions = []
# brilliance_distributions = []
#
# for i in range(t):
#     print(i)
#     pa_graph = make_pa_graph(n,m)
#     degree_distributions += [normalized_degree_distribution(pa_graph, n)]
#     brilliance_distributions += [normalized_vertex_brilliance_distribution(pa_graph, n)]
# pa_graph_dd = average_normalized_distribution(degree_distributions)
# pa_graph_vb = average_normalized_distribution(brilliance_distributions)
# create_degree_distribution_plot(pa_graph_dd, 'pa_graph-' + str(n) + '-' + str(m) + '-degree', 'PA Graph')
# create_vertex_brilliance_distribution_plot(pa_graph_vb, 'pa_graph-' + str(n) + '-' + str(m) + '-brilliance', 'PA Graph')
#
# for i in range(t):
#     print(i)
#     pa_graph_nx = nx.to_dict_of_lists(nx.barabasi_albert_graph(n, m))
#     degree_distributions += [normalized_degree_distribution(pa_graph_nx, n)]
#     brilliance_distributions += [normalized_vertex_brilliance_distribution(pa_graph_nx, n)]
# pa_graph_dd_nx = average_normalized_distribution(degree_distributions)
# pa_graph_vb_nx = average_normalized_distribution(brilliance_distributions)
# create_degree_distribution_plot(pa_graph_dd_nx, 'nx-pa_graph-' + str(n) + '-' + str(m) + '-degree', 'PA Graph nx')
# create_vertex_brilliance_distribution_plot(pa_graph_vb_nx, 'nx-pa_graph-' + str(n) + '-' + str(m) +
#                                                 '-brilliance', 'PA Graph nx')

# TODO find the right m and k so that the number of nodes is 1559 or close and p and q so that the number of edges
# is close to 43617

# (60,26,0.22,0.03), (39,40,0.17,0.03), (120,13,0.26,0.04) give almost the same numbers of nodes and edges as coauthorhip graph

# m = 120
# k = 13
# p = 0.26
# q = 0.04
#
# t = 100
#
# degree_distributions = []
# brilliance_distributions = []
#
# for i in range(t):
#     print(i)
#     ring_group_graph = make_ring_group_graph(m, k, p, q)
#     degree_distributions += [normalized_degree_distribution(ring_group_graph, m * k)]
#     brilliance_distributions += [normalized_vertex_brilliance_distribution(ring_group_graph, m * k)]
# ring_group_graph_dd = average_normalized_distribution(degree_distributions)
# create_degree_distribution_plot(ring_group_graph_dd,
#                                            'ring_group_graph-' + str(m) + '-' + str(k) + '-'
#                                            + str(p) + '-' + str(q) + '-degree',
#                                            'Ring Group Graph')
# ring_group_graph_vb = average_normalized_distribution(brilliance_distributions)
# create_vertex_brilliance_distribution_plot(ring_group_graph_vb,
#                                            'ring_group_graph-' + str(m) + '-' + str(k) + '-'
#                                            + str(p) + '-' + str(q) + '-brilliance',
#                                            'Ring Group Graph')

# ring_group_graph = make_ring_group_graph(m, k, p, q)
# ring_group_graph = make_ring_group_graph(k, m, p, q)

# ring_group_graph_dd = normalized_degree_distribution(ring_group_graph, m * k)
# create_degree_distribution_plot(ring_group_graph_dd,
#                                            'ring_group_graph-' + str(k) + '-' + str(m) + '-'
#                                            + str(p) + '-' + str(q) + '-degree',
#                                            'Ring Group Graph')
# ring_group_graph_vb = normalized_vertex_brilliance_distribution(ring_group_graph, k * m)
# create_vertex_brilliance_distribution_plot(ring_group_graph_vb,
#                                            'ring_group_graph-' + str(k) + '-' + str(m) + '-'
#                                            + str(p) + '-' + str(q) + '-brilliance',
#                                            'Ring Group Graph')
