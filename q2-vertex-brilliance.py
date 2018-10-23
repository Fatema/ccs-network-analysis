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

    for line in graph:
        if read_vertices and line != '' and 'Edges' not in line:
            vertex_key, vertex_val = line.strip().split(' ',1)
            vertices_dict[vertex_key] = vertex_val

        if read_edges and line != '':
            node_v, node_u, weight = line.strip().split(' ')
            if node_v == node_u:
                continue
            if node_v in answer_graph:
                answer_graph[node_v].add(node_u)
            else:
                answer_graph[node_v] = {node_u}
            if node_u in answer_graph:
                answer_graph[node_u].add(node_v)
            else:
                answer_graph[node_u] = {node_v}

        if 'Vertices' in line:
            read_vertices = True
            print('starting to read vertices')

        if 'Edges' in line:
            read_vertices = False
            read_edges = True
            print('starting to read edges')

    print("Loaded graph with", len(answer_graph.keys()), "nodes")

    return vertices_dict, answer_graph


def vertex_brilliance(graph, source):
    """finds the distance (the length of the shortest path) from the source to
    every other vertex in the same component using breadth-first search, and
    returns the value of the largest distance found"""
    kstar = {}
    for v in graph[source]:
        kstar[v] = 0
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
    plt.xlabel('vertex brilliance')
    plt.ylabel('Normalized Rate')
    plt.title('Vertex Brilliance Distribution of ' + plot_name)
    plt.plot(xdata, ydata, marker='.', linestyle='None', color='b')
    plt.savefig('distributions/q2/' + plot_file_name + '.png')


# graph = {1: {3,2}, 2: {3,1}, 3: {1,2,4,6}, 4: {3,5}, 5: {4,7,6}, 6: {3,7,5,8}, 7: {6,5,8,9}, 8: {6,7}, 9: {7}}
#
# print(vertex_brilliance_distribution(graph))

vertices_dict, coauthorship_graph = load_graph("coauthorship.txt")
vertex_brilliance_distribution_coauthorship = normalized_vertex_brilliance_distribution(coauthorship_graph, len(coauthorship_graph.keys()))
create_vertex_brilliance_distribution_plot(vertex_brilliance_distribution_coauthorship,'coauthorship-vertex-brilliance','coauthorship')
