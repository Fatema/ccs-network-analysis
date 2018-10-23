import queue


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
        if  minkstar == 0: break
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


graph = {1: {3,2}, 2: {3,1}, 3: {1,2,4,6}, 4: {3,5}, 5: {4,7,6}, 6: {3,7,5,8}, 7: {6,5,8,9}, 8: {6,7}, 9: {7}}

print(vertex_brilliance_distribution(graph))