import random


def make_random_graph(num_nodes, prob):
    """Returns a random undirected graph with the specified number of nodes and edge probability.
    The nodes of the graph are numbered 0 to num_nodes - 1.  For every pair of nodes, i and j, the pair is connected with probability prob. """
    #initialize empty graph
    random_graph = {}
    for vertex in range(num_nodes): random_graph[vertex] = []
    #consider each vertex
    for vertex in range(num_nodes):
        #consider each neighbour with greater value
        for neighbour in range(vertex + 1, num_nodes):
            random_number = random.random()
            #add edge from vertex to neighbour with probability prob
            if random_number < prob:
                #maybe this method for set union is deprecated in python 3
                random_graph[vertex] = set(random_graph[vertex]) | {neighbour}
                random_graph[neighbour] = set(random_graph[neighbour]) | {vertex}
    return random_graph


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


def search_random_graph(graph, total_nodes, prob, start_vertex, target_vertex):
    search_time = 0
    average_neighbour_num = total_nodes * prob # Todo look into how you can possibly use this to improve the search time
    found_target = False
    random_move = 0
    while not found_target:
        moved = False
        print('current search time',search_time,'for', start_vertex,target_vertex)
        if start_vertex == target_vertex: break
        neighbours_num = len(graph[start_vertex])
        neighbours = list(graph[start_vertex])
        random.shuffle(neighbours)
        for i in range(neighbours_num):
            neighbour = neighbours[i]
            search_time += 1
            if target_vertex == neighbour:
                found_target = True
                break
            if random.random() < 0.05:
                start_vertex = neighbour
                moved = True
                random_move += 1
                print('random move', start_vertex, target_vertex, 'at', search_time)
                break
        if not moved and not found_target:
            start_vertex = random.choice(neighbours)
    print('total time taken to find target', search_time)
    print('number of random moves', random_move)
    return search_time

n = 1000
p = 0.1

graph = make_random_graph(n,p)

search_time = {}
avg_search_time = 0
start_v = random.randint(0,n - 1)
target_u = random.randint(0,n - 1)

for i in range(100):
    search_time[(start_v,target_u, i)]= search_random_graph(graph,n,p,start_v,target_u)
    avg_search_time += search_random_graph(graph,n,p,start_v,target_u)

print(graph)
print(search_time)
print(avg_search_time / 100)



