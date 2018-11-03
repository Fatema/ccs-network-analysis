import random
import networkx as nx


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


def search_random_graph(graph, total_nodes, prob, s, t):
    search_time = 0
    average_neighbour_num = total_nodes * prob
    random_move = 0

    while True:
        moved = False

        if s == t: break

        neighbours_num = len(graph[s])
        neighbours = list(graph[s])
        random.shuffle(neighbours)

        # at max iterate until the average number of neighbours
        num_iterations = max(average_neighbour_num, neighbours_num)

        for i in range(num_iterations):
            u = neighbours[i]
            search_time += 1
            if t == u:
                s = u
                moved = True
                break
            # probability that an edge might exist between u and t
            if random.random() < p:
                s = u
                moved = True
                random_move += 1
                break

        if not moved:
            s = random.choice(neighbours)
    # print('total time taken to find target', search_time)
    # print('number of random moves', random_move)
    return search_time


def search_random_graph_v2(graph, total_nodes, prob, start_vertex, target_vertex):
    search_time = 0
    average_neighbour_num = total_nodes * prob # Todo look into how you can possibly use this to improve the search time
    found_target = False
    while not found_target:
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
            else:
                start_vertex = neighbour
                break
    # print('total time taken to find target', search_time)
    return search_time


def search_ring_group_graph(graph, m, k, p, q, v, t):
    search_time = 0

    t_group = t // k

    max_num_edges_p = k**2 * m + k*(k-1)*m / 2
    max_num_edges = m*k * (m*k - 1) / 2

    avg_num_edges_p = max_num_edges_p * p
    avg_num_edges_q = (max_num_edges - max_num_edges_p) * q

    while True:
        # find the group label of the vertex
        v_group = v // k

        # reset moved flag
        moved = False

        # if the vertex id equals the target id stop
        if v == t: break

        # set d(v) for the vertex (number of neighbours)
        neighbours_num = len(graph[v])

        # set the neighbours array and shuffle it
        neighbours = list(graph[v])
        random.shuffle(neighbours)

        adjacent_groups = False

        # set the number of iterations and the adjacent flag
        if v_group == t_group or abs(v_group - t_group) == 0 or abs(v_group-t_group) == m - 1:
            num_iterations = max(int(round(neighbours_num * avg_num_edges_p / (avg_num_edges_p + avg_num_edges_q))),1)
            adjacent_groups = True
        else:
            num_iterations = max(int(round(neighbours_num * avg_num_edges_q / (avg_num_edges_p + avg_num_edges_q))),1)

        # keep track of adjacent neighbours
        adjacent_neighbour = -1
        closest_not_adjacent_v = -1

        for i in range(num_iterations):
            u = neighbours[i]
            u_group = u // k

            search_time += 1

            if t == u:
                v = u
                moved = True
                break

            # check if t is adjacent to one of the neighbours
            elif u_group == t_group or abs(u_group - t_group) == 0 or abs(u_group-t_group) == m - 1:
                # if v is not in adjacent group, u could be connected to t with probability p
                if not adjacent_groups and random.random() < p:
                    v = u
                    moved = True
                    break
                else:
                    # v is adjacent so there is a chance that t is connected to v
                    # so record u and look at other neighbours
                    adjacent_neighbour = u

            # if t is not adjacent to u and v, see if u is closer to t than v in terms of groups
            else:
                abs_u_t = abs(t_group - u_group)
                abs_v_t = abs(t_group - v_group)

                diff_u_t = min(abs_u_t, m - abs_u_t)
                diff_v_t = min(abs_v_t, m - abs_v_t)

                # as this is a ring both sides must be considered
                if diff_u_t < diff_v_t:
                    closest_not_adjacent_v = u
                    if random.random() < q:
                        v = u
                        moved = True
                        break
                # u is not closer to t than v and there hasn't been a vertex found to be close
                elif closest_not_adjacent_v == -1:
                    # if u and v are in the same group
                    if diff_v_t == diff_u_t:
                        closest_not_adjacent_v = u
                    # if non of the neighbours are close to t or in the same group as v
                    elif adjacent_neighbour == -1:
                        closest_not_adjacent_v = u
        if not moved:
            # no adjacent neighbour to t is found in terms of group
            if adjacent_neighbour == -1:
                v = closest_not_adjacent_v
            else:
                v = adjacent_neighbour
    return search_time


def search_ring_group_graph_v2(graph, m, k, p, q, v, t):
    search_time = 0

    t_group = t // k

    max_num_edges_p = k**2 * m + k*(k-1)*m / 2
    max_num_edges = m*k * (m*k - 1) / 2

    avg_num_edges_p = max_num_edges_p * p
    avg_num_edges_q = (max_num_edges - max_num_edges_p) * q

    while True:
        # find the group label of the vertex
        v_group = v // k

        # reset moved flag
        moved = False

        # if the vertex id equals the target id stop
        if v == t: break

        # set d(v) for the vertex (number of neighbours)
        neighbours_num = len(graph[v])

        # set the neighbours array and shuffle it
        neighbours = list(graph[v])
        random.shuffle(neighbours)

        adjacent_groups = False

        # set the number of iterations and the adjacent flag
        if v_group == t_group or abs(v_group - t_group) == 0 or abs(v_group-t_group) == m - 1:
            num_iterations = max(int(round(neighbours_num * avg_num_edges_p / (avg_num_edges_p + avg_num_edges_q))),1)
            adjacent_groups = True
        else:
            num_iterations = max(int(round(neighbours_num * avg_num_edges_q / (avg_num_edges_p + avg_num_edges_q))),1)

        # keep track of adjacent neighbours
        adjacent_neighbour = -1
        closest_not_adjacent_v = -1

        for i in range(num_iterations):
            u = neighbours[i]
            u_group = u // k

            search_time += 1

            if t == u:
                v = u
                moved = True
                break

            # check if t is adjacent to one of the neighbours
            elif u_group == t_group or abs(u_group - t_group) == 0 or abs(u_group-t_group) == m - 1:
                # if v is not in adjacent group, u could be connected to t with probability p
                if not adjacent_groups:
                    v = u
                    moved = True
                    break
                else:
                    # v is adjacent so there is a chance that t is connected to v
                    # so record u and look at other neighbours
                    adjacent_neighbour = u
            # if t is not adjacent to u and v, see if u is closer to t than v in terms of groups
            else:
                abs_u_t = abs(t_group - u_group)
                abs_v_t = abs(t_group - v_group)

                diff_u_t = min(abs_u_t, m - abs_u_t)
                diff_v_t = min(abs_v_t, m - abs_v_t)

                # as this is a ring both sides must be considered
                if diff_u_t < diff_v_t:
                    closest_not_adjacent_v = u
                    if random.random() < q:
                        v = u
                        moved = True
                        break
                # u is not closer to t than v and there hasn't been a vertex found to be close
                elif closest_not_adjacent_v == -1:
                    # if u and v are in the same group
                    if diff_v_t == diff_u_t:
                        closest_not_adjacent_v = u
                    # if non of the neighbours are close to t or in the same group as v
                    elif adjacent_neighbour == -1:
                        closest_not_adjacent_v = u
        if not moved:
            # no adjacent neighbour to t is found in terms of group
            if adjacent_neighbour == -1:
                v = closest_not_adjacent_v
            else:
                v = adjacent_neighbour
    return search_time


def search_ring_group_graph_v3(graph, m, k, p, q, v, t):
    search_time = 0

    t_group = t // k

    while True:
        # reset moved flag
        moved = False

        # if the vertex id equals the target id stop
        if v == t: break

        # set d(v) for the vertex (number of neighbours)
        neighbours_num = len(graph[v])

        # set the neighbours array and shuffle it
        neighbours = list(graph[v])
        random.shuffle(neighbours)

        adjacent_groups = False
        adjacent_neighbour = -1

        for i in range(neighbours_num):
            u = neighbours[i]
            u_group = u // k

            search_time += 1

            rand = random.random()

            if t == u:
                v = u
                moved = True
                break

            # check if t is adjacent to one of the neighbours
            elif u_group == t_group or abs(u_group - t_group) == 0 or abs(u_group-t_group) == m - 1:
                # if v is not in adjacent group
                if not adjacent_groups and rand < p:
                    v = u
                    moved = True
                    break
                else:
                    # v is adjacent so there is a chance that t is connected to v
                    # so record u and look at other neighbours
                    adjacent_neighbour = u
            # if t is not adjacent to u and v, see if u is closer to t than v in terms of groups
            else:
                if rand < q:
                    v = u
                    moved = True
                    break

        if not moved:
            # no adjacent neighbour to t is found in terms of group
            if adjacent_neighbour == -1:
                v = neighbours[neighbours_num - 1]
            else:
                v = adjacent_neighbour
    return search_time


def run_search_ring_group_graph():
    m = 10
    k = 50
    p = 0.25
    q = 0.05

    # sample of 100 pair of vertices
    t1 = 100
    # number of iterations for searching in one pair of vertices
    t2 = 1000
    # number of graph instances
    t3 = 20

    avg_search_time = {}

    for l in range(t3):
        graph = make_nx_ring_group_graph(m,k,p,q)

        if not nx.is_connected(graph):
            print('not connected')
            continue

        graph = nx.to_dict_of_lists(graph)

        print(graph)

        # search_time = {}
        search_time_v3 = {}

        for i in range(t1):
            start_v = random.randint(0, m * k - 1)
            target_u = random.randint(0, m * k - 1)
            print(l,i, start_v, target_u)
            # avg_search_time = 0
            avg_search_time_v3 = 0

            for j in range(t2):
                # avg_search_time += search_ring_group_graph(graph,m,k,p,q,start_v,target_u)
                avg_search_time_v3 += search_ring_group_graph_v3(graph,m,k,p,q,start_v,target_u)

            # search_time[(start_v,target_u,i)] = avg_search_time / t2
            search_time_v3[(start_v,target_u,i)] = avg_search_time_v3 / t2

        avg_search_time[l] = int(round(sum(search_time_v3.values())/ t1))

    return avg_search_time


run_search_ring_group_graph()

# n = 1000
# p = 0.1
#
# graph = make_random_graph(n,p)
#
# search_time = {}
# avg_search_time = 0
# search_time_v2 = {}
# avg_search_time_v2 = 0
# start_v = random.randint(0,n - 1)
# target_u = random.choice(list(graph[start_v]))
#
# for i in range(1000):
#     print(i)
#     search_time[(start_v,target_u, i)] = search_random_graph(graph,n,p,start_v,target_u)
#     search_time_v2[(start_v,target_u, i)] = search_random_graph_v2(graph,n,p,start_v,target_u)
#     avg_search_time += search_time[(start_v,target_u, i)]
#     avg_search_time_v2 += search_time_v2[(start_v,target_u, i)]
#
# print(graph)
# print(search_time)
# print(avg_search_time / 1000)
# print(search_time_v2)
# print(avg_search_time_v2 / 1000)






