import random
import networkx as nx
import matplotlib.pyplot as plt


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


def make_nx_random_graph(num_vertices, p):
    # initialize the graph
    random_graph = nx.Graph()

    for v in range(num_vertices):
        for u in range(v + 1, num_vertices):
            random_number = random.random()
            if random_number < p:
                random_graph.add_edge(v, u)

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
    average_neighbour_num = int(round(total_nodes * prob))
    random_move = 0

    while True:
        moved = False

        if s == t: break

        neighbours = list(graph[s])
        neighbours_num = len(neighbours)

        random.shuffle(neighbours)

        # at max iterate until the average number of neighbours
        num_iterations = min(average_neighbour_num, neighbours_num)

        for i in range(num_iterations):
            u = neighbours[i]
            search_time += 1
            if t == u:
                s = u
                moved = True
                break
            # probability that an edge might exist between u and t
            if random.random() < prob:
                s = u
                moved = True
                random_move += 1
                break

        if not moved:
            s = neighbours[neighbours_num - 1]
    # print('v1 total time taken to find target', search_time)
    # print('number of random moves', random_move)
    return search_time


def search_random_graph_v2(graph, total_nodes, prob, s, t):

    search_time = 0
    random_move = 0

    while True:
        moved = False

        if s == t: break

        neighbours = list(graph[s])
        neighbours_num = len(neighbours)

        random.shuffle(neighbours)

        # at max iterate until the average number of neighbours
        num_iterations = neighbours_num

        for i in range(num_iterations):
            u = neighbours[i]
            search_time += 1
            if t == u:
                s = u
                moved = True
                break
            # probability that an edge might exist between u and t
            if random.random() < prob:
                s = u
                moved = True
                random_move += 1
                break

        if not moved:
            s = neighbours[neighbours_num - 1]
    # print('v2 total time taken to find target', search_time)
    # print('number of random moves', random_move)
    return search_time


def search_ring_group_graph(graph, m, k, p, q, v, t):
    search_time = 0

    t_group = t // k

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

        # keep track of adjacent neighbours
        adjacent_neighbour = -1
        closest_not_adjacent_v = -1

        for i in range(neighbours_num):
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

        if not moved:
            if adjacent_neighbour == -1:
                if closest_not_adjacent_v == -1:
                    v = neighbours[neighbours_num - 1]
                else:
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

def search_time_distribution(search_time_instances):
    search_distribution = {}
    # consider each vertex
    for instance in search_time_instances:
        # update degree_distribution
        search_time = search_time_instances[instance]
        if search_time in search_distribution:
            search_distribution[search_time] += 1
        else:
            search_distribution[search_time] = 1
    return search_distribution


def normalized_search_time_distribution(search_time_instances, num_instances):
    """Takes a graph and computes the normalized distribution of the degrees of the graph.
    Returns a dictionary whose keys correspond to in-degrees of nodes in the graph and values are the
    fraction of nodes with that in-degree.
    Degrees with no corresponding nodes in the graph are not included in the dictionary."""
    unnormalized_dist = search_time_distribution(search_time_instances)
    normalized_dist = {}
    for degree in unnormalized_dist:
        normalized_dist[degree] = 1.0 * unnormalized_dist[degree] / num_instances
    return normalized_dist


def plot_search_time_num_instances(search_time_ins, plot_file_name, plot_name):
    # clears plot
    plt.clf()

    xdata = []
    ydata = []

    for search_time in search_time_ins:
        xdata += [search_time]
        ydata += [search_time_ins[search_time]]

    # plot degree distribution
    plt.xlabel('search time')
    plt.ylabel('number of instances')
    plt.title('Search time for ' + plot_name)
    plt.plot(xdata, ydata, marker='.', linestyle='None', color='g')
    plt.savefig('distributions/q3/' + plot_file_name + '.png')


def plot_search_time_q(search_time_ins, plot_file_name, plot_name):
    # clears plot
    plt.clf()

    xdata = []
    ydata = []

    for search_time in search_time_ins:
        xdata += [search_time]
        ydata += [search_time_ins[search_time]]

    # plot degree distribution
    plt.xlabel('q')
    plt.ylabel('search time')
    plt.title('Search time vs q for ' + plot_name)
    plt.plot(xdata, ydata, marker='.', linestyle='None', color='g')
    plt.savefig('distributions/q3/' + plot_file_name + '.png')


def plot_n_search_time_p(p_search_time_dists, plot_file_name, plot_name):
    # clears plot
    plt.clf()

    colours = ['b', 'g', 'r', 'c', 'm', 'y', 'k']

    i = 0
    num_colours = len(colours)

    # plot degree distribution
    plt.xlabel('p')
    plt.ylabel('search time')
    plt.title('Search time vs p for ' + plot_name)
    for dist in p_search_time_dists:
        q, p_search_time = dist
        xdata = []
        ydata = []
        for p in p_search_time:
            xdata += [p]
            ydata += [p_search_time[p]]
        plt.plot(xdata, ydata, marker='.', linestyle='-', color=colours[i % num_colours], label='q = ' + str(q))
        i += 1
    plt.legend(loc='upper right')
    plt.savefig('distributions/q3/' + plot_file_name + '.png')


def plot_n_search_time_num_instances(all_times, plot_file_name, plot_name):
    colours = ['b', 'g', 'r', 'c', 'm', 'y', 'k']

    i = 0
    num_colours = len(colours)

    # clears plot
    plt.clf()

    # plot degree distribution
    plt.xlabel('search time')
    plt.ylabel('number of instances')
    plt.title('Search time for ' + plot_name)
    for time in all_times:
        param, search_time_ins = time
        xdata = []
        ydata = []
        for search_time in search_time_ins:
            xdata += [search_time]
            ydata += [search_time_ins[search_time]]
        plt.plot(xdata, ydata, marker='.', linestyle='None', color=colours[i % num_colours], label= 'params = ' + str(param))
        i += 1
    plt.legend(loc='upper right')
    plt.savefig('distributions/q3/' + plot_file_name + '.png')


def run_search_ring_group_graph():
    params = [(20,50,0.45,0.01), (20,50,0.45,0.001), (20,50,0.45,0.0001)]

    # sample of 50 pair of vertices
    t1 = 50
    # number of iterations for searching in one pair of vertices
    t2 = 50
    # number of graph instances
    t3 = 50

    all_times = []
    search_time_q = {}

    for param in params:
        m, k, p, q = param

        avg_search_time = {}

        for l in range(t3):
            graph = make_nx_ring_group_graph(m,k,p,q)

            if not nx.is_connected(graph):
                print('not connected')
                continue

            graph = nx.to_dict_of_lists(graph)

            print(graph)

            search_time_v1 = {}
            # search_time_v3 = {}

            for i in range(t1):
                start_v = random.randint(0, m * k - 1)
                target_u = random.randint(0, m * k - 1)
                print(l,i, start_v, target_u)
                avg_search_time_v1 = 0
                # avg_search_time_v3 = 0

                for j in range(t2):
                    avg_search_time_v1 += search_ring_group_graph(graph,m,k,p,q,start_v,target_u)
                    # avg_search_time_v3 += search_ring_group_graph_v3(graph,m,k,p,q,start_v,target_u)

                search_time_v1[(start_v,target_u,i)] = avg_search_time_v1 / t2
                # search_time_v3[(start_v,target_u,i)] = avg_search_time_v3 / t2

            # avg_search_time[l] = int(round(sum(search_time_v3.values())/ t1))
            avg_search_time[l] = int(round(sum(search_time_v1.values())/ t1))

        search_time_dist = search_time_distribution(avg_search_time)

        plot_search_time_num_instances(search_time_dist, 'ring_group_graph-' + str(m) + '-'
                                   + str(k) + '-' + str(p)+ '-' + str(q) + '-search_time_v1', 'Ring Group Graph')

        all_times += [(param, search_time_dist)]
        search_time_q[q] = sum(avg_search_time.values()) / t3

    plot_n_search_time_num_instances(all_times, 'ring_group_graph-search_time_v1', 'Ring Group Graph')
    plot_search_time_q(search_time_q, 'ring_group_graph-search_time-q_v1', 'Ring Group Graph')
    return all_times


def run_search_ring_group_graph_p():
    params = [(20,50,0.01), (20,50,0.001), (20,50,0.0001)]
    # params = [(50,20,0.01), (50,20,0.001), (50,20,0.0001)]

    # sample of 50 pair of vertices
    t1 = 50
    # number of iterations for searching in one pair of vertices
    t2 = 10
    # number of graph instances
    t3 = 20

    all_times = []
    search_time_p_dists = []

    for param in params:
        m, k, q = param

        search_time_p = {}

        p = round(q + 1 / 10, 7)

        print(m, k, p, q)

        while p < 1:

            avg_search_time = {}

            for l in range(t3):
                graph = make_nx_ring_group_graph(m,k,p,q)

                if not nx.is_connected(graph):
                    print('not connected')
                    continue

                graph = nx.to_dict_of_lists(graph)

                search_time_v1 = {}

                for i in range(t1):
                    start_v = random.randint(0, m * k - 1)
                    target_u = random.randint(0, m * k - 1)

                    avg_search_time_v1 = 0

                    for j in range(t2):
                        avg_search_time_v1 += search_ring_group_graph(graph,m,k,p,q,start_v,target_u)

                    search_time_v1[(start_v,target_u,i)] = avg_search_time_v1 / t2

                avg_search_time[l] = int(round(sum(search_time_v1.values())/ t1))

            if avg_search_time.values() is not {}:
                search_time_p[p] = sum(avg_search_time.values()) / t3
                print(p, search_time_p[p])

            p = round(p * 1.2, 7)

        search_time_p_dists += [(q,search_time_p)]

    plot_n_search_time_p(search_time_p_dists, 'ring_group_graph-search_time-p-kGreater', 'Ring Group Graph')

    return all_times


def run_search_random_graph():
    n = 1000
    p = 0.1

    # sample of 50 pair of vertices
    t1 = 50
    # number of iterations for searching in one pair of vertices
    t2 = 10
    # number of graph instances
    t3 = 20

    ins_avg_search_time = {}
    ins_avg_search_time_v3 = {}

    for l in range(t3):
        graph = make_nx_random_graph(n, p)

        if not nx.is_connected(graph):
            print('not connected')
            continue

        graph = nx.to_dict_of_lists(graph)

        search_time = {}
        search_time_v3 = {}

        for i in range(t1):
            start_v = random.randint(0, n - 1)
            target_u = random.randint(0, n - 1)
            print(l, i, start_v, target_u)
            avg_search_time_v1 = 0
            avg_search_time_v3 = 0

            for j in range(t2):
                avg_search_time_v1 += search_random_graph(graph,n,p,start_v,target_u)
                avg_search_time_v3 += search_random_graph_v2(graph,n,p,start_v,target_u)

            search_time[(start_v,target_u,i)] = avg_search_time_v1 / t2
            search_time_v3[(start_v, target_u, i)] = avg_search_time_v3 / t2

        ins_avg_search_time[l] = int(round(sum(search_time.values()) / t1))
        ins_avg_search_time_v3[l] = int(round(sum(search_time_v3.values()) / t1))

    search_time_dist = normalized_search_time_distribution(ins_avg_search_time, t3)
    search_time_dist_v3 = normalized_search_time_distribution(ins_avg_search_time, t3)

    plot_search_time_num_instances(search_time_dist,'random_graph-' + str(n) + '-' + str(p) + '-search_time',
                                   'Random Graph')
    plot_search_time_num_instances(search_time_dist_v3,'random_graph-' + str(n) + '-' + str(p) + '-search_time_v3',
                                   'Random Graph')

    return ins_avg_search_time, ins_avg_search_time_v3


run_search_ring_group_graph_p()








