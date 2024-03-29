import math
import random
import networkx as nx
import matplotlib.pyplot as plt
from tqdm import tqdm


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

    # make it less likely to more to a node if p is more 0.5 as the graph is more connected
    if prob > 0.5:
        prob = 1 - prob

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
                break

        if not moved:
            s = neighbours[neighbours_num - 1]
    return search_time


def search_random_graph_v2(graph, total_nodes, prob, s, t):

    search_time = 0

    # make it less likely to more to a node if p is more 0.5 as the graph is more connected
    if prob > 0.5:
        prob = 1 - prob

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
                break

        if not moved:
            s = neighbours[neighbours_num - 1]
    return search_time


def search_ring_group_graph(graph, m, k, p, q, v, t):
    search_time = 0

    t_group = t // k

    # make it less likely to more to a node if p is more 0.5 as the graph is more connected
    if p > 0.5:
        p = 1 - p

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


def search_ring_group_graph_v2(graph, m, k, p, q, v, t):
    search_time = 0

    t_group = t // k

    prob_close = p

    # make it less likely to more to a node if p is more 0.5 as the graph is more connected
    if p > 0.5:
        p = 1 - p


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

            rand_num = random.random()

            search_time += 1

            if t == u:
                v = u
                moved = True
                break

            # check if t is adjacent to one of the neighbours
            elif u_group == t_group or abs(u_group - t_group) == 0 or abs(u_group-t_group) == m - 1:
                # if v is not in adjacent group, u could be connected to t with probability p
                if not adjacent_groups and rand_num < p:
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
                    if rand_num < prob_close and adjacent_neighbour == -1:
                        v = u
                        moved = True
                        break

                if rand_num < q:
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
    plt.ylabel('normalized number of instances')
    plt.title('Search time for ' + plot_name)
    plt.plot(xdata, ydata, marker='.', linestyle='None', color='g')
    plt.savefig('distributions/q3/search-time/' + plot_file_name + '.png')


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

    plt.ylim(0,1000)

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


def box_plot_search_time(data, xaxis, plot_file_name, plot_name):
    # clears plot
    plt.clf()

    fig, ax = plt.subplots()

    plt.ylabel('search time')
    plt.title('Search time for ' + plot_name)

    red_square = dict(markerfacecolor='r', marker='.')
    ax.boxplot(data, flierprops=red_square, whis=0.75)

    ax.set_xticklabels(xaxis)

    plt.savefig('distributions/q3/box/' + plot_file_name + '.png')


def plot_n_search_time_num_instances(all_times, plot_file_name, plot_name):
    colours = ['b', 'g', 'r', 'c', 'm', 'y', 'k']

    i = 0
    num_colours = len(colours)

    # clears plot
    plt.clf()

    # plot degree distribution
    plt.xlabel('search time')
    plt.ylabel('normalized number of instances')
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
    # params = [(20,50,0.45,0.01), (20,50,0.45,0.001), (20,50,0.45,0.0001)]
    # params = [(50,20,0.45,0.01), (50,20,0.45,0.001), (50,20,0.45,0.0001)]
    # params = [(50,20,0.45,0.01), (20,20,0.45,0.01), (20,50,0.45,0.01)]
    # params = [(50,20,0.45,0.01), (50,20,0.45,0.1)]
    # params = [(20,50,0.45,0.01), (20,50,0.45,0.1)]
    params = [(50,20,0.45,0.01), (20,20,0.45,0.01), (10,20,0.45,0.01)]
    # params = [(20,40,0.45,0.01), (20,20,0.45,0.01), (20,10,0.45,0.01)]
    # params = [(50,20,0.35,0.01), (50,20,0.45,0.01), (50,20,0.55,0.01)]
    # params = [(6,50,0.35,0.01), (6,50,0.45,0.01), (6,50,0.55,0.01)]
    # params = [(6,100,0.45,0.01), (6,100,0.45,0.1)]

    # sample of 30000 pair of vertices
    t1 = 30000
    # number of iterations for searching in one pair of vertices
    t2 = 1
    # number of graph instances
    t3 = 10

    all_times = []
    all_times_v2 = []

    for param in params:
        m, k, p, q = param

        cumulative_dist = {}

        for l in range(t3):
            graph = make_nx_ring_group_graph(m,k,p,q)

            if not nx.is_connected(graph):
                print('not connected')
                continue

            graph = nx.to_dict_of_lists(graph)

            avg_search_time_v1 = {}

            for i in tqdm(range(t1)):
                start_v = random.randint(0, m * k - 1)
                target_u = random.randint(0, m * k - 1)

                avg_search_time_v1[i] = 0

                for j in range(t2):
                    avg_search_time_v1[i] += search_ring_group_graph(graph,m,k,p,q,start_v,target_u)

                avg_search_time_v1[i] = avg_search_time_v1[i] / t2

            search_time_dist = normalized_search_time_distribution(avg_search_time_v1,t1)
            plot_search_time_num_instances(search_time_dist, 'ring_group_graph-' + str(m) + '-'
                                               + str(k) + '-' + str(p) + '-' + str(q) + '-search_time_v1-' + str(l),
                                               'Ring Group Graph')
            box_plot_search_time(list(avg_search_time_v1.values()), [str(param)], 'box-ring_group_graph-' + str(m) + '-'
                                               + str(k) + '-' + str(p) + '-' + str(q) + '-search_time_v1-' + str(l),
                                               'Ring Group Graph')

            print(sum(avg_search_time_v1.values())/t1)

            for search_time in search_time_dist:
                if search_time in cumulative_dist:
                    cumulative_dist[search_time] += search_time_dist[search_time]
                else:
                    cumulative_dist[search_time] = search_time_dist[search_time]

        average_dist = {}

        for deg in cumulative_dist:
            average_dist[deg] = cumulative_dist[deg] / t3

        all_times += [(param, average_dist)]

        avg_dist_box = [k for k in average_dist for i in range(math.floor(average_dist[k] * t1))]

        all_times_v2 += [avg_dist_box]

        plot_search_time_num_instances(average_dist, 'ring_group_graph-' + str(m) + '-'
                                       + str(k) + '-' + str(p) + '-' + str(q) + '-search_time_v1-avg',
                                       'Ring Group Graph')
        box_plot_search_time(avg_dist_box, [str(param)], 'box-ring_group_graph-' + str(m) + '-'
                             + str(k) + '-' + str(p) + '-' + str(q) + '-search_time_v1-avg',
                             'Ring Group Graph')

    plot_n_search_time_num_instances(all_times, 'ring_group_graph-50-20-10-search_time_v1-avg', 'Ring Group Graph')
    box_plot_search_time(all_times_v2, [str(param) for param in params],
                         'box-ring_group_graph-50-20-10-search_time_v1-avg',
                         'Ring Group Graph')

    return all_times


def run_search_random_graph():
    # params = [(100,0.1), (250,0.1), (500,0.1)]
    params = [(250,0.1), (250,0.2), (250,0.3)]

    # sample of 10000 pair of vertices
    t1 = 10000
    # number of iterations for searching in one pair of vertices
    t2 = 1
    # number of graph instances
    t3 = 10

    all_times = []
    all_times_v2 = []

    for param in params:
        n, p = param

        cumulative_dist = {}

        for l in range(t3):
            graph = make_nx_random_graph(n, p)

            if not nx.is_connected(graph):
                print('not connected')
                continue

            graph = nx.to_dict_of_lists(graph)

            search_time = {}

            for i in tqdm(range(t1)):
                start_v = random.randint(0, n - 1)
                target_u = random.randint(0, n - 1)
                avg_search_time_v1 = 0

                for j in range(t2):
                    avg_search_time_v1 += search_random_graph(graph,n,p,start_v,target_u)

                search_time[i] = avg_search_time_v1 / t2

            search_time_dist = normalized_search_time_distribution(search_time, t1)

            plot_search_time_num_instances(search_time_dist,'random_graph-' + str(n) + '-' + str(p) + '-search_time_v1',
                                       'Random Graph')
            box_plot_search_time(list(search_time.values()), [str(param)], 'box-random_graph-' + str(n) + '-'
                                 + '-' + str(p) + '-search_time_v1-' + str(l),
                                 'Random Graph')

            for search_time in search_time_dist:
                if search_time in cumulative_dist:
                    cumulative_dist[search_time] += search_time_dist[search_time]
                else:
                    cumulative_dist[search_time] = search_time_dist[search_time]

        average_dist = {}

        for deg in cumulative_dist:
                average_dist[deg] = cumulative_dist[deg] / t3

        all_times += [(param, average_dist)]

        avg_dist_box = [k for k in average_dist for i in range(math.floor(average_dist[k] * t1))]

        all_times_v2 += [avg_dist_box]

        plot_search_time_num_instances(average_dist, 'random_graph-' + str(n) + '-'
                                        + str(p) + '-search_time_v1-avg',
                                       'Random Graph')
        box_plot_search_time(avg_dist_box, [str(param)], 'box-random_graph-' + str(n) + '-'
                              + str(p) + '-search_time_v1-avg',
                             'Random Graph')

    plot_n_search_time_num_instances(all_times, 'random_graph-0.1-search_time_v1', 'Random Graph')
    box_plot_search_time(all_times_v2, [str(param) for param in params], 'box-random_graph-n-search_time_v1-avg',
                         'Random Graph')

    return all_times


run_search_ring_group_graph()








