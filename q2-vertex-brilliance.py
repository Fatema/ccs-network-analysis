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
        print('kstar for ',source, ' ', kstar)
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


graph = {'a': {3,7,1}, 1: {'a',3,2}, 2: {3,1}, 3: {'a', 1,2,4,6}, 4: {3,5}, 5: {4,7,6}, 6: {3,7,5,8}, 7: {'a',6,5,8,9}, 8: {6,7}, 9: {7}}

for v in graph:
    brilliance = vertex_brilliance(graph, v)
    print('b(',v,')=',brilliance)
