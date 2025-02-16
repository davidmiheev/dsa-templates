print('top_sort')


def top_sort(graph):
    '''
    Topological sort of a directed acyclic graph (DAG)
    Topological sort is useful when you need to perform a series of tasks
    where each task has some dependencies on other tasks.

    Time complexity: :math:`O(V + E)`, where V is the number of vertices and E is the number of edges

    Space complexity: :math:`O(V)`

    :param graph: dict
    :return: list
    '''
    seen = set()
    current_path = set()
    order = []
    def dfs(vert):
        nonlocal order

        seen.add(vert)
        current_path.add(vert)

        for adj in graph[vert]:
            
            if adj not in seen:
                dfs(adj)
            elif adj in current_path: # check if graph has a cycle
                print('Cycle detected')
                order = ['Cycle detected']
                return # if graph has a cycle, then there is no topological sort

        current_path.discard(vert)
        order.append(vert)


    for vert in graph:
        if vert in seen: 
            continue
        dfs(vert)

    return order[::-1]

graph = {0: [1, 2], 1: [3], 2: [3], 3: [4], 4: [0]}
print(top_sort(graph))
