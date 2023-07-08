# Topological sort is a special ordering of vertices in a directed acyclic graph (DAG)
# where for each directed edge from vertex A to vertex B, vertex A comes before vertex B in the ordering.
#
# Topological sort is useful when you need to perform a series of tasks
# where each task has some dependencies on other tasks.
#
# For example, if you have a list of tasks to perform,
# and some of those tasks depend on other tasks being completed first,
# you could represent this as a graph with a directed edge from task A to task B
# if task A must be completed before task B.
#
# Then, a topological sort of the graph would give you an ordering of the tasks
# that respects the dependencies.
#
# Topological sort can be used to detect cycles in a graph.
# If the graph has a cycle, then there is no topological sort.
#
# Topological sort can be used to find the shortest path in a weighted graph (if there are no cycles).
#
# Topological sort can be used to find the strongly connected components of a graph.
#
# Topological sort can be used to order the compilation tasks for a program


# Path: top_sort.py

def top_sort(graph):
    seen = set()
    current_path = set()
    order = []
    def dfs(vert):
        seen.add(vert)
        current_path.add(vert)
        for adj in graph[vert]:
            if adj not in seen:
                dfs(adj)
            elif adj in current_path: # check if graph has a cycle
                raise ValueError('Cycle detected') # if graph has a cycle, then there is no topological sort

        current_path.discard(vert)
        order.append(vert)


    for vert in graph:
        if vert in seen: continue
        dfs(vert)

    return order[::-1]

graph = {0: [1, 2], 1: [3, 0], 2: [3], 3: [4], 4: []}
print(top_sort(graph))
