Algos on graph
================

Construction of a graph from a list of edges
--------------------------------------------
.. code-block:: python

    from collections import defaultdict
    from collections import deque

    edges = [[1,2], [2,3], [5,2], [1, 5]]
    graph = defaultdict(list)

    for a, b in edges:
        graph[a] += [b]
        graph[b] += [a] # this line for undirected graphs only

BFS
---
Breadth-first search is a particular search on a graph.
Breadth-first search is useful when you need to find the length of a shortest path from one vertex to another (in unweighted graph).
Also, with BFS you can traverse a graph
Breadth-first search implementation is non-recursive.

.. autofunction:: bfs.bfs.bfs

Maze problem:
    For instance, you have a matrix `(m x n)` with zeros and ones,
    you can move from the current cell to one of four adjacent cells whenever that adjacent cell is filled by zero
    initial position is `(0, 0)`, find the length of a shortest path to the right-bottom cell `(m - 1, n - 1)`,
    if you can't reach the right-bottom cell, return `-1`

0 1 0 0 0

0 1 0 1 0

0 1 0 1 0

0 0 0 1 0

answer: 13

.. autofunction:: bfs.bfs.bfs_matrix



DFS
---
Depth-First Search (DFS) is an algorithm for traversing or searching tree or graph data structures.
One starts at the root (selecting some arbitrary node as the root in the case of a graph)
and explores as far as possible along each branch before backtracking.
DFS is typically implemented recursively, but can be implemented iteratively as well.

Time complexity: :math:`O(V + E)`

Space complexity: :math:`O(V)`,
where V is the number of vertices and E is the number of edges in the graph.

Applications:
    * Topological sorting

    * Finding connected components

    * Finding bridges and articulation points

    * Finding strongly connected components

    * Solving puzzles with only one solution, such as mazes

.. automodule:: dfs.dfs
    :members:
    :undoc-members:
    :show-inheritance:

Dijkstra
--------

.. automodule:: dijkstra.dijkstra
    :members:
    :undoc-members:
    :show-inheritance:
