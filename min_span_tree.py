from collections import defaultdict
from heapq import heappop, heappush, heapify
from typing import List
from union_find import UnionFind

print('min_span_tree')

def minimumCost(n: int, connections: List[List[int]]) -> int:
    '''
    Find the minimum cost to connect all the cities.
    Kruskal's algorithm is a greedy algorithm that finds a minimum spanning tree for a connected weighted undirected graph,
    using Union-Find data structure

    Time complexity: :math:`O(E \log E)`, where E is the number of edges in the graph

    Space complexity: :math:`O(E)`

    :type n: int
    :type connections: List[List[int]]
    :rtype: int
    '''
    uf = UnionFind(set(range(1, n+1)))
    pq = []
    for x, y, cost in connections:
        pq.append((cost, x, y))

    heapify(pq)
    ans = 0
    while pq and n > 1:
        cost, x, y = heappop(pq)
        if not uf.are_connected(x, y):
            uf.join(x, y)
            n -= 1
            ans += cost

    return ans if n == 1 else -1

connections = [[1,2,5],[1,3,6],[2,3,1]]
n = 3
print(minimumCost(n, connections))

connections = [[1,2,3],[3,4,4]]
n = 4
print(minimumCost(n, connections))

connections = [[1,2,1],[2,3,2],[1,3,2]]
n = 3
print(minimumCost(n, connections))
