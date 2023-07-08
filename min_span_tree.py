# Minimum Spanning Tree (MST) is a subset of the edges of a connected, edge-weighted undirected graph that connects all the vertices together,
# without any cycles and with the minimum possible total edge weight.
#
#
# # Path: min_span_tree.py

from collections import defaultdict
from heapq import heappop, heappush, heapify
from typing import List


# Kruskal's algorithm
# Time complexity: O(E log E)
# Space complexity: O(E)
# where E is the number of edges in the graph
#
# Kruskal's algorithm is a greedy algorithm that finds a minimum spanning tree for a connected weighted undirected graph.
# Using Union-Find data structure.

class UnionFind:
    def __init__(self, vertices: int):
        self.group = defaultdict(int)
        self.rank = defaultdict(int)
        for i in range(1, vertices+1):
            self.group[i] = i

    def find(self, node: int) -> int:
        if self.group[node] != node:
            self.group[node] = self.find(self.group[node])
        return self.group[node]

    def join(self, node1: int, node2: int):
        group1 = self.find(node1)
        group2 = self.find(node2)

        if group1 == group2: return

        if self.rank[group1] > self.rank[group2]:
            self.group[group2] = group1
        elif self.rank[group1] < self.rank[group2]:
            self.group[group1] = group2
        else:
            self.group[group1] = group2
            self.rank[group2] += 1

    def are_connected(self, node1: int, node2: int) -> bool:
        return self.find(node1) == self.find(node2)


def minimumCost(n: int, connections: List[List[int]]) -> int:
    uf = UnionFind(n)
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
