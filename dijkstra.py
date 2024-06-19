# Dijkstra algo
# Finding the shortest path in a weighted graph

from collections import defaultdict
from heapq import heappop, heappush
from typing import List

def shortestPath(edges: List[List[int]], costs: List[int], start: int, end: int) -> int:
    '''
    Find the shortest path in a weighted graph from the start node to the end node.
    Dijkstra's algorithm is a greedy algorithm that finds the shortest path between nodes in a graph.

    Time complexity: :math:`O(E \log V)`
    Space complexity: :math:`O(E + V)`

    :type edges: List[List[int]]
    :type costs: List[int]
    :type start: int
    :type end: int
    :rtype: int
    '''
    graph = defaultdict(list)
    for j, (a, b) in enumerate(edges):
        graph[a] += [(b, costs[j])]
        graph[b] += [(a, costs[j])]

    seen = {start: 0}
    pq = [(0, start)]
    while pq:
        dist, vert = heappop(pq)
        if vert == end: return dist
        for neighbor, cost in graph[vert]:
            if neighbor in seen and seen[neighbor] <= dist+cost: continue
            seen[neighbor] = dist+cost
            heappush(pq, (dist+cost, neighbor))

    return -1
