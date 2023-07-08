# Dijkstra algo
# Finding the shortest path in a weighted graph

from collections import defaultdict
from heapq import heappop, heappush
from typing import List

def shortestPath(edges: List[List[int]], costs: List[int], start: int, end: int) -> int:
    graph = defaultdict(list)
    for j, (a, b) in enumerate(edges):
        graph[a] += [(b, costs[j])]
        graph[b] += [(a, costs[j])]

    seen = {}
    pq = [(0, start)]
    while pq:
        dist, vert = heappop(pq)
        seen[vert] = dist
        if vert == end: return dist
        for neighbor, cost in graph[vert]:
            if neighbor in seen and seen[neighbor] <= dist+cost: continue
            heappush(pq, (dist+cost, neighbor))

    return -1
