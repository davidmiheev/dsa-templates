# Dijkstra algo
# Finding the shortest path in a weighted graph

def shortestPath(self, node1: int, node2: int) -> int:
    seen = {}
    pq = [(0, node1)]
    while pq:
        dist, vert = heappop(pq)
        seen[vert] = dist
        if vert == node2: return dist
        for neighbor, cost in self.graph[vert]:
            if neighbor in seen and seen[neighbor] <= dist+cost: continue
            heappush(pq, (dist+cost, neighbor))
        
    return -1
