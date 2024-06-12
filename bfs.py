# Breadth-first search is a particular search on a graph.
# Breadth-first search is useful when you need to find the length of a shortest path from one vertex to another (in unweighted graph).
# Also, with BFS you can traverse a graph
# Time complexity O(n + m), where n is the number of vertices and m is thr number of edges
# Breadth-first search implementation is non-recursive.
# Here is an example:
from collections import defaultdict
from collections import deque

edges = [[1,2], [2,3], [5,2], [1, 5]]
target = 5
graph = defaultdict(list)
for a, b in edges:
    graph[a] += [b]
    graph[b] += [a] # this line for undirected graphs only

def bfs(start):
    seen = set()
    q = deque([(start, 0)])
    while q:
        vert, dist = q.popleft()
        seen.add(vert)
        if vert == target: return dist
        for adj in graph[vert]:
            if adj in seen: continue
            q.append((adj, dist+1))

    return -1

print(bfs(1), bfs(3))

# If you use BFS for search on matrices, you need modification:
# Labirint problem:
# For instance, you have a matrix (m x n) with zeros and ones,
# you can move from the current cell to one of four adjacent cells whenever that adjacent cell is filled by zero
# initial position is (0, 0), find the length of a shortest path to the right-bottom cell (m-1, n-1),
# if you can't reach the right-bottom cell, return -1

# 0 1 0 0 0
# 0 1 0 1 0
# 0 1 0 1 0
# 0 0 0 1 0

# answer: 13

labirint = [[0, 1, 0, 0, 0], [0, 1, 0, 1, 0], [0, 0, 0, 1, 0], [0, 0, 0, 1, 0]]
m, n = len(labirint), len(labirint[0])
dirs = [(1, 0), (0, 1), (-1, 0), (0, -1)]
start = (0, 0)
def bfs_matrix(start):
    seen = set([start])
    q = deque([(start[0], start[1], 0)])
    j = 0
    while q:
        j += 1
        r, c, dist = q.popleft()
        # print(j, (r, c))
        # seen.add((r, c))
        if (r, c) == (m-1, n-1): return dist
        for dr, dc in dirs:
            if r+dr >= m or r+dr < 0 or c+dc >= n or c+dc < 0: continue
            if (r+dr, c+dc) in seen or labirint[r+dr][c+dc]: continue
            q.append((r+dr, c+dc, dist+1))
            seen.add((r+dr, c+dc)) # to escape repetitions in queue in matrix case

    return -1

# print(seen)
print(bfs_matrix(start))
