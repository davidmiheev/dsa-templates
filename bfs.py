from collections import defaultdict
from collections import deque

print('BFS')

edges = [[1,2], [2,3], [5,2], [1, 5]]
graph = defaultdict(list)

for a, b in edges:
    graph[a] += [b]
    graph[b] += [a] # this line for undirected graphs only

def bfs(start, target):
    '''
    Find the length of a shortest path from the start node to the target node in an unweighted graph.

    Time complexity: :math:`O(n + m)`, where n is the number of vertices and m is the number of edges
    '''
    seen = set()
    q = deque([(start, 0)])

    while q:

        vert, dist = q.popleft()
        seen.add(vert)

        if vert == target: 
            return dist
        
        for adj in graph[vert]:

            if adj in seen: 
                continue

            q.append((adj, dist + 1))

    return -1

print(bfs(1, 5), bfs(3, 5))

# If you use BFS for search on matrices, you need modification:


walls = [[0, 1, 0, 0, 0], [0, 1, 0, 1, 0], [0, 0, 0, 1, 0], [0, 0, 0, 1, 0]]
m, n = len(walls), len(walls[0])
dirs = [(1, 0), (0, 1), (-1, 0), (0, -1)]
start = (0, 0)
def bfs_matrix(start):
    '''
    Find the length of a shortest path to the right-bottom cell (m-1, n-1) in a matrix.

    Time complexity: :math:`O(mn)`

    Space complexity: :math:`O(mn)`
    where m is the number of rows and n is the number of columns in the matrix

    :param start: tuple
    :return: int
    '''
    seen = set([start])
    q = deque([(start[0], start[1], 0)])

    while q:

        r, c, dist = q.popleft()
        
        if (r, c) == (m - 1, n - 1):
            return dist
        
        for dr, dc in dirs:

            if r + dr >= m or r + dr < 0 or c + dc >= n or c + dc < 0: 
                continue

            if (r + dr, c + dc) not in seen and not walls[r + dr][c + dc]:
                q.append((r+dr, c+dc, dist+1))
                seen.add((r+dr, c+dc)) # to escape repetitions in queue in matrix case

    return -1

print(bfs_matrix(start))
