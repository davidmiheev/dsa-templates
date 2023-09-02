# Depth-First Search (DFS) is an algorithm for traversing or searching tree or graph data structures.
# One starts at the root (selecting some arbitrary node as the root in the case of a graph)
# and explores as far as possible along each branch before backtracking.
#
# DFS is typically implemented recursively, but can be implemented iteratively as well.
#
# Time complexity: O(V + E)
#
# Space complexity: O(V)
#
# where V is the number of vertices and E is the number of edges in the graph.
#
# Applications:
# Topological sorting
# Finding connected components
# Finding bridges and articulation points
# Finding strongly connected components
# Solving puzzles with only one solution, such as mazes
#
# Example problem: https://leetcode.com/problems/number-of-enclaves/
#
# You are given an m x n binary matrix grid, where 0 represents a sea cell and 1 represents a land cell.

# A move consists of walking from one land cell to another adjacent (4-directionally) land cell or walking off the boundary of the grid.

# Return the number of land cells in grid for which we cannot walk off the boundary of the grid in any number of moves.
#
from typing import List

class Solution:
    def numEnclaves(self, grid: List[List[int]]) -> int:
        m, n = len(grid), len(grid[0])
        dirs = [(1,0), (0,1), (-1,0), (0,-1)]
        visited = set()
        def dfs(cell):
            if cell in visited: return False
            i, j = cell
            if not grid[i][j]: return False
            visited.add(cell)
            if i == m-1 or j == n-1 or not i or not j:
                return False
            ret = True
            for di, dj in dirs:
                if (i+di, j+dj) in visited: continue
                if not grid[i+di][j+dj]: continue
                ret &= dfs((i+di, j+dj))

            return ret

        ret = 0
        for i in range(1, m-1):
            for j in range(1, n-1):
                temp = len(visited)
                if dfs((i, j)): ret += len(visited) - temp

        return ret

sol = Solution()
grid = [[0,0,0,0],[1,0,1,0],[0,1,1,0],[0,0,0,0]]
print(sol.numEnclaves(grid))

# Example problem: https://leetcode.com/problems/number-of-provinces/
#
# There are n cities. Some of them are connected, while some are not. If city a is connected directly with city b,
# and city b is connected directly with city c, then city a is connected indirectly with city c.
#
# A province is a group of directly or indirectly connected cities and no other cities outside of the group.
#
# You are given an n x n matrix isConnected where isConnected[i][j] = 1 if the ith city and the jth city are directly connected,
# and isConnected[i][j] = 0 otherwise.
#
# Return the total number of provinces.
#

from functools import cache

class Solution:
    def findCircleNum(self, isConnected: List[List[int]]) -> int:
        n = len(isConnected)
        seen = set()
        @cache
        def dfs(i):
            seen.add(i)
            for j in range(n):
                if not isConnected[i][j]: continue
                if j in seen: continue
                dfs(j)

        n_provinces = 0
        for i in range(n):
            m = len(seen)
            dfs(i)
            if len(seen) > m:
                n_provinces += 1

        return n_provinces

sol = Solution()
isConnected = [[1,1,0],[1,1,0],[0,0,1]]
print(sol.findCircleNum(isConnected))
