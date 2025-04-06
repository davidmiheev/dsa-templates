from typing import List
from collections import defaultdict

print('Union-Find')

class UnionFind:
    '''
    Union-Find data structure

    Variant with size of component attribute (for optimization)
    and number of components attribute
    '''
    def __init__(self, n: int, vertices = None):
        '''
        Initialize the Union-Find data structure with the given number of vertices.

        Time complexity: :math:`O(n)`, where n is the number of vertices.

        Space complexity: :math:`O(n)`, where n is the number of vertices.
        '''
        if vertices:
            self.root = {v: v for v in vertices}
            self.size = {v: 1 for v in vertices}
            self.components = len(vertices)
        else:
            self.root = list(range(n + 1))
            self.size = [1] * (n + 1)
            self.components = n

    def find(self, x: int) -> int:
        '''
        Find the parent of the given node.

        Time complexity: :math:`O(a(n))`

        Space complexity: :math:`O(1)`
        '''
        if self.root[x] == x:
            return x
        self.root[x] = self.find(self.root[x])
        return self.root[x]

    def union(self, x: int, y: int) -> int:
        '''
        Union the components that the given nodes belong to.

        Time complexity: :math:`O(a(n))`, where a(n) is the inverse Ackermann function (very slow growing function).

        Space complexity: :math:`O(1)`
        '''
        x = self.find(x)
        y = self.find(y)
        if x == y:
            return 0

        if self.size[x] > self.size[y]:
            self.size[x] += self.size[y]
            self.root[y] = x
        else:
            self.size[y] += self.size[x]
            self.root[x] = y

        self.components -= 1
        return 1

    def are_connected(self, x: int, y: int) -> bool:
        '''
        Check if the given nodes are connected

        Time complexity: :math:`O(a(n))`, where a(n) is the inverse Ackermann function (very slow growing function).

        Space complexity: :math:`O(1)`
        '''
        return self.find(x) == self.find(y)



class UnionFindSimple:
    '''
    Union-Find

    Variant without rank attribute (for simplicity)
    '''
    def __init__(self, size):
        self.root = [i for i in range(size)]

    def find(self, x):
        if x == self.root[x]: 
            return x
        self.root[x] = self.find(self.root[x])
        return self.root[x]

    def union(self, x, y):
        rootX = self.find(x)
        rootY = self.find(y)
        if rootX != rootY:
            self.root[rootY] = rootX

    def are_connected(self, x, y):
        return self.find(x) == self.find(y)


# Example usage:
vertices = [1, 2, 3, 4, 5, 6, 7]
uf = UnionFind(7, vertices)
graph = [(1, 2), (2, 3), (1, 3), (1, 4), (2, 5), (6, 7)]
connected_components = len(vertices)
cycles = 0
for x, y in graph:
    if uf.are_connected(x, y):
        cycles += 1
    else:
        connected_components -= 1
        uf.union(x, y)

print(f'N of connected components: {connected_components}, N of cycles: {cycles}')
