from collections import defaultdict

print('Union-Find')

class UnionFind:
    '''
    Union-Find data structure

    Union-Find data structure is a disjoint-set data structure
    that keeps track of a set of elements partitioned into a number of disjoint (non-overlapping) subsets.

    Variant with rank attribute (for optimization)

    '''
    def __init__(self, vertices: set):
        '''
        Initialize the Union-Find data structure with the given set of vertices.

        Time complexity: :math:`O(n)`, where n is the number of vertices.

        :type vertices: set
        :rtype: None
        '''
        self.group = defaultdict(int)
        self.rank = defaultdict(int)
        for i in vertices:
            self.group[i] = i

    def find(self, node: int) -> int:
        '''
        Find the group that the given node belongs to.


        :type node: int
        :rtype: int
        '''
        if self.group[node] != node:
            self.group[node] = self.find(self.group[node])
        return self.group[node]

    def join(self, node1: int, node2: int):
        '''
        Join the groups that the given nodes belong to.

        :type node1: int
        :type node2: int
        :rtype: None
        '''
        group1 = self.find(node1)
        group2 = self.find(node2)

        # node1 and node2 already belong to same group.
        if group1 == group2: return

        if self.rank[group1] > self.rank[group2]:
            self.group[group2] = group1
        elif self.rank[group1] < self.rank[group2]:
            self.group[group1] = group2
        else:
            self.group[group1] = group2
            self.rank[group2] += 1

    def are_connected(self, node1: int, node2: int) -> bool:
        '''
        Check if the given nodes belong to the same group.

        :type node1: int
        :type node2: int
        :rtype: bool
        '''
        return self.find(node1) == self.find(node2)


# :

# UnionFind class
class UnionFindSimple:
    '''
    Union-Find

    Variant without rank attribute (for simplicity)
    '''
    def __init__(self, size):
        '''

        '''
        self.root = [i for i in range(size)]

    def find(self, x):

        if x == self.root[x]: return x
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
vertices = {1, 2, 3, 4, 5, 6, 7}
uf = UnionFind(vertices)
graph = [(1, 2), (2, 3), (1, 3), (1, 4), (2, 5), (6, 7)]
connected_components = len(vertices)
cycles = 0
for x, y in graph:
    if uf.are_connected(x, y):
        cycles += 1
    else:
        connected_components -= 1
        uf.join(x, y)

print(f'N of connected components: {connected_components}, N of cycles: {cycles}')
