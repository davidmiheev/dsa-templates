# Fenwick Tree (Binary Indexed Tree)
#
# Fenwick Tree is a data structure that can efficiently update elements and calculate prefix sums in a table of numbers.
#
# Given an array of n numbers, we define prefix sum up to index i as the sum of the first i elements in the array.
#
# Fenwick Tree allows us to do both operations (update and prefix sum) in O(log n) time.
#
# The idea is to represent the numbers as a tree, where the value of each node is the sum of the numbers in that subtree.
#
# The tree structure allows us to efficiently update elements and calculate prefix sums.
#

from typing import List

class FenwickTree:
    '''
    Fenwick Tree (Binary Indexed Tree)

    Fenwick Tree is a data structure that can efficiently update elements
    and calculate prefix sums in a table of numbers.


    '''

    def __init__(self, nums: List[int]):
        '''
        Initialize the Fenwick Tree with the given array of numbers.

        :type nums: List[int]
        :rtype: None
        '''
        self.size = len(nums)
        self.tree = [0] * (self.size + 1)
        self.data = [0] * (self.size)

        for i in range(self.size):
            self.update(i, nums[i])

    def update(self, index: int, val: int) -> None:
        '''
        Update the value of the element at index i to be val.

        :type index: int
        :type val: int
        :rtype: None
        '''
        # Calculate the delta
        delta = val - self.data[index]
        self.data[index] = val
        # Increment the index to match the 1-based index of the Fenwick Tree
        index += 1
        while index <= self.size:
            self.tree[index] += delta
            index += index & -index


    def pref(self, index: int) -> int:
        '''
        Calculate the prefix sum up to index i.

        :type index: int
        :rtype: int
        '''
        index += 1
        result = 0
        while index > 0:
            result += self.tree[index]
            index -= index & -index

        return result

    def sumRange(self, left: int, right: int) -> int:
        '''
        Calculate the sum of elements from index l to r (with 0-based indexing).

        :type left: int
        :type right: int
        :rtype: int
        '''
        return self.pref(right) - self.pref(left-1)


nums = [1, 3, 5]
fenwick = FenwickTree(nums)
print(fenwick.sumRange(0, 2))
fenwick.update(1, 2)
print(fenwick.sumRange(0, 2))

class Fenwick:
    def __init__(self, n, nums=None):
        self.n = n
        self.arr = [0 for _ in range(n + 1)]
        if nums:
            self.arr[1:] = nums
            self.__build()

    def __build(self):
        for i in range(1, self.n):
            if i + (i & -i) <= self.n:
                self.arr[i + (i & -i)] += self.arr[i]

    def add(self, idx, dlt):
        idx += 1
        while idx <= self.n:
            self.arr[idx] += dlt
            idx += idx & -idx

    def query(self, ql, qr):
        # [ql, qr)
        return self.pref(qr) - self.pref(ql)

    def pref(self, qr):
        # [0, qr)
        ans = 0
        while qr:
            ans += self.arr[qr]
            qr -= qr & -qr
        return ans

    def suff(self, ql):
        # [ql, n)
        return self.pref(self.n) - self.pref(ql)



fenwick = Fenwick(10, [5,4,6,7,2,3,1,2,3,4])
print(fenwick.arr)
print(fenwick.query(0, 10))
print(fenwick.query(0, 5))
print(fenwick.query(5, 10))
print(fenwick.query(2, 8))
print(fenwick.add(3, 7))
print(fenwick.arr)
print(fenwick.query(0, 10))
print(fenwick.query(0, 5))
print(fenwick.query(5, 10))
print(fenwick.query(2, 8))

# Applications of Fenwick Tree
#
# 1. Keeping track of positions problem
#
# Given an array of n numbers, we need to efficiently answer q queries of the form:
    # What is the position of the first element that is greater than or equal to x?
    # What is the position of the last element that is less than or equal to x?
    # What is the position of the first element that is greater than x?
    # What is the position of the last element that is less than x?
    # What is the sum of the first k elements?
    # What is the sum of the elements from l to r?
    # What is the smallest element in the first k elements?
    # What is the greatest element in the first k elements?
    # What is the smallest element in the elements from l to r?
    # What is the greatest element in the elements from l to r?
    # What is the number of elements in the first k elements that are greater than or equal to x?
    # What is the number of elements in the elements from l to r that are greater than or equal to x?
