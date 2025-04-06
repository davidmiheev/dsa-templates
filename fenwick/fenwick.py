from typing import List

print('fenwick')

class FenwickTree:
    '''
    Fenwick Tree (Binary Indexed Tree)

    Fenwick Tree is a data structure that can efficiently update elements
    and calculate prefix sums in a table of numbers.


    '''

    def __init__(self, nums: List[int]):
        '''
        Initialize the Fenwick Tree with the given array of numbers.

        Time complexity: :math:`O(n)`, where n is the number of elements in the array.

        Space complexity: :math:`O(n)`, where n is the number of elements in the array.

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

        Time complexity: :math:`O(\log n)`, where n is the number of elements in the array.

        Space complexity: :math:`O(1)`

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

        Time complexity: :math:`O(\log n)`, where n is the number of elements in the array.

        Space complexity: :math:`O(1)`

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

        Time complexity: :math:`O(\log n)`, where n is the number of elements in the array.

        Space complexity: :math:`O(1)`

        :type left: int
        :type right: int
        :rtype: int
        '''
        return self.pref(right) - self.pref(left-1)


nums = [1, 3, 5, 7, 9, 11]
fenwick = FenwickTree(nums)
print(fenwick.sumRange(0, 2))
fenwick.update(1, 2)
print(fenwick.sumRange(0, 2))
print(fenwick.tree)
print(fenwick.data)


# class Fenwick:
#     def __init__(self, n, nums=None):
#         self.n = n
#         self.arr = [0 for _ in range(n + 1)]
#         if nums:
#             self.arr[1:] = nums
#             self.__build()

#     def __build(self):
#         for i in range(1, self.n):
#             if i + (i & -i) <= self.n:
#                 self.arr[i + (i & -i)] += self.arr[i]

#     def add(self, idx, dlt):
#         idx += 1
#         while idx <= self.n:
#             self.arr[idx] += dlt
#             idx += idx & -idx

#     def query(self, ql, qr):
#         # [ql, qr)
#         return self.pref(qr) - self.pref(ql)

#     def pref(self, qr):
#         # [0, qr)
#         ans = 0
#         while qr:
#             ans += self.arr[qr]
#             qr -= qr & -qr
#         return ans

#     def suff(self, ql):
#         # [ql, n)
#         return self.pref(self.n) - self.pref(ql)
