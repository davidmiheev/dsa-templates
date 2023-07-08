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
# The following picture shows the idea of Fenwick Tree:
    # https://www.geeksforgeeks.org/wp-content/uploads/BITSum.png
    # https://www.geeksforgeeks.org/wp-content/uploads/BITUpdate.png
    # https://www.geeksforgeeks.org/wp-content/uploads/BITUpdate2.png
    # https://www.geeksforgeeks.org/wp-content/uploads/BITSum2.png
    # https://www.geeksforgeeks.org/wp-content/uploads/BITSum3.png
    # https://www.geeksforgeeks.org/wp-content/uploads/BITSum4.png
    # https://www.geeksforgeeks.org/wp-content/uploads/BITSum5.png
    # https://www.geeksforgeeks.org/wp-content/uploads/BITSum6.png
    # https://www.geeksforgeeks.org/wp-content/uploads/BITSum7.png

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
