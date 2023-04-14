# Classic application: search element in sorted array with O(log n) time complexity
# In python there is built-in library for binary search: bisect
# Example with this library:

import bisect

arr = [1,3,4,5,6,7]
print(bisect.bisect_left(arr, 5))

# More advanced application is solving of optimisation problems.
# Optimisation problems are problems which ask to find the maximum or minimum for some value
# This technique is called "binary search on answer"
# Example:

# You are given a 0-indexed integer array nums and an integer p. Find p pairs of indices of nums such that the maximum difference amongst all the pairs is minimized.
# Also, ensure no index appears more than once amongst the p pairs.
# Note that for a pair of elements at the index i and j, the difference of this pair is |nums[i] - nums[j]|, where |x| represents the absolute value of x.
# Return the minimum maximum difference among all p pairs. We define the maximum of an empty set to be zero.

# Solution:

class Solution:
    def minimizeMax(self, nums: list[int], p: int) -> int:
        if not p: return 0
        nums.sort()
        diffs = [a-b for a,b in zip(nums[1:], nums[:-1])]
        left, right = min(diffs), max(diffs)
        while left <= right:
            mid = (left+right)//2
            cnt, j = 0, 0
            while j < len(diffs):
                if diffs[j] <= mid:
                    cnt += 1
                    j += 1
                j += 1
            if cnt >= p: right = mid - 1
            else: left = mid + 1

        return left

sol = Solution()
print(sol.minimizeMax([10,1,2,7,1,3], p = 2))

# Template
def optimize(nums, p):
    left, right = min(nums), max(nums)
    def feasible(val):
        # ...
        return val >= p

    while left <= right:
        mid = (left+right)//2
        if feasible(mid): right = mid - 1
        else: left = mid + 1

    return left
