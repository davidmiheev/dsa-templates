from functools import cache

print('Dynamic Programming')

def kadane_recursive(arr):
    '''
    Given an integer array `nums`, find the contiguous subarray (containing at least one number)
    which has the largest sum and return its sum.

    It is recursive version of Kadane's Algorithm.
    Example of Recursive DP (I prefer recursive DP since it's much more comfortable to me)

    Time complexity: :math:`O(n)`

    Space complexity: :math:`O(n)`


    :type arr: List[int]
    :rtype: int
    '''

    @cache
    def dp(i):
        if i == 0: 
            return arr[0], arr[0]
        max_ending_here, max_so_far = dp(i-1)
        max_ending_here = max(arr[i], max_ending_here + arr[i])
        max_so_far = max(max_so_far, max_ending_here)
        return max_ending_here, max_so_far

    _, ans = dp(len(arr)-1)
    return ans

print("max_subarray_sum (recursive_dp):", kadane_recursive([1, 2, 3, -7, 5])) # 6

def max_subarray_sum(arr):
    '''
    Given an integer array nums, find the contiguous subarray (containing at least one number)
    which has the largest sum and return its sum.

    Classic Kadane's Algorithm

    Time complexity: :math:`O(n)`

    Space complexity: :math:`O(1)`

    :type arr: List[int]
    :rtype: int
    '''
    max_so_far = arr[0]
    max_ending_here = arr[0]
    for i in range(1, len(arr)):
        max_ending_here = max(arr[i], max_ending_here + arr[i])
        max_so_far = max(max_so_far, max_ending_here)

    return max_so_far

print("max_subarray_sum:", max_subarray_sum([1, 2, 3, -7, 5])) # 6


def max_subarray_product(arr):
    '''
    Given an integer array nums, find the contiguous subarray (containing at least one number)
    which has the largest product and return its product.

    Time complexity: :math:`O(n)`

    Space complexity: :math:`O(1)`

    :type arr: List[int]
    :rtype: int
    '''
    max_so_far = arr[0]
    max_ending_here = arr[0]
    min_ending_here = arr[0]
    for i in range(1, len(arr)):
        temp = max_ending_here
        max_ending_here = max(arr[i], max_ending_here * arr[i], min_ending_here * arr[i])
        min_ending_here = min(arr[i], temp * arr[i], min_ending_here * arr[i])
        max_so_far = max(max_so_far, max_ending_here)

    return max_so_far

print("maximum_subarray_product:", max_subarray_product([2, 3, -2, -3, -5]))


def max_subarray_sum_circular(arr):
    '''
    Given a circular integer array nums of length n,
    return the maximum possible sum of a non-empty subarray of nums.
    A circular array means the end of the array connects to the beginning of the array.

    Time complexity: :math:`O(n)`

    Space complexity: :math:`O(1)`

    :type arr: List[int]
    :rtype: int
    '''
    max_so_far = arr[0]
    max_ending_here = arr[0]
    min_ending_here = arr[0]
    min_so_far = arr[0]
    total = arr[0]
    for i in range(1, len(arr)):
        total += arr[i]
        temp = max_ending_here
        max_ending_here = max(arr[i], max_ending_here + arr[i])
        min_ending_here = min(arr[i], min_ending_here + arr[i])
        max_so_far = max(max_so_far, max_ending_here)
        min_so_far = min(min_so_far, min_ending_here)

    if max_so_far < 0:
        return max_so_far
    else:
        return max(max_so_far, total - min_so_far)

print("max_subarray_sum_circular:", max_subarray_sum_circular([1, -2, 3, -2])) # 3




# DP with bitmask
# TODO


#####################
# Range Sum Queries #
#####################
def range_sum_queries(arr):
    '''
    Given an array of n numbers, we need to efficiently answer q queries of the form:
        What is the sum of elements from index l to r (with 0-based indexing)?

    Time complexity: :math:`O(n)`

    Space complexity: :math:`O(n)`

    :type arr: List[int]
    :rtype: Callable[[int, int], int]
    '''
    n = len(arr)
    dp = [0] * (n + 1)
    for i in range(1, n + 1):
        dp[i] = dp[i - 1] + arr[i - 1]

    def query(l, r):
        return dp[r + 1] - dp[l]

    return query

print("range_sum_queries:", range_sum_queries([1, 2, 3, 4, 5])(1, 3)) # 9

# What is the minimum/maximum element from index l to r (with 0-based indexing)?
def range_min_queries(arr):
    '''
    Given an array of n numbers, we need to efficiently answer q queries of the form:
        What is the minimum/maximum element from index l to r (with 0-based indexing)?

    Time complexity: :math:`O(n^2)`

    Space complexity: :math:`O(n^2)`

    :type arr: List[int]
    :rtype: Tuple[Callable[[List[int]], List[int]], Callable[[List[int], int, int], int]]
    '''
    def build(arr):
        n = len(arr)
        dp = [[0] * n for _ in range(n)]
        for i in range(n):
            dp[i][i] = arr[i]
        for i in range(n):
            for j in range(i + 1, n):
                dp[i][j] = min(dp[i][j - 1], arr[j]) # change to max for max queries
        return dp


    return build, lambda dp, l, r: dp[l][r]

build = range_min_queries([1, 2, 3, 4, 5])[0]
print("range_min_queries:", range_min_queries([1, 2, 3, 4, 5])[1](build([1, 2, 3, 4, 5]), 1, 3)) # 2

# What is the greatest common divisor of elements from index l to r (with 0-based indexing)?
# What is the least common multiple of elements from index l to r (with 0-based indexing)?
from math import gcd, lcm
def range_gcd_queries(arr):
    '''
    Given an array of n numbers, we need to efficiently answer q queries of the form:
        What is the greatest common divisor of elements from index l to r (with 0-based indexing)?

    Time complexity: :math:`O(n^2)`

    Space complexity: :math:`O(n^2)`

    :type arr: List[int]
    :rtype: Tuple[Callable[[List[int]], List[int]], Callable[[List[int], int, int], int]]
    '''

    def build(arr):
        n = len(arr)
        dp = [[0] * n for _ in range(n)]
        for i in range(n):
            dp[i][i] = arr[i]
        for i in range(n):
            for j in range(i + 1, n):
                dp[i][j] = gcd(dp[i][j - 1], arr[j]) # change to lcm for lcm queries
        return dp

    def query(dp, l, r):
        return dp[l][r]

    return build, query

dp = range_gcd_queries([2, 4, 6, 8, 10])[0]
print("range_gcd_queries:", range_gcd_queries([2, 4, 6, 8, 10])[1](dp([2, 4, 6, 8, 10]), 1, 3)) # 2
