# Kadane's Algorithm
# Maximum Subarray Sum
# Maximum Subarray Product
# Maximum Subarray Sum Circular

def max_subarray_sum(arr):
    max_so_far = arr[0]
    max_ending_here = arr[0]
    for i in range(1, len(arr)):
        max_ending_here = max(arr[i], max_ending_here + arr[i])
        max_so_far = max(max_so_far, max_ending_here)

    return max_so_far

def max_subarray_product(arr):
    max_so_far = arr[0]
    max_ending_here = arr[0]
    min_ending_here = arr[0]
    for i in range(1, len(arr)):
        temp = max_ending_here
        max_ending_here = max(arr[i], max_ending_here * arr[i], min_ending_here * arr[i])
        min_ending_here = min(arr[i], temp * arr[i], min_ending_here * arr[i])
        max_so_far = max(max_so_far, max_ending_here)

    return max_so_far

def max_subarray_sum_circular(arr):
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

# Range Sum Queries
#
# Given an array of n numbers, we need to efficiently answer q queries of the form:
#
# What is the sum of elements from index l to r (with 0-based indexing)?
def range_sum_queries(arr):
    n = len(arr)
    dp = [0] * (n + 1)
    for i in range(1, n + 1):
        dp[i] = dp[i - 1] + arr[i - 1]

    def query(l, r):
        return dp[r + 1] - dp[l]

    return query

# What is the minimum/maximum element from index l to r (with 0-based indexing)?
def range_min_queries(arr):

    def build(arr):
        n = len(arr)
        dp = [[0] * n for _ in range(n)]
        for i in range(n):
            dp[i][i] = arr[i]
        for i in range(n):
            for j in range(i + 1, n):
                dp[i][j] = min(dp[i][j - 1], arr[j]) # change to max for max queries
        return dp

    def query(dp, l, r):
        return dp[l][r]

    return build, query

# What is the greatest common divisor of elements from index l to r (with 0-based indexing)?
# What is the least common multiple of elements from index l to r (with 0-based indexing)?
from math import gcd, lcm
def range_gcd_queries(arr):

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

# What is the bitwise AND/OR/XOR of elements from index l to r (with 0-based indexing)?
def range_and_queries(arr):

    def build(arr):
        n = len(arr)
        dp = [[0] * n for _ in range(n)]
        for i in range(n):
            dp[i][i] = arr[i]
        for i in range(n):
            for j in range(i + 1, n):
                dp[i][j] = dp[i][j - 1] & arr[j] # change to | for OR queries, ^ for XOR queries
        return dp

    def query(dp, l, r):
        return dp[l][r]

    return build, query

# What is the number of elements from index l to r (with 0-based indexing) that are greater than or equal to x?

def range_count_queries(arr, x):

    def build(arr):
        n = len(arr)
        dp = [[0] * n for _ in range(n)]
        for i in range(n):
            dp[i][i] = 1 if arr[i] >= x else 0
        for i in range(n):
            for j in range(i + 1, n):
                dp[i][j] = dp[i][j - 1] + (1 if arr[j] >= x else 0)
        return dp

    def query(dp, l, r):
        return dp[l][r]

    return build, query

# Recursive DP

from functools import cache

def top_down_memoization(arr):

    @cache
    def dp(i):
        if i < 0: return 0
        return max(dp(i - 1), dp(i - 1) + arr[i])

    return dp(len(arr)-1)

# Bottom-up with tabulation

def bottom_up_tabulation(arr):

    n = len(arr)
    dp = [0] * (n + 1)
    dp[1] = arr[0]
    for i in range(2, n + 1):
        dp[i] = max(dp[i - 1], dp[i - 2] + arr[i - 1])
    return dp[n]

# Bottom-up with tabulation + space optimization

def bottom_up_tabulation_space_optimization(arr):

    n = len(arr)
    dp = [0] * 3
    dp[1] = arr[0]
    for i in range(2, n + 1):
        dp[i % 3] = max(dp[(i - 1) % 3], dp[(i - 2) % 3] + arr[i - 1])
    return dp[n % 3]

# State machine

# 1. State machine with transition table

from math import inf

def state_machine_with_transition_table(arr):

    def build(arr):
        n = len(arr)
        dp = [[0] * 2 for _ in range(n + 1)]
        dp[0][0] = 0
        dp[0][1] = -inf
        for i in range(1, n + 1):
            dp[i][0] = max(dp[i - 1][0], dp[i - 1][1] + arr[i - 1])
            dp[i][1] = max(dp[i - 1][1], dp[i - 1][0] - arr[i - 1])

        return dp[n][0]

    return build(arr)

# 2. State machine with transition function

def state_machine_with_transition_function(arr):

    def build(arr):
        n = len(arr)
        dp = [0] * 2
        dp[0] = 0
        dp[1] = -inf
        for i in range(1, n + 1):
            dp[0], dp[1] = max(dp[0], dp[1] + arr[i - 1]), max(dp[1], dp[0] - arr[i - 1])

        return dp[0]

    return build(arr)

# DP with bitmask
#
