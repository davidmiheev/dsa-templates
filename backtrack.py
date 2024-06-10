# Here is a template for backtracking approach
#
#
"""
def backtrack(curr, OTHER_ARGUMENTS...):
    if (BASE_CASE):
        # modify the answer
        return

    ans = 0
    for (ITERATE_OVER_INPUT):
        # modify the current state
        ans += backtrack(curr, OTHER_ARGUMENTS...)
        # undo the modification of the current state

    return ans
"""

# Example of backtracking approach
# Problem: Find all possible combinations of k numbers that add up to a number n, given that only numbers from 1 to 9 can be used and each combination should be a unique set of numbers.
#

def combinationSum3(k, n):

    def backtrack(start, k, n, path, res):
        if k == 0 and n == 0:
            res.append(path)
            return
        if k == 0 or n == 0:
            return
        for i in range(start, 10):
            backtrack(i + 1, k - 1, n - i, path + [i], res)

    res = []
    backtrack(1, k, n, [], res)
    return res

# Example: Find all possible castles positions in a chessboard
#
# Problem: Given a chessboard of size n x n, find all possible positions to place n castles such that no two castles attack each other.
#
def nCastles(n):

    def backtrack(row, n, path, res):
        if row == n:
            res.append(path)
            return
        for col in range(n):
            if is_safe(path, row, col):
                backtrack(row + 1, n, path + [(row, col)], res)

    def is_safe(path, row, col):
        for r, c in path:
            if c == col or abs(row - r) == abs(col - c):
                return False
        return True

    res = []
    backtrack(0, n, [], res)
    return res

# The n-queens puzzle is the problem of placing n queens on an n x n chessboard
# such that no two queens attack each other.
# Given an integer n, return the number of distinct solutions to the n-queens puzzle.
#

def totalNQueens(n: int) -> int:
    arrangement = []
    ans = 0
    def backtrack(j):
        nonlocal ans
        if j in arrangement: return
        for r, c in enumerate(arrangement):
            if abs(j-c) == len(arrangement) - r:
                return

        arrangement.append(j)
        if len(arrangement) == n:
            ans += 1
            arrangement.pop()
            return

        for i in range(n):
            backtrack(i)

        arrangement.pop()

    for i in range(n):
        backtrack(i)

    return ans

n = 8
print(totalNQueens(n))
