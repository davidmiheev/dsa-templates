from typing import List
import math
# Description: Backtracking algorithms
print('backtrack')

def permutations(nums: List[int]) -> List[List[int]]:
    '''
    Given a collection of distinct integers, return all possible permutations.

    Time complexity: :math:`O(n!)`

    Space complexity: :math:`O(n!)`
    '''
    def backtrack(path, nums, res):
        if not nums:
            res.append(path)
            return
        for i in range(len(nums)):
            backtrack(path + [nums[i]], nums[:i] + nums[i+1:], res)

    res = []
    backtrack([], nums, res)
    return res

def combinationSum(k, n):
    '''
    Problem: Find all possible combinations of k numbers that add up to a number n,
    given that only numbers from 1 to 9 can be used and each combination should be a unique set of numbers.

    Time complexity: :math:`O(9^k)`

    Space complexity: :math:`O(9^k)`
    '''

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


def nCastles(n):
    '''
    Problem: Given a chessboard of size n x n,
    find all possible positions to place n castles such that no two castles attack each other.

    Time complexity: :math:`O(n!)`

    Space complexity: :math:`O(n!)`
    '''
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


#

def nQueens(n: int) -> List[List[str]]:
    '''
    The n-queens puzzle is the problem of placing n queens on an n x n chessboard
    such that no two queens attack each other.
    Problem: Given an integer n, return all distinct solutions to the n-queens puzzle.
    You may return the answer in any order.

    Each solution contains a distinct board configuration of the n-queens' placement, where 'Q' and '.' both indicate a queen and an empty space, respectively.

    Time complexity: :math:`O(n!)`

    Space complexity: :math:`O(n!)`
    '''
    arrangement = []
    ans = []
    template = [['.']*n for _ in range(n)]
    def backtrack(j):
        nonlocal template
        if j in arrangement: return
        for r, c in enumerate(arrangement):
            if abs(j-c) == len(arrangement) - r:
                return

        arrangement.append(j)
        if len(arrangement) == n:
            for r, c in enumerate(arrangement):
                template[r][c] = 'Q'
                template[r] = ''.join(template[r])

            ans.append(template)
            template = [['.']*n for _ in range(n)]
            arrangement.pop()
            return

        for i in range(n):
            backtrack(i)

        arrangement.pop()

    for i in range(n):
        backtrack(i)

    return ans

n = 8
print(permutations([1, 2, 3]))
print(combinationSum(3, 12))
print(len(nCastles(n)), math.factorial(n))
print(nQueens(n))
