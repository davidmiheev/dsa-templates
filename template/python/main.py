from math import sqrt, isqrt, gcd, inf, ceil, log2, comb, factorial, atan2, pi
from collections import deque, defaultdict, Counter
from itertools import permutations, combinations, product
from bisect import bisect_left, bisect_right
from heapq import heappush, heappop, heapify
from functools import cache, reduce, cmp_to_key
from random import getrandbits
import sys


mod = 998244353


def power(base, exp):
  result = 1
    
  base %= mod
  while exp > 0:
    if exp % 2 == 1:
      result = (result * base) % mod
            
    exp //= 2
    base = (base * base) % mod
  
       
  return result
      

def sieve(limit):
  s = list(range(limit + 1))
  s[0] = s[1] = -1
  
  for i in range(2, limit + 1):
    if s[i] == i:
      for j in range(i * i, limit + 1, i):
        if s[j] == j:
          s[j] = i
  
  return s


def prime_factors(n, s):
  factors = set()
  while n > 1:
    p = s[n]
    if p == -1:
      break
    factors.add(p)
    while n % p == 0:
      n //= p
  
  return factors 
  
  
class FenwickTree:
  def __init__(self, n):
    self.size = n
    self.tree = [0] * (self.size + 1)

  def update(self, index: int, delta: int) -> None:
    index += 1
    while index <= self.size:
      self.tree[index] += delta
      index += index & -index

  def pref(self, index: int) -> int:
    index += 1
    result = 0
    while index > 0:
      result += self.tree[index]
      index -= index & -index

    return result


  
class UnionFind:

  def __init__(self, n: int):
    self.root = list(range(n + 1))
    self.size = [1] * (n + 1)
    self.components = n
    self.max_size = 1

  def find(self, x: int) -> int:
    if self.root[x] == x:
        return x
    self.root[x] = self.find(self.root[x])
    return self.root[x]

  def union(self, x: int, y: int) -> int:
    x = self.find(x)
    y = self.find(y)
    if x == y:
        return 0

    if self.size[x] > self.size[y]:
        self.size[x] += self.size[y]
        self.max_size = max(self.max_size, self.size[x])
        self.root[y] = x
    else:
        self.size[y] += self.size[x]
        self.max_size = max(self.max_size, self.size[y])
        self.root[x] = y

    self.components -= 1
    return 1

  def are_connected(self, x: int, y: int) -> bool:
    return self.find(x) == self.find(y)


        
input = sys.stdin.readline


def solve():


    s = input().strip()
    n = int(input())

    print(s[n:-n])

   



if __name__ == '__main__':
  # run = input()
  t = 1 # int(input())
  for j in range(t):
    solve()
    
    
