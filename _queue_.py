from typing import List
from collections import deque

print('Queue')

class Queue:

    def __init__(self):
        self.front = []
        self.back = []


    def push(self, x: int) -> None:
        self.back += [x]


    def pop(self) -> int:
        if not self.front:
            self.front, self.back = self.back, []
            self.front.reverse()

        return self.front.pop()


    def peek(self) -> int:
        if not self.front:
            self.front, self.back = self.back, []
            self.front.reverse()

        return self.front[-1]


    def empty(self) -> bool:
        if not self.back and not self.front:
            return True

        return False


def maxSlidingWindow(nums: List[int], k: int) -> List[int]:
    '''
    Given an array of integers and a window size k, return the maximum value in each window of size k.
    '''
    n = len(nums)
    deq = deque()

    ans = []

    for i in range(n):
        while deq and nums[i] > nums[deq[-1]]:
            deq.pop()

        if deq and i-k+1 > deq[0]:
            deq.popleft()

        deq.append(i)

        if i >= k-1:
            ans.append(nums[deq[0]])

    return ans


nums = [1,3,-1,-3,5,3,6,7]
k = 3
print(maxSlidingWindow(nums, k)) # [3,3,5,5,6,7]

q = Queue()
q.push(1)
q.push(2)
print(q.peek()) # 1
print(q.pop()) # 1
print(q.empty()) # False
