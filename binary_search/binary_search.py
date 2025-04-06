print('bisect')

def bisect_1(nums, target):
    '''
    Given a sorted (in ascending order) integer array `nums` of n elements and a target value,
    find the target value in the array.

    Time complexity: :math:`O(\log n)`

    Space complexity: :math:`O(1)`


    :type nums: List[int]
    :type target: int
    :rtype: int
    '''
    left, right = 0, len(nums) - 1
    while left <= right:
        mid = (left + right) // 2
        if nums[mid] == target:
            return mid
        elif nums[mid] < target:
            left = mid + 1
        else:
            right = mid - 1

    # End Condition: left > right
    return -1


def bisect_2(nums, target):
    '''
    Given a sorted (in ascending order) integer array nums of n elements and a target value,
    find the target value in the array.
    The array may contain duplicates.

    Time complexity: :math:`O(\log n)`

    Space complexity: :math:`O(1)`

    :type nums: List[int]
    :type target: int
    '''
    left, right = 0, len(nums) - 1
    while left < right:
        mid = (left + right) // 2
        if nums[mid] == target:
            return mid
        elif nums[mid] < target:
            left = mid + 1
        else:
            right = mid

    # End Condition: left == right
    # Post-processing:
    if nums[left] == target:
        return left

    return -1


def bisect_3(nums, target):
    '''
    Given a sorted (in ascending order) integer array nums of n elements and a target value,
    find the target value in the array.
    The array may contain duplicates.

    Time complexity: :math:`O(\log n)`

    Space complexity: :math:`O(1)`

    :type nums: List[int]
    :type target: int
    '''
    left, right = 0, len(nums) - 1
    while left + 1 < right:
        mid = (left + right) // 2
        if nums[mid] == target:
            return mid
        elif nums[mid] < target:
            left = mid
        else:
            right = mid

    # End Condition: left + 1 == right
    # Post-processing:
    if nums[left] == target: 
        return left
    if nums[right] == target:
        return right

    return -1


def minimizeMax(nums: list[int], p: int) -> int:
    '''
    Example problem: https://leetcode.com/problems/minimize-the-maximum-difference-of-pairs/
    You are given a 0-indexed integer array nums and an integer `p`.
    Find p pairs of indices of nums such that the maximum difference amongst all the pairs is minimized.
    Also, ensure no index appears more than once amongst the `p` pairs.
    Note that for a pair of elements at the index `i` and `j`, the difference of this pair is `|nums[i] - nums[j]|`,
    where `|x|` represents the absolute value of `x`.
    Return the minimum maximum difference among all p pairs. We define the maximum of an empty set to be zero.

    Time complexity: :math:`O(n \log n)`

    Space complexity: :math:`O(n)`
    '''
    if not p:
        return 0
    nums.sort()
    diffs = [a - b for a, b in zip(nums[1:], nums[:-1])]
    left, right = min(diffs), max(diffs)

    while left <= right:

        mid = (left + right) // 2
        cnt, j = 0, 0

        while j < len(diffs):

            if diffs[j] <= mid:
                cnt += 1
                j += 1
            
            j += 1
        
        if cnt >= p:
            right = mid - 1
        else: 
            left = mid + 1

    return left

print(minimizeMax([10, 1, 2, 7, 1, 3], p = 2))
