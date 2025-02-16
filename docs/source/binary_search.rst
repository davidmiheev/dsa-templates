
Binary search
=======================

Binary search is a technique used to quickly locate a specific item in a sorted list.
The algorithm works by comparing the target value to the middle element of the list.
If the target value is equal to the middle element, the search is successful.
If the target value is less than the middle element, the search continues on the lower half of the list.
If the target value is greater than the middle element, the search continues on the upper half of the list.
This process continues until the target value is found or the search is unsuccessful.

In python there is built-in library for binary search: bisect

Example with this library:

.. code-block:: python

    import bisect

    arr = [1,3,4,5,6,7]
    print(bisect.bisect_left(arr, 5))

Basic Template
--------------

Basic Template of binary search without use of built-in library:

.. code-block:: python

    def fn(arr, target):
        left = 0
        right = len(arr) - 1
        while left <= right:
            mid = (left + right) // 2
            if arr[mid] == target:
                # do something
                return
            if arr[mid] > target:
                right = mid - 1
            else:
                left = mid + 1

        # left is the insertion point
        return left

Binare search: Version 1
------------------------

Most basic and elementary form of Binary Search.
Search Condition can be determined without comparing to the element's neighbors (or use specific elements around it).
No post-processing required because at each step, you are checking to see if the element has been found.
If you reach the end, then you know the element is not found

.. autofunction:: binary_search.bisect_1

Binary search: Version 2
------------------------
An advanced way to implement Binary Search.
Use the element's right neighbor to determine if the condition is met and decide whether to go left or right.
Guarantees Search Space is at least 2 in size at each step.
Post-processing required. Loop/Recursion ends when you have 1 element left. Need to assess if the remaining element meets the condition.

.. autofunction:: binary_search.bisect_2

Binary search: Version 3
------------------------

An alternative way to implement Binary Search.
Use element's neighbors to determine if condition is met and decide whether to go left or right
Gurantees Search Space is at least 3 in size at each step.
Post-processing required. Loop/Recursion ends when you have 2 elements left. Need to assess if the remaining elements meet the condition.

.. autofunction:: binary_search.bisect_3

Binary Search on Answer
-----------------------

More advanced application is solving of optimisation problems.
Optimisation problems are problems which ask to find the maximum or minimum for some value
This technique is called **"binary search on answer"**

Template for optimisation problems:

.. code-block:: python

   def optimize(nums, p):
       left, right = min(nums), max(nums)
       def feasible(val):
           # ...
           return val >= p

       while left <= right:
           mid = (left + right) // 2
           if feasible(mid):
               right = mid - 1
           else: left = mid + 1

       return left


.. autofunction:: binary_search.minimizeMax
