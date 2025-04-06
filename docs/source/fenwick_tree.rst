Fenwick Tree
================

Fenwick Tree (Binary Indexed Tree)

Fenwick Tree is a data structure that can efficiently update elements and calculate prefix sums in a table of numbers.

Given an array of n numbers, we define prefix sum up to index i as the sum of the first i elements in the array.

Fenwick Tree allows us to do both operations (update and prefix sum) in O(log n) time.

The idea is to represent the numbers as a tree, where the value of each node is the sum of the numbers in that subtree.

The tree structure allows us to efficiently update elements and calculate prefix sums.

The tree is constructed in such a way that each node is responsible for a range of elements in the array.

The range of elements that a node is responsible for is determined by the least significant bit of the node's index.

For example, if the index of a node is 6, then the range of elements that the node is responsible for is [5, 6].

.. automodule:: fenwick.fenwick
    :members:
    :undoc-members:
    :show-inheritance:
    :inherited-members:

Applications
------------

Applications of Fenwick Tree

1. Keeping track of positions problem

Given an array of n numbers, we need to efficiently answer q queries of the form:
    What is the position of the first element that is greater than or equal to x?

    What is the position of the last element that is less than or equal to x?

    What is the position of the first element that is greater than x?

    What is the position of the last element that is less than x?

    What is the sum of the first k elements?

    What is the sum of the elements from l to r?

    What is the smallest element in the first k elements?

    What is the greatest element in the first k elements?

    What is the smallest element in the elements from l to r?

    What is the greatest element in the elements from l to r?

    What is the number of elements in the first k elements that are greater than or equal to x?

    What is the number of elements in the elements from l to r that are greater than or equal to x?
