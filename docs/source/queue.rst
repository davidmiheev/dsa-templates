Queue and Deque
===============

Queue is data structure that follows the FIFO (First In First Out) principle.
The first element that is added to the queue is the first one to be removed.

The queue is a linear data structure that can be implemented using an array or a linked list.

The queue has two main operations:

- Enqueue: add an element to the end of the queue
- Dequeue: remove the element at the front of the queue

Dequeue is a double-ended queue that allows adding and removing elements from both ends.

Queue is an important data structure that is used in many algorithms and problems.
For example, the breadth-first search algorithm uses a queue to traverse a graph.

Here are collected problems that can be solved using the queue (deque) data structure alongwith implementation of queue.

Queue Implementation using Array
--------------------------------
.. autoclass:: _queue_._queue_.Queue
    :members:

Max Sliding Window
------------------
.. autofunction:: _queue_._queue_.maxSlidingWindow
