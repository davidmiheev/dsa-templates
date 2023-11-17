# Interval tree is a tree data structure to hold intervals. Specifically, it allows one to efficiently find all intervals that overlap with any given interval or point.
#


# Path: interval_tree.py
#

class IntervalTree:

    def __init__(self, intervals):
        self.root = None
        self.intervals = intervals
        self.__build()

    def __build(self):
        self.root = self.__build_helper(self.intervals)

    def __build_helper(self, intervals):
        if not intervals:
            return None

        mid = len(intervals) // 2
        root = TreeNode(intervals[mid])
        root.left = self.__build_helper(intervals[:mid])
        root.right = self.__build_helper(intervals[mid+1:])
        return root

    def query(self, interval):
        return self.__query_helper(self.root, interval)

    def __query_helper(self, root, interval):
        if not root:
            return []

        if root.interval.overlap(interval):
            return [root.interval] + self.__query_helper(root.left, interval) + self.__query_helper(root.right, interval)

        if root.left and root.left.max >= interval.start:
            return self.__query_helper(root.left, interval)

        return self.__query_helper(root.right, interval)

    def add(self, interval):
        self.intervals.append(interval)
        self.__build()

    def remove(self, interval):
        self.intervals.remove(interval)
        self.__build()

    def __repr__(self):
        return str(self.intervals)


class TreeNode:

    def __init__(self, interval):
        self.interval = interval
        self.max = interval.end
        self.left = None
        self.right = None

    def __repr__(self):
        return str(self.interval)


class Interval:

    def __init__(self, start, end):
        self.start = start
        self.end = end

    def overlap(self, other):
        return (self.start <= other.end and self.end >= other.start)

    def __repr__(self):
        return f'[{self.start}, {self.end}]'


intervals = [Interval(15, 20), Interval(10, 30), Interval(17, 19), Interval(5, 20), Interval(12, 15), Interval(30, 40)]
interval_tree = IntervalTree(intervals)
print(interval_tree)
print(interval_tree.query(Interval(14, 16)))
print(interval_tree.query(Interval(21, 23)))
print(interval_tree.query(Interval(29, 31)))
print(interval_tree.query(Interval(39, 41)))
print(interval_tree.query(Interval(9, 11)))
print(interval_tree.query(Interval(19, 21)))
