Topological Sort
========================
Topological sort is a special ordering of vertices in a directed acyclic graph (DAG)
where for each directed edge from vertex A to vertex B, vertex A comes before vertex B in the ordering.

Topological sort is useful when you need to perform a series of tasks
where each task has some dependencies on other tasks.

For example, if you have a list of tasks to perform,
and some of those tasks depend on other tasks being completed first,
you could represent this as a graph with a directed edge from task A to task B
if task A must be completed before task B.

Then, a topological sort of the graph would give you an ordering of the tasks
that respects the dependencies.

Topological sort can be used to detect cycles in a graph.
If the graph has a cycle, then there is no topological sort.

Topological sort can be used to find the shortest path in a weighted graph (if there are no cycles).

Topological sort can be used to find the strongly connected components of a graph.

Topological sort can be used to order the compilation tasks for a program

Topological sort can be used to schedule tasks in a multi-core processor

.. automodule:: top_sort.top_sort
    :members:
    :undoc-members:
    :show-inheritance:
