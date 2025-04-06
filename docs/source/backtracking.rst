Backtracking
=============
Backtracking is a general algorithm for finding all (or some) solutions to some computational problems,
notably constraint satisfaction problems, that incrementally builds candidates to the solutions,
and abandons a candidate ("backtracks") as soon as it determines that the candidate cannot possibly be completed to a valid solution.

Backtracking can be applied only for problems which admit the concept of a "partial candidate solution" and a relatively quick test of whether it can possibly be completed to a valid solution.

Backtracking is often much faster than brute force enumeration of all complete candidates, since it can eliminate a large number of candidates with a single test.

Backtracking is an important tool for solving constraint satisfaction problems, such as crosswords, verbal arithmetic, Sudoku, and many other puzzles.

Here is a template for backtracking approach:

.. code-block:: python

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


.. autofunction:: backtrack.backtrack.permutations

.. autofunction:: backtrack.backtrack.combinationSum

.. autofunction:: backtrack.backtrack.nCastles

.. autofunction:: backtrack.backtrack.nQueens
