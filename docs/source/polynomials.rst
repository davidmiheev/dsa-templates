Polynomials, Convolution and Linear Algebra
===========================================

NTT-based convolutions, power series, matrices and recurrences.

NTT Convolution
---------------

.. rust-api:: convolution
   :description: Convolutions by the number-theoretic transform.

Polynomials and Power Series
----------------------------

.. rust-api:: polynomial
   :description: Polynomials and formal power series over an NTT-friendly prime field.

Polynomials mod p
-----------------

.. rust-api:: poly_mod_p
   :description: Polynomials over ``F_p`` for a runtime prime.

Set Power Series
----------------

.. rust-api:: set_power_series
   :description: Subset transforms, subset convolution, exp, log and composition.

Multiplicative Convolution
--------------------------

.. rust-api:: multiplicative_convolution
   :description: Convolutions where the indices multiply.

Multivariate Convolution
------------------------

.. rust-api:: multivariate_convolution
   :description: Products of multivariate polynomials stored as dense arrays.

Min-Plus Convolution
--------------------

.. rust-api:: min_plus
   :description: Min-plus convolution of convex, concave and arbitrary sequences.

Linear Recurrences
------------------

.. rust-api:: linear_recurrence
   :description: Linear recurrences over ``Z/MZ``.

Matrices mod p
--------------

.. rust-api:: modmatrix
   :description: Dense matrices over ``Z/MZ``.

Bit Matrix over F2
------------------

.. rust-api:: bit_matrix
   :description: Matrices over ``F_2`` packed into ``u64`` words.

Spanning Tree Counting
----------------------

.. rust-api:: spanning_count
   :description: Counting spanning structures with determinants.
