print('all_subs')

def generate_subsequences(arr):
    '''
    Given an array of integers, arr, return all possible subsequences of the array.

    Time complexity: :math:`O(2^n)`

    :type arr: List[int]
    :rtype: List[List[int]]
    '''
    n = len(arr)

    for i in range(2 ** n):
        subseq = []
        for j in range(n):
            if i & (1 << j):
                subseq.append(arr[j])

        yield subseq

arr = [1, 2, 3]
for sub in generate_subsequences(arr):
    print(sub)
