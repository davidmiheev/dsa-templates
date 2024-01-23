def generate_subsequences(arr):
    n = len(arr)
    result = []

    for i in range(2**n):
        subseq = []
        for j in range(n):
            if i & (1 << j):
                subseq.append(arr[j])

        result.append(subseq)

    return result

arr = [1, 2, 3]
print(generate_subsequences(arr))
