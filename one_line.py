from itertools import accumulate, chain, combinations
import operator

print('one-line')

arr = [1,2,3,4,5,6,7,8]
def prefix_sums(arr: list[int]) -> list[int]:
    '''
    calculate array of prefix sums
    '''
    return [0] + list(accumulate(arr, operator.add))

print(f'init arr = {arr}, prefix sums = {prefix_sums(arr)}')

matrix = [[1, 2, 3], [4, 5, 6]]
def transpose(matrix: list[list[int]]) -> list[int]:
    '''
    transpose matrix given as python nested list
    '''
    return list(map(list, zip(*matrix)))

print(f'init matrix is {matrix}, transposed matrix is {transpose(matrix)}')

dict1 = {'a': 1, 'b': 2}
dict2 = {'b': 3, 'c': 4}
def merge_dicts(dict1: dict, dict2: dict) -> dict:
    '''
    merge two python dictionaries
    '''
    return {**dict1, **dict2}

print(f'init dict1 = {dict1}, init dict2 = {dict2}, merged dict = {merge_dicts(dict1, dict2)}')

def power_set(arr: list[int]) -> list[list[int]]:
    '''
    generate power set of given array
    '''
    return list(chain.from_iterable(combinations(arr, r) for r in range(len(arr) + 1)))

arr = [1, 2, 3]
print(f'power set of {arr}: {power_set(arr)}')

def flatten(lst: list[list[int]]) -> list[int]:
    '''
    flatten nested list
    '''
    return [item for sublist in lst for item in sublist]

arr = [[1, 2], [3, 4], [5]]
print(f'init arr is {arr}, flattened arr is {flatten(arr)}')
