from itertools import accumulate 
import operator

print('one-line')

arr = [1,2,3,4,5,6,7,8]
def prefix_sums(arr: list[int]) -> list[int]:
  '''
  calculate array of prefix sums
  '''
  return [0] + list(accumulate(arr, operator.add))

print(prefix_sums(arr))

matrix = [[1, 2, 3], [4, 5, 6]]
def transpose(matrix: list[list[int]]) -> list[int]:
  '''
  transpose matrix given as python nested list
  '''
  return list(map(list, zip(*matrix)))

print(transpose(matrix))

dict1 = {'a': 1, 'b': 2}
dict2 = {'b': 3, 'c': 4}
def merge_dicts(dict1: dict, dict2: dict) -> dict:
  '''
  merge two python dictionaries
  '''
  return {**dict1, **dict2}
