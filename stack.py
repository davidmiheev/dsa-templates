print('stack')

def is_valid(s: str) -> bool:
    '''
    Given a string containing just the characters '(', ')', '{', '}', '[' and ']',
    determine if the input string is valid.

    An input string is valid if:
        Open brackets must be closed by the same type of brackets.


    Time complexity: :math:`O(n)`

    Space complexity: :math:`O(n)`
    '''
    stack = []
    for ch in s:
        if ch == ')' and stack and stack[-1] == '(':
            stack.pop()
        elif ch == ']' and stack and stack[-1] == '[':
            stack.pop()
        elif ch == '}' and stack and stack[-1] == '{':
            stack.pop()
        else:
            stack.append(ch)

    return not stack


s = "()[()]{}"
print(is_valid(s))
s = "([)]"
print(is_valid(s))
