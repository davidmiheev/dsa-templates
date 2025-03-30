from collections import defaultdict, deque

print('tree')
class TreeNode:
    def __init__(self, val=0, left=None, right=None):
        self.val = val
        self.left = left
        self.right = right
    
    def build_tree(self, values) -> 'TreeNode':
        '''
        Build a binary tree from a list of values.
        The values are given in level order.
        '''
        if not values:
            self.val = None
            return self
        self.val = values[0]
        queue = deque([self])
        values = values[1:][::-1]
        while values:
            node = queue.popleft()
            val = values.pop()
            if val is not None:
                node.left = TreeNode(val)
                queue.append(node.left)
            if not values:
                break
            val = values.pop()
            if val is not None:
                node.right = TreeNode(val)
                queue.append(node.right)

        return self

    def level_order(self) -> list[int]:
        '''
        Traverse the tree in level order and return the values of the nodes.
        '''
        
        queue = deque([self])
        levels = []
        
        while queue:
            node = queue.popleft()
            levels.append(node.val)
            
            if node.left:
                queue.append(node.left)
            if node.right:
                queue.append(node.right)
        
        return levels
    
    def inorder(self) -> list[int]:
        '''
        Traverse the tree in inorder and return the values of the nodes.
        '''
        l = self.left.inorder() if self.left else []
        r = self.right.inorder() if self.right else []
        return l + [self.val] + r
    
    def preorder(self) -> list[int]:
        '''
        Traverse the tree in preorder and return the values of the nodes.
        '''
        l = self.left.preorder() if self.left else []
        r = self.right.preorder() if self.right else []
        return [self.val] + l + r
    
    def postorder(self) -> list[int]:
        '''
        Traverse the tree in postorder and return the values of the nodes.
        '''
        l = self.left.postorder() if self.left else []
        r = self.right.postorder() if self.right else []
        return l + r + [self.val]

root = TreeNode()
root = root.build_tree([])
print(f'level order: {root.level_order()}')
print(f'inorder: {root.inorder()}')
print(f'preorder: {root.preorder()}')
print(f'postorder: {root.postorder()}')
