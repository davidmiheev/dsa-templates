print('trie')

class TrieNode:
    '''
    Trie Node
    '''

    def __init__(self, isValid=False, val=''):
        self.children = [None for _ in range(26)]
        self.valid = isValid

class Trie:
    '''
    Trie (prefix tree) implementation.
    Trie is a tree-like data structure whose nodes store the letters of an alphabet
    by structuring the nodes in a particular way, words and strings
    can be retrieved from the structure by traversing down a branch path of the tree.

    No node in the tree stores the key associated with that node;
    instead, its position in the tree defines the key with which it is associated.

    All the descendants of a node have a common prefix of the string associated with that node,
    and the root is associated with the empty string.

    Values are not necessarily associated with every node.

    Rather, values tend only to be associated with leaves, and with some inner nodes that correspond to keys of interest.
    '''

    def __init__(self):
        '''
        Initialize the Trie data structure.

        '''
        self.root = TrieNode()

    def insert(self, word: str) -> None:
        '''
        Inserts a word into the trie.

        Time complexity: :math:`O(n)`, where n is the length of the word.

        :type word: str
        :rtype: None
        '''
        cur = self.root
        for ch in word:
            if cur.children[ord(ch)-97] is None:
                cur.children[ord(ch)-97] = TrieNode()
            cur = cur.children[ord(ch)-97]

        cur.valid = True


    def search(self, word: str) -> bool:
        '''
        Returns if the word is in the trie.

        Time complexity: :math:`O(n)`, where n is the length of the word.

        :type word: str
        :rtype: bool
        '''
        cur = self.root
        for ch in word:
            if cur.children[ord(ch)-97] is None:
                return False
            cur = cur.children[ord(ch)-97]

        return cur.valid


    def startsWith(self, prefix: str) -> bool:
        '''
        Returns if there is any word in the trie that starts with the given prefix.

        Time complexity: :math:`O(n)`, where n is the length of the prefix.

        :type prefix: str
        :rtype: bool
        '''
        cur = self.root
        for ch in prefix:
            if cur.children[ord(ch)-97] is None:
                return False
            cur = cur.children[ord(ch)-97]

        return True


# Your Trie object will be instantiated and called as such:
obj = Trie()
obj.insert("asdf")
obj.insert("as")
obj.insert("asdd")
obj.insert("asdd")
obj.insert("apple")
print(obj.search("as"))
print(obj.search("asdd"))
print(obj.search("asddd"))
print(obj.search("asdddd"))
print(obj.search("apple"))
print(obj.startsWith("app"))
print(obj.startsWith("a"))
print(obj.startsWith("as"))
# param_3 = obj.startsWith(prefix)
