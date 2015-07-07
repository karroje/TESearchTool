import sys
import re

class tree:
    """Implementation of a (directed) tree.
       Each element of the children field is
       a (tree_node, weight) tuple"""

    def __init__(self, label, children, sequence):
        """Constructor: create a root node with
        the specified label and child set.
        Assumes children is a tuple of tree objects."""
        self.label = label
        self.children = children
        self.sequence = sequence    # each node corresponds to a sequence

    def getLabel(self):
        """Get the label of the node"""
        return self.label

    def child(self, i):
        """Get the ith child.   Assumes
        the tree has at least i+1 children."""
        return self.children[i]

    def set_seq(self, s):   # new
        """set the sequence attribute to s"""
        self.sequence = s

    def numChildren(self):
        """Returns the number of children."""
        return len(self.children)


    def isLeaf(self):
        """Indicates if the tree has children."""
        return self.numChildren()==0

    def countLeaves(self):
        if self.isLeaf():
            return 1
        else:
            counter = 0
            for i in range(self.numChildren()):
                counter += self.child(i)[0].countLeaves()
            return counter

    def totalWeight(self):
        weight = sum([t.totalWeight() + int(w) for t,w in self.children])
        return weight

    def leafMap(self, f):
        leaf_list = []
        self.leafMapHelper(f, leaf_list)
        return leaf_list

    def leafMapHelper(self, f, leaf_list):
        if self.isLeaf():
            leaf_list.append(f(self.label))
            return leaf_list
        else:
            for i in range(self.numChildren()):
                self.child(i)[0].leafMapHelper(f, leaf_list)
            return leaf_list

    def leafList(self):
        return self.leafMap(lambda x: x)

    def leafSeqs(self):   # new
        leaf_seq_list = []
        self.leafSeqsHelper(leaf_seq_list)
        return leaf_seq_list

    def leafSeqsHelper(self, leaf_seq_list):  # new
        if self.isLeaf():
            leaf_seq_list.append(self.sequence)
            return leaf_seq_list
        else:
            for i in range(self.numChildren()):
                self.child(i)[0].leafSeqsHelper(leaf_seq_list)
            return leaf_seq_list

    def nodeDepths(self):
        dic = {}     # init an empty dictionary, how to define the dic['root'] = 0?
        self.nodeDepthsHelper(dic, 0)
        return dic

    def nodeDepthsHelper(self, dic, cur_depth):
        if self.isLeaf() == False:      # while it has children
            dic[self.label] = cur_depth
            for i in range(self.numChildren()):
                cur_depth += int(self.child(i)[1])
                if self.child(i)[0].label != None:     # if this node has a label
                    dic[self.child(i)[0].label] = cur_depth # child(self, i)[1] is the weight of this edge
                self.child(i)[0].nodeDepthsHelper(dic, cur_depth)
                cur_depth -= int(self.child(i)[1])
            return dic
        else:
            return dic

    def printNewick(self):
        if self.isLeaf():
            return self.label
        else:
            return "(" + ",".join([t.printNewick() + ":" + str(w) for t,w in self.children]) + ")" + self.label
            
        

clauseRE = re.compile("(.*):([^:]*)")
labelRE = re.compile("(\(.*\))?([^\(\)]*)")

def parseNewick(s):
    """Create a tree from a passed string"""
    r = re.search("^(.*):(\d+(.\d*)?)\s+$", s)
    if r:
        return tree(label = "", children = [(parseNewick(r.group(1)), float(r.group(2)))], sequence = "")        
    s2,l = labelRE.search(s).groups()
    if s2 == None:
        return tree(label=l, children = [], sequence = "")    
    return tree(label = l, children = [(parseNewick(t),w) for t,w in splitTopLevel(s)], sequence = "")


def splitTopLevel(s):
    """Given a string in Newick Format, split it at the top level ---
    returning a list of the strings for each subtree"""
    if s[0] != "(":
        return clauseRE.search(s).groups()
    breakPoints = [0]
    p = 1
    level = 1
    while p < len(s):
        if level == 1 and (s[p] in (",", ")")):
            breakPoints.append(p)
        if s[p] == "(":
            level = level + 1
        elif s[p] == ")":
            level = level - 1
        p = p + 1
    results = []
    for i in range(len(breakPoints) - 1):
        results.append(s[breakPoints[i]+1:breakPoints[i+1]].rstrip().strip())
    return [clauseRE.search(tree).groups() for tree in results]
