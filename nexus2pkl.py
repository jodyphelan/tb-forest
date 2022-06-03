import sys
from Bio import Phylo
import ete3
import os
from uuid import uuid4
from collections import namedtuple
# M = namedtuple("mutation","pos nuc")
from copy import copy
import random
import pathogenprofiler as pp
import os
import pickle


def parse_mutations(muts):
    results = set()
    for m in muts:
        anc = m[0]
        mut = m[-1]
        pos = int(m[1:-1])
        results.add((pos,mut))
    return results

def nexus2ete3(filename):
    tree = list(Phylo.parse(filename, "nexus"))[0]
    mutations = {}
    for c in tree.find_clades():
        if c.name==None:
            c.name = copy(c.confidence)
            c.confidence=None

        id =  c.name
        if c.comment!=None:
            mutations[id] = parse_mutations(c.comment.replace('[&mutations=\"','').replace('"]','').split(","))
        c.comment = None

    tmp = str(uuid4())
    Phylo.write(tree,tmp,"newick")
    import ete3
    t = ete3.Tree(tmp,1)
    for n in t.traverse():
        n.add_features(mutations = mutations.get(n.name,None))
        n.dist = len(mutations.get(n.name,[]))
    os.remove(tmp)
    return t


t = nexus2ete3(sys.argv[1])
pickle.dump(t,open(sys.argv[2],"wb"))