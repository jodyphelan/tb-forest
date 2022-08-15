from collections import defaultdict
import sys
import argparse
from Bio import Phylo
import ete3
import os
from uuid import uuid4
import pathogenprofiler as pp
from copy import copy 
import pickle

def nexus2ete3(filename):
    def parse_mutations(muts):
        results = set()
        for m in muts:
            anc = m[0]
            mut = m[-1]
            pos = int(m[1:-1])
            results.add((pos,mut))
        return results
        
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
    t = ete3.Tree(tmp,1)
    for n in t.traverse():
        n.add_features(mutations = mutations.get(n.name,None))
        n.dist = len(mutations.get(n.name,[]))
    os.remove(tmp)
    return t


def load_missing_positions(vcf):
    positions = defaultdict(set)
    for l in pp.cmd_out(f"bcftools query -f '[%POS\t%SAMPLE\t%GT\n]' {vcf} | awk '$3==\".\"'"):
        row = l.strip().split()
        positions[row[1]].add(int(row[0]))
    return positions

def main(args):
    missing_positions = load_missing_positions(args.vcf)
    t = nexus2ete3(args.nexus)
    for leaf in t.get_leaf_names():
        mutations = (t & leaf).mutations
        if mutations==None: continue
        for mut in list(mutations):
            if mut[0] in missing_positions[leaf]:
                mutations.remove(mut)
    pickle.dump(t,open(args.out+".pkl","wb"))

parser = argparse.ArgumentParser(description='tbprofiler script',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--nexus',type=str,help='Nexus file')
parser.add_argument('--vcf',type=str,help='Nexus file')
parser.add_argument('--out',type=str,help='Nexus file')
parser.set_defaults(func=main)
args = parser.parse_args()
args.func(args)