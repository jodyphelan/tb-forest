import sys
import argparse
import pickle
import ete3


def czb(t,cut=0):
    sys.stderr.write("Collapsing zero length branches\n")
    change_made = True
    leaves = set(t.get_leaf_names())
    while change_made:
        change_made = False
        for n in t.traverse():
            if n.is_root():
                continue
            if n.dist<=cut:
                p = n.get_ancestors()[0]
                for c in list(n.children):
                    p.add_child(c.detach())
                if n.is_leaf() and n.name not in leaves:
                    n.detach()
                    change_made = True
                    break
        
    return t

def main(args):
    t = pickle.load(open(args.tree,'rb'))
    positions = set([int(x.strip()) for x in open(args.positions)])
    for n in t.traverse():
        if n.mutations:
            n.mutations = [m for m in n.mutations if m[0] not in positions]
            n.dist = len(n.mutations)
    t = czb(t)
    pickle.dump(t,open(args.out+".pkl",'wb'))
    t.write(format=1,outfile=args.out+".newick")


parser = argparse.ArgumentParser(description='tbprofiler script',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--tree',type=str,help='Tree file in pickle format')
parser.add_argument('--positions',type=str,help='File with positions to remove')
parser.add_argument('--out',type=str,help='Output file')
parser.set_defaults(func=main)
args = parser.parse_args()
args.func(args)