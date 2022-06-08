import argparse
import ete3

def is_monophyletic(t,leaves):
    if set(t.get_common_ancestor(leaves).get_leaf_names())==set(leaves):
        return True
    else:
        return False

def main(args):
    t1 = ete3.Tree(args.tree1,format=1)
    t2 = ete3.Tree(args.tree2,format=1)

    i = 0
    for n in t1.traverse():
        if n.is_leaf():
            continue
        leaf_names = n.get_leaf_names()
        if not is_monophyletic(t2,leaf_names):
            print(n)
            quit("Found different node")
            
        i+=1

    print(f"{i} nodes test... all are the same!")
parser = argparse.ArgumentParser(description='tbprofiler script',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--tree1',type=str,help='Pickle file containing the master tree',required=True)
parser.add_argument('--tree2',type=str,help='New sample to add to the tree',required=True)
parser.set_defaults(func=main)
args = parser.parse_args()
args.func(args)