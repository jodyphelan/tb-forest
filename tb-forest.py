from dataclasses import replace
from Bio import Phylo
import sys
import ete3
import os
from uuid import uuid4
from copy import copy
import pathogenprofiler as pp
import os
import pickle
import argparse



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

def load_vcf_mutations(vcf,ref):
    fa = pp.fasta(ref).fa_dict['Chromosome']
    positions = set()
    sample_mutations = set()
    for l in open(vcf):
        if "#" in l:continue
        row = l.strip().split()
        sample_mutations.add((int(row[1]),row[4]))
        positions.add(int(row[1]))
    for i in range(len(fa)):
        if i+1 in positions:
            continue
        sample_mutations.add((i+1,fa[i:i+1]))
    return sample_mutations

def compare_sample_to_branch(sample_mutations,branch_mutations):
    br_common = len(branch_mutations.intersection(sample_mutations))
    br_uniq = len(branch_mutations.difference(sample_mutations))
    return br_common/(br_common+br_uniq)




def local_phylo_reconstruct(samples,ref,snippy_dir,outgroup,exclude_bed,working_dir="/tmp",clean=True):
    current_dir = os.getcwd()
    tmpdir = f"{working_dir}/{str(uuid4())}"
    os.mkdir(tmpdir)
    os.chdir(tmpdir)
    print(tmpdir)
    snippy_dirs = " ".join([f"{snippy_dir}/{s}" for s in samples])
    pp.run_cmd(f"snippy-core --mask {exclude_bed} -r {ref} {snippy_dirs}")
    pp.run_cmd(f"fasta_remove_seq.py --fasta core.full.aln --seqs Reference > core.noref.aln")
    pp.run_cmd(f"iqtree -s core.noref.aln -m GTR+F+I -czb")
    pp.run_cmd(f"tree_root_on_outgroup.py  --tree core.noref.aln.treefile --outfile core.noref.aln.rooted.treefile --outgroup {outgroup}")
    pp.run_cmd(f"treetime ancestral --aln core.vcf --vcf-reference {ref} --tree core.noref.aln.rooted.treefile  --outdir ancestral")
    t = nexus2ete3("ancestral/annotated_tree{}.nexus")
    os.chdir(current_dir)
    if clean:
        pp.run_cmd(f"rm -rf {tmpdir}")
    return t
    
    
def replace_node(phy,samples_to_replace,replacement_node):
    if len(samples_to_replace)==1:
        node_to_replace = phy & samples_to_replace[0]
    else:
        node_to_replace = phy.get_common_ancestor(samples_to_replace)
    pp.warninglog("*"*40)
    pp.debug("Replacing:")
    print(node_to_replace)
    pp.debug("With:")
    print(replacement_node)
    pp.warninglog("*"*40)
    node_parent = node_to_replace.get_ancestors()[0]
    node_to_replace.detach()
    node_parent.add_child(replacement_node)

def replace_node_children(phy,samples_to_replace,replacement_node):
    if len(samples_to_replace)==1:
        node_to_replace = phy & samples_to_replace[0]
    else:
        node_to_replace = phy.get_common_ancestor(samples_to_replace)
    
    pp.warninglog("*"*40)
    pp.debug("Replacing:")
    print(node_to_replace)
    pp.debug("With:")
    print(replacement_node)
    pp.warninglog("*"*40)
    node_to_replace.write(format=1,outfile="test.newick")
    for n in list(node_to_replace.children):
        pp.debug("Detaching:")
        print(n)
        n.detach()
    
    pp.debug("Node now looks like this:")
    print(node_to_replace)
    for n in replacement_node.children:
        pp.debug("Attaching:")
        print(n)
        node_to_replace.add_child(n)
    pp.debug("Node now looks like this:")
    print(node_to_replace)

def is_monophyletic(t,leaves):
    if set(t.get_common_ancestor(leaves).get_leaf_names())==set(leaves):
        return True
    else:
        return False
def czb(t,cut=0):
    sys.stderr.write("Collapsing zero length branches\n")
    change_made = True
    while change_made:
        change_made = False
        for n in t.traverse():
            if n.is_root():
                continue
            if n.dist<=cut:
                p = n.get_ancestors()[0]
                for c in list(n.children):
                    p.add_child(c.detach())
                if n.is_leaf() and "DRR" not in n.name and "SRR" not in n.name and "ERR" not in n.name:
                    n.detach()
                    change_made = True
                    break
        
    return t

def get_close_samples(t,sample_name,cutoff=10):
    samples = []
    for s in t.get_leaf_names():
        if s==sample_name: 
            continue
        if t.get_distance(sample_name,s)<=cutoff:
            samples.append(s)
    return samples
    
def flatten(l):
    return [item for sublist in l for item in sublist]

def print_marked_branch(t,branch_name):
    for n in t.traverse():
        if n.is_leaf():
            n.add_features(tmp=n.name)

    (t & branch_name).add_features(tmp="******")

    print(t.get_ascii(attributes=["tmp"],show_internal=True))

def main(args):
    args.snippy_dir = os.path.abspath(args.snippy_dir)
    t = pickle.load(open(args.master_tree,"rb"))

    from collections import Counter
    if Counter(t.get_leaf_names())["ERR4553785"]>1:
        quit("ERR4553785 is not unique")
    t = czb(t,cut=1)
    # new_sample = "SRR8651557"
    input_vcf = f"{args.snippy_dir}/{args.new_sample}/snps.vcf"
    # outgroup = "ERR4553785"
    sample_mutations = load_vcf_mutations(input_vcf,args.ref)

    def dont_go_further(node):
        if node.mutations==None: 
            return False
        if compare_sample_to_branch(sample_mutations,node.mutations)<0.1:
            return True
        else:
            return False

    results = []
    for n in t.traverse(is_leaf_fn=dont_go_further):
        if n.is_root():
            continue
        if n.mutations==None:
            continue
        
        results.append((n.name,compare_sample_to_branch(sample_mutations,n.mutations),len(n.get_ancestors())))

    print(results)
    filtered_results = [r for r in results if r[1]>0.5]
    print(filtered_results)
    nodeA = t & filtered_results[-1][0]
    print(nodeA.mutations)
    if nodeA.is_leaf():
        nodeA = nodeA.get_ancestors()[0]

    # nodeA = nodeA.get_ancestors()[0]
    while True:
        nodeD = nodeA.get_ancestors()[0]
        outclade_nodes = [n for n in nodeD.children if n.name!=nodeA.name]
        outclade = flatten([n.get_leaf_names()[:3] for n in outclade_nodes])
        pp.debug("Outclade: %s" % str(outclade))



        print_marked_branch(nodeA.get_ancestors()[0],nodeA.name)
        children_nodes = nodeA.children
        children_node_reps = [set(n.get_leaf_names()[:2] + n.get_leaf_names()[-2:]) for n in children_nodes]
        representative_children = flatten(children_node_reps)
        tmp_samples =  set(representative_children + outclade + [args.new_sample, args.outgroup])
        



        x = local_phylo_reconstruct(
            samples = tmp_samples,
            ref = args.ref,
            snippy_dir = args.snippy_dir,
            outgroup = args.outgroup,
            working_dir=args.working_dir,
            exclude_bed=args.exclude_bed,
            clean=args.no_clean
        )
        x = czb(x,cut=1)
        if not is_monophyletic(x,representative_children + [args.new_sample]):
            nodeA = nodeA.get_ancestors()[0]
            pp.warninglog("Not monophyletic, going back to a higher node")
        else:
            break
    print("*"*40)
    print(x.get_ascii(attributes=["name","dist"],show_internal=True))
    print("*"*40)
    for n in x.traverse():

        if n.is_leaf(): continue
        if n.is_root(): continue
        n.name = str(uuid4())

    
    pp.debug("Outclade samples: %s" % outclade)
    pp.debug("Representative children samples: %s" % representative_children)
    
    print(children_node_reps)
    for i in range(len(children_nodes)):
        print(children_node_reps[i])
        if len(children_node_reps[i])==1:
            pp.warninglog("Node only has one sample... doing nothing")
            children_nodes.pop(0)
        elif is_monophyletic(x,children_node_reps[i]):
            pp.infolog("Node is monophyletic and will be replaced")
            replace_node(x,children_node_reps[i],children_nodes[0].detach())
        else:
            
            if len(t.get_common_ancestor(children_node_reps[i]))==len(children_node_reps[i]):
                pp.warninglog("Clade has been broken up but is complete... skipping")
            else:
                pp.errlog("Don't know what do do here",True)
            children_nodes.pop(0)
            

    reconstructed_node = x.get_common_ancestor(representative_children + [args.new_sample])
    nodeA.detach()
    nodeD.add_child(reconstructed_node)



    print(t)
    if args.outgroup in reconstructed_node.get_leaf_names():
        print("Outroup is in the reconstructed tree")
        quit()


    for n in reconstructed_node.traverse():
        if n.is_root(): continue
        if n.mutations==None: continue
        if n.get_ancestors()[0].mutations==None: continue
        for m in list(n.mutations):
            if m in n.get_ancestors()[0].mutations:
                n.mutations.remove(m)
        n.dist = len(n.mutations)

    
    
    pickle.dump(t,open(args.out+".pkl","wb"))
    t.write(format=1,outfile=args.out+".newick")
    


parser = argparse.ArgumentParser(description='tbprofiler script',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--master-tree',type=str,help='Pickle file containing the master tree',required=True)
parser.add_argument('--new-sample',type=str,help='New sample to add to the tree',required=True)
parser.add_argument('--out',type=str,help='Output file prefix',required=True)
parser.add_argument('--snippy-dir',type=str,help='Directory containing snippy output')
parser.add_argument('--outgroup',type=str,help='Outgroup to use for rooting')
parser.add_argument('--working-dir',type=str,help='Working directory',default="/tmp")
parser.add_argument('--ref',type=str,help='Reference genome',default="/home/jody/refgenome/MTB-h37rv_asm19595v2-eg18.fa")
parser.add_argument('--exclude-bed',required = True)
parser.add_argument('--no-clean',action="store_false",help='Do not clean up working directory')
parser.set_defaults(func=main)
args = parser.parse_args()
args.func(args)