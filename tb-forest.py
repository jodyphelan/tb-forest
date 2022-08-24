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
from tqdm import tqdm
from collections import Counter, defaultdict

__version__ = "0.0.5"

def nexus2ete3(filename):
    def parse_mutations(muts):
        results = set()
        for m in muts:
            if m=="": continue
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

def walk_paths(n, sample_mutations, paths=None, current_path=None,skips=0):
    if paths is None:
        paths = []
    if current_path is None:
        current_path = []

    
 
    if n.mutations:
        overlap = compare_sample_to_branch(sample_mutations,n.mutations)
        current_path.append((n.name,overlap,len(n.mutations)))
    else:
        overlap = 0
        current_path.append((n.name,overlap,0))
    
    
        
    
    if len(n.children) == 0:
        paths.append(current_path)
    else:
        if n.mutations:
            if len(n.mutations)<2:
                if compare_sample_to_branch(sample_mutations,n.mutations)<0.1:
                    skips+=1
                if skips>1:
                    paths.append(current_path)
                else:
                    for child in n.children:
                        walk_paths(child, sample_mutations, paths, list(current_path),skips)
            else:
                if compare_sample_to_branch(sample_mutations,n.mutations)<0.1:
                    paths.append(current_path)
                else:
                    for child in n.children:
                        walk_paths(child, sample_mutations, paths, list(current_path),skips)
        else:
            for child in n.children:
                walk_paths(child, sample_mutations, paths, list(current_path),skips)
    return paths


def get_paths(n,sample_mutations):
    paths = walk_paths(n,sample_mutations)
    for p in paths:
        while p[-1][1]==None or p[-1][1]<0.5:
            p.pop()
            if len(p)==0:
                break
    return sorted(paths,key=lambda x:sum([n[1]*n[2] for n in x]))[-1]



def load_vcf_mutations(vcf,ref):
    fa = pp.fasta(ref).fa_dict['Chromosome']
    positions = set()
    sample_mutations = set()
    for l in tqdm(pp.cmd_out(f"bcftools norm -a {vcf}")):
        if "#" in l:continue
        row = l.strip().split()
        sample_mutations.add((int(row[1]),row[4]))
        positions.add(int(row[1]))
    for i in tqdm(range(len(fa))):
        if i+1 in positions:
            continue
        sample_mutations.add((i+1,fa[i:i+1]))
    return sample_mutations

# def load_mixed_positions(fasta):
#     x  = set()
#     fa = pp.fasta(fasta).fa_dict['Chromosome']
#     positions = set()
#     for i in tqdm(range(len(fa))):
#         x.add(fa[i])
#         if fa[i]=="N":
#             positions.add(i+1)
#     print(len(x))
#     return positions

def load_mixed_positions(vcf):
    positions = set()
    for l in open(vcf):
        if "#" in l: continue
        if "TYPE=ins" in l: continue
        if "TYPE=del" in l: continue
        row = l.strip().split("\t")
        p = int(row[1])
        genotype = row[9].split(":")[0]
        x = genotype.split("/")
        if x[0]==x[1]:
            continue
        for r,a in zip(row[3],row[4]):
            if r!=a:
                positions.add(p)
            p+=1     
    return positions

def compare_sample_to_branch(sample_mutations,branch_mutations):
    br_common = len(branch_mutations.intersection(sample_mutations))
    br_uniq = len(branch_mutations.difference(sample_mutations))
    return br_common/(br_common+br_uniq)


def load_missing_positions(vcf):
    positions = defaultdict(set)
    for l in pp.cmd_out(f"bcftools query -f '[%POS\t%SAMPLE\t%GT\n]' {vcf} | awk '$3==\".\"'"):
        row = l.strip().split()
        positions[row[1]].add(int(row[0]))
    return positions

def local_phylo_reconstruct(samples,ref,snippy_dir,outgroup,exclude_bed,working_dir="/tmp",clean=True,id=None,new_sample=None):
    current_dir = os.getcwd()
    if id:
        os.chdir(f"{working_dir}/{id}")
        t = nexus2ete3("ancestral/annotated_tree{}.nexus")
        os.chdir(current_dir)
        return t
    tmpdir = f"{working_dir}/{str(uuid4())}"
    os.mkdir(tmpdir)
    os.chdir(tmpdir)
    print(tmpdir)
    snippy_dirs = " ".join([f"{snippy_dir}/{s}" for s in samples])
    # pp.run_cmd(f"snippy-core --mask {exclude_bed} -r {ref} {snippy_dirs}")
    # pp.run_cmd(f"fasta_remove_seq.py --fasta core.full.aln --seqs Reference > core.noref.aln")
    # pp.run_cmd(f"mv core.noref.aln core.full.aln")
    pp.run_cmd(f"snippy-cloud.py --mask {exclude_bed} -r {ref} {snippy_dirs} ")
    pp.run_cmd(f"iqtree -s core.full.aln -m GTR+F+I -czb")
    pp.run_cmd(f"tree_root_on_outgroup.py  --tree core.full.aln.treefile --outfile core.full.aln.rooted.treefile --outgroup {outgroup}")
    pp.run_cmd(f"treetime ancestral --aln core.vcf --vcf-reference {ref} --tree core.full.aln.rooted.treefile  --outdir ancestral")
    t = nexus2ete3("ancestral/annotated_tree{}.nexus")
    # pp.run_cmd("python /home/jody/github/tb-forest/correct_missing_positions_asr.py --nexus ancestral/annotated_tree{}.nexus --vcf core.vcf  --out tree")
    # t = pickle.load(open("tree.pkl","rb"))
    
    
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

def get_close_samples(t,sample_name,cutoff=10):
    samples = []
    for s in t.get_leaf_names():
        if s==sample_name: 
            continue
        d = t.get_distance(sample_name,s)
        if d<=cutoff:
            samples.append((s,d))
    return samples
    
def flatten(l):
    return [item for sublist in l for item in sublist]

def print_marked_branch(t,branch_name):
    for n in t.traverse():
        if n.is_leaf():
            n.add_features(tmp=n.name)

    (t & branch_name).add_features(tmp="******")

    print(t.get_ascii(attributes=["tmp"],show_internal=True))

def main_cut(args):
    t = pickle.load(open(args.master_tree,"rb"))
    samples = get_close_samples(t,args.sample,args.cutoff)
    print(samples)

def rename_internal_nodes(t):
    for n in t.traverse():
        if n.is_leaf(): continue
        if n.is_root(): continue
        n.name = str(uuid4())

def check_if_clusters_with_one_leaf(x,new_sample):
    y = (x & new_sample).get_ancestors()[0]

    if len(y.children)==2 and y.children[0].is_leaf() and y.children[1].is_leaf():
        l = [l for l in y.get_leaf_names() if l!=args.new_sample][0]
        return l
    else:
        return False

    


def check_for_mix(mixed_positions,t):
    high_mixed_branches = []
    for n in t.traverse():
        if n.mutations:
            branch_positions = set([x[0] for x in n.mutations if n.mutations])
            f = len(branch_positions.intersection(mixed_positions))/len(branch_positions)
            if f>0.5:
                high_mixed_branches.append((n.name,f,f*len(branch_positions)))
    pp.debug(high_mixed_branches)
    if len(high_mixed_branches)>0:
        return True
    else:
        return False
    

def main_add(args):
    sys.stderr.write(f"Using version: {__version__}\n")
    args.snippy_dir = os.path.abspath(args.snippy_dir)
    print("Loading tree")
    t = pickle.load(open(args.master_tree,"rb"))
    tcopy = t.copy("deepcopy")
    original_samples = set(t.get_leaf_names())
    print("Loaded tree")
    t = czb(t,cut=0)
    print("finished")
    # new_sample = "SRR8651557"
    input_vcf = f"{args.snippy_dir}/{args.new_sample}/snps.vcf"
    raw_vcf = f"{args.snippy_dir}/{args.new_sample}/snps.raw.vcf"
    input_consensus = f"{args.snippy_dir}/{args.new_sample}/snps.aligned.fa"
    # outgroup = "ERR4553785"
    sample_mutations = load_vcf_mutations(input_vcf,input_consensus)
    mixed_positions = load_mixed_positions(raw_vcf)

    path = get_paths(t,sample_mutations)
    
    nodeA = t & path[-1][0]
    print(nodeA.mutations)
    if nodeA.is_leaf():
        nodeA = nodeA.get_ancestors()[0]


    #### check for mixed infection
    if check_for_mix(mixed_positions,t):
        pp.warninglog("Mixed infection detected")
        pickle.dump(t,open(args.out+".pkl","wb"))
        t.write(format=1,outfile=args.out+".newick")
        quit()
        
    
    # nodeA = nodeA.get_ancestors()[0]
    
    nodeD = nodeA.get_ancestors()[0]
    nodeO = nodeD
    outclade = []
    children_nodes = nodeA.children
    children_node_reps = [set(n.get_leaf_names()[:2] + n.get_leaf_names()[-2:]) for n in children_nodes]
    representative_children = flatten(children_node_reps)

    while True:
        outclade_nodes = [n for n in nodeO.children if n.name!=nodeA.name]
        outclade = outclade + flatten([n.get_leaf_names()[:3] for n in outclade_nodes])
        pp.debug("Outclade: %s" % str(outclade))

        tmp_samples =  set(representative_children + outclade + [args.new_sample, args.outgroup])
        



        x = local_phylo_reconstruct(
            samples = tmp_samples,
            ref = args.ref,
            snippy_dir = args.snippy_dir,
            outgroup = args.outgroup,
            working_dir=args.working_dir,
            exclude_bed=args.exclude_bed,
            clean=args.no_clean,
            new_sample = args.new_sample,
            # id="6ab5bd82-7d8c-4c11-9a3e-c52655d46801"
        )
        x = czb(x,cut=0)
        print("*"*40)
        print(x.get_ascii(attributes=["name","dist"],show_internal=True))
        print("*"*40)
        if not is_monophyletic(x,representative_children + [args.new_sample]):
            if check_if_clusters_with_one_leaf(x,args.new_sample):
                break
            else:
                # print_marked_branch(nodeA.get_ancestors()[0],nodeA.name)
                # print(x)
                nodeO = nodeO.get_ancestors()[0]
                pp.warninglog("Not monophyletic, going back to a higher node")
        else:
            pp.successlog("Is monophyletic, continuing...")
            break

    
    rename_internal_nodes(x)

    
    pp.debug("Outclade samples: %s" % outclade)
    pp.debug("Representative children samples: %s" % representative_children)
    
    
    err = None
    print(children_node_reps)
    if check_if_clusters_with_one_leaf(x,args.new_sample):
        clustered_sample = check_if_clusters_with_one_leaf(x,args.new_sample)
        replacement_node = (x & args.new_sample).get_ancestors()[0].detach()
        replace_node(t,[clustered_sample],replacement_node)
    else:
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
                    children_nodes.pop(0)
                else:
                    err = "Don't know what to do"
                    break

        if err:
            pp.warninglog("Doing something complicated  ")
            print(nodeA)
            tmp_samples = set(nodeA.get_leaf_names() + outclade + [args.new_sample, args.outgroup])
            x = local_phylo_reconstruct(
                samples = tmp_samples,
                ref = args.ref,
                snippy_dir = args.snippy_dir,
                outgroup = args.outgroup,
                working_dir=args.working_dir,
                exclude_bed=args.exclude_bed,
                clean=args.no_clean,
                new_sample = args.new_sample,
                # id="6ab5bd82-7d8c-4c11-9a3e-c52655d46801"
            )
            rename_internal_nodes(x)
            nodeA_leaves = nodeA.get_leaf_names()
            replacement_node = x.get_common_ancestor(nodeA_leaves + [args.new_sample])
            replace_node(tcopy,nodeA.get_leaf_names(),replacement_node)
            t = tcopy
        else:
            print(x)
            print(representative_children)
            reconstructed_node = x.get_common_ancestor(representative_children + [args.new_sample])
            nodeA.detach()
            nodeD.add_child(reconstructed_node)

            for n in reconstructed_node.traverse():
                if n.is_root(): continue
                if n.mutations==None: continue
                if n.get_ancestors()[0].mutations==None: continue
                for m in list(n.mutations):
                    if m in n.get_ancestors()[0].mutations:
                        n.mutations.remove(m)
                n.dist = len(n.mutations)
            
            if args.outgroup in reconstructed_node.get_leaf_names():
                print("Outroup is in the reconstructed tree")
                quit()

    new_tree_samples = set(t.get_leaf_names()) 
    if len(original_samples) +  1 != len(new_tree_samples):
        pp.errlog(new_tree_samples - original_samples)
        pp.errlog(original_samples - new_tree_samples)
        pp.errlog("Lost some samples",True)
    
    print(t)
    

    pickle.dump(t,open(args.out+".pkl","wb"))
    t.write(format=1,outfile=args.out+".newick")
    
    close_samples = get_close_samples(t,args.new_sample)
    print(close_samples)
    x = t.get_common_ancestor([x[0] for x in close_samples] + [args.new_sample])

    positions = []
    for n in t.traverse():
        if n.mutations==None: continue
        for m in n.mutations:
            positions.append(m[0])
    print([m for m in Counter(positions).most_common() if m[1]>1])

parser = argparse.ArgumentParser(description='TB-Forest pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
subparsers = parser.add_subparsers(help="Task to perform")

parser_sub = subparsers.add_parser('add', help='Run whole profiling pipeline', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_sub.add_argument('--master-tree',type=str,help='Pickle file containing the master tree',required=True)
parser_sub.add_argument('--new-sample',type=str,help='New sample to add to the tree',required=True)
parser_sub.add_argument('--out',type=str,help='Output file prefix',required=True)
parser_sub.add_argument('--snippy-dir',type=str,help='Directory containing snippy output')
parser_sub.add_argument('--outgroup',type=str,help='Outgroup to use for rooting')
parser_sub.add_argument('--working-dir',type=str,help='Working directory',default="/tmp")
parser_sub.add_argument('--ref',type=str,help='Reference genome',default="/home/jody/refgenome/MTB-h37rv_asm19595v2-eg18.fa")
parser_sub.add_argument('--exclude-bed',required = True)
parser_sub.add_argument('--no-clean',action="store_false",help='Do not clean up working directory')
parser_sub.set_defaults(func=main_add)

parser_sub = subparsers.add_parser('cut', help='Run whole profiling pipeline', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_sub.add_argument('--master-tree',type=str,help='Pickle file containing the master tree',required=True)
parser_sub.add_argument('--sample',type=str,help='Pickle file containing the master tree',required=True)
parser_sub.add_argument('--cutoff',type=int,help='New sample to add to the tree',required=True)
parser_sub.set_defaults(func=main_cut)

args = parser.parse_args()
if hasattr(args, 'func'):
    args.func(args)
else:
    parser.print_help(sys.stderr)






# parser = argparse.ArgumentParser(description='NTM-Profiler pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
# subparsers = parser.add_subparsers(help="Task to perform")




# # Update database #
# parser_sub = subparsers.add_parser('update_db', help='Update all databases', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
# parser_sub.add_argument('--branch','-b',default="main",help='Storage directory')

# parser_sub.add_argument('--dir','-d',default=".",help='Storage directory')
# parser_sub.set_defaults(func=main_cut)







# args = parser.parse_args()
# if hasattr(args, 'func'):
#     args.func(args)
# else:
#     parser.print_help(sys.stderr)