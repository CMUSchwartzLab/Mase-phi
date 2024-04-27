
import copy
import gzip
from zipfile import ZipFile
import json
import numpy as np
from visualize import render_tumor_tree

def aggregation2summary(tree_distribution_aggregation):
    tree_distribution_summary = {'cp_tree': [], 'node_dict': [], 'node_dict_name': [], 'node_dict_re': [], 'tree_structure': [], 'freq': [], 'clonal_freq': [], 'vaf_frac': []}
    node_dict, node_dict_name, node_dict_re, tree_structure, final_tree_cp, clonal_freq, vaf_frac = (tree_distribution_aggregation[key] for key in [
        'node_dict', 'node_dict_name', 'node_dict_re', 'tree_structure', 'cp_tree', 'clonal_freq', 'vaf_frac'
    ])
    print(tree_structure)
    for tree_idx in range(len(tree_distribution_aggregation["tree_structure"])):
        combine_tree(node_dict[tree_idx], node_dict_name[tree_idx], node_dict_re[tree_idx], tree_structure[tree_idx], final_tree_cp[tree_idx], clonal_freq[tree_idx], vaf_frac[tree_idx], 'phylowgs',
                     tree_distribution_summary)
    return tree_distribution_summary

def subtree(tree_structure, node_dict, submut_list):
    node_dict_sub = {}
    tree_unit = TreeUnit(tree_structure)
    for node, mut_list in node_dict.items():
        node_keep = False
        for mut in mut_list:
            if mut in submut_list:
                node_keep = True
                node_dict_sub.setdefault(node, [])
                node_dict_sub[node].append(mut)
        if not node_keep:
            if not tree_unit.is_root(node):
                print('Remove ', node)
                tree_unit.delete_node(node)
    return tree_unit.tree, node_dict_sub


class TreeUnit:
    def __init__(self, tree_structure):
        self.tree = {parent: [child for child in children] for parent, children in tree_structure.items()}
        self.cp_tree = generate_cp(self.tree)

    def delete_node(self, idx):
        if self.is_root(idx):
            if self.num_children(idx) > 1:
                raise('Cannot delete root node with more than one child!')
            elif self.num_children(idx) == 1:
                child = self.tree[idx][0]
                del self.cp_tree[child]
                del self.tree[idx]
        elif self.is_leaf(idx):
            if idx not in self.cp_tree.keys():
                idx = int(idx)
            parent = self.cp_tree[idx]
            del self.cp_tree[idx]
            if self.num_children(parent) == 1:
                del self.tree[parent]
            else:
                self.tree[parent].remove(idx)
        else:
            parent = self.cp_tree[idx]
            children = self.tree[idx]
            self.tree[parent].remove(idx)
            del self.cp_tree[idx]
            for child in children:
                self.cp_tree[child] = parent
                self.tree[parent].append(child)
            del self.tree[idx]


    def is_leaf(self, idx):
        if idx not in self.tree.keys():
            return True
        else:
            return False

    def is_root(self, idx):
        if idx in self.tree.keys() and idx not in self.cp_tree.keys():
            return True
        else:
            return False

    def num_children(self, idx):
        if self.is_leaf(idx):
            return 0
        else:
            return len(self.tree[idx])
def combine_tree(node_dict, node_dict_name, node_dict_re, tree_structure, tree_cp, clonal_freq, vaf_frac, method, tree_distribution):
    tree_distribution_deepcopy = copy.deepcopy(tree_distribution)
    found_flag = False
    for idx in range(len(tree_distribution_deepcopy['cp_tree'])):
        if tree_cp == tree_distribution_deepcopy['cp_tree'][idx]:
            tree_distribution_deepcopy['freq'][idx] += 1
            if method == 'phylowgs':
                match_dict = match_trees_phylowgs(node_dict_re, tree_distribution_deepcopy['node_dict_re'][idx])
            elif method == 'deconv':
                match_dict = match_trees_deconv(node_dict_re, tree_structure, tree_distribution_deepcopy['node_dict_re'][idx], tree_distribution_deepcopy['tree_structure'][idx])

            for node, freq in clonal_freq.items():
                tree_distribution_deepcopy['clonal_freq'][idx][match_dict[node]].append(freq[0])
            for node, frac in vaf_frac.items():
                tree_distribution_deepcopy['vaf_frac'][idx][match_dict[node]].append(frac[0])
            found_flag = True
            break
    if not found_flag:
        tree_distribution_deepcopy['cp_tree'].append(tree_cp)
        tree_distribution_deepcopy['node_dict'].append(node_dict)
        tree_distribution_deepcopy['node_dict_name'].append(node_dict_name)
        tree_distribution_deepcopy['node_dict_re'].append(node_dict_re)
        tree_distribution_deepcopy['tree_structure'].append(tree_structure)
        tree_distribution_deepcopy['freq'].append(1)
        tree_distribution_deepcopy['clonal_freq'].append(clonal_freq)
        tree_distribution_deepcopy['vaf_frac'].append(vaf_frac)
    return tree_distribution_deepcopy

def match_trees_phylowgs(node_dict_re, target_node_dict_re):
    match_dict = {0: 0}
    for muts, node in target_node_dict_re.items():
        #print(node_dict_re)
        match_dict[node_dict_re[muts]] = node
    #print(match_dict, 'm')
    return match_dict

def match_trees_deconv(node_dict_re, tree_structure, target_node_dict_re, target_tree_structure):
    match_dict = {}
    for muts, node in target_node_dict_re.items():
        match_dict[node_dict_re[muts]] = node
    for node in tree_structure.keys():
        if node not in match_dict.keys():
            match_dict[node] = node
    #print(match_dict, 'm')
    return match_dict


def generate_cp(tree):
    cp_tree = {c: p for p in tree.keys() for c in tree[p]}
    return cp_tree # child: parent


def process_phylowgs_output(summ_file, muts_file, mutass_file):
    with gzip.open(summ_file, "r") as f:
        #with open(summ_file, "r") as f:
        j_summ = json.loads(f.read().decode('utf-8'))
        best_tree = None
        best_tree_llh = -np.infty
        for tree in j_summ['trees'].keys():
            if j_summ['trees'][tree]['llh'] > best_tree_llh:
                best_tree_llh = j_summ['trees'][tree]['llh']
                best_tree = tree
        #print('Total number of trees: ', int(tree)+1, 'Index of best tree: ', best_tree)
    with ZipFile(mutass_file, 'r') as zip:
        with zip.open(best_tree + ".json", ) as g:
            tree_detail = json.load(g)
    mut_assignments = tree_detail['mut_assignments']
    with gzip.open(muts_file, "r") as k:
        j_muts = json.loads(k.read().decode('utf-8'))
    final_tree_cp = {}
    for p, c_list in j_summ['trees'][best_tree]['structure'].items():
        p_list = ['normal'] if p == '0' else [j_muts['ssms'][m]['name'] for m in tree_detail['mut_assignments'][p]['ssms']]
        p_final = tuple(sorted(p_list))
        p_list_idx = [] if p == '0' else [m for m in tree_detail['mut_assignments'][p]['ssms']]
        p_idx_final = tuple(sorted(p_list))
        for c in c_list:
            c_final = tuple(sorted([j_muts['ssms'][m]['name'] for m in tree_detail['mut_assignments'][str(c)]['ssms']]))
            final_tree_cp[c_final] = p_final
    tree_structure_orig = j_summ['trees'][best_tree]['structure']
    tree_structure = {int(i): tree_structure_orig[i] for i in tree_structure_orig.keys()}
    node_dict = {}
    node_dict_name = {}
    node_dict_re = {}
    for p in tree_detail['mut_assignments'].keys():
        node_dict_name[int(p)] = list(sorted([j_muts['ssms'][m]['name'] for m in tree_detail['mut_assignments'][p]['ssms']]))
        node_dict[int(p)] = list(sorted([m for m in tree_detail['mut_assignments'][p]['ssms']]))
        node_dict_re[tuple(sorted([m for m in tree_detail['mut_assignments'][p]['ssms']]))] = int(p)
    prev_mat = []
    population_dict = j_summ['trees'][best_tree]['populations']
    #print(population_dict, tree_structure)
    clonal_freq = {}
    vaf_frac = {}
    for node, populations in population_dict.items():
        clonal_freq[int(node)] = []
        vaf_frac[int(node)] = []
        for sample_idx in range(len(populations['cellular_prevalence'])):
            prev_frac = populations['cellular_prevalence'][sample_idx] - sum([population_dict[str(child)]['cellular_prevalence'][sample_idx] for child in tree_structure[int(node)]] if int(node) in tree_structure.keys() else [0])
            clonal_freq[int(node)] += [prev_frac]
            vaf_frac[int(node)] += [populations['cellular_prevalence'][sample_idx]]
            prev_mat.append({'fraction': prev_frac, 'sample': sample_idx, 'clone': node})
        clonal_freq[int(node)] = [clonal_freq[int(node)]]
        vaf_frac[int(node)] = [vaf_frac[int(node)]]
    return tree_structure, node_dict, node_dict_name, node_dict_re, final_tree_cp, prev_mat, clonal_freq, vaf_frac


def find_root(tree_structure):
    children_set = []
    for parent, children in tree_structure.items():
        children_set += children
    for parent in tree_structure.keys():
        if int(parent) not in children_set:
            return str(parent)

def node_dict_name2node_dict(node_dict_name, idx2name):
    node_dict = {}
    name2idx = {value: key for key, value in idx2name.items()}
    for node, muts_idx in node_dict_name.items():
        node_dict[node] = []
        for mut_idx in muts_idx:
            node_dict[node].append(name2idx[mut_idx])
    return node_dict


if __name__ == '__main__':
    simulation_num = 1
    directory = f"/data/liquid_biopsy/simulations/mask/multi-tumor/boot_simulation_7_800_3000_50_50_50/2_2/Simulation_{simulation_num}/1_2/common"
    chain_num = 5
    result = directory + f"/result_{chain_num}"
    summ_file = result + '.summ.json.gz'
    muts_file = result + '.muts.json.gz'
    mutass_file = result + '.mutass.zip'

    tree_structure, node_dict, node_dict_name, node_dict_re, final_tree_cp, prev_mat, clonal_freq, vaf_frac = process_phylowgs_output(summ_file, muts_file, mutass_file)
    g = render_tumor_tree(tree_structure, node_dict_name)
    g.render(filename=directory + f"/phylowgs_{simulation_num}_inter")