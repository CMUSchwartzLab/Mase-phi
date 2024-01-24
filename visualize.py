from graphviz import Digraph
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd


def collapse_nodes(U, C, E, A, W, threshold=0.0, only_leaf=False):
    print("Loading collapse nodes")
    # generate the tree
    tree = ModifyTree(E)
    if not only_leaf:
        # collapse the branches with 0 length
        branch_remove_idx = []
        for i in range(tree.N-1, -1, -1):
            for j in range(tree.N-1, -1, -1):
                if int(E[i, j]) == 1 and sum(W[j, :]) == 0:
                    branch_remove_idx.append(j)
        for node in branch_remove_idx:
            target = tree.cp_tree[node]
            U[:, target] += U[:, node]
            tree.delete_node(node)

        # collapse the nodes with 0 frequency
        freq_remove_idx = []
        freq_leaf_remove_idx = []
        for i in range(tree.N-1, -1, -1):
            if i in branch_remove_idx:
                continue
            if tree.is_root(i):
                continue
            if np.mean(U[:, i]) <= threshold:
                if tree.num_children(i) == 1:
                    freq_remove_idx.append(i)
                elif tree.is_leaf(i):
                    freq_leaf_remove_idx.append(i)
        print(freq_remove_idx)
        for node in freq_remove_idx:
            target = tree.tree[node][0]
            parent = tree.cp_tree[node]
            tree.delete_node(node)
            W[target, :] += W[node, :]
        for node in freq_leaf_remove_idx:
            tree.delete_node(node)
    else:
        # collapse the branches with 0 length and the child of the branch doesn't belong to leaf nodes
        branch_remove_idx = []
        print("Collapsing branch with length 0")
        for i in range(tree.N - 1, -1, -1):
            for j in range(tree.N - 1, -1, -1):
                if int(E[i, j]) == 1 and sum(W[j, :]) == 0 and not tree.is_leaf(j):
                    branch_remove_idx.append(j)
        for node in branch_remove_idx:
            target = tree.cp_tree[node]
            U[:, target] += U[:, node]
            tree.delete_node(node)

        # collapse the leaf nodes with 0 frequency
        freq_remove_idx = []
        freq_leaf_remove_idx = []
        print('Collapsing leaf nodes with frequency 0')
        for i in range(tree.N - 1, -1, -1):
            if i in branch_remove_idx:
                continue
            if np.mean(U[:, i]) <= threshold and tree.is_leaf(i):
                freq_leaf_remove_idx.append(i)
        for node in freq_leaf_remove_idx:
            tree.delete_node(node)
        for i in range(tree.N - 1, -1, -1):
            if tree.num_children(i) == 1:
                freq_remove_idx.append(i)
        for node in freq_remove_idx:
            target = tree.tree[node][0]
            parent = tree.cp_tree[node]
            tree.delete_node(node)
            W[target, :] += W[node, :]

    # delete those nodes
    remove_idx = branch_remove_idx + freq_remove_idx + freq_leaf_remove_idx
    print('Nodes ', remove_idx, 'will be collapsed.')
    U_new = np.delete(U, remove_idx, axis=1)
    C_new = np.delete(C, remove_idx, axis=0)
    A_new = np.delete(A, remove_idx, axis=0)
    A_new = np.delete(A_new, remove_idx, axis=1)
    E_new = np.delete(tree.E, remove_idx, axis=0)
    E_new = np.delete(E_new, remove_idx, axis=1)
    W_new = np.delete(W, remove_idx, axis=0)
    print("collapse", U_new.shape, C_new.shape)
    return U_new, C_new, E_new, A_new, W_new

class ModifyTree:
    def __init__(self, E):
        self.cp_tree = {}
        self.tree = {}
        self.E = E
        self.N = len(E)
        for i in range(self.N - 1, -1, -1):
            for j in range(self.N - 1, -1, -1):
                if int(E[i, j]) == 1:
                    self.cp_tree[j] = i
                    if i not in self.tree.keys():
                        self.tree[i] = [j]
                    else:
                        self.tree[i].append(j)

    def delete_node(self, idx):
        if self.is_root(idx):
            if self.num_children(idx) > 1:
                raise('Cannot delete root node with more than one child!')
            child = self.tree[idx][0]
            del self.cp_tree[child]
            del self.tree[idx]
        elif self.is_leaf(idx):
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
                self.E[parent, child] = 1
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


def add_prefix_tree(mutation):
    mutation_prefixed = {}
    for node, mut_list in mutation.items():
        mutation_prefixed[node] = ["mut_" + str(mut) for mut in mut_list]
    return mutation_prefixed

def render_tumor_tree(tree_structure, node_dict):
    w = Digraph(format='png')
    #w._repr_svg_(True)
    edge_idx = 0
    root = root_searching(tree_structure)
    for p, c_list in tree_structure.items():
        p_list = ['normal'] if p == root else node_dict[p]
        p_final = tuple(set(p_list))
        ne = len(p_final)
        for c in c_list:
            edge_idx += 1
            print(c)
            if c not in node_dict:
                c_type = str(c)
            else:
                c_type = c
            c_final = tuple(set(node_dict[c_type]))
            ne_ = len(c_final)
            if ne >= 10:
                p_str = str(p) + ' ' + ' '.join((str(e) for e in p_final[:ne//3])) + '\n' + ' '.join((str(e) for e in p_final[ne//3:ne//3*2])) + '\n' + ' '.join((str(e) for e in p_final[ne//3*2:]))
                if ne_ >= 10:
                    c_str = str(c) + ' ' + ' '.join((str(e) for e in c_final[:ne_//3])) + '\n' + ' '.join((str(e) for e in c_final[ne_//3:ne_//3*2])) + '\n' + ' '.join((str(e) for e in c_final[ne_//3*2:]))
                else:
                    c_str = str(c) + ' ' + ' '.join((str(e) for e in c_final))
            else:
                p_str = str(p) + ' ' + ' '.join((str(e) for e in p_final))
                if ne_ >= 10:
                    c_str = str(c) + ' ' + ' '.join((str(e) for e in c_final[:ne_//3])) + '\n' + ' '.join((str(e) for e in c_final[ne_//3:ne_//3*2])) + '\n' + ' '.join((str(e) for e in c_final[ne_//3*2:]))
                else:
                    c_str = str(c) + ' ' + ' '.join((str(e) for e in c_final))
            w.edge(p_str, c_str, "b" + str(edge_idx))
    #display(w)
    return w

def df2dict(df):
    idx2name = {}
    for i in range(len(df)):
        if isinstance(df.loc[i]["Gene"], str):
            name = df.loc[i]["Gene"]
        else:
            name = str(df.loc[i]["Chromosome"]) + '_' + str(df.loc[i]["Genomic Position"])
        idx2name.setdefault(i, name)
    return idx2name

def W2node_dict(W_node, idx2name=None):
    node_dict = {}
    N, m = W_node.shape
    for i in range(N):
        node_dict.setdefault(i, [])
        for j in range(m):
            if W_node[i, j] == 1:
                if idx2name is not None:
                    node_dict[i].append(idx2name[j])
                else:
                    node_dict[i].append(j)
    return node_dict


def generate_tree(cp_tree):
    tree = {}
    for child, parent in cp_tree.items():
        if parent in tree.keys():
            tree[parent].append(child)
        else:
            tree[parent] = [child]
    return tree

def generate_cp(tree):
    return {c: p  for p in tree.keys() for c in tree[p]} # child: parent

def root_searching(tree):  # O(depth of tree) <= O(k)
    tree_cp = generate_cp(tree)
    start_node = list(tree_cp.keys())[0]
    iter_count = 0
    while True:
        iter_count += 1
        start_node = tree_cp[start_node]
        if start_node not in tree_cp.keys():
            break
        if iter_count >= 100:
            print("The directed tree exists self-loop.")
            return None
    return start_node
def E2tree(E):
    tree = {}
    N = E.shape[0]
    for i in range(N):
        for j in range(N):
            if E[i, j] == 1:
                tree.setdefault(i, [])
                tree[i].append(j)
    return tree

def analyze_tree_distribution(tree_distribution, directory, patient_num, type, fig=False):
    for idx in range(len(tree_distribution['freq'])):
        tree_structure = tree_distribution['tree_structure'][idx]
        cp_tree = tree_distribution['cp_tree'][idx]
        node_dict = tree_distribution['node_dict'][idx]
        node_dict_name = tree_distribution['node_dict_name'][idx]
        freq = tree_distribution['freq'][idx]
        clonal_freq = tree_distribution['clonal_freq'][idx]

        ## clonal prevalence
        prev_mat = []
        for node, freqs in clonal_freq.items():
            for freq_sample in freqs:
                prev_blood = freq_sample[0]
                prev_tissue = freq_sample[1]
                prev_mat.append({'fraction': prev_blood, 'sample': 'blood', 'clone': node})
                prev_mat.append({'fraction': prev_tissue, 'sample': 'tissue', 'clone': node})
        df_prev = pd.DataFrame(prev_mat)
        # visualization
        #print(tree_structure, node_dict_name)
        if fig:
            g = render_tumor_tree(tree_structure, node_dict_name)
            g.render(filename=directory / f'{patient_num}_tree_dist{idx}_{type}')

            plt.figure()
            sns.barplot(data=df_prev, x='clone',y='fraction', hue='sample')
            plt.title(str(patient_num) + '_' + str(idx) + '_freq' + str(freq))
            plt.savefig(directory / f'{patient_num}_freq_dist{idx}_{type}.png')