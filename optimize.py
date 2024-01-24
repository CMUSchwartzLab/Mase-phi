import numpy as np
import matplotlib.pyplot as plt
import graphviz
from graphviz import Digraph
from scipy.special import comb, perm
import os
import sys
from itertools import combinations, permutations
from collections import deque
from functools import reduce
import copy
from ete3 import Tree
import re
import random
from pathlib import Path
from scipy.stats import norm
from itertools import combinations
import copy
import json
import os






def wald_test(freq_hat_1, freq_hat_2, correction_rate, relation='ancestor', depth=100, alpha=0.05):
    '''
    return True if reject
    '''
    #print(freq_hat_1, freq_hat_2)
    assert relation in ['ancestor', 'descendant', 'same']
    if relation == 'ancestor':
        W = (freq_hat_2 - freq_hat_1) / np.sqrt((freq_hat_1 * (1 - freq_hat_1) + freq_hat_2 * (1 - freq_hat_2)) / depth)
        z = norm.ppf(alpha / correction_rate)
    elif relation == 'descendant':
        W = (freq_hat_1 - freq_hat_2) / np.sqrt((freq_hat_1 * (1 - freq_hat_1) + freq_hat_2 * (1 - freq_hat_2)) / depth)
        z = norm.ppf(alpha / correction_rate)
    elif relation == 'same':
        W = np.abs((freq_hat_1 - freq_hat_2) / np.sqrt(
            (freq_hat_1 * (1 - freq_hat_1) + freq_hat_2 * (1 - freq_hat_2)) / depth))
        z = norm.ppf(alpha / correction_rate / 2)
    #print(W)
    if W > - z:
        return True, W, - z
    else:
        return False, W, z


def simulate_freq(tree, k, alpha=0.3, beta=0.3):
    ### need to rewrite
    freq_true = np.random.beta(alpha, beta, k)
    freq_obs = np.random.dirichlet(freq_true)
    freq_sum = {}
    bfs_order = bfs_structure(tree)
    for node_idx in range(len(bfs_order) - 1, -1, -1):
        node = bfs_order[node_idx]
        if node not in simulate_tree_template.keys():
            freq_sum[node] = freq_obs[node]
        else:
            freq_sum[node] = sum([freq_sum[child_node] for child_node in tree[node]]) + freq_obs[node]
    return freq_sum


# test block
def calculate_tree_entropy(tree_freq_list, rejected_tree_indices):
    tree_freq = np.array(tree_freq_list)
    if len(rejected_tree_indices) != 0:
        tree_freq[rejected_tree_indices] = 0
    tree_freq = tree_freq[tree_freq != 0]
    tree_freq = (tree_freq).astype('float')/np.sum(tree_freq)
    entropy = - np.sum(tree_freq * np.log(tree_freq))
    return entropy


def calculate_square_sum(tree_freq_list, rejected_tree_indices):
    tree_freq = np.array(tree_freq_list)
    if len(rejected_tree_indices) != 0:
        tree_freq[rejected_tree_indices] = 0
    tree_freq = tree_freq[tree_freq != 0]
    tree_freq = (tree_freq).astype('float')/np.sum(tree_freq)
    sq_sum = 1 - np.sum(np.square(tree_freq))
    return sq_sum


# utility functions

def mut2node(node_dict):
    mut2node_dict = {}
    for node, mut_list in node_dict.items():
        for mut in mut_list:
            mut2node_dict[mut] = int(node)
    return mut2node_dict


def bfs_structure(tree):  # O(k)
    order = []
    root = find_root(tree)
    q = deque([root])
    while len(q) != 0:
        node = q.popleft()
        order.append(node)
        if str(node) in tree.keys():
            for child in tree[str(node)]:
                q.append(child)
        elif node in tree.keys():
            for child in tree[node]:
                q.append(child)
    return order

def bfs(root, tree):  #O(k)
    order = []
    q = deque([root])
    while len(q) != 0:
        node = q.popleft()
        order.append(node)
        if str(node) in tree.keys():
            for child in tree[str(node)]:
                q.append(child)
        elif node in tree.keys():
            for child in tree[node]:
                q.append(child)
    return order

def find_root(tree):
    non_root = []
    for item in tree.values():
        non_root += list(item)
    for node in tree.keys():
        if int(node) not in non_root:
            return int(node)


def ancestor2descendant(tree):
    order = bfs_structure(tree)
    a2d = np.zeros((len(order), len(order)))
    for node in order[::-1]:
        if node in tree.keys():
            for child in tree[node]:
                a2d[int(node)][int(child)] = 1
                a2d[int(node)] += a2d[int(child)]
    return a2d


def generate_cp(tree):
    return {int(c): int(p) for p in tree.keys() for c in tree[p]}  # child: parent


def generate_tree(cp_tree):
    tree = {}
    for child, parent in cp_tree.items():
        if parent in tree.keys():
            tree[parent].append(child)
        else:
            tree[parent] = [child]
    return tree


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


### count ancestor-descendant relationships of all pairs of mutations
def find_all_ancestors(tree, node_dict):
    root = root_searching(tree)
    cp_tree = generate_cp(tree)
    order = bfs(root, tree)
    ancestors_dict = {}
    ancestors_node_dict = {}
    for node in order:
        if node is root:
            ancestors_node_dict.setdefault(root, [])
            continue
        parent = cp_tree[node]
        ancestors_node_dict.setdefault(node, ancestors_node_dict[parent] + [parent])  # inherit the ancestors of the parent
        mut_anc = []
        for n in ancestors_node_dict[node]:
            if n != root:
                mut_anc += node_dict[n]
        for mut in node_dict[node]:
            ancestors_dict.setdefault(mut, mut_anc)
    return ancestors_dict, ancestors_node_dict


def create_ancestor_descendant_matrix(tree, node_dict, gene2idx):
    ancestors_dict, ancestors_node_dict = find_all_ancestors(tree, node_dict)
    num_muts = len(ancestors_dict.keys())
    anc_des_matrix = np.zeros((num_muts, num_muts))
    print(ancestors_dict)
    for mut, ancestors in ancestors_dict.items():
        if len(ancestors) >= 1:
            index = (np.array([gene2idx[anc] for anc in ancestors]), np.repeat(gene2idx[mut], len(ancestors)))
            anc_des_matrix[index] += 1
    return anc_des_matrix


def create_same_clone_matrix(tree, node_dict, gene2idx):
    root = root_searching(tree)
    order = bfs(root, tree)
    sam_clo_matrix = np.zeros((len(gene2idx.keys()), len(gene2idx.keys())))

    for node in order:
        if node != root:
            muts = [gene2idx[mut] for mut in node_dict[node]]
            if len(muts) >= 2:
                indices = list(permutations(muts, r=2))
                for idx in indices:
                    sam_clo_matrix[idx] = 1
    return sam_clo_matrix


def calculate_entropy(matrix):
    return -np.sum(matrix * np.log(matrix, out=np.zeros_like(matrix), where=(matrix != 0)), axis=0)

def calculate_relation_matrix(tree_list, node_list, gene2idx,  tree_freq_list=None):
    num_tree = len(tree_list)
    if tree_freq_list is None:
        tree_freq_list = np.ones(num_tree)
    sum_tree_freq = np.sum(tree_freq_list)
    for i in range(num_tree):
        tree = tree_list[i]
        node_dict = node_list[i]
        anc_des_matrix = create_ancestor_descendant_matrix(tree, node_dict, gene2idx)
        sam_clo_matrix = create_same_clone_matrix(tree, node_dict, gene2idx)
        if i == 0:
            anc_des_matrix_sum = np.zeros(anc_des_matrix.shape)
            sam_clo_matrix_sum = np.zeros(anc_des_matrix.shape)
        anc_des_matrix_sum += anc_des_matrix * tree_freq_list[i]
        sam_clo_matrix_sum += sam_clo_matrix * tree_freq_list[i]
    des_anc_matrix_sum = anc_des_matrix_sum.T
    no_rel_matrix_sum = np.ones(anc_des_matrix.shape) * sum_tree_freq - anc_des_matrix_sum - des_anc_matrix_sum - sam_clo_matrix_sum
    full_matrix_sum = np.concatenate((anc_des_matrix_sum[np.newaxis, :] / sum_tree_freq,
                                      des_anc_matrix_sum[np.newaxis, :] / sum_tree_freq,
                                      sam_clo_matrix_sum[np.newaxis, :] / sum_tree_freq,
                                      no_rel_matrix_sum[np.newaxis, :] / sum_tree_freq), axis=0)
    return full_matrix_sum

def calculate_pairwise_uncertainty(full_matrix_sum, method='entropy'):
    if method == 'entropy':
        return calculate_entropy(full_matrix_sum)


def tree2E(tree, k):
    E = np.zeros((k, k))
    for parent, children_list in tree.items():
        p = int(parent)
        for child in children_list:
            c = int(child)
            E[p, c] = 1
    return E


def tree2E_list(tree_list, k_list):
    E_list = []
    for idx, (tree, k) in enumerate(zip(tree_list, k_list)):
        E = tree2E(tree, k)
        E_list.append(E)
    return E_list


def create_M(node_dict, gene2idx, N):
    num_node = len(node_dict.keys())+1
    M = np.zeros((num_node, N))
    for node, mut_list in node_dict.items():
        n = int(node)
        for mut in mut_list:
            m = gene2idx[mut]
            M[n, m] = 1
    return M

def create_M_list(node_list, gene2idx, N):
    M_list = []
    for idx, node_dict in enumerate(node_list):
        num_node = len(node_dict.keys())+1
        M = create_M(node_dict, gene2idx, N)
        M_list.append(M)
    return M_list

def create_k_list(node_list):
    k_list = []
    for node_dict in node_list:
        num_node = len(node_dict.keys()) + 1
        k_list.append(num_node)
    return k_list


def calculate_F_Fhat(clonal_freq, tree_dict, sample_idx=0):
    node_num = len(clonal_freq.keys())
    F = np.zeros((node_num))
    F_hat = np.zeros((node_num))
    for node in clonal_freq.keys():
        F[node] = clonal_freq[node][sample_idx]
        prev_blood = clonal_freq[node][sample_idx] - sum(
            [clonal_freq[child][sample_idx] for child in tree_dict[node]] if node in tree_dict.keys() else [0])
        F_hat[node] = prev_blood
    return F, F_hat

def calculate_F_Fhat_from_pcr(clonal_freq_actual, tree_dict, total_node_num):
    F = np.zeros((total_node_num))
    F_hat = np.zeros((total_node_num))
    for node in clonal_freq_actual.keys():
        F[node] = clonal_freq_actual[node]
        if node in tree_dict.keys():
            prev_blood = clonal_freq_actual[node]
            for child in tree_dict[node]:
                if child in clonal_freq_actual:
                    prev_blood -= clonal_freq_actual[child]
                else:
                    prev_blood = 0
                    break
        else:
            prev_blood = clonal_freq_actual[node]
        F_hat[node] = prev_blood
    return F, F_hat

def create_F_F_hat_list(clonal_freq_list, tree_list, sample_idx=0):
    F_list, F_hat_list = [], []
    for idx, (clonal_freq, tree_dict) in enumerate(zip(clonal_freq_list, tree_list)):
        F, F_hat = calculate_F_Fhat(clonal_freq, tree_dict, sample_idx)
        F_list.append(F)
        F_hat_list.append(F_hat)
    return F_list, F_hat_list

def create_F_F_hat_list_from_pcr(clonal_freq_list_actual, tree_list, clonal_freq_list_origin):
    F_list, F_hat_list = [], []
    for idx, (clonal_freq, tree_dict) in enumerate(zip(clonal_freq_list_actual, tree_list)):
        total_node_num = len(clonal_freq_list_origin[idx].keys())
        F, F_hat = calculate_F_Fhat_from_pcr(clonal_freq, tree_dict, total_node_num)
        F_list.append(F)
        F_hat_list.append(F_hat)
    return F_list, F_hat_list