import matplotlib.pyplot as plt
import numpy as np
from itertools import combinations, permutations
from collections import deque

import pandas as pd
from scipy.stats import norm, chi2

from pathlib import Path
import pickle
from analyze import *
import math
from scipy.integrate import quad
import random

def adjust_tree_distribution_struct(tree_list, node_dict_list, read_depth, ddpcr_marker_counts, marker_idx2gene, alpha):
    accepted_tree_indices = []
    for tree_idx in range(len(tree_list)):
        tree_structure = tree_list[tree_idx]
        node_dict = node_dict_list[tree_idx]
        bool_reject = test_single_tree_struct(tree_structure, node_dict, read_depth, ddpcr_marker_counts, marker_idx2gene, alpha)
        if not bool_reject:
            accepted_tree_indices.append(tree_idx)
    return accepted_tree_indices

def adjust_tree_distribution_frac(clonal_freq_list,node_dict_list, read_depth, ddpcr_marker_counts, marker_idx2gene, alpha):
    accepted_tree_indices = []
    for tree_idx in range(len(clonal_freq_list)):
        clonal_freq_dict = clonal_freq_list[tree_idx]
        node_dict = node_dict_list[tree_idx]
        bool_reject = test_single_tree_frac(clonal_freq_dict, node_dict, read_depth, ddpcr_marker_counts, marker_idx2gene, alpha)
        if not bool_reject:
            accepted_tree_indices.append(tree_idx)
    return accepted_tree_indices

def adjust_tree_distribution_struct_bayesian(tree_list, node_dict_list, tree_freq_list, read_depth, ddpcr_marker_counts, marker_idx2gene):
    updated_tree_freq_list = []
    marker_idx_list = random.choice(list(combinations(list(range(len(ddpcr_marker_counts))), 2)))
    print(marker_idx_list)
    for tree_idx in range(len(tree_list)):
        tree_structure = tree_list[tree_idx]
        node_dict = node_dict_list[tree_idx]
        tree_freq = tree_freq_list[tree_idx]
        updated_tree_freq = update_single_tree_fractions(tree_structure, node_dict, tree_freq, read_depth, ddpcr_marker_counts, marker_idx2gene, marker_idx_list)
        updated_tree_freq_list.append(updated_tree_freq)
    updated_tree_freq_list_np = np.array(updated_tree_freq_list)
    updated_tree_freq_list_np /= np.sum(updated_tree_freq_list_np) / 100
    return updated_tree_freq_list_np.tolist()

def update_single_tree_fractions(tree_structure, node_dict, tree_freq, read_depth_list, ddpcr_marker_counts, marker_idx2gene, marker_idx_list, lower_bound=0, upper_bound=1):
    n_markers = len(ddpcr_marker_counts)
    a2d_matrix = ancestor2descendant(tree_structure)
    mut2node_dict = mut2node(node_dict)
    [marker_idx1, marker_idx2] = marker_idx_list
    node1, node2 = mut2node_dict[marker_idx2gene[marker_idx1]], mut2node_dict[marker_idx2gene[marker_idx2]]
    read_depth_sublist = [int(np.round(read_depth_list[marker_idx1])), int(np.round(read_depth_list[marker_idx2]))]
    read_counts_sublist = [int(np.round(ddpcr_marker_counts[marker_idx1])), int(np.round(ddpcr_marker_counts[marker_idx2]))]
    print(read_depth_sublist, read_counts_sublist)
    if node1 == node2:
        relation = 'same'
    elif a2d_matrix[node1, node2] == 1:
        relation = 'ancestor'
    elif a2d_matrix[node1, node2] == 0 and a2d_matrix[node2, node1] == 1:
        relation = 'descendant'
    else:
        relation = 'null'
    if relation == 'same':
        conditional_prob = outer_integral_single(read_depth_sublist[0], read_depth_sublist[1], read_counts_sublist[0], read_counts_sublist[1], lower_bound, upper_bound)
    else:
        conditional_prob = outer_integral(read_depth_sublist[0], read_depth_sublist[1], read_counts_sublist[0], read_counts_sublist[1], relation, lower_bound, upper_bound)
    updated_tree_freq = tree_freq * conditional_prob
    return updated_tree_freq

def ncr(n, r): ##for python 3.7
    return math.factorial(n) // math.factorial(r) // math.factorial(n-r)
def integrand(f2, f1, d1, d2, r1, r2):
    try:
        # return math.comb(d1, r1) * (f1 ** r1) * ((1 - f1) ** (d1 - r1)) * math.comb(d2, r2) * (f2 ** r2) * ((1 - f2) ** (d2 - r2))
        return math.exp(math.log(math.comb(d1, r1)) + r1 * math.log(f1) + (d1 - r1) * math.log(1 - f1) +
               math.log(math.comb(d2, r2)) + r2 * math.log(f2) + (d2 - r2) * math.log(1 - f2))
    except AttributeError:
        return ncr(d1, r1) * (f1 ** r1) *((1 - f1)**(d1-r1)) * ncr(d2, r2) * (f2 ** r2) * ((1-f2) ** (d2-r2))

def recursive_integral(f1, d1, d2, r1, r2, relation,lower_bound, upper_bound):
    integral, error = quad(integrand, f2_start(f1, relation, lower_bound), f2_end(f1, relation, upper_bound), args=(f1, d1, d2, r1, r2))
    return integral

def outer_integral(d1, d2, r1, r2, relation, lower_bound, upper_bound):
    integral, error = quad(recursive_integral, f1_start(lower_bound), f1_end(upper_bound), args=(d1, d2, r1, r2,relation,lower_bound, upper_bound))
    return integral

def integrand_single(f, d1, d2, r1, r2):
    try:
        return math.comb(d1, r1) * (f ** r1) * ((1 - f) ** (d1 - r1)) * math.comb(d2, r2) * (f ** r2) * ((1 - f) ** (d2 - r2))
    except AttributeError:
        return ncr(d1, r1) * (f ** r1) *((1 - f)**(d1-r1)) * ncr(d2, r2) * (f ** r2) * ((1-f) ** (d2-r2))

def outer_integral_single(d1, d2, r1, r2, lower_bound, upper_bound):
    integral, error = quad(integrand_single, f1_start(lower_bound), f1_end(upper_bound), args=(d1, d2, r1, r2))
    return integral

def f2_start(f1, relation, lower_bound):
    if relation == 'descendant':
        return f1
    elif relation == 'ancestor' or 'null':
        return lower_bound
    elif relation == 'same':
        raise AttributeError
def f2_end(f1, relation, upper_bound):
    if relation == 'descendant':
        return upper_bound
    elif relation == 'ancestor':
        return f1
    elif relation == 'null':
        return 1 - f1
    elif relation == 'same':
        raise AttributeError


def f1_start(lower_bound=0):
    return lower_bound

def f1_end(upper_bound=1):
    return upper_bound


def update_tree_distribution(tree_distribution, accepted_tree_indices):
    updated_tree_distribution = {}
    for key, item_list in tree_distribution.items():
        print(key, len(item_list))
        updated_item_list = [item_list[idx] for idx in accepted_tree_indices]
        updated_tree_distribution[key] = updated_item_list
    return updated_tree_distribution
def update_tree_distribution_bayesian(tree_distribution, update_tree_freq_list):
    updated_tree_distribution = {}
    for key, item_list in tree_distribution.items():
        print(key, len(item_list))
        if key != 'freq':
            updated_tree_distribution[key] = item_list
        else:
            updated_tree_distribution['freq'] = update_tree_freq_list
    return updated_tree_distribution

def mut2node(node_dict):
    mut2node_dict = {}
    for node, mut_list in node_dict.items():
        for mut in mut_list:
            mut2node_dict[mut] = int(node)
    return mut2node_dict

def find_root(tree):
    non_root = []
    for item in tree.values():
        non_root += list(item)
    for node in tree.keys():
        if node not in non_root:
            return node

def bfs_structure(tree):  # O(k)
    order = []
    root = find_root(tree)
    q = deque([root])
    while len(q) != 0:
        node = q.popleft()
        order.append(node)
        if node in tree.keys():
            for child in tree[node]:
                q.append(child)
    return order

def ancestor2descendant(tree):
    order = bfs_structure(tree)
    a2d = np.zeros((len(order), len(order)))
    for node in order[::-1]:
        if node in tree.keys():
            for child in tree[node]:
                a2d[int(node)][int(child)] = 1
                a2d[int(node)] += a2d[int(child)]
    return a2d

def test_single_tree_struct(tree_structure, node_dict, read_depth_list, ddpcr_marker_counts, marker_idx2gene, alpha):
    n_markers = len(ddpcr_marker_counts)
    a2d_matrix = ancestor2descendant(tree_structure)
    mut2node_dict = mut2node(node_dict)
    for (marker_idx1, marker_idx2) in list(combinations(list(range(len(ddpcr_marker_counts))), 2)):
        node1, node2 = mut2node_dict[marker_idx2gene[marker_idx1]], mut2node_dict[marker_idx2gene[marker_idx2]]
        read_depth_sublist = [read_depth_list[marker_idx1], read_depth_list[marker_idx1]]

        freq_hat_1, freq_hat_2 = ddpcr_marker_counts[marker_idx1]/read_depth_sublist[0], ddpcr_marker_counts[marker_idx2]/read_depth_sublist[1]
        if node1 == node2:
            relation = 'same'
        elif a2d_matrix[node1, node2] == 1:
            relation = 'ancestor'
        elif a2d_matrix[node1, node2] == 0 and a2d_matrix[node2, node1] == 1:
            relation = 'descendant'
        else:
            relation = 'null'
        if relation == 'null':
            bool_reject = False
            continue

        # node_correct - 1 is due to the fact that the normal node 0 contains no mutations, thus being omitted
        bool_reject, W, z = wald_test(freq_hat_1, freq_hat_2, n_markers * (n_markers - 1) / 2, relation, read_depth_sublist, alpha=alpha)
        if bool_reject:
            print(relation, bool_reject)
            return bool_reject
        print(relation, bool_reject, W, -z)
    return bool_reject


def test_single_tree_frac(clonal_freq_dict, node_dict, read_depth_list, ddpcr_marker_counts, marker_idx2gene, alpha):
    n_markers = len(ddpcr_marker_counts)
    mut2node_dict = mut2node(node_dict)
    for marker_idx in range(n_markers):
        marker_count = ddpcr_marker_counts[marker_idx]
        node = mut2node_dict[marker_idx2gene[marker_idx]]
        read_depth = read_depth_list[marker_idx]
        clonal_freq = clonal_freq_dict[node][0][0]
        # node_correct - 1 is due to the fact that the normal node 0 contains no mutations, thus being omitted
        bool_reject, T, chi_square = chi_square_test(clonal_freq, read_depth, marker_count, n_markers, alpha=alpha)
        if bool_reject:
            print(bool_reject, T, chi_square)
            return bool_reject
        print(bool_reject, T, chi_square)
    return bool_reject

def chi_square_test(clonal_freq, read_depth, marker_count, correction_rate,alpha=0.05):
    expected_num = [clonal_freq * read_depth, read_depth *(1-clonal_freq)]
    T = (marker_count - expected_num[0])**2 / expected_num[0] + (read_depth - marker_count - expected_num[1])**2 / expected_num[1]
    chi_square = chi2.ppf(1-alpha/correction_rate, 1)
    if T > chi_square:
        return True, T, chi_square
    else:
        return False, T, chi_square

def wald_test(freq_hat_1, freq_hat_2, correction_rate, relation='ancestor', depth=[100,100], alpha=0.05):
    '''
    return True if reject
    '''
    print(freq_hat_1, freq_hat_2)
    assert relation in ['ancestor', 'descendant', 'same']
    if relation == 'ancestor':
        W = (freq_hat_2 - freq_hat_1) / (np.sqrt(
            (freq_hat_1 * (1 - freq_hat_1))/depth[0]) + np.sqrt(freq_hat_2 * (1 - freq_hat_2)/ depth[1]))
        z = norm.ppf(alpha / correction_rate)
    elif relation == 'descendant':
        W = (freq_hat_1 - freq_hat_2) / (np.sqrt(
            (freq_hat_1 * (1 - freq_hat_1))/depth[0]) + np.sqrt(freq_hat_2 * (1 - freq_hat_2) / depth[1]))
        z = norm.ppf(alpha / correction_rate)
    elif relation == 'same':
        W = np.abs((freq_hat_1 - freq_hat_2) / (np.sqrt(
            (freq_hat_1 * (1 - freq_hat_1))/depth[0]) + np.sqrt(freq_hat_2 * (1 - freq_hat_2) / depth[1])))
        z = norm.ppf(alpha / correction_rate / 2)
    #print(W)
    if W > - z:
        return True, W, - z
    else:
        return False, W, z


if __name__ == '__main__':
    ddpcr = [{
        'gene': 'MLH1', 'mut': 100.21, 'WT': 1080.75, 'liquid_biopsy_sample': 'i', 'patient':256
    },
        {
        'gene': 'TP53', 'mut': 204.6, 'WT': 837.65, 'liquid_biopsy_sample': 'i', 'patient':256
        }]
    n_markers = 2
    blood_sample_idx = 0
    df_ddpcr = pd.DataFrame(ddpcr)
    marker_idx2gene = {i: df_ddpcr["gene"][i] for i in range(len(df_ddpcr)) }
    directory = Path('/data/liquid_biopsy/MOONSHOT2/bootstrap')
    patient_num = 256
    method = 'phylowgs'
    tree_distribution_file = directory / f'{patient_num}_bootstrap_{method}/{method}_bootstrap_summary.pkl'
    with open(directory / tree_distribution_file, 'rb') as f:
        tree_distribution_summary = pickle.load(f)

    gene_list = []
    gene2idx = {}
    inter = pd.read_excel(directory / f"xlsx/Patient_{patient_num}_bootstrap.xlsx",
                          sheet_name='common_blood_tissue_no_germline', index_col=0)
    inter = inter[inter["Allele Frequency_x"] < 0.9]
    inter = inter[inter["Allele Frequency_y"] < 0.9]
    calls = inter
    gene2idx = {'s' + str(i): i for i in range(len(inter))}
    genename2idx = {inter["Gene"][i] : i for i in range(len(inter))}
    gene_list = list(gene2idx.keys())
    gene_name_list = []
    gene_count = {}
    for i in range(inter.shape[0]):
        gene = calls.iloc[i, 0]
        if gene in gene_name_list:
            gene_count[gene] += 1
            gene = gene + '_' + str(gene_count[gene])
        else:
            gene_count[gene] = 1
        if not isinstance(gene, str):
            gene = calls.iloc[i, 1] + '_' + str(calls.iloc[i, 2])
        gene_name_list.append(gene)

    tree_list, node_list, node_name_list, tree_freq_list = tree_distribution_summary['tree_structure'], tree_distribution_summary[
        'node_dict'], tree_distribution_summary['node_dict_name'], tree_distribution_summary['freq']

    # scrub node_list
    node_list_scrub = []
    for node_dict in node_list:
        temp = {}
        for key, values in node_dict.items():
            temp.setdefault(int(key), values)
        node_list_scrub.append(temp)

    adjust_algo = 'bayesian'
    ddpcr_marker_counts = list(df_ddpcr["mut"])
    read_depth_list = list(df_ddpcr["mut"] + df_ddpcr["WT"])
    updated_tree_freq_list = adjust_tree_distribution_struct_bayesian(tree_list, node_name_list,tree_freq_list, read_depth_list, ddpcr_marker_counts, marker_idx2gene)
    updated_tree_distribution_summary = update_tree_distribution_bayesian(tree_distribution_summary,
                                                                          updated_tree_freq_list)

    tree_distribution_file_summary_updated = directory / f'{patient_num}_bootstrap_{method}/{method}_bootstrap_summary_updated_{adjust_algo}_{n_markers}_{blood_sample_idx}_bayesian.pkl'
    with open(tree_distribution_file_summary_updated, 'wb') as f:
        pickle.dump(updated_tree_distribution_summary, f)
    print(tree_freq_list, updated_tree_freq_list)

    plt.figure(figsize=(6, 18))
    timepoint = [0, 1]
    for i in range(len(tree_freq_list)):
        if tree_freq_list[i] > updated_tree_freq_list[i]:
            color = "tab:blue"
        else:
            print(i, tree_freq_list[i], updated_tree_freq_list[i])
            color = "tab:orange"
        plt.plot(timepoint, [tree_freq_list[i], updated_tree_freq_list[i]], marker='o', linestyle='-', color=color, ms=5)

    # Add labels and a legend
    plt.xlabel('time points')
    plt.ylabel('tree weights')
    plt.xticks([0,1])
    plt.savefig(f"Patient_{patient_num}_weights_change.eps",bbox_inches='tight')