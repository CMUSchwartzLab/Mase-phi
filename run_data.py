from optimize import *
from optimize_fraction import *
import pandas as pd
from zipfile import ZipFile
import json
import gzip
import pickle


directory=Path('/data/liquid_biopsy/MOONSHOT2/bootstrap')
# apply the real data to method
# step1. generate the input format
'''     gene_list: a list of total sequenced genes
        tree_list: a list of dict for tree structure
        node_list: a list of dict for tree node - mutation_list
        clonal_freq_list: a list of dicts of the frequency for each node
        tree_freq_list: a list of the frequencies of each possible tree
        n_markers: number for the gene markers
'''
patient_num = 256
method = 'phylowgs'
tree_distribution_file = directory / f'{patient_num}_bootstrap_{method}/{method}_bootstrap_aggregation.pkl'
with open(directory / tree_distribution_file, 'rb') as f:
    tree_distribution = pickle.load(f)

gene_list = []
gene2idx = {}
inter = pd.read_excel(directory / f"xlsx/Patient_{patient_num}_bootstrap.xlsx", sheet_name='common_blood_tissue_no_germline', index_col=0)
inter = inter[inter["Allele Frequency_x"] < 0.9]
inter = inter[inter["Allele Frequency_y"] < 0.9]
calls = inter
gene2idx = {'s' + str(i): i for i in range(len(inter))}
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



tree_list, node_list, clonal_freq_list, tree_freq_list = tree_distribution['tree_structure'], tree_distribution['node_dict'],tree_distribution['vaf_frac'],tree_distribution['freq']

#scrub node_list
node_list_scrub = []
for node_dict in node_list:
    temp = {}
    for key, values in node_dict.items():
        temp.setdefault(int(key), values)
    node_list_scrub.append(temp)

clonal_freq_list_scrub = []
for clonal_freq_dict in clonal_freq_list:
    temp = {}
    for key, values in clonal_freq_dict.items():
        temp.setdefault(int(key), values[0])
    clonal_freq_list_scrub.append(temp)

# tree_list_sub = [tree_list[i] for i in range(5) ]
# node_list_sub = [node_list_scrub[i] for i in range(5) ]
# clonal_freq_list_sub = [clonal_freq_list_scrub[i] for i in range(5) ]
# tree_freq_list_sub = [tree_freq_list[i] for i in range(5)]


selected_markers1_genename_ordered = []
obj1_ordered = []

for n_markers in range(1, len(gene_name_list) + 1):
    #selected_markers1, obj = select_markers_fractions_weighted_single(gene_list, n_markers, tree_list, node_list_scrub, clonal_freq_list_scrub, gene2idx, tree_freq_list, 0)
    selected_markers1, obj = select_markers_fractions_weighted_overall(gene_list, n_markers, tree_list, node_list_scrub, clonal_freq_list_scrub, gene2idx, tree_freq_list)
    selected_markers1_genename = [gene_name_list[int(i[1:])] for i in selected_markers1]
    obj1_ordered.append(obj)
    if len(selected_markers1_genename) == 1:
        selected_markers1_genename_ordered.append(selected_markers1_genename[0])
    else:
        diff_set = set(selected_markers1_genename).difference(set(selected_markers1_genename_ordered))
        selected_markers1_genename_ordered.append(list(diff_set)[0])
print(selected_markers1_genename_ordered, obj1_ordered)

position1 = list(range(len(obj1_ordered)))

plt.figure(figsize=(8, 5))
plt.plot(position1, obj1_ordered, 'o-', label='tracing-fractions')
plt.xticks(position1, selected_markers1_genename_ordered, rotation=30)
plt.legend()
plt.savefig(f'{patient_num}_tracing_subclones.eps', format='eps')
plt.show()

read_depth=1500
lam1=1  # lam1 is tree fractions
lam2=0  # lam2 is tree distributions
selected_markers2_genename_ordered = []
obj2_ordered = []
for n_markers in range(1, len(gene_name_list) + 1):
    selected_markers2, obj_frac, obj_struct = select_markers_tree_gp(gene_list, n_markers, tree_list, node_list_scrub, clonal_freq_list_scrub, gene2idx, tree_freq_list, read_depth=read_depth, lam1=lam1, lam2=lam2)
    selected_markers2_genename = [gene_name_list[int(i[1:])] for i in selected_markers2]
    obj2_ordered.append((obj_frac, obj_struct))
    if len(selected_markers2_genename) == 1:
        selected_markers2_genename_ordered.append(selected_markers2_genename[0])
    else:
        selected_markers2_genename_ordered.append(
            list(set(selected_markers2_genename).difference(set(selected_markers2_genename_ordered)))[0])

print(selected_markers2_genename_ordered, obj2_ordered)

import seaborn as sns

obj2_frac_ordered = [obj2_ordered[i][0] for i in range(len(obj2_ordered))]
obj2_struct_ordered = [obj2_ordered[i][1] for i in range(len(obj2_ordered))]
position2 = list(range(len(obj2_ordered)))



plt.figure(figsize=(8, 5))
plt.plot(position2, obj2_frac_ordered, 'o-', color='tab:orange',label='trees-fractions')
plt.xticks(position2, selected_markers2_genename_ordered, rotation=30)
plt.legend()
plt.savefig(f'{patient_num}_trees_fractions_{lam1}_{lam2}.eps', format='eps')
plt.show()


plt.figure(figsize=(8, 5))
plt.plot(position2, obj2_struct_ordered, 'o-', color='tab:green',label='trees-structure')
plt.xticks(position2, selected_markers2_genename_ordered, rotation=30)
plt.legend()
plt.savefig(f'{patient_num}_trees_structures_{lam1}_{lam2}.eps', format='eps')
plt.show()

read_depth=1500
lam1=0  # lam1 is tree fractions
lam2=1  # lam2 is tree distributions
selected_markers2_genename_ordered = []
obj2_ordered = []
for n_markers in range(1, len(gene_name_list) + 1):
    selected_markers2, obj_frac, obj_struct = select_markers_tree_gp(gene_list, n_markers, tree_list, node_list_scrub, clonal_freq_list_scrub, gene2idx, tree_freq_list, read_depth=read_depth, lam1=lam1, lam2=lam2)
    selected_markers2_genename = [gene_name_list[int(i[1:])] for i in selected_markers2]
    obj2_ordered.append((obj_frac, obj_struct))
    if len(selected_markers2_genename) == 1:
        selected_markers2_genename_ordered.append(selected_markers2_genename[0])
    else:
        selected_markers2_genename_ordered.append(
            list(set(selected_markers2_genename).difference(set(selected_markers2_genename_ordered)))[0])

print(selected_markers2_genename_ordered, obj2_ordered)

import seaborn as sns

obj2_frac_ordered = [obj2_ordered[i][0] for i in range(len(obj2_ordered))]
obj2_struct_ordered = [obj2_ordered[i][1] for i in range(len(obj2_ordered))]
position2 = list(range(len(obj2_ordered)))



plt.figure(figsize=(8, 5))
plt.plot(position2, obj2_frac_ordered, 'o-', color='tab:orange',label='trees-fractions')
plt.xticks(position2, selected_markers2_genename_ordered, rotation=30)
plt.legend()
plt.savefig(f'{patient_num}_trees_fractions_{lam1}_{lam2}_{read_depth}.eps', format='eps')
plt.show()


plt.figure(figsize=(8, 5))
plt.plot(position2, obj2_struct_ordered, 'o-', color='tab:green',label='trees-structure')
plt.xticks(position2, selected_markers2_genename_ordered, rotation=30)
plt.legend()
plt.savefig(f'{patient_num}_trees_structures_{lam1}_{lam2}_{read_depth}.eps', format='eps')
plt.show()

lam1=1  # lam1 is tree fractions
lam2=0  # lam2 is tree distributions
selected_markers2_genename_ordered = []
obj2_ordered = []
for n_markers in range(1, len(gene_name_list) + 1):
    selected_markers2, obj_frac, obj_struct = select_markers_tree_gp(gene_list, n_markers, tree_list, node_list_scrub, clonal_freq_list_scrub, gene2idx, tree_freq_list, read_depth=read_depth, lam1=lam1, lam2=lam2)
    selected_markers2_genename = [gene_name_list[int(i[1:])] for i in selected_markers2]
    obj2_ordered.append((obj_frac, obj_struct))
    if len(selected_markers2_genename) == 1:
        selected_markers2_genename_ordered.append(selected_markers2_genename[0])
    else:
        selected_markers2_genename_ordered.append(
            list(set(selected_markers2_genename).difference(set(selected_markers2_genename_ordered)))[0])

print(selected_markers2_genename_ordered, obj2_ordered)

import seaborn as sns

obj2_frac_ordered = [obj2_ordered[i][0] for i in range(len(obj2_ordered))]
obj2_struct_ordered = [obj2_ordered[i][1] for i in range(len(obj2_ordered))]
position2 = list(range(len(obj2_ordered)))



plt.figure(figsize=(8, 5))
plt.plot(position2, obj2_frac_ordered, 'o-', color='tab:orange',label='trees-fractions')
plt.xticks(position2, selected_markers2_genename_ordered, rotation=30)
plt.legend()
plt.savefig(f'{patient_num}_trees_fractions_{lam1}_{lam2}_{read_depth}.eps', format='eps')
plt.show()


plt.figure(figsize=(8, 5))
plt.plot(position2, obj2_struct_ordered, 'o-', color='tab:green',label='trees-structure')
plt.xticks(position2, selected_markers2_genename_ordered, rotation=30)
plt.legend()
plt.savefig(f'{patient_num}_trees_structures_{lam1}_{lam2}_{read_depth}.eps', format='eps')
plt.show()
