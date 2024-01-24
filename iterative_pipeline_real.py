from optimize import *
from optimize_fraction import *
import pandas as pd
from zipfile import ZipFile
import json
import gzip
import pickle
from analyze import *
from optimize import *
from adjust_tree_distribution import *
import matplotlib.pyplot as plt


patient_name='CRUK0041'
used_num_blood=0
used_num_tissue=5
type='common'
directory = Path(f'/home/xuecongf/liquid_biopsy/abbosh')
liquid_biopsy_directory = Path("/home/xuecongf/liquid_biopsy/abbosh/DateSample")

#directory = Path(f'Abbosh_et_al/data/bootstrap')
#liquid_biopsy_directory = Path("Abbosh_et_al/DateSample")
n_markers = 2
algo = "struct"
file = directory / f'bootstrap/{patient_name}/{used_num_blood}_{used_num_tissue}/{type}'
bootstrap_num = 100
method = 'phylowgs'
num_chain=5
directory_xlsx = directory / f'{patient_name}_subset.xlsx'
inter = pd.read_excel(directory_xlsx, sheet_name='tissue_no_germline_0', index_col=0)
filtered = inter.copy(deep=True)

calls = filtered
gene2idx = {'s' + str(i): i for i in range(len(inter))}
print(gene2idx)
gene_list = list(gene2idx.keys())
gene_name_list = []
gene_count = {}
for i in range(inter.shape[0]):
    gene = calls["Gene"].loc[i]
    if gene in gene_name_list:
        gene_count[gene] += 1
        gene = gene + '_' + str(gene_count[gene])
    else:
        gene_count[gene] = 1
    if not isinstance(gene, str):
        gene = str(calls["Chromosome"][i]) + '_' + str(calls["Genomic Position"][i])
    gene_name_list.append(gene)

num_marker = len(gene_list)

directory_ddpcr = Path("/home/xuecongf/liquid_biopsy/abbosh/DateSample")
file_ddpcr = directory_ddpcr / f"{patient_name}__DateSample.xlsx"
ddpcr_raw = pd.read_excel(file_ddpcr, sheet_name=None, index_col=0)
date_list = []
for date in ddpcr_raw.keys():
    date_list.append(pd.to_datetime(date, format="%Y-%m-%d"))
date_keys_sorted = [ts.strftime("%Y-%m-%d") for ts in sorted(date_list)]

for order_idx in range(0, len(date_keys_sorted)):
    if order_idx == 0:
        tree_distribution_file_summary = file / f'{method}_bootstrap_summary.pkl'
        with open(tree_distribution_file_summary, 'rb') as f:
            tree_distribution_summary = pickle.load(f)
        tree_list_summary, node_list_summary, node_name_list_summary, tree_freq_list_summary = \
            tree_distribution_summary['tree_structure'], tree_distribution_summary[
                'node_dict'], tree_distribution_summary['node_dict_name'], tree_distribution_summary['freq']
        tree_distribution_file = file / f'{method}_bootstrap_aggregation.pkl'
        with open(file / tree_distribution_file, 'rb') as f:
            tree_distribution = pickle.load(f)
        tree_list, node_list, clonal_freq_list, tree_freq_list = tree_distribution['tree_structure'], tree_distribution[
            'node_dict'], tree_distribution['vaf_frac'], tree_distribution['freq']

    else:
        tree_distribution_file_summary = file / f'{method}_bootstrap_summary_updated_{algo}_{n_markers}_{order_idx-1}_bayesian.pkl'
        with open(tree_distribution_file_summary, 'rb') as f:
            tree_distribution_summary = pickle.load(f)
        tree_list_summary, node_list_summary, node_name_list_summary, tree_freq_list_summary = \
            tree_distribution_summary['tree_structure'], tree_distribution_summary[
                'node_dict'], tree_distribution_summary['node_dict_name'], tree_distribution_summary['freq']
        tree_list, node_list, tree_freq_list = tree_distribution_summary['tree_structure'], tree_distribution_summary[
            'node_dict'], tree_distribution_summary['freq']
        clonal_freq_list = []
        for idx in range(len(tree_distribution_summary['vaf_frac'])):
            clonal_freq_dict = tree_distribution_summary['vaf_frac'][idx]
            clonal_freq_dict_new = {}
            for node, freqs in clonal_freq_dict.items():
                clonal_freq_dict_new[node] = [list(np.array(freqs).mean(axis=0))]
            clonal_freq_list.append(clonal_freq_dict_new)

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

    sample_df = pd.read_excel(liquid_biopsy_directory / f"{patient_name}__DateSample.xlsx", sheet_name=None, index_col=0)
    subset_markers = list(sample_df[list(sample_df.keys())[0]].index)
    subset_list = list(inter[inter.Gene.isin(subset_markers)].index)
    subset_markers_s = list([f"s{i}" for i in subset_list])
    gene2idx_sub = {subset_markers_s[i]: i for i in range(len(subset_markers_s))}

    read_depth=90000
    lam1 = 0  # lam1 is tree fractions
    lam2 = 1  # lam2 is tree distributions

    selected_markers2, obj_frac, obj_struct = select_markers_tree_gp(gene_list, n_markers, tree_list, node_list, clonal_freq_list, gene2idx, tree_freq_list, read_depth=read_depth, lam1=lam1, lam2=lam2, focus_sample_idx=0)
    selected_markers2_genename = [gene_name_list[int(i[1:])] for i in selected_markers2]
    ddpcr = []

    blood_sample_idx = date_keys_sorted[order_idx]
    ddpcr_raw_sample = ddpcr_raw[blood_sample_idx]
    for gene in selected_markers2_genename:
        ddpcr.append({
            'gene': gene, 'mut': ddpcr_raw_sample["MutDOR"].loc[gene],
            'WT': ddpcr_raw_sample["DOR"].loc[gene], 'liquid_biopsy_sample': blood_sample_idx
        })

    df_ddpcr_2 = pd.DataFrame(ddpcr)
    marker_idx2gene = {i: df_ddpcr_2["gene"][i] for i in range(len(df_ddpcr_2))}


    tree_list_summary, node_list_summary, node_name_list_summary, tree_freq_list_summary = tree_distribution_summary[
        'tree_structure'], tree_distribution_summary[
        'node_dict'], tree_distribution_summary['node_dict_name'], tree_distribution_summary['freq']
    ddpcr_marker_counts = list(df_ddpcr_2["mut"])
    read_depth_list = list(df_ddpcr_2["mut"] + df_ddpcr_2["WT"])

    updated_tree_freq_list = adjust_tree_distribution_struct_bayesian(tree_list_summary, node_name_list_summary,
                                                                      tree_freq_list_summary, read_depth_list,
                                                                      ddpcr_marker_counts, marker_idx2gene)
    updated_tree_distribution_summary = update_tree_distribution_bayesian(tree_distribution_summary,
                                                                          updated_tree_freq_list)

    tree_distribution_file_summary_updated = file / f'{method}_bootstrap_summary_updated_{algo}_{n_markers}_{order_idx}_bayesian.pkl'
    with open(tree_distribution_file_summary_updated, 'wb') as f:
        pickle.dump(updated_tree_distribution_summary, f)