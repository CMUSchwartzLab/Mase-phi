import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import graphviz
import json
import os
from pathlib import Path
import warnings
import argparse
from analyze import *
from visualize import *
import pickle
from optimize_fraction import *
from adjust_tree_distribution import *
from optimize import *

def get_args(argv):
    parser = argparse.ArgumentParser(prog='run_simu_select_adjust.py',)
    parser.add_argument('-n', '--num_node', type=int, dest='num_node', default=7)
    parser.add_argument('-dt', '--depth_tissue', type=int, dest='depth_tissue', default=800)
    parser.add_argument('-db', '--depth_blood', type=int, dest='depth_blood', default=3000)
    parser.add_argument('-mr', '--mut_rate', type=int, dest='mut_rate', default=50)
    parser.add_argument('-rr', '--recover_rate', type=float, dest='recover_rate', default=0.5)
    parser.add_argument('-mp', '--mask_proportion', type=float, dest='mask_proportion', default=0.5)
    parser.add_argument('-dir', '--directory', type=str, dest='directory', default='/data/liquid_biopsy/simulations/mask/multi-tumor')
    parser.add_argument('-d', '--index', type=int, dest='simulation_idx', default=1)
    parser.add_argument('-nt', '--num_tissue', type=int, dest='num_tissue', default=2)
    parser.add_argument('-nb', '--num_blood', type=int, dest='num_blood', default=2)
    parser.add_argument('-unt', '--used_num_tissue', type=int, dest='used_num_tissue', default=2)
    parser.add_argument('-unb', '--used_num_blood', type=int, dest='used_num_blood', default=1)
    parser.add_argument('-b', '--num_bootstrap', type=int, dest='num_bootstrap', default=100)
    parser.add_argument('-m', '--method', type=str, dest='method', default="phylowgs")
    parser.add_argument('-t', '--type', type=str, dest='type', default="common")
    parser.add_argument('-mk', '--num_marker', type=int, dest='num_marker', default=3)
    parser.add_argument('-bi', '--blood_sample_idx', type=int, dest='blood_sample_idx', default=0)
    return vars(parser.parse_args(argv))

def main_select_adjust(argv):
    args = get_args(argv)
    num_node = args['num_node']
    main_dir = Path(args['directory'])
    depth_tissue = args['depth_tissue']
    depth_blood = args['depth_blood']
    mut_rate = args['mut_rate']
    recover_rate = args['recover_rate']
    simulation_idx = args['simulation_idx']
    mask_proportion = args['mask_proportion']
    num_tissue = args['num_tissue']
    num_blood = args['num_blood']
    used_num_tissue = args['used_num_tissue']
    used_num_blood = args['used_num_blood']
    num_bootstrap = args['num_bootstrap']
    method = args['method']
    type = args['type']
    n_markers = args['num_marker']
    blood_sample_idx = args['blood_sample_idx']
    directory_path = main_dir / f'boot_simulation_{num_node}_{depth_tissue}_{depth_blood}_{mut_rate}_{int(recover_rate*100)}_{int(mask_proportion*100)}/{num_blood}_{num_tissue}'
    num_chains=5
    #for file in directory_path.glob(f'Simulation*'):
    file = directory_path / f'Simulation_{simulation_idx}'
    patient_num = int(file.stem.split('_')[-1])
    print(patient_num)
    node_dict_true_file = file / "intermediate_value.json"
    with open(node_dict_true_file, 'rb') as f:
        node_dict_true = json.load(f)['Mutation']
    mut2node_dict = mut2node(node_dict_true)

    directory_xlsx = file / f'Bootstrap/simulate_bootstrap.xlsx'
    df = pd.read_excel(directory_xlsx, index_col=0, sheet_name='common_blood_tissue_no_germline')
    tree_distribution_summary = {'cp_tree': [], 'node_dict': [], 'node_dict_name': [], 'node_dict_re': [], 'tree_structure': [], 'freq': [], 'clonal_freq': [], 'vaf_frac': []}
    tree_distribution_aggre = {'cp_tree': [], 'node_dict': [], 'node_dict_name': [], 'node_dict_re': [], 'tree_structure': [], 'freq': [], 'clonal_freq': [], 'vaf_frac': []}
    for bootstrap_idx in range(1, num_bootstrap + 1, 1):
        summ_file = file / f"Bootstrap/{used_num_blood}_{used_num_tissue}/{type}/bootstrap{bootstrap_idx}/result_{num_chains}.summ.json.gz"
        muts_file = file / f"Bootstrap/{used_num_blood}_{used_num_tissue}/{type}/bootstrap{bootstrap_idx}/result_{num_chains}.muts.json.gz"
        mutass_file = file / f"Bootstrap/{used_num_blood}_{used_num_tissue}/{type}/bootstrap{bootstrap_idx}/result_{num_chains}.mutass.zip"
        try:
            tree_structure, node_dict, node_dict_name, node_dict_re, final_tree_cp, prev_mat, clonal_freq, vaf_frac = process_phylowgs_output(summ_file, muts_file, mutass_file)
            #print(vaf_frac)
        except:
            print(file, bootstrap_idx)
            continue
        tree_distribution_aggre['tree_structure'].append(tree_structure)
        tree_distribution_aggre['cp_tree'].append(final_tree_cp)
        tree_distribution_aggre['node_dict'].append(node_dict)
        tree_distribution_aggre['node_dict_re'].append(node_dict_re)
        tree_distribution_aggre['node_dict_name'].append(node_dict_name)
        tree_distribution_aggre['freq'].append(1)
        tree_distribution_aggre['clonal_freq'].append(clonal_freq)
        tree_distribution_aggre['vaf_frac'].append(vaf_frac)
        tree_distribution_summary = combine_tree(node_dict, node_dict_name, node_dict_re, tree_structure, final_tree_cp, clonal_freq, vaf_frac, 'phylowgs', tree_distribution_summary)


        #print(tree_distribution_aggre['vaf_frac'])

    # f.close()
    analyze_tree_distribution(tree_distribution_summary, file, patient_num, type)
    best_bootstrap_idx = np.argmax(tree_distribution_summary['freq'])
    tree_structure = tree_distribution_summary['tree_structure'][best_bootstrap_idx]
    node_dict = tree_distribution_summary['node_dict'][best_bootstrap_idx]
    results_dict = {'node_dict_name': node_dict,'tree_structure': tree_structure}
    with open(file / f"Bootstrap/{used_num_blood}_{used_num_tissue}/common/{patient_num}_results_bootstrap_{type}_best.json", 'w') as f:
        json.dump(results_dict, f)
    # g = render_tumor_tree(tree_structure, node_dict)
    # g.render(filename=directory / f"{patient_num}_results_bootstrap_{type}_best_inter")
    tree_distribution_file_summary = file / f'Bootstrap/{used_num_blood}_{used_num_tissue}/common/{method}_bootstrap_summary.pkl'
    # with open(tree_distribution_file_summary, 'wb') as g:
    #     pickle.dump(tree_distribution_summary, g)
    #
    # with open(file / f'Bootstrap/{used_num_blood}_{used_num_tissue}/common/{method}_bootstrap_aggregation.pkl', 'wb') as g:
    #     pickle.dump(tree_distribution_aggre, g)



    simu_idx = file.stem.split('_')[-1]
    print("Simulation", simu_idx)
    inter = pd.read_excel(file / f'Bootstrap/simulate_bootstrap.xlsx', sheet_name='common_blood_tissue_no_germline',
                          index_col=0)
    calls = inter
    gene2idx = {'s' + str(i): i for i in range(len(inter))}
    gene_list = list(gene2idx.keys())
    gene_name_list = []
    for i in range(inter.shape[0]):
        gene = calls["Gene"][i]
        if not isinstance(gene, str):
            gene = str(calls["Chromosome"][i]) + '_' + str(calls["Genomic Position"][i])
        gene_name_list.append(gene)
    t_idx = used_num_tissue
    tree_list_aggre, node_list_aggre, clonal_freq_list_aggre, tree_freq_list_aggre = tree_distribution_aggre['tree_structure'], \
        tree_distribution_aggre['node_dict'], tree_distribution_aggre['vaf_frac'], tree_distribution_aggre['freq']

    node_list_aggre_scrub = []
    for node_dict in node_list_aggre:
        temp = {}
        for key, values in node_dict.items():
            temp.setdefault(int(key), values)
        node_list_aggre_scrub.append(temp)

    clonal_freq_list_aggre_scrub = []
    for clonal_freq_dict in clonal_freq_list_aggre:
        temp = {}
        for key, values in clonal_freq_dict.items():
            temp.setdefault(int(key), values[0])
        clonal_freq_list_aggre_scrub.append(temp)




    ddpcr_raw = pd.read_excel(f"{file}/ddPCR/ddPCR_data.xlsx", sheet_name="Sheet1", index_col=0)




    alpha_list = [0.01, 0.05, 0.1, 0.2]
    read_depth = 1500
    algo_list = ['struct', 'frac']
    algo = 'frac'
    lam1 = 1  # lam1 is tree fractions
    lam2 = 0  # lam2 is tree structure
    # selected_markers2_genename_ordered = []
    selected_markers2, obj_frac, obj_struct = select_markers_tree_gp(gene_list, n_markers, tree_list_aggre,
                                                                     node_list_aggre, clonal_freq_list_aggre,
                                                                     gene2idx, tree_freq_list_aggre,
                                                                     read_depth=read_depth, lam1=lam1, lam2=lam2)
    selected_markers2_genename_frac = [gene_name_list[int(i[1:])] for i in selected_markers2]
    selected_markers2_file = file / f'Bootstrap/{used_num_blood}_{used_num_tissue}/{type}/{method}_bootstrap_frac_{n_markers}.pkl'
    with open(selected_markers2_file, 'wb') as f:
        pickle.dump(selected_markers2_genename_frac, f)
    ddpcr = []
    blood_sample_idx = 0
    for gene in selected_markers2_genename_frac:
        ddpcr.append({
            'gene': gene, 'mut': ddpcr_raw[blood_sample_idx * 2].loc[gene],
            'WT': ddpcr_raw[blood_sample_idx * 2 + 1].loc[gene], 'liquid_biopsy_sample': blood_sample_idx
        })

    df_ddpcr_2 = pd.DataFrame(ddpcr)
    marker_idx2gene = {i: df_ddpcr_2["gene"][i] for i in range(len(df_ddpcr_2))}
    tree_list_summary, node_list_summary, node_name_list_summary, tree_freq_list_summary = tree_distribution_summary[
        'tree_structure'], tree_distribution_summary[
        'node_dict'], tree_distribution_summary['node_dict_name'], tree_distribution_summary['freq']
    ddpcr_marker_counts = list(df_ddpcr_2["mut"])
    read_depth_list = list(df_ddpcr_2["mut"] + df_ddpcr_2["WT"])
    for alpha in alpha_list:
        for adjust_algo in ['struct', 'frac']:
            if adjust_algo == 'struct':
                accept_tree_indices_struct = adjust_tree_distribution_struct(tree_list_summary, node_name_list_summary, read_depth_list,
                                                                             ddpcr_marker_counts, marker_idx2gene, alpha)
                updated_tree_distribution_summary = update_tree_distribution(tree_distribution_summary,
                                                                             accept_tree_indices_struct)
            else:
                clonal_freq_list_aggre, node_dict_list_aggre = tree_distribution_aggre["clonal_freq"], \
                tree_distribution_aggre["node_dict_name"]
                accept_tree_indices_frac = adjust_tree_distribution_frac(clonal_freq_list_aggre, node_dict_list_aggre,
                                                                         read_depth_list, ddpcr_marker_counts,
                                                                         marker_idx2gene, alpha)
                updated_tree_distribution_aggre = update_tree_distribution(tree_distribution_aggre,
                                                                           accept_tree_indices_frac)
                updated_tree_distribution_summary = aggregation2summary(updated_tree_distribution_aggre)
            tree_distribution_file_summary_updated = file / f'Bootstrap/{used_num_blood}_{used_num_tissue}/common/{method}_bootstrap_summary_updated_{algo}_{n_markers}_{blood_sample_idx}_{alpha}_{adjust_algo}.pkl'
            with open(tree_distribution_file_summary_updated, 'wb') as f:
                pickle.dump(updated_tree_distribution_summary, f)


    algo = 'struct'
    lam1 = 0  # lam1 is tree fractions
    lam2 = 1  # lam2 is tree distributions
    selected_markers2, obj_frac, obj_struct = select_markers_tree_gp(gene_list, n_markers, tree_list_aggre,
                                                                     node_list_aggre, clonal_freq_list_aggre,
                                                                     gene2idx, tree_freq_list_aggre,
                                                                     read_depth=read_depth, lam1=lam1, lam2=lam2)
    selected_markers2_genename_struc = [gene_name_list[int(i[1:])] for i in selected_markers2]
    selected_markers2_file = file / f'Bootstrap/{used_num_blood}_{used_num_tissue}/{type}/{method}_bootstrap_struct_{n_markers}.pkl'

    with open(selected_markers2_file, 'wb') as f:
        pickle.dump(selected_markers2_genename_struc, f)
    ddpcr = []
    for gene in selected_markers2_genename_struc:
        ddpcr.append({
            'gene': gene, 'mut': ddpcr_raw[blood_sample_idx * 2].loc[gene],
            'WT': ddpcr_raw[blood_sample_idx * 2 + 1].loc[gene], 'liquid_biopsy_sample': blood_sample_idx
        })

    df_ddpcr_2 = pd.DataFrame(ddpcr)
    marker_idx2gene = {i: df_ddpcr_2["gene"][i] for i in range(len(df_ddpcr_2))}
    tree_list_summary, node_list_summary, node_name_list_summary, tree_freq_list_summary = tree_distribution_summary[
        'tree_structure'], tree_distribution_summary[
        'node_dict'], tree_distribution_summary['node_dict_name'], tree_distribution_summary['freq']

    ddpcr_marker_counts = list(df_ddpcr_2["mut"])
    read_depth_list = list(df_ddpcr_2["mut"] + df_ddpcr_2["WT"])
    for alpha in alpha_list:
        for adjust_algo in algo_list:
            if adjust_algo == 'struct':
                accept_tree_indices_struct = adjust_tree_distribution_struct(tree_list_summary, node_name_list_summary, read_depth_list,
                                                                             ddpcr_marker_counts, marker_idx2gene, alpha)
                updated_tree_distribution_summary = update_tree_distribution(tree_distribution_summary,
                                                                             accept_tree_indices_struct)
            else:
                clonal_freq_list_aggre, node_dict_list_aggre = tree_distribution_aggre["clonal_freq"], \
                tree_distribution_aggre["node_dict_name"]
                accept_tree_indices_frac = adjust_tree_distribution_frac(clonal_freq_list_aggre, node_dict_list_aggre,
                                                                         read_depth_list, ddpcr_marker_counts,
                                                                         marker_idx2gene, alpha)
                updated_tree_distribution_aggre = update_tree_distribution(tree_distribution_aggre,
                                                                           accept_tree_indices_frac)
                updated_tree_distribution_summary = aggregation2summary(updated_tree_distribution_aggre)
            tree_distribution_file_summary_updated = file / f'Bootstrap/{used_num_blood}_{used_num_tissue}/common/{method}_bootstrap_summary_updated_{algo}_{n_markers}_{blood_sample_idx}_{alpha}_{adjust_algo}.pkl'
            with open(tree_distribution_file_summary_updated, 'wb') as f:
                pickle.dump(updated_tree_distribution_summary, f)

def main_select_adjust_bayesian(argv):
    args = get_args(argv)
    num_node = args['num_node']
    main_dir = Path(args['directory'])
    depth_tissue = args['depth_tissue']
    depth_blood = args['depth_blood']
    mut_rate = args['mut_rate']
    recover_rate = args['recover_rate']
    simulation_idx = args['simulation_idx']
    mask_proportion = args['mask_proportion']
    num_tissue = args['num_tissue']
    num_blood = args['num_blood']
    used_num_tissue = args['used_num_tissue']
    used_num_blood = args['used_num_blood']
    num_bootstrap = args['num_bootstrap']
    method = args['method']
    type = args['type']
    n_markers = args['num_marker']
    blood_sample_idx = args['blood_sample_idx']

    directory_path = main_dir / f'boot_simulation_{num_node}_{depth_tissue}_{depth_blood}_{mut_rate}_{int(recover_rate*100)}_{int(mask_proportion*100)}/{num_blood}_{num_tissue}'
    num_chains=5
    #for file in directory_path.glob(f'Simulation*'):
    file = directory_path / f'Simulation_{simulation_idx}'
    patient_num = int(file.stem.split('_')[-1])
    print(patient_num)
    node_dict_true_file = file / "intermediate_value.json"
    with open(node_dict_true_file, 'rb') as f:
        node_dict_true = json.load(f)['Mutation']
    mut2node_dict = mut2node(node_dict_true)

    directory_xlsx = file / f'Bootstrap/simulate_bootstrap.xlsx'
    df = pd.read_excel(directory_xlsx, index_col=0, sheet_name='common_blood_tissue_no_germline')
    tree_distribution_summary = {'cp_tree': [], 'node_dict': [], 'node_dict_name': [], 'node_dict_re': [], 'tree_structure': [], 'freq': [], 'clonal_freq': [], 'vaf_frac': []}
    tree_distribution_aggre = {'cp_tree': [], 'node_dict': [], 'node_dict_name': [], 'node_dict_re': [], 'tree_structure': [], 'freq': [], 'clonal_freq': [], 'vaf_frac': []}
    for bootstrap_idx in range(1, num_bootstrap + 1, 1):
        summ_file = file / f"Bootstrap/{used_num_blood}_{used_num_tissue}/{type}/bootstrap{bootstrap_idx}/result_{num_chains}.summ.json.gz"
        muts_file = file / f"Bootstrap/{used_num_blood}_{used_num_tissue}/{type}/bootstrap{bootstrap_idx}/result_{num_chains}.muts.json.gz"
        mutass_file = file / f"Bootstrap/{used_num_blood}_{used_num_tissue}/{type}/bootstrap{bootstrap_idx}/result_{num_chains}.mutass.zip"
        try:
            tree_structure, node_dict, node_dict_name, node_dict_re, final_tree_cp, prev_mat, clonal_freq, vaf_frac = process_phylowgs_output(summ_file, muts_file, mutass_file)
            #print(vaf_frac)
        except:
            print(file, bootstrap_idx)
            continue
        tree_distribution_aggre['tree_structure'].append(tree_structure)
        tree_distribution_aggre['cp_tree'].append(final_tree_cp)
        tree_distribution_aggre['node_dict'].append(node_dict)
        tree_distribution_aggre['node_dict_re'].append(node_dict_re)
        tree_distribution_aggre['node_dict_name'].append(node_dict_name)
        tree_distribution_aggre['freq'].append(1)
        tree_distribution_aggre['clonal_freq'].append(clonal_freq)
        tree_distribution_aggre['vaf_frac'].append(vaf_frac)
        tree_distribution_summary = combine_tree(node_dict, node_dict_name, node_dict_re, tree_structure, final_tree_cp, clonal_freq, vaf_frac, 'phylowgs', tree_distribution_summary)


        #print(tree_distribution_aggre['vaf_frac'])

    # f.close()
    analyze_tree_distribution(tree_distribution_summary, file, patient_num, type)
    best_bootstrap_idx = np.argmax(tree_distribution_summary['freq'])
    tree_structure = tree_distribution_summary['tree_structure'][best_bootstrap_idx]
    node_dict = tree_distribution_summary['node_dict'][best_bootstrap_idx]
    results_dict = {'node_dict_name': node_dict,'tree_structure': tree_structure}
    with open(file / f"Bootstrap/{used_num_blood}_{used_num_tissue}/common/{patient_num}_results_bootstrap_{type}_best.json", 'w') as f:
        json.dump(results_dict, f)
    # g = render_tumor_tree(tree_structure, node_dict)
    # g.render(filename=directory / f"{patient_num}_results_bootstrap_{type}_best_inter")
    tree_distribution_file_summary = file / f'Bootstrap/{used_num_blood}_{used_num_tissue}/common/{method}_bootstrap_summary.pkl'
    # with open(tree_distribution_file_summary, 'wb') as g:
    #     pickle.dump(tree_distribution_summary, g)
    #
    # with open(file / f'Bootstrap/{used_num_blood}_{used_num_tissue}/common/{method}_bootstrap_aggregation.pkl', 'wb') as g:
    #     pickle.dump(tree_distribution_aggre, g)



    simu_idx = file.stem.split('_')[-1]
    print("Simulation", simu_idx)
    inter = pd.read_excel(file / f'Bootstrap/simulate_bootstrap.xlsx', sheet_name='common_blood_tissue_no_germline',
                          index_col=0)
    calls = inter
    gene2idx = {'s' + str(i): i for i in range(len(inter))}
    gene_list = list(gene2idx.keys())
    gene_name_list = []
    for i in range(inter.shape[0]):
        gene = calls["Gene"][i]
        if not isinstance(gene, str):
            gene = str(calls["Chromosome"][i]) + '_' + str(calls["Genomic Position"][i])
        gene_name_list.append(gene)
    t_idx = used_num_tissue
    tree_list_aggre, node_list_aggre, clonal_freq_list_aggre, tree_freq_list_aggre = tree_distribution_aggre['tree_structure'], \
        tree_distribution_aggre['node_dict'], tree_distribution_aggre['vaf_frac'], tree_distribution_aggre['freq']

    node_list_aggre_scrub = []
    for node_dict in node_list_aggre:
        temp = {}
        for key, values in node_dict.items():
            temp.setdefault(int(key), values)
        node_list_aggre_scrub.append(temp)

    clonal_freq_list_aggre_scrub = []
    for clonal_freq_dict in clonal_freq_list_aggre:
        temp = {}
        for key, values in clonal_freq_dict.items():
            temp.setdefault(int(key), values[0])
        clonal_freq_list_aggre_scrub.append(temp)




    ddpcr_raw = pd.read_excel(f"{file}/ddPCR/ddPCR_data.xlsx", sheet_name="Sheet1", index_col=0)

    read_depth = 1500
    algo = 'frac'
    lam1 = 1  # lam1 is tree fractions
    lam2 = 0  # lam2 is tree structure
    selected_markers2_genename_ordered = []
    selected_markers2, obj_frac, obj_struct = select_markers_tree_gp(gene_list, n_markers, tree_list_aggre,
                                                                     node_list_aggre, clonal_freq_list_aggre,
                                                                     gene2idx, tree_freq_list_aggre,
                                                                     read_depth=read_depth, lam1=lam1, lam2=lam2)
    selected_markers2_genename_struc = [gene_name_list[int(i[1:])] for i in selected_markers2]
    selected_markers2_file = file / f'Bootstrap/{used_num_blood}_{used_num_tissue}/{type}/{method}_bootstrap_{algo}_{n_markers}.pkl'

    # with open(selected_markers2_file, 'wb') as f:
    #     pickle.dump(selected_markers2_genename_struc, f)
    with open(selected_markers2_file, 'rb') as f:
        selected_markers2_genename_frac = pickle.load(f)
    ddpcr = []
    for gene in selected_markers2_genename_frac:
        ddpcr.append({
            'gene': gene, 'mut': ddpcr_raw[blood_sample_idx * 2].loc[gene],
            'WT': ddpcr_raw[blood_sample_idx * 2 + 1].loc[gene], 'liquid_biopsy_sample': blood_sample_idx
        })

    df_ddpcr_2 = pd.DataFrame(ddpcr)
    marker_idx2gene = {i: df_ddpcr_2["gene"][i] for i in range(len(df_ddpcr_2))}
    tree_list_summary, node_list_summary, node_name_list_summary, tree_freq_list_summary = tree_distribution_summary[
        'tree_structure'], tree_distribution_summary[
        'node_dict'], tree_distribution_summary['node_dict_name'], tree_distribution_summary['freq']
    ddpcr_marker_counts = list(df_ddpcr_2["mut"])
    read_depth_list = list(df_ddpcr_2["mut"] + df_ddpcr_2["WT"])

    updated_tree_freq_list = adjust_tree_distribution_struct_bayesian(tree_list_summary, node_name_list_summary,
                                tree_freq_list_summary, read_depth_list, ddpcr_marker_counts, marker_idx2gene)
    updated_tree_distribution_summary = update_tree_distribution_bayesian(tree_distribution_summary, updated_tree_freq_list)

    tree_distribution_file_summary_updated = file / f'Bootstrap/{used_num_blood}_{used_num_tissue}/common/{method}_bootstrap_summary_updated_{algo}_{n_markers}_{blood_sample_idx}_bayesian.pkl'
    with open(tree_distribution_file_summary_updated, 'wb') as f:
        pickle.dump(updated_tree_distribution_summary, f)


    algo = 'struct'
    lam1 = 0  # lam1 is tree fractions
    lam2 = 1  # lam2 is tree distributions
    # selected_markers2, obj_frac, obj_struct = select_markers_tree_gp(gene_list, n_markers, tree_list_aggre,
    #                                                                  node_list_aggre, clonal_freq_list_aggre,
    #                                                                  gene2idx, tree_freq_list_aggre,
    #                                                                  read_depth=read_depth, lam1=lam1, lam2=lam2)
    #selected_markers2_genename_struc = [gene_name_list[int(i[1:])] for i in selected_markers2]
    selected_markers2_file = file / f'Bootstrap/{used_num_blood}_{used_num_tissue}/{type}/{method}_bootstrap_{algo}_{n_markers}.pkl'

    # with open(selected_markers2_file, 'wb') as f:
    #     pickle.dump(selected_markers2_genename_struc, f)
    with open(selected_markers2_file, 'rb') as f:
        selected_markers2_genename_frac = pickle.load(f)
    ddpcr = []
    blood_sample_idx = 0
    for gene in selected_markers2_genename_frac:
        ddpcr.append({
            'gene': gene, 'mut': ddpcr_raw[blood_sample_idx * 2].loc[gene],
            'WT': ddpcr_raw[blood_sample_idx * 2 + 1].loc[gene], 'liquid_biopsy_sample': blood_sample_idx
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

    tree_distribution_file_summary_updated = file / f'Bootstrap/{used_num_blood}_{used_num_tissue}/common/{method}_bootstrap_summary_updated_{algo}_{n_markers}_{blood_sample_idx}_bayesian.pkl'
    with open(tree_distribution_file_summary_updated, 'wb') as f:
        pickle.dump(updated_tree_distribution_summary, f)


def random_select_adjust_bayesian(argv):
    args = get_args(argv)
    num_node = args['num_node']
    main_dir = Path(args['directory'])
    depth_tissue = args['depth_tissue']
    depth_blood = args['depth_blood']
    mut_rate = args['mut_rate']
    recover_rate = args['recover_rate']
    simulation_idx = args['simulation_idx']
    mask_proportion = args['mask_proportion']
    num_tissue = args['num_tissue']
    num_blood = args['num_blood']
    used_num_tissue = args['used_num_tissue']
    used_num_blood = args['used_num_blood']
    num_bootstrap = args['num_bootstrap']
    method = args['method']
    type = args['type']
    n_markers = args['num_marker']
    blood_sample_idx = args['blood_sample_idx']
    directory_path = main_dir / f'boot_simulation_{num_node}_{depth_tissue}_{depth_blood}_{mut_rate}_{int(recover_rate*100)}_{int(mask_proportion*100)}/{num_blood}_{num_tissue}'
    num_chains=5
    #for file in directory_path.glob(f'Simulation*'):
    file = directory_path / f'Simulation_{simulation_idx}'
    patient_num = int(file.stem.split('_')[-1])
    print(patient_num)
    node_dict_true_file = file / "intermediate_value.json"
    with open(node_dict_true_file, 'rb') as f:
        node_dict_true = json.load(f)['Mutation']
    mut2node_dict = mut2node(node_dict_true)

    directory_xlsx = file / f'Bootstrap/simulate_bootstrap.xlsx'
    df = pd.read_excel(directory_xlsx, index_col=0, sheet_name='common_blood_tissue_no_germline')
    tree_distribution_summary = {'cp_tree': [], 'node_dict': [], 'node_dict_name': [], 'node_dict_re': [], 'tree_structure': [], 'freq': [], 'clonal_freq': [], 'vaf_frac': []}
    tree_distribution_aggre = {'cp_tree': [], 'node_dict': [], 'node_dict_name': [], 'node_dict_re': [], 'tree_structure': [], 'freq': [], 'clonal_freq': [], 'vaf_frac': []}
    for bootstrap_idx in range(1, num_bootstrap + 1, 1):
        summ_file = file / f"Bootstrap/{used_num_blood}_{used_num_tissue}/{type}/bootstrap{bootstrap_idx}/result_{num_chains}.summ.json.gz"
        muts_file = file / f"Bootstrap/{used_num_blood}_{used_num_tissue}/{type}/bootstrap{bootstrap_idx}/result_{num_chains}.muts.json.gz"
        mutass_file = file / f"Bootstrap/{used_num_blood}_{used_num_tissue}/{type}/bootstrap{bootstrap_idx}/result_{num_chains}.mutass.zip"
        try:
            tree_structure, node_dict, node_dict_name, node_dict_re, final_tree_cp, prev_mat, clonal_freq, vaf_frac = process_phylowgs_output(summ_file, muts_file, mutass_file)
        except:
            print(file, bootstrap_idx)
            continue
        tree_distribution_aggre['tree_structure'].append(tree_structure)
        tree_distribution_aggre['cp_tree'].append(final_tree_cp)
        tree_distribution_aggre['node_dict'].append(node_dict)
        tree_distribution_aggre['node_dict_re'].append(node_dict_re)
        tree_distribution_aggre['node_dict_name'].append(node_dict_name)
        tree_distribution_aggre['freq'].append(1)
        tree_distribution_aggre['clonal_freq'].append(clonal_freq)
        tree_distribution_aggre['vaf_frac'].append(vaf_frac)
        tree_distribution_summary = combine_tree(node_dict, node_dict_name, node_dict_re, tree_structure, final_tree_cp, clonal_freq, vaf_frac, 'phylowgs', tree_distribution_summary)


        #print(tree_distribution_aggre['vaf_frac'])

    # f.close()
    analyze_tree_distribution(tree_distribution_summary, file, patient_num, type)
    best_bootstrap_idx = np.argmax(tree_distribution_summary['freq'])
    tree_structure = tree_distribution_summary['tree_structure'][best_bootstrap_idx]
    node_dict = tree_distribution_summary['node_dict'][best_bootstrap_idx]
    results_dict = {'node_dict_name': node_dict,'tree_structure': tree_structure}
    with open(file / f"Bootstrap/{used_num_blood}_{used_num_tissue}/common/{patient_num}_results_bootstrap_{type}_best.json", 'w') as f:
        json.dump(results_dict, f)
    # g = render_tumor_tree(tree_structure, node_dict)
    # g.render(filename=directory / f"{patient_num}_results_bootstrap_{type}_best_inter")
    tree_distribution_file_summary = file / f'Bootstrap/{used_num_blood}_{used_num_tissue}/common/{method}_bootstrap_summary.pkl'
    # with open(tree_distribution_file_summary, 'wb') as g:
    #     pickle.dump(tree_distribution_summary, g)
    #
    # with open(file / f'Bootstrap/{used_num_blood}_{used_num_tissue}/common/{method}_bootstrap_aggregation.pkl', 'wb') as g:
    #     pickle.dump(tree_distribution_aggre, g)



    simu_idx = file.stem.split('_')[-1]
    print("Simulation", simu_idx)
    inter = pd.read_excel(file / f'Bootstrap/simulate_bootstrap.xlsx', sheet_name='common_blood_tissue_no_germline',
                          index_col=0)
    calls = inter
    gene2idx = {'s' + str(i): i for i in range(len(inter))}
    gene_list = list(gene2idx.keys())
    gene_name_list = []
    for i in range(inter.shape[0]):
        gene = calls["Gene"][i]
        if not isinstance(gene, str):
            gene = str(calls["Chromosome"][i]) + '_' + str(calls["Genomic Position"][i])
        gene_name_list.append(gene)
    t_idx = used_num_tissue
    tree_list_aggre, node_list_aggre, clonal_freq_list_aggre, tree_freq_list_aggre = tree_distribution_aggre['tree_structure'], \
        tree_distribution_aggre['node_dict'], tree_distribution_aggre['vaf_frac'], tree_distribution_aggre['freq']

    node_list_aggre_scrub = []
    for node_dict in node_list_aggre:
        temp = {}
        for key, values in node_dict.items():
            temp.setdefault(int(key), values)
        node_list_aggre_scrub.append(temp)

    clonal_freq_list_aggre_scrub = []
    for clonal_freq_dict in clonal_freq_list_aggre:
        temp = {}
        for key, values in clonal_freq_dict.items():
            temp.setdefault(int(key), values[0])
        clonal_freq_list_aggre_scrub.append(temp)




    ddpcr_raw = pd.read_excel(f"{file}/ddPCR/ddPCR_data.xlsx", sheet_name="Sheet1", index_col=0)

    read_depth = 1500
    algo = 'random'
    selected_markers2_genename_ordered = []
    selected_markers2 = random.sample(gene_list, n_markers)
    selected_markers2_genename_random = [gene_name_list[int(i[1:])] for i in selected_markers2]
    selected_markers2_file = file / f'Bootstrap/{used_num_blood}_{used_num_tissue}/{type}/{method}_bootstrap_{algo}_{n_markers}.pkl'

    with open(selected_markers2_file, 'wb') as f:
        pickle.dump(selected_markers2_genename_random, f)
    with open(selected_markers2_file, 'rb') as f:
        selected_markers2_genename_random = pickle.load(f)
    ddpcr = []
    for gene in selected_markers2_genename_random:
        ddpcr.append({
            'gene': gene, 'mut': ddpcr_raw[blood_sample_idx * 2].loc[gene],
            'WT': ddpcr_raw[blood_sample_idx * 2 + 1].loc[gene], 'liquid_biopsy_sample': blood_sample_idx
        })

    df_ddpcr_2 = pd.DataFrame(ddpcr)
    marker_idx2gene = {i: df_ddpcr_2["gene"][i] for i in range(len(df_ddpcr_2))}
    tree_list_summary, node_list_summary, node_name_list_summary, tree_freq_list_summary = tree_distribution_summary[
        'tree_structure'], tree_distribution_summary[
        'node_dict'], tree_distribution_summary['node_dict_name'], tree_distribution_summary['freq']
    ddpcr_marker_counts = list(df_ddpcr_2["mut"])
    read_depth_list = list(df_ddpcr_2["mut"] + df_ddpcr_2["WT"])

    updated_tree_freq_list = adjust_tree_distribution_struct_bayesian(tree_list_summary, node_name_list_summary,
                                tree_freq_list_summary, read_depth_list, ddpcr_marker_counts, marker_idx2gene)
    updated_tree_distribution_summary = update_tree_distribution_bayesian(tree_distribution_summary, updated_tree_freq_list)

    tree_distribution_file_summary_updated = file / f'Bootstrap/{used_num_blood}_{used_num_tissue}/common/{method}_bootstrap_summary_updated_{algo}_{n_markers}_{blood_sample_idx}_bayesian.pkl'
    with open(tree_distribution_file_summary_updated, 'wb') as f:
        pickle.dump(updated_tree_distribution_summary, f)



if __name__ == '__main__':
    #main_select_adjust(sys.argv[1:])
    #main_select_adjust_bayesian(sys.argv[1:])
    #random_select_adjust_bayesian(sys.argv[1:])