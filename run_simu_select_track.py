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
    parser = argparse.ArgumentParser(prog='run_simu_select_track.py',)
    parser.add_argument('-n', '--num_node', type=int, dest='num_node', default=10)
    parser.add_argument('-dt', '--depth_tissue', type=int, dest='depth_tissue', default=800)
    parser.add_argument('-db', '--depth_blood', type=int, dest='depth_blood', default=3000)
    parser.add_argument('-mr', '--mut_rate', type=int, dest='mut_rate', default=100)
    parser.add_argument('-rr', '--recover_rate', type=float, dest='recover_rate', default=0.5)
    parser.add_argument('-mp', '--mask_proportion', type=float, dest='mask_proportion', default=0.5)
    parser.add_argument('-dir', '--directory', type=str, dest='directory', default='/data/liquid_biopsy/simulations/mask/multi-tumor')
    parser.add_argument('-nt', '--num_tissue', type=int, dest='num_tissue', default=3)
    parser.add_argument('-nb', '--num_blood', type=int, dest='num_blood', default=3)
    parser.add_argument('-unt', '--used_num_tissue', type=int, dest='used_num_tissue', default=3)
    parser.add_argument('-unb', '--used_num_blood', type=int, dest='used_num_blood', default=1)
    parser.add_argument('-b', '--num_bootstrap', type=int, dest='num_bootstrap', default=100)
    parser.add_argument('-m', '--method', type=str, dest='method', default="phylowgs")
    parser.add_argument('-t', '--type', type=str, dest='type', default="common")
    parser.add_argument('-mk', '--num_marker', type=int, dest='num_marker', default=3)
    return vars(parser.parse_args(argv))

def input_argv(argv):
    args = get_args(argv)
    num_node = args['num_node']
    main_dir = args['directory']
    depth_tissue = args['depth_tissue']
    depth_blood = args['depth_blood']
    mut_rate = args['mut_rate']
    recover_rate = args['recover_rate']
    mask_proportion = args['mask_proportion']
    num_tissue = args['num_tissue']
    num_blood = args['num_blood']
    used_num_tissue = args['used_num_tissue']
    used_num_blood = args['used_num_blood']
    num_bootstrap = args['num_bootstrap']
    method = args['method']
    type = args['type']
    n_markers = args['num_marker']
    return num_node, main_dir, depth_tissue, depth_blood, mut_rate, recover_rate, mask_proportion, num_tissue, num_blood, \
        used_num_tissue, used_num_blood, num_bootstrap, method, type, n_markers

def main_select_track(num_node, main_dir, depth_tissue, depth_blood, mut_rate, recover_rate, mask_proportion, num_tissue, num_blood, \
        used_num_tissue, used_num_blood, num_bootstrap, method, type, n_markers):
    main_dir = Path(main_dir)
    directory_path = main_dir / f'boot_simulation_{num_node}_{depth_tissue}_{depth_blood}_{mut_rate}_{int(recover_rate*100)}_{int(mask_proportion*100)}/{num_blood}_{num_tissue}'
    num_chains=5
    for file in directory_path.glob(f'Simulation*'):
        simulation_idx = int(file.stem.split("_")[-1])
    #file = directory_path / f'Simulation_{simulation_idx}'
        patient_num = simulation_idx
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

        selected_markers1_genename_ordered = []
        obj1_ordered = []

        num_marker = len(gene_list)
        print(num_marker)
        df_fractions = []
        for n_markers in range(1, num_marker + 1):
            print(clonal_freq_list_aggre)
            selected_markers1_random = random.sample(gene_list, n_markers)
            print(selected_markers1_random)
            selected_markers1_genename_random = [gene_name_list[int(i[1:])] for i in selected_markers1_random]
            #selected_markers1, obj = select_markers_fractions_weighted_single(gene_list, n_markers, tree_list, node_list_scrub, clonal_freq_list_scrub, gene2idx, tree_freq_list, 0)
            try:
                selected_markers1, obj = select_markers_fractions_weighted_overall(gene_list, n_markers, tree_list_aggre, node_list_aggre, clonal_freq_list_aggre_scrub, gene2idx, tree_freq_list_aggre, None)
            except AttributeError:
                break
            selected_markers1_genename = [gene_name_list[int(i[1:])] for i in selected_markers1]
            ddpcr = []
            blood_sample_idx = 1
            for gene in selected_markers1_genename:
                ddpcr.append({
                    'gene': gene, 'mut': ddpcr_raw[blood_sample_idx * 2].loc[gene],
                    'WT': ddpcr_raw[blood_sample_idx * 2 + 1].loc[gene], 'liquid_biopsy_sample': blood_sample_idx
                })

            df_ddpcr_1 = pd.DataFrame(ddpcr)
            ddpcr_random = []
            blood_sample_idx = 1
            for gene in selected_markers1_genename_random:
                ddpcr_random.append({
                    'gene': gene, 'mut': ddpcr_raw[blood_sample_idx * 2].loc[gene],
                    'WT': ddpcr_raw[blood_sample_idx * 2 + 1].loc[gene], 'liquid_biopsy_sample': blood_sample_idx
                })

            df_ddpcr_1_random = pd.DataFrame(ddpcr_random)

            obj1_ordered.append(obj)
            if len(selected_markers1_genename) == 1:
                selected_markers1_genename_ordered.append(selected_markers1_genename[0])
            else:
                diff_set = set(selected_markers1_genename).difference(set(selected_markers1_genename_ordered))
                selected_markers1_genename_ordered.append(list(diff_set)[0])
            clonal_freq_list_mean = []
            clonal_freq_list_mean_random = []
            for node_dict in node_list_aggre:
                clonal_freq_new = {}
                mut2node_bootstrap = mut2node(node_dict)
                for marker in selected_markers1:
                    node = mut2node_bootstrap[marker]
                    if node not in clonal_freq_new.keys():
                        clonal_freq_new.setdefault(node, [])
                    marker_name = gene_name_list[int(marker[1:])]
                    clonal_freq_new[node].append(df_ddpcr_1["mut"].loc[df_ddpcr_1["gene"] == marker_name].item() / (df_ddpcr_1["mut"].loc[df_ddpcr_1["gene"] == marker_name].item() + df_ddpcr_1["WT"].loc[df_ddpcr_1["gene"] == marker_name].item()))
                clonal_freq_mean = {}
                for key, value in clonal_freq_new.items():
                    clonal_freq_mean[key] = np.mean(np.array(value))
                clonal_freq_list_mean.append(clonal_freq_mean)

                clonal_freq_new_random = {}
                mut2node_bootstrap = mut2node(node_dict)
                for marker in selected_markers1_random:
                    node = mut2node_bootstrap[marker]
                    if node not in clonal_freq_new_random.keys():
                        clonal_freq_new_random.setdefault(node, [])
                    marker_name = gene_name_list[int(marker[1:])]
                    clonal_freq_new_random[node].append(df_ddpcr_1_random["mut"].loc[df_ddpcr_1_random["gene"] == marker_name].item() / (
                                df_ddpcr_1_random["mut"].loc[df_ddpcr_1_random["gene"] == marker_name].item() + df_ddpcr_1_random["WT"].loc[
                            df_ddpcr_1_random["gene"] == marker_name].item()))
                clonal_freq_mean_random = {}
                for key, value in clonal_freq_new_random.items():
                    clonal_freq_mean_random[key] = np.mean(np.array(value))
                clonal_freq_list_mean_random.append(clonal_freq_mean_random)
            F_list_pcr, F_hat_list_pcr = create_F_F_hat_list_from_pcr(clonal_freq_list_mean, tree_list_aggre, clonal_freq_list_aggre)
            F_list_pcr_random, F_hat_list_pcr_random = create_F_F_hat_list_from_pcr(clonal_freq_list_mean_random, tree_list_aggre, clonal_freq_list_aggre_scrub)
            trackable_fraction_list = [sum(F_hat_list_pcr[i]) for i in range(len(F_hat_list_pcr))]
            trackable_fraction = np.mean(trackable_fraction_list)
            trackable_fraction_list_random = [sum(F_hat_list_pcr_random[i]) for i in range(len(F_hat_list_pcr_random))]
            trackable_fraction_random = np.mean(trackable_fraction_list_random)
            df_fractions.append({'gene_rank': n_markers, 'trackable_fractions': trackable_fraction, 'simu_idx': simulation_idx,
                                 'trackable_fractions_random': trackable_fraction_random})
        df_fractions = pd.DataFrame(df_fractions)
        print('df', df_fractions)
        df_fraction_file = file / f'Bootstrap/{used_num_blood}_{used_num_tissue}/common/{method}_fraction_dataframe.csv'
        df_fractions.to_csv(df_fraction_file, sep ='\t')
        print(selected_markers1_genename_ordered, obj1_ordered)




if __name__ == '__main__':
    #num_node, main_dir, depth_tissue, depth_blood, mut_rate, recover_rate, mask_proportion, num_tissue, num_blood, \
    #    used_num_tissue, used_num_blood, num_bootstrap, method, type, n_markers = input(sys.argv[1:])
    main_dir = '/data/liquid_biopsy/simulations/mask/multi-tumor'
    parameter_set_list = [(10, 800, 3000, 100, 0.5, 0.5, 3, 3, 3, 1),
                          (10, 800, 3000, 100, 0.5, 0.3, 3, 3, 3, 1)]
    num_bootstrap = 100
    method = 'phylowgs'
    type = "common"
    for (num_node, depth_tissue, depth_blood, mut_rate, recover_rate, mask_proportion, num_tissue, num_blood,
        used_num_tissue, used_num_blood) in parameter_set_list:
        main_select_track(num_node, main_dir, depth_tissue, depth_blood, mut_rate, recover_rate, mask_proportion, num_tissue, num_blood,
        used_num_tissue, used_num_blood, num_bootstrap, method, type, n_markers=None)