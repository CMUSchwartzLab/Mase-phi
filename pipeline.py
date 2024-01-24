import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import graphviz
import json
import os
from pathlib import Path
import warnings
from simulate import simulate_clonal_tree, simulate_variant_reads, simulation_excel, mask_recent_clones_tissue
from bootstrap import phyloWGS_output, output_bootstrap, merge_phyloWGS_input_bootstrap, merge_phyloWGS_input
import argparse
sys.path.insert(0, '..')
from visualize import add_prefix_tree, render_tumor_tree

######################################

def simulation(num_nodes, depth_tissue, depth_blood, mut_rate, directory_path, num_bootstraps, surgery=False,
               recover_rate=0.5, mask_proportion=0.0, num_blood=2, num_tissue=2, used_num_blood=1, used_num_tissue=2):
    '''
    This function intended to integrate clonal tree & variant mutation simulation
    with bootstrap in pre-surgery scenario & post_surgery scenario:

        1) Simulate a clonal tree and its variant reads;
        2) Make the xlsx file and save them into a json file for later accessing;
        3) Based on the excel, use the bootstrap to generate bootstrapped data;
        4) Output the bootstrapped excel with the plyloWGS for each round of bootstrapping;

    Input:
    num_nodes: number of nodes for the simulated tree;
    depth: a poisson mean of sequencing depth;
    mut_rate: a poisson mean of number of mutations;
    directory_path: the path for directing the output simulation data;
    num_bootstraps: the number of bootstrap used;
    surgery: a boolean variable for whether we need the post-surgery simulation or not;
    recover_rate: the parameter for representing the extent of recovery after surgery;

    '''

    '''
    ###### 1) Simulate clonal trees and its variant reads; 
    '''
    while True:
        try:
            tree, freq, freq_T_list, freq_B_list = simulate_clonal_tree(k=num_nodes, recover_rate=recover_rate, num_blood=num_blood, num_tissue=num_tissue)  # parameter k
            print("freq_T" , freq_T_list, "freq_B", freq_B_list)
            break
        except:
            continue
    if mask_proportion != 0:
        print(freq_T_list)
        freq_T_list = mask_recent_clones_tissue(tree, freq_T_list, 0, mask_proportion)
        print('freq_T_masked', freq_T_list)

    mutation, variant_T_list, reference_T_list, variant_B_list, reference_B_list = simulate_variant_reads(
        tree, freq_T_list, freq_B_list, depth_tissue, depth_blood, [mut_rate / (num_nodes - 1)])

    '''
    ###### 2) Make the xlsx file and save each of them into a separated json file;
            - The path of the simulated_variant_read output can be changed;
    '''

    simulation_excel(variant_B_list, reference_B_list, variant_T_list, reference_T_list, Path(directory_path), used_num_blood, used_num_tissue)
    for t_idx in range(used_num_tissue):
        for b_idx in range(used_num_blood):
            phyloWGS_output(Path(os.path.join(directory_path, "simulated_mut.xlsx")), Path(directory_path), 0, surgery=False,
                  sheet_type="all", num_blood=1, num_tissue=t_idx+1)

    # if surgery:
    #     sample_type = '_post'
    #     simulation_excel(variant_B_post, reference_B_post, variant_T_list, reference_T_list, Path(directory_path), True)
    #
    #     phyloWGS_output(Path(os.path.join(directory_path, "simulated_mut_post.xlsx")), Path(directory_path), 0,
    #                   surgery=True, sheet_type="all")
    #     merge_phyloWGS_input(directory_path / 'union/ssm_data.txt', directory_path / 'union/_post/ssm_data_post.txt')
    #     #merge_phyloWGS_input(directory_path / 'common/ssm_data.txt', directory_path / 'common/_post/ssm_data_post.txt')
    sample_type=''
    mutation_prefix = add_prefix_tree(mutation)
    g = render_tumor_tree(tree_structure=tree, node_dict=mutation_prefix)
    g.render(filename=directory_path / 'true_tree')



    data_list = {"tree": tree,
                 "freq": freq,
                 "freq_T": freq_T_list,
                 "freq_B": freq_B_list,
                 "Mutation": mutation_prefix,
                 "Variant_T": variant_T_list,
                 "Reference_T": reference_T_list,
                 "Variant_B": variant_B_list,
                 "Reference_B": reference_B_list,
                 }
    with open(f"{directory_path}/intermediate_value.json", "w") as f:
        json_data = json.dumps(data_list, indent=4, separators=(',', ':'))
        f.write(json_data)

    # if surgery==True:

    ''' 
    ###### 3) Based on the excel, use the bootstrap to generate bootstrapped data;
    '''
    if num_bootstraps > 0:
        directory = Path(directory_path)  # Make the input in Path format;

        # Perform the bootstrapping for the input number of iteration;
        directory_bootstrap = Path(os.path.join(directory, "Bootstrap"))
        os.makedirs(directory_bootstrap, exist_ok=True)
        # Make a new directory for each iteration of bootstrapping;

        if surgery == True:
            directory_bootstrap_post = Path(os.path.join(directory, "Bootstrap_post"))
            os.makedirs(directory_bootstrap_post, exist_ok=True)

        if surgery:
            files_list = ['simulated_mut.xlsx', 'simulated_mut_post.xlsx']
        else:
            files_list = ['simulated_mut.xlsx']
        for file_name in files_list:

            file = directory / file_name
            # The number of the patient
            # Read the blood_no_germline and tissue_no germline sheet;
            # Store the whole file in a dictionary of dataframe;
            dfs = pd.read_excel(file, sheet_name=None, index_col=0, engine="openpyxl")

            dfs = output_bootstrap(dfs, num_bootstraps=num_bootstraps, used_num_blood=num_blood, used_num_tissue=num_tissue)

            if "_post" in file_name:
                path_new_bootstrap = Path(directory_bootstrap_post / f'simulate{sample_type}_bootstrap.xlsx')
            else:
                path_new_bootstrap = Path(directory_bootstrap / f'simulate{sample_type}_bootstrap.xlsx')

            with pd.ExcelWriter(path_new_bootstrap) as writer:

            # Write the bootstrap AF of each patient in this bootstrap directory
                for sheet_name, data in dfs.items():
                    data.to_excel(writer, sheet_name=sheet_name)


            '''
            ###### 4) Based on the bootstrap excel to generate the plyloWGS output;
            '''
            if "_post" in file.stem:
                phyloWGS_output(path_new_bootstrap, directory_bootstrap_post, num_bootstraps, True, sheet_type="common", num_tissue=2, num_blood=1)
                #merge_phyloWGS_input_bootstrap(directory_bootstrap / 'union', directory_bootstrap_post / 'union', num_bootstrap)
                #merge_phyloWGS_input_bootstrap(directory_bootstrap / 'common', directory_bootstrap_post / 'common', num_bootstrap)
            else:
                phyloWGS_output(path_new_bootstrap, directory_bootstrap, num_bootstraps, False, sheet_type="common",
                                        num_tissue=num_tissue, num_blood=1)
                if num_blood != 0:
                    phyloWGS_output(path_new_bootstrap, directory_bootstrap, num_bootstraps, False, sheet_type="common",
                                        num_tissue=num_tissue, num_blood=0)
                phyloWGS_output(path_new_bootstrap, directory_bootstrap, num_bootstraps, False, sheet_type="common",
                                num_tissue=1, num_blood=1)


def get_args(argv):
    parser = argparse.ArgumentParser(prog='pipeline.py',)
    parser.add_argument('-n', '--num_node', type=int, dest='num_node', default=10)
    parser.add_argument('-dt', '--depth_tissue', type=int, dest='depth_tissue', default=800)
    parser.add_argument('-db', '--depth_blood', type=int, dest='depth_blood', default=3000)
    parser.add_argument('-mr', '--mut_rate', type=int, dest='mut_rate', default=50)
    parser.add_argument('-rr', '--recover_rate', type=float, dest='recover_rate', default=0.5)
    parser.add_argument('-s', '--surgery', type=bool, dest='surgery', default=False)
    parser.add_argument('-mp', '--mask_proportion', type=float, dest='mask_proportion', default=0.5)
    parser.add_argument('-d', '--index', type=int, dest='simulation_idx', default=1)
    parser.add_argument('-dir', '--directory', type=str, dest='directory', default='MOONSHOT2')
    parser.add_argument('-nt', '--num_tissue', type=int, dest='num_tissue', default=2)
    parser.add_argument('-nb', '--num_blood', type=int, dest='num_blood', default=2)
    parser.add_argument('-unt', '--used_num_tissue', type=int, dest='used_num_tissue', default=2)
    parser.add_argument('-unb', '--used_num_blood', type=int, dest='used_num_blood', default=2)
    parser.add_argument('-b', '--num_bootstrap', type=int, dest='num_bootstrap', default=100)
    return vars(parser.parse_args(argv))

def run_simulation(argv):
    args = get_args(argv)
    num_node = args['num_node']
    main_dir = args['directory']
    depth_tissue = args['depth_tissue']
    depth_blood = args['depth_blood']
    mut_rate = args['mut_rate']
    recover_rate = args['recover_rate']
    mask_proportion = args['mask_proportion']
    simulation_idx = args['simulation_idx']
    num_tissue = args['num_tissue']
    num_blood = args['num_blood']
    used_num_tissue = args['used_num_tissue']
    used_num_blood = args['used_num_blood']
    num_bootstrap = args['num_bootstrap']
    bool_surgery = args['surgery']
    directory_path = main_dir + f'/boot_simulation_{num_node}_{depth_tissue}_{depth_blood}_{mut_rate}_{int(recover_rate*100)}_{int(mask_proportion*100)}'
    track = 1
    print(track, simulation_idx, directory_path)
    while track <= 1:
        directory_sample = Path(os.path.join(directory_path, f"{num_blood}_{num_tissue}"))
        # if os.path.exists(directory_sample):
        #     track += 1
        #     print(track)
        #     continue
        # else:
        directory_simulation = Path(os.path.join(directory_sample, f"Simulation_{simulation_idx}"))
        os.makedirs(directory_simulation)
        print("start simulation")
        simulation(num_node, depth_tissue, depth_blood, mut_rate, directory_simulation, num_bootstrap, bool_surgery, recover_rate, mask_proportion, num_blood, num_tissue, used_num_blood, used_num_tissue)
        track += 1


################   Main Function Implementation  #################
##################################################################
if __name__ == "__main__":
    run_simulation(sys.argv[1:])
    # num_node = 7
    # depth_tissue = 800
    # depth_blood = 3000
    # mut_rate = 50
    # recover_rate = 0.5
    # mask_proportion = 0.5
    # directory_path = "/data/liquid_biopsy/simulations/mask/multi-tumor/boot_simulation_" + str(num_node) + "_" + str(depth_tissue) + "_" + str(depth_blood) + "_" + str(mut_rate) + "_" + str(int(recover_rate*100)) + "_" + str(int(mask_proportion*100))
    # # This section of lines are for running the main pipeline for 10 simulations;
    # num_tissue = 2
    # num_blood = 2
    # num_bootstrap = 100
    # bool_surgery = False
    #
    # simulation_idx = 1
    # track = 1
    # while track <= 1:
    #     print(track)
    #
    #     directory_simulation = Path(os.path.join(directory_path, f"{num_blood}_{num_tissue}/Simulation_{simulation_idx}"))
    #     if os.path.exists(directory_simulation):
    #         track += 1
    #         continue
    #     else:
    #         os.makedirs(directory_simulation, exist_ok=True)
    #         simulation(num_node, depth_tissue, depth_blood, mut_rate, directory_simulation, num_bootstrap, bool_surgery, recover_rate, mask_proportion, num_blood, num_tissue)
    #         track += 1
    #         #except:
    #             #track = track
