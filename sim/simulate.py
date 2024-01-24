### Simulation file

from collections import deque
from copy import deepcopy
import numpy as np
import pandas as pd
from pathlib import Path
import json
from scipy.stats import truncnorm
import random

from bootstrap import phyloWGS_output


# simulation, read count
def ddPCR(freq_blood_dict_list, volume=20, cycle=40, amp_eff=0.95):
    """
        This function is used for simulating the ddPCR process
        on blood samples with controls over the volume and 
        the number of cycle for the sample;

        Volume: the sample volume (default will be 20ul);
        Cycle: the number of PCR cycle this sample went through;
        amp_eff: 

    """
    ddPCR_list = []

    for blood_dict in freq_blood_dict_list:
        new_dict = {key: 0 for key in blood_dict.keys()}
        ddPCR_list.append(new_dict)

        accum_freq = [] # make a list of accumulated freq for each time point of blood sample;
        freq_tracker = 0
        for key, value in blood_dict.items():
            freq_tracker += value
            accum_freq.append(freq_tracker)

        for droplet in range(volume*1000):
            sample_dice = random.random()
            tracker = 0
            while accum_freq[tracker] < sample_dice: # check which range does the random number drop into;
                tracker += 1
            
            exp_num = 0
            amp_dice = random.random()
            for cycle_num in range(cycle):
                
                if amp_eff > amp_dice: # check whether this cycle of amplification is successful;
                    exp_num += 1 # if successful, expoential number +1, otherwise stop;
            
            ddPCR_list[-1][tracker] += 2^exp_num
    
    return ddPCR_list


def simulate_clonal_tree(d=3, k=5, a=6, b=6, recover_rate=0.5, growth_rate=0.1, num_blood=1, num_tissue=2):
    """
    input:
    d: maximum degree of each node
    k: number of nodes
    a, b: the parameter for the beta distribution
    recover_rate: the parameter for the change of VAF after surgery
    growth_rate: the parameter for the change of VAF due to regrowth of tumour 

    output:
    tree: dict represents tree structure
    freq_obs: frequencies for each node
    """

    tree = {}
    root = 0
    nodes = [root]
    nodes_deque = deque([root])
    child_node = root
    freq_true = np.random.beta(a, b, k)
    freq_true = freq_true.tolist()
    freq_obs_tissue = np.random.dirichlet(freq_true, num_tissue)
    freq_obs_tissue_dict_list = []
    for tissue_id in range(num_tissue):
        freq_obs_tissue_dict = {i: freq_obs_tissue[tissue_id, i] for i in range(freq_obs_tissue.shape[1])}
        freq_obs_tissue_dict_list.append(freq_obs_tissue_dict)
    freq_obs_blood = np.random.dirichlet(freq_true)
    freq_obs_blood[0] = 0
    freq_obs_blood_dict_list = []
    new_freq_obs_blood_normal = np.random.rand() / 10 + 0.9
    freq_obs_blood_dict = {i: freq_obs_blood[i] / np.sum(freq_obs_blood) * (
                1 - new_freq_obs_blood_normal) if i != 0 else new_freq_obs_blood_normal for i in
                           range(len(freq_obs_blood))}
    freq_obs_blood_dict_list.append(freq_obs_blood_dict)

    ## blood draw after surgery
    freq_obs_blood_dict_post = deepcopy(freq_obs_blood_dict)

    # After surgery, all nodes freq exceptnormal cell will decrease randomly;
    # The normal node freq will increase accordingly by the sum of all decrease in the rest of nodes;
    tracker = k - 1
    total_increment = 0
    while tracker != 0:
        # freq_blood_decrease = np.random.binomial(n=100, p=recover_rate)/100*freq_obs_blood_dict_post[tracker]
        mean = freq_obs_blood_dict[tracker]
        sd = mean * recover_rate
        freq_blood_decrease = truncnorm.rvs(-mean / sd, 0, loc=mean, scale=sd, size=None)
        freq_obs_blood_dict_post[tracker] = freq_obs_blood_dict_post[tracker] - freq_blood_decrease
        total_increment += freq_blood_decrease
        tracker -= 1

    freq_obs_blood_dict_post[0] += total_increment

    if num_blood == 1:
        freq_obs_blood_dict_list = [freq_obs_blood_dict]
    elif num_blood == 2:
        freq_obs_blood_dict_list = [freq_obs_blood_dict, freq_obs_blood_dict_post]
    else:
        freq_obs_blood_dict_list = [freq_obs_blood_dict, freq_obs_blood_dict_post]
        for tp in range(0, (num_blood-2)):
            # each time point after sugery, append a dict in the dict_list
            freq_obs_blood_dict_list.append({})
            total_freq = 1
            for node in range(1,k):
                # iterate over every mutant nodes:
                #new_tumor_growth = (freq_obs_blood_dict_list[0][node]-freq_obs_blood_dict_list[1][node]) * growth_rate
                new_tumor_growth = freq_obs_blood_dict_list[1][node] * growth_rate
                print("before: ", new_tumor_growth)

                new_tumor_growth = np.random.normal(loc=new_tumor_growth, scale=new_tumor_growth/3, size=None)
                # Add randomization to the new tumor growth on a normalization curve;
                print("after: ", new_tumor_growth)

                freq_obs_blood_dict_list[-1][node] = freq_obs_blood_dict_list[-2][node] + new_tumor_growth
                total_freq = freq_obs_blood_dict_list[-2][0] - new_tumor_growth
                freq_obs_blood_dict_list[-1][0] =  total_freq
                freq_obs_blood_dict_list[-1] = {k: freq_obs_blood_dict_list[-1][k] for k in sorted(freq_obs_blood_dict_list[-1])}
        

    print("True Freq:\n", freq_true, "\n")
    ### simulate the ground truth tumor tree
    end_loop = False
    while not end_loop:
        current_node = nodes_deque.popleft()
        while True:
            if current_node == root:
                child_num = np.random.choice(range(1, d + 1))
            else:
                child_num = np.random.choice(d + 1)
            if child_num + len(nodes) < k:
                break
            elif child_num + len(nodes) == k:
                end_loop = True
                break
        if current_node not in tree.keys() and child_num != 0:
            tree.setdefault(current_node, [])
        for child in range(child_num):
            child_node += 1
            nodes.append(child_node)
            nodes_deque.append(child_node)
            tree[current_node].append(child_node)
    print("Tree:\n", tree, "\n")
    return tree, freq_true, freq_obs_tissue_dict_list, freq_obs_blood_dict_list
####TODO: fix the dict vs dict_list for blood

def bfs(root, tree):
    """
        This function performs the breadth first search on the root and its given tree; 
    """
    order = []
    q = deque([root])
    while len(q) != 0:
        node = q.popleft()
        order.append(node)
        if node in tree.keys():
            for child in tree[node]:
                q.append(child)
    return order


## simulate the frequency of each mutation by summing up the frequencies of itself and children
def cumulative_Freq(root, tree, freq_list):
    """
        This function used bfs func and given freq-list for both blood and tissue to calculate
        the given culmutative freq of the child nodes for each parent node;
    """
    bfs_order = bfs(root, tree)
    freq_sum = {}
    for node_idx in range(len(bfs_order) - 1, -1, -1):
        node = bfs_order[node_idx]
        if node not in tree.keys():
            freq_sum[node] = freq_list[node]
        else:
            freq_sum[node] = sum([freq_sum[child_node] for child_node in tree[node]]) + freq_list[node]
    return {i: freq_sum[i] for i in freq_sum.keys()}


def mask_recent_clones_tissue(tree, freq_obs_tissue_dict_list, root, mask_proportion=0.5):
    total_nodes_list = bfs(root, tree)
    num_mask = np.round(mask_proportion * len(total_nodes_list))
    mask_nodes_list = []
    while True:
        mask_node = np.random.choice(total_nodes_list)
        clade = bfs(mask_node, tree)
        num_mask -= len(clade)
        if num_mask == 0:
            mask_nodes_list += clade
            total_nodes_list = list(set(total_nodes_list).difference(set(clade)))
            break
        elif num_mask > 0:
            mask_nodes_list += clade
            total_nodes_list = list(set(total_nodes_list).difference(set(clade)))
        else:
            num_mask += len(clade)
            continue
    masked_freq_total = 0
    freq_obs_tissue_dict_new = {}
    freq_obs_tissue_dict_new_list = [{} for idx in range(len(freq_obs_tissue_dict_list))]
    masked_freq_total = np.zeros((len(freq_obs_tissue_dict_list)))
    for node in freq_obs_tissue_dict_list[0].keys():
        if node in mask_nodes_list:
            for freq_obs_tissue_dict_idx in range(len(freq_obs_tissue_dict_list)):
                masked_freq_total[freq_obs_tissue_dict_idx] += freq_obs_tissue_dict_list[freq_obs_tissue_dict_idx][node]
                freq_obs_tissue_dict_new_list[freq_obs_tissue_dict_idx][node] = 0.0
    for node in freq_obs_tissue_dict_list[0]:
        if node in total_nodes_list:
            for freq_obs_tissue_dict_idx in range(len(freq_obs_tissue_dict_list)):
                freq_obs_tissue_dict_new_list[freq_obs_tissue_dict_idx][node] = freq_obs_tissue_dict_list[freq_obs_tissue_dict_idx][node] / (1 - masked_freq_total[freq_obs_tissue_dict_idx])
    return freq_obs_tissue_dict_new_list


def simulate_variant_reads(tree, freq_tissue_list, freq_blood_list, depth_tissue=50, depth_blood=50,
                           mutation_rate=10):
    '''
        Input:
        tree: dict represents tree structure
        freq: frequencies for each node
        depth_SAMPLE: a poisson mean of sequencing depth
        mutation_rate: a poisson mean of number of mutations

        simulate a set of mutations for each node, the number of mutations holds a poisson distribution
        simulate the reference read counts which holds a poisson distribution
        and variant read counts for each mutation which holds a binomial distribution with the parameter of frequencies

        output:
        node_mutations: dict for tree node - mutation_list pair, i.e. {0: [0,1,2], 1:[3,4], 2:[5,6,7]}
        mutation_refer_count: dict for mutation - reference read counts
        mutation_variant_count: dict for mutation - variant read counts
    '''

    # Use bfs and cumulative_Freq to calculate out the cumulative freq for each node:
    freq_tissue_list = [cumulative_Freq(0, tree, freq_tissue_list[i]) for i in range(len(freq_tissue_list))]
    freq_blood_list = [cumulative_Freq(0, tree, freq_blood_list[i]) for i in range(len(freq_blood_list))]

    # Because genes have two copies of allele, it is important to divide all freq by 2:
    freq_tissue_list = [{k: v / 2 for k, v in freq_tissue.items()} for freq_tissue in freq_tissue_list]
    freq_blood_list = [{k: v / 2 for k, v in freq_blood.items()} for freq_blood in freq_blood_list]

    print("Tissue AFreq:\n", freq_tissue_list, "\n")
    print("Blood AFreq:\n", freq_blood_list, "\n")


    mutation_list = np.random.poisson(mutation_rate, (len(freq_tissue_list[0]) - 1))
    # Set a minimum number of mutations for each node to be 1:

    for i in range(len(mutation_list)):
        if mutation_list[i] == 0:
            mutation_list[i] = 1

    mutation_list = np.insert(mutation_list, 0, 0)
    # Number of Mutations in each node;
    # Eg: [12,11,13,8,9,10,11]

    depth_list_T_list = [np.random.poisson(depth_tissue, mutation_list.sum()) for i in range(len(freq_tissue_list))]
    depth_list_B_list = [np.random.poisson(depth_blood, mutation_list.sum()) for i in range(len(freq_blood_list))]
    # Create Tissue and Blood list of depth with length equals to the sum of all mutations;

    variant_list_T_list = [np.zeros(mutation_list.sum(), dtype=int) for i in range(len(freq_tissue_list))]
    reference_list_T_list = [np.zeros(mutation_list.sum(), dtype=int) for i in range(len(freq_tissue_list))]

    variant_list_B_list = [np.zeros(mutation_list.sum(), dtype=int) for i in range(len(freq_blood_list))]
    reference_list_B_list = [np.zeros(mutation_list.sum(), dtype=int) for i in range(len(freq_blood_list))]


    tracker = 0
    node_mutation = {}
    mutation_variant_T_list = [{} for i in range(len(freq_tissue_list))]
    mutation_reference_T_list = [{} for i in range(len(freq_tissue_list))]
    mutation_variant_B_list = [{} for i in range(len(freq_blood_list))]
    mutation_reference_B_list = [{} for i in range(len(freq_blood_list))]

    # iterate over the mutation list for each node and choose the depth
    # and perform binomial distribution based on the value and number of mutations
    for node in freq_tissue_list[0].keys():
        num_mutation = mutation_list[node]
        start = tracker
        tracker += num_mutation
        end = tracker
        node_mutation[node] = [i for i in range(start, end)]

        for tissue_idx in range(len(freq_tissue_list)):
            variant_list_T_list[tissue_idx][start:end] = np.random.binomial(depth_list_T_list[tissue_idx][start:end], freq_tissue_list[tissue_idx][node])
            reference_list_T_list[tissue_idx][start:end] = np.subtract(depth_list_T_list[tissue_idx][start:end], variant_list_T_list[tissue_idx][start:end])
        for blood_idx in range(len(freq_blood_list)):
            variant_list_B_list[blood_idx][start:end] = np.random.binomial(depth_list_B_list[blood_idx][start:end], freq_blood_list[blood_idx][node])
            reference_list_B_list[blood_idx][start:end] = np.subtract(depth_list_B_list[blood_idx][start:end], variant_list_B_list[blood_idx][start:end])


        # Depth :  [50,51,53,46,47,49...]
        # Variant: [ 2, 3, 4, 5, 1, 2...]

    for i in range(mutation_list.sum()):
        # assign TISSUE and BLOOD variant and refer reads for each mutation dictionary;
        for tissue_idx in range(len(freq_tissue_list)):
            mutation_variant_T_list[tissue_idx][i] = variant_list_T_list[tissue_idx][i]
            mutation_reference_T_list[tissue_idx][i] = reference_list_T_list[tissue_idx][i]

        for blood_idx in range(len(freq_blood_list)):
            mutation_variant_B_list[blood_idx][i] = variant_list_B_list[blood_idx][i]
            mutation_reference_B_list[blood_idx][i] = reference_list_B_list[blood_idx][i]


    # Converting all dictionary's value into JSON serializable;
    mutation_variant_T_list = [{k: int(v) for k, v in mutation_variant_T.items()} for mutation_variant_T in mutation_variant_T_list]
    mutation_reference_T_list = [{k: int(v) for k, v in mutation_reference_T.items()} for mutation_reference_T in mutation_reference_T_list]
    mutation_variant_B_list = [{k: int(v) for k, v in mutation_variant_B.items()} for mutation_variant_B in mutation_variant_B_list]
    mutation_reference_B_list = [{k: int(v) for k, v in mutation_reference_B.items()} for mutation_reference_B in mutation_reference_B_list]

    return node_mutation, mutation_variant_T_list, mutation_reference_T_list, mutation_variant_B_list, mutation_reference_B_list


def simulation_phyloWGS(v_blood, a_blood, v_tissue, a_tissue, output_path):
    """
        This function used the simulated output data of mutations
        to write them back to the plyWGS format can be ready to 
        run in the plyWGS program;
    """
    d_blood = dict([(k, v_blood[k] + a_blood[k]) for k in set(v_blood) & set(a_blood)])
    d_tissue = dict([(k, v_tissue[k] + a_tissue[k]) for k in set(v_tissue) & set(a_tissue)])
    # The total number of read depth for blood and tissue samples;

    plyWGS_format = []
    count = 0
    for i in range(len(a_blood)):
        # Iterate over all simulated mutations
        plyWGS_format.append({'id': 's' + str(count), 'gene': f"mut_{i}", 'a': str(a_blood[i]) + ',' + str(a_tissue[i]),
                              'd': str(d_blood[i]) + ',' + str(d_tissue[i]), 'mu_r': 0.999, 'mu_v': 0.499})
        count += 1

    df_plyWGS = pd.DataFrame(plyWGS_format)
    df_plyWGS.to_csv(output_path / f'ssm_data.txt', index=False, sep='\t')


def simulation_excel(v_blood_list, a_blood_list, v_tissue_list, a_tissue_list, output_path, used_num_blood=1, used_num_tissue=1, surgery=False):
    '''
        This function use the simulated output data of mutations
        to write them back to an excel that can be ready to perform
        further implementations;
    '''
    if surgery:
        sample_type = "_post"
    else:
        sample_type = ""

    d_tissue_list = []
    freq_tissue_list = []
    for i, (v_tissue, a_tissue) in enumerate(zip(v_tissue_list, a_tissue_list)):
        d_tissue = dict([(k, v_tissue[k] + a_tissue[k]) for k in set(v_tissue) & set(a_tissue)])
        freq_tissue = dict([(k, v_tissue[k] / d_tissue[k]) for k in set(v_tissue) & set(d_tissue)])
        d_tissue_list.append(d_tissue)
        freq_tissue_list.append(freq_tissue)
    d_blood_list = []
    freq_blood_list = []
    for i, (v_blood, a_blood) in enumerate(zip(v_blood_list, a_blood_list)):
        d_blood = dict([(k, v_blood[k] + a_blood[k]) for k in set(v_blood) & set(a_blood)])
        freq_blood = dict([(k, v_blood[k] / a_blood[k]) for k in set(v_blood) & set(d_blood)])
        d_blood_list.append(d_blood)
        freq_blood_list.append(freq_blood)
    excel_format_union = []
    excel_format_blood_list = [[] for i in range(len(freq_blood_list))]
    excel_format_tissue_list = [[] for i in range(len(freq_tissue_list))]
    count = 0
    for i in range(len(a_blood_list[0])):
        # Iterate over all simulated mutations
        for blood_idx in range(len(freq_blood_list)):
            freq_blood = freq_blood_list[blood_idx]
            d_blood = d_blood_list[blood_idx]
            if freq_blood[i] > 0.0:
                excel_format_blood_list[blood_idx].append({
                                 'Gene': f"mut_{i}",
                                 'Chromosome': '',
                                 'Genomic Position': '',
                                 'Reference Call': '',
                                 'Alternative Call': '',
                                 'Allele Frequency_b{}'.format(blood_idx): str(freq_blood[i]),
                                 'Depth_b{}'.format(blood_idx): str(d_blood[i]),
                                 'P-Dot Notation': '',
                                 'C-Dot Notation': '',
                                 'Consequence(s)': '',
                                 })
        for tissue_idx in range(len(freq_tissue_list)):
            freq_tissue = freq_tissue_list[tissue_idx]
            d_tissue = d_tissue_list[tissue_idx]
            if freq_tissue[i] > 0.0:
                excel_format_tissue_list[tissue_idx].append({
                                 'Gene': f"mut_{i}",
                                 'Chromosome': '',
                                 'Genomic Position': '',
                                 'Reference Call': '',
                                 'Alternative Call': '',
                                 'Allele Frequency_t{}'.format(tissue_idx): str(freq_tissue[i]),
                                 'Depth_t{}'.format(tissue_idx): str(d_tissue[i]),
                                 'P-Dot Notation': '',
                                 'C-Dot Notation': '',
                                 'Consequence(s)': '',})
            count += 1


    df_excel_blood_list = [pd.DataFrame(excel_format_blood) for excel_format_blood in excel_format_blood_list]
    df_excel_tissue_list = [pd.DataFrame(excel_format_tissue) for excel_format_tissue in excel_format_tissue_list]
    df_excel_tissue_common = df_excel_tissue_list[0].copy(deep=True)
    #df_excel_tissue_common.rename(columns={'Allele Frequency': 'Allele Frequency_t0', 'Depth': 'Depth_t0'}, inplace=True)
    for i in range(1, used_num_tissue):
        df_excel_tissue_common = pd.merge(df_excel_tissue_common, df_excel_tissue_list[i], how='inner', on=['Gene', 'Chromosome', 'Genomic Position'])
        #df_excel_tissue_common.rename(columns={'Allele Frequency': 'Allele Frequency_t' + str(i), 'Depth': 'Depth_t' + str(i)}, inplace=True)
    df_excel_tissue_union = df_excel_tissue_list[0].copy(deep=True)
    df_excel_tissue_union.rename(columns={'Allele Frequency': 'Allele Frequency_t0', 'Depth':'Depth_t0'}, inplace=True)
    for i in range(1, used_num_tissue):
        df_excel_tissue_union = pd.merge(df_excel_tissue_union, df_excel_tissue_list[i], how='outer',
                                                              on=['Gene', 'Chromosome', 'Genomic Position'])
        #df_excel_tissue_union.rename(columns={'Allele Frequency': 'Allele Frequency_t' + str(i), 'Depth': 'Depth_t' + str(i)},inplace=True)
    print(df_excel_tissue_union.columns)

    df_excel_blood_common = df_excel_blood_list[0].copy(deep=True)
    #df_excel_blood_common.rename(columns={'Allele Frequency': 'Allele Frequency_b0', 'Depth': 'Depth_b0'}, inplace=True)
    for i in range(1, used_num_blood):
        df_excel_blood_common = pd.merge(df_excel_blood_common, df_excel_blood_list[i], how='inner', on=['Gene', 'Chromosome', 'Genomic Position'])
        df_excel_blood_common.rename(columns={'Allele Frequency': 'Allele Frequency_b' + str(i), 'Depth': 'Depth_b' + str(i)}, inplace=True)
    df_excel_blood_union = df_excel_blood_list[0].copy(deep=True)
    df_excel_blood_union.rename(columns={'Allele Frequency': 'Allele Frequency_b0', 'Depth':'Depth_b0'}, inplace=True)
    for i in range(1, used_num_blood):
        df_excel_blood_union = pd.merge(df_excel_blood_union, df_excel_blood_list[i], how='outer',
                                                              on=['Gene', 'Chromosome', 'Genomic Position'])
        df_excel_blood_union.rename(columns={'Allele Frequency': 'Allele Frequency_b' + str(i), 'Depth': 'Depth_b' + str(i)},inplace=True)
    print(df_excel_blood_union.columns)
    df_excel_union = pd.merge(df_excel_blood_union, df_excel_tissue_union, how='outer', on=['Gene', 'Chromosome', 'Genomic Position'])
    df_excel_common = pd.merge(df_excel_blood_common, df_excel_tissue_common, how='inner', on=['Gene', 'Chromosome', 'Genomic Position'])

    # writer = pd.ExcelWriter(output_path/f"simulated_mut{surgery}.xlsx", engine='xlsxwriter')
    with pd.ExcelWriter(output_path / f"simulated_mut{sample_type}.xlsx") as writer:
        for i in range(len(freq_tissue_list)):
            df_excel_tissue_list[i].to_excel(writer, sheet_name='tissue_no_germline_' + str(i))
        df_excel_tissue_common.to_excel(writer, sheet_name='common_tissue_no_germline')
        for i in range(len(freq_blood_list)):
            df_excel_blood_list[i].to_excel(writer, sheet_name='blood_no_germline_' + str(i))
        df_excel_blood_common.to_excel(writer, sheet_name='common_blood_no_germline')
        df_excel_common.to_excel(writer, sheet_name='common_blood_tissue_no_germline')
        df_excel_union.to_excel(writer, sheet_name='union_blood_tissue_no_germline')

    # writer.close()
    # The total number of read depth for blood and tissue samples;


##################################################################
################   Main Function Implementation  #################
##################################################################
if __name__ == "__main__":
    tree, freq, freq_T, freq_B = simulate_clonal_tree(k=5, growth_rate=0.1, num_blood=2)
    freq_T_masked=mask_recent_clones_tissue(tree, freq_T, 0)
    print(freq_T, freq_B)
    for dict in freq_B:
        print(dict,"\n")

    print(ddPCR(freq_B))
