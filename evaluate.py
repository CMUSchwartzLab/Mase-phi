import numpy as np
from ete3 import Tree
#import CASet, DISC
from visualize import *
import glob
import json
import pandas as pd
import seaborn as sns
import sys
from pathlib import Path
sys.path.append("/home/cindyfu/Documents/softwares/stereodist")
sys.path.append("/Users/cindyfu/Documents/softwares/stereodist")
sys.path.append("/home/xuecongf/stereodist")
import CASet, DISC
import matplotlib.pyplot as plt
import gzip
import json
import glob
from zipfile import ZipFile
from analyze import subtree, process_phylowgs_output

def tree2nwk(tree_structure, node_dict):
    cp_dict = tree2cp_dict(tree_structure, node_dict)
    pc_table = cp_dict2pc_table(cp_dict)
    tree = Tree.from_parent_child_table(pc_table)
    tree_nwk = tree.write(format=8, format_root_node=True).replace("..", ",")
    return tree_nwk


def tree2cp_dict(tree_structure, node_dict):
    cp_dict = {}
    for parent, children in tree_structure.items():
        if parent in tree_structure.keys():
            if parent not in node_dict.keys():
                parent_node = ('root',)
            else:
                parent_node = node_dict[parent]
        else:
            parent_node = node_dict[parent]

        for child in children:
            if isinstance(child, int):
                #cp_dict[tuple(node_dict[str(child)])] = tuple(parent_node)
            #else:
                cp_dict[tuple(node_dict[child])] = tuple(parent_node)
    return cp_dict


def co_clustering_matrix(node_list):
    num_mut = 0
    for node, mut_list in node_list.items():
        num_mut += len(mut_list)
    m = np.zeros((num_mut, num_mut))
    for node, mut_list in node_list.items():
        for mut1 in mut_list:
            for mut2 in mut_list:
                m[mut1, mut2] = 1
    return m



def phylowgs2W(mut_assignment, df_ssm_data):
    num_snv = len(df_ssm_data)
    W = np.zeros((len(mut_assignment.keys())+1, num_snv))
    for node in mut_assignment.keys():
        for id_snv in mut_assignment[node]['ssms']:
            W[int(node), int(id_snv[1:])] = 1
    return W

def E2pctable(E):
    edge_list = np.where(E == 1)
    pctable = []
    for i in range(len(edge_list[0])):
        pctable.append((edge_list[0][i], edge_list[1][i]))
    return pctable




def cp_dict2pc_table(cp_dict):
    pc_table = []
    for child, parent in cp_dict.items():
        if len(parent) >1 and len(child) > 1:
            pc_table.append(['{' + '..'.join(str(e) for e in parent) + '}', '{' + '..'.join(str(e) for e in child) + '}', 1])
        elif len(parent) >1:
            pc_table.append(['{' + '..'.join(str(e) for e in parent) + '}', '..'.join(str(e) for e in child), 1])
        elif len(child) > 1:
            pc_table.append([ '..'.join(str(e) for e in parent), '{' + '..'.join(str(e) for e in child) + '}', 1])
        else:
             pc_table.append([ '..'.join(str(e) for e in parent), '..'.join(str(e) for e in child), 1])
    return np.array(pc_table)


class EvaluateTreeDistance:
    def __init__(self, tree_structure_true, node_dict_true, tree_structure_est, node_dict_est):
        self.tree_structure_true = tree_structure_true
        self.tree_structure_est = tree_structure_est
        self.node_dict_true = node_dict_true
        self.node_dict_est = node_dict_est

    def evaluate(self, metric):
        distance = getattr(self, metric)()
        print(distance)
        return distance

    def caset(self):
        metric = CASet.caset_union
        self.tree_nwk_true = tree2nwk(self.tree_structure_true, self.node_dict_true)
        self.tree_nwk_est = tree2nwk(self.tree_structure_est, self.node_dict_est)
        return metric(self.tree_nwk_true, self.tree_nwk_est)

    def disc(self):
        metric = DISC.disc_union
        self.tree_nwk_true = tree2nwk(self.tree_structure_true, self.node_dict_true)
        self.tree_nwk_est = tree2nwk(self.tree_structure_est, self.node_dict_est)
        return metric(self.tree_nwk_true, self.tree_nwk_est)


def generate_mutation_tree(pctable, W, if_start0, l):
    mutation_tree = {}
    for parent, child in pctable:
        if if_start0:
            parent_true = parent
            child_true = child
        else:
            parent_true = parent - 1
            child_true = child - 1
        if len(np.where(W[parent_true] == 1)[0]) == 0:
            parent_mut = ('root',)
        else:
            parent_mut = tuple('sv' + str(i) if i < l else 'snv' + str(i) for i in np.where(W[parent_true] == 1)[0])
        child_mut = tuple('sv' + str(i) if i < l else 'snv' + str(i) for i in np.where(W[child_true] == 1)[0])
        if parent_mut not in mutation_tree.keys():
            mutation_tree[parent_mut] = []
        mutation_tree[parent_mut].append(child_mut)
    return mutation_tree


if __name__ == "__main__":
    df_dist = []
    num_node = 10
    depth_tissue = 800
    depth_blood = 3000
    mut_rate = 50
    recover_rate = 0.5
    mask_proportion = 0.5
    # deconv parameter
    num_node_deconv = 5
    num_iter = 20
    time_limit = 1000
    num_restarts = 5
    deconv_subdir = f"{num_node_deconv}_{num_iter}_{time_limit}_{num_restarts}_0"
    '''
    directory_path = "sim/simulations/mask/simulation_" + str(num_node) + "_" + str(depth_tissue) + "_" + str(depth_blood) + "_" + str(mut_rate) + "_" + str(int(recover_rate*100)) + "_" + str(int(mask_proportion*100))

    for simulation in glob.glob(directory_path + "/*"):
        simu_idx = int(simulation.split("_")[-1])
        excel_link = simulation + "/simulated_mut.xlsx"
        inter = pd.read_excel(excel_link, sheet_name='common_blood_tissue_no_germline', index_col=0)
        with open(simulation + "/subtree.json", 'r') as f:
             data = json.load(f)
        subtree_structure_true = data["subtree_structure"]
        subnode_dict_true = data["subnode_dict"]
        with open(simulation + "/common/result_phylowgs.json", 'r') as f:
             data = json.load(f)
        subtree_structure_phylowgs = data["tree_structure"]
        subnode_dict_phylowgs = data["node_dict"]

        with open(simulation + "/deconv_inter/" + deconv_subdir + "/results.json", 'r') as f:
             data = json.load(f)
        subtree_structure_deconv = data["tree_structure"]
        subnode_dict_deconv = data["node_dict"]

        dist_phylowgs = EvaluateTreeDistance(subtree_structure_true, subnode_dict_true, subtree_structure_phylowgs, subnode_dict_phylowgs)
        dist_deconv = EvaluateTreeDistance(subtree_structure_true, subnode_dict_true, subtree_structure_deconv, subnode_dict_deconv)

        df_dist.append((dist_phylowgs.caset(), "PhyloWGS", "CASet", num_node, depth_tissue, depth_blood, mut_rate, recover_rate, mask_proportion))
        df_dist.append((dist_deconv.caset(), "deconv", "CASet", num_node, depth_tissue, depth_blood, mut_rate, recover_rate, mask_proportion))
        df_dist.append((dist_phylowgs.disc(), "PhyloWGS", "DISC", num_node, depth_tissue, depth_blood, mut_rate, recover_rate, mask_proportion))
        df_dist.append((dist_deconv.disc(), "deconv", "DISC", num_node, depth_tissue, depth_blood, mut_rate, recover_rate, mask_proportion))

    df_dist = pd.DataFrame(df_dist, columns=['distance', 'method', 'metric', 'num node', 'depth tissue', 'depth blood', 'mutation rate',
                                                 'recover rate', 'mask proportion'])
    ax = sns.catplot(data=df_dist, x='mutation rate', y='distance', col='metric', kind='boxen', hue='method')
    plt.show()
    # # This section of lines are for running the main pipeline for 10 simulations;
    # for simulation in glob.glob(directory_path + "/*"):
    #     excel_link = simulation + "/simulated_mut.xlsx"
    #     inter = pd.read_excel(excel_link, sheet_name='common_blood_tissue_no_germline', index_col=0)
    #     submut_list = list(inter.Gene)
    #     with open(simulation + "/intermediate_value.json", 'r') as f:
    #          data = json.load(f)
    #     tree_structure = data["tree"]
    #     node_dict = data["Mutation"]
    #     subtree_structure, subnode_dict = subtree(tree_structure, node_dict, submut_list)
    #     subdata = {"subtree_structure": subtree_structure, "subnode_dict": subnode_dict}
    #     with open(simulation + "/subtree.json", 'w') as ff:
    #         json.dump(subdata, ff)
    #     g = render_tumor_tree(subtree_structure, subnode_dict)
    #     g.render(filename=simulation + '/true_tree_common')
    '''
    num_tissue=3
    num_blood=1
    num_chains=5
    directory_path = "sim/simulations/mask/multi-tumor/simulation_" + str(num_node) + "_" + str(depth_tissue) + "_" + str(depth_blood) + "_" + str(mut_rate) + "_" + str(int(recover_rate*100)) + "_" + str(int(mask_proportion*100)) + "/" + str(num_blood) + "_" + str(num_tissue)
    name='common'
    directory_path = Path(directory_path)
    for simulation in directory_path.glob("Simulation*"):
        print(simulation)
        simu_idx = simulation.stem.split("_")[-1]
        excel_link = simulation / f"simulated_mut.xlsx"
        inter = pd.read_excel(excel_link, sheet_name='common_blood_tissue_no_germline', index_col=0)
        with open(simulation / f"intermediate_value.json", 'r') as f:
            data = json.load(f)
        tree_structure = data["tree"]
        node_dict = data["Mutation"]
        print(node_dict)
        submut_list = list(inter["Gene"])
        subtree_structure, subnode_dict = subtree(tree_structure, node_dict, submut_list)
        subdata = {"subtree_structure": subtree_structure, "subnode_dict": subnode_dict}
        with open(simulation / f"subtree.json", 'w') as ff:
            json.dump(subdata, ff)
        with open(simulation / f"subtree.json", 'r') as f:
            data = json.load(f)
        subtree_structure_true = {key:[int(it) for it in value] for key, value in data["subtree_structure"].items()}
        subnode_dict_true = data["subnode_dict"]

        ### generate result for phylowgs outputs:
        for t_idx in range(1, num_tissue+1):
            file1=simulation
            print(file1)
            patient_num = int(file1.stem.split('_')[-1])
            print(patient_num)
            file = file1 / f"{num_blood}_{t_idx}"
            summ_file = file / f"{name}/result_{num_chains}.summ.json.gz"
            muts_file = file / f"{name}/result_{num_chains}.muts.json.gz"
            mutass_file = file / f"{name}/result_{num_chains}.mutass.zip"
            tree_structure, node_dict, node_dict_name, node_dict_re, final_tree_cp, prev_mat, clonal_freq, vaf_frac = process_phylowgs_output(summ_file, muts_file, mutass_file)


            phylowgs_dict = {"tree_structure": tree_structure, "node_dict": node_dict, "node_dict_name": node_dict_name}
            with open(file / f"{name}/result_phylowgs_{num_chains}.json", 'w') as ff:
                json.dump(phylowgs_dict, ff)


            with open(file / f"{name}/result_phylowgs_{num_chains}.json", 'r') as f:
                data = json.load(f)
            subtree_structure_phylowgs = data["tree_structure"]
            subnode_dict_phylowgs = data["node_dict_name"]
            g = render_tumor_tree(subtree_structure_phylowgs, subnode_dict_phylowgs)
            g.render(filename=file / f"{name}/tree_common_phylowgs")

            with open(file / f"deconv_inter/{deconv_subdir}/results.json", 'r') as f:
                 data = json.load(f)
            subtree_structure_deconv = data["tree_structure"]
            subnode_dict_deconv = data["node_dict_name"]
            g = render_tumor_tree(subtree_structure_deconv, subnode_dict_deconv)
            g.render(filename=file / f"deconv_inter/{deconv_subdir}/tree_common_phylowgs")


            dist_phylowgs = EvaluateTreeDistance(subtree_structure_true, subnode_dict_true, subtree_structure_phylowgs, subnode_dict_phylowgs)
            dist_deconv = EvaluateTreeDistance(subtree_structure_true, subnode_dict_true, subtree_structure_deconv, subnode_dict_deconv)

            df_dist.append((dist_phylowgs.caset(), "PhyloWGS", "CASet", num_node, depth_tissue, depth_blood, mut_rate, recover_rate, mask_proportion, t_idx))
            df_dist.append((dist_deconv.caset(), "deconv", "CASet", num_node, depth_tissue, depth_blood, mut_rate, recover_rate, mask_proportion, t_idx))
            df_dist.append((dist_phylowgs.disc(), "PhyloWGS", "DISC", num_node, depth_tissue, depth_blood, mut_rate, recover_rate, mask_proportion, t_idx))
            df_dist.append((dist_deconv.disc(), "deconv", "DISC", num_node, depth_tissue, depth_blood, mut_rate, recover_rate, mask_proportion, t_idx))

    df_dist = pd.DataFrame(df_dist, columns=['distance', 'method', 'metric', 'num node', 'depth tissue', 'depth blood',
                                             'mutation rate', 'recover rate', 'mask proportion', 'num tissue'])
    ax = sns.catplot(data=df_dist, x='method', y='distance', col='metric', kind='boxen', hue='num tissue')
    plt.savefig("tree_accuracy.eps")
    plt.show()
    # # This section of lines are for running the main pipeline for 10 simulations;
    # for simulation in glob.glob(directory_path + "/*"):
    #     excel_link = simulation + "/simulated_mut.xlsx"
    #     inter = pd.read_excel(excel_link, sheet_name='common_blood_tissue_no_germline', index_col=0)
    #     submut_list = list(inter.Gene)
    #     with open(simulation + "/intermediate_value.json", 'r') as f:
    #          data = json.load(f)
    #     tree_structure = data["tree"]
    #     node_dict = data["Mutation"]
    #     subtree_structure, subnode_dict = subtree(tree_structure, node_dict, submut_list)
    #     subdata = {"subtree_structure": subtree_structure, "subnode_dict": subnode_dict}
    #     with open(simulation + "/subtree.json", 'w') as ff:
    #         json.dump(subdata, ff)
    #     g = render_tumor_tree(subtree_structure, subnode_dict)
    #     g.render(filename=simulation + '/true_tree_common')
