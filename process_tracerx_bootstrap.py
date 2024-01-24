from pathlib import Path
import pandas as pd
import pickle
from visualize import *
from analyze import *
from optimize import *

patient='CRUK0041'
num_blood=0
num_tissue=5
type='common'
directory = Path(f'/home/xuecongf/liquid_biopsy/abbosh/bootstrap')
file = directory / f'{patient}/{num_blood}_{num_tissue}/{type}'
bootstrap_num = 100
method = 'phylowgs'
num_chain=5
#patient_num = int(file.stem.split('_')[0])
#type = file.stem.split('_')[-1]
directory_xlsx = directory / f'{patient}_bootstrap.xlsx'
df = pd.read_excel(directory_xlsx, index_col=0, sheet_name='common_blood_tissue_no_germline')
gene2idx = {'s' + str(i): i for i in range(len(df))}
idx2gene = {i: 's' + str(i) for i in range(len(df))}
#mastro_format = []
tree_distribution = {'cp_tree': [], 'node_dict': [], 'node_dict_name': [], 'node_dict_re': [], 'tree_structure': [], 'freq': [], 'clonal_freq': [], 'vaf_frac': []}
tree_aggregation = {'cp_tree': [], 'node_dict': [], 'node_dict_name': [], 'node_dict_re': [], 'tree_structure': [], 'freq': [], 'clonal_freq': [], 'vaf_frac': []}
# f = open(f'mastro_input_{patient_num}.txt', 'w'):
for bootstrap_idx in range(1, bootstrap_num + 1, 1):
    print(bootstrap_idx)
    if bootstrap_idx == 16:
        continue
    summ_file = file / f"bootstrap{bootstrap_idx}/result_{num_chain}.summ.json.gz"
    muts_file = file / f"bootstrap{bootstrap_idx}/result_{num_chain}.muts.json.gz"
    mutass_file = file / f"bootstrap{bootstrap_idx}/result_{num_chain}.mutass.zip"
    tree_structure, node_dict, node_dict_name, node_dict_re, final_tree_cp, prev_mat, clonal_freq, vaf_frac = process_phylowgs_output(summ_file, muts_file, mutass_file)
    #print(tree_structure)
    #sc_mat = create_same_clone_matrix(tree_structure, node_dict, gene2idx)
    #ad_mat = create_ancestor_descendant_matrix(tree_structure, node_dict, gene2idx)
    #print(sc_mat, ad_mat)
    #relation_list = []
    # write mastro input
    # for i in range(sc_mat.shape[0]):
    #     for j in range(i+1, sc_mat.shape[1]):
    #         if sc_mat[i, j] == 0 and ad_mat[i, j] == 1:
    #             relation_list.append(idx2gene[i] + '->-' + idx2gene[j])
    #         if sc_mat[i, j] == 1 and ad_mat[i, j] == 0:
    #             relation_list.append(idx2gene[i] + '-?-' + idx2gene[j])
    #         if sc_mat[i, j] == 0 and ad_mat[j, i] == 1:
    #             relation_list.append(idx2gene[j] + '->-' + idx2gene[i])
    #         if sc_mat[i, j] == 0 and ad_mat[i, j] == 0:
    #             relation_list.append(idx2gene[j] + '-/-' + idx2gene[i])
    # print('.'.join(relation_list), file=f)
    tree_distribution = combine_tree(node_dict, node_dict_name, node_dict_re, tree_structure, final_tree_cp, clonal_freq, vaf_frac, 'phylowgs', tree_distribution)
    tree_aggregation['tree_structure'].append(tree_structure)
    tree_aggregation['cp_tree'].append(final_tree_cp)
    tree_aggregation['node_dict'].append(node_dict)
    tree_aggregation['node_dict_re'].append(node_dict_re)
    tree_aggregation['node_dict_name'].append(node_dict_name)
    tree_aggregation['freq'].append(1)
    tree_aggregation['clonal_freq'].append(clonal_freq)
    tree_aggregation['vaf_frac'].append(vaf_frac)

# f.close()
analyze_tree_distribution(tree_distribution, directory, patient, type, fig=True)
best_bootstrap_idx = np.argmax(tree_distribution['freq'])
tree_structure = tree_distribution['tree_structure'][best_bootstrap_idx]
node_dict = tree_distribution['node_dict'][best_bootstrap_idx]
results_dict = {'node_dict_name': node_dict,'tree_structure': tree_structure}
with open(directory / f"{patient}_results_bootstrap_{type}_best.json", 'w') as f:
    json.dump(results_dict, f)
# g = render_tumor_tree(tree_structure, node_dict)
# g.render(filename=directory / f"{patient_num}_results_bootstrap_{type}_best_inter")
with open(file / f'{method}_bootstrap_summary.pkl', 'wb') as g:
    pickle.dump(tree_distribution, g)
with open(file / f'{method}_bootstrap_aggregation.pkl', 'wb') as g:
    pickle.dump(tree_aggregation, g)