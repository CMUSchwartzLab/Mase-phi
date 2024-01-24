import copy
import pickle
from visualize import *
from analyze import *
from optimize import *
import pandas as pd
from evaluate import *
import seaborn as sns

directory=Path('/data/liquid_biopsy/simulations/mask/multi-tumor')


num_node = 10
depth_tissue = 800
depth_blood = 3000
mut_rate = 100
recover_rate = 0.5
mask_proportion = 0.5
num_blood = 3
num_tissue = 3
used_num_blood = 1
used_num_tissue = 3
dir_simu = directory / f'boot_simulation_{num_node}_{depth_tissue}_{depth_blood}_{mut_rate}_{int(recover_rate*100)}_{int(mask_proportion*100)}/{num_blood}_{num_tissue}'

bootstrap_num = 100
method = 'phylowgs'
type="common"
num_chains=5
df_dist = []
t_idx = used_num_tissue
df_fractions_total = None
for file in dir_simu.glob(f'Simulation*'):
    tree_fraction_file = file / f'Bootstrap/{used_num_blood}_{used_num_tissue}/common/{method}_fraction_dataframe.csv'
    df_fractions = pd.read_csv(tree_fraction_file, sep='\t', index_col=0)
    if df_fractions_total is None:
        df_fractions_total = copy.deepcopy(df_fractions)
    else:
        df_fractions_total = pd.concat([df_fractions_total, df_fractions])

fig, ax = plt.subplots(1,1, figsize=(8,5))
simu_num = 10
last_end = 0
xtick_label = []
for simulation_idx in range(1, simu_num + 1):
    xtick_label.append(last_end)
    df_fractions_selected = df_fractions_total[df_fractions_total["simu_idx"] == simulation_idx]
    df_fractions_selected["loc"] = last_end + 0.2*(df_fractions_selected["gene_rank"])
    last_end = max(df_fractions_selected["loc"])

    ax1 = sns.lineplot(data=df_fractions_selected, x='loc', y='trackable_fractions', legend=None, ax=ax, label='optim', c='tab:blue')
    ax2 = sns.lineplot(data=df_fractions_selected, x='loc', y='trackable_fractions_random', legend=None, ax=ax, label='random', c='tab:orange')

ax.legend(["optim", 'random'])

#ax2.legend(bbox_to_anchor=(1,1))
ax.set_xticks(xtick_label)
ax.set_xticklabels(range(1, simu_num + 1))
ax.set_xlabel("simu_idx")
ax.set_ylabel("trackable fractions")
plt.savefig(f"evaluate_tree_fraction_{num_node}_{depth_tissue}_{depth_blood}_{mut_rate}_{int(recover_rate*100)}_{int(mask_proportion*100)}_{used_num_blood}_{used_num_tissue}.eps", bbox_inches='tight')
