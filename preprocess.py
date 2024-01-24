import glob
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

directory=Path('MOONSHOT2')
tissue_bc_dir = directory / 'BC-0011-486-i-2_CombinedVariantOutput.tsv'
tissue_dir = directory / 'ST-190-486_CombinedVariantOutput.tsv'
blood_bc_dir = directory / 'CF-0241-486-BC_CombinedVariantOutput.tsv'
blood_dir = directory / 'CF-0019-486-i_CombinedVariantOutput.tsv'

# scrub all the raw files
files = [tissue_bc_dir, tissue_dir, blood_bc_dir, blood_dir]
annot = ['tissue_bc', 'tissue', 'blood_bc', 'blood']
for idx in range(len(files)):
    file = files[idx]
    print_bool = False
    with open(file, 'r') as f:
        for line in f:
            if print_bool and len(line.strip()) > 5:
                print(line, file=f_out, end='')
            if line.strip() == '[Small Variants]':
                print_bool = True
                annot_file = annot[idx] + '.tsv'
                file_new = directory / annot_file
                f_out = open(file_new, 'w')
        f_out.close()

blood = pd.read_csv(directory / 'blood.tsv', sep='\t')
tissue = pd.read_csv(directory /'tissue.tsv', sep='\t')
blood_bc = pd.read_csv(directory /'blood_bc.tsv', sep='\t')
tissue_bc = pd.read_csv(directory /'tissue_bc.tsv', sep='\t')

blood_bc_inter = pd.merge(blood, blood_bc, how='inner', on=['Gene', 'Chromosome', 'Genomic Position'])
blood_bc_union = pd.merge(blood, blood_bc, how='outer', on=['Gene', 'Chromosome', 'Genomic Position'])
blood_no_germline = blood_bc_union[blood_bc_union.iloc[:,[13]].isnull().any(axis=1)].iloc[:, :10]
tissue_bc_inter = pd.merge(tissue, tissue_bc, how='inner', on=['Gene', 'Chromosome', 'Genomic Position'])
tissue_bc_union = pd.merge(tissue, tissue_bc, how='outer', on=['Gene', 'Chromosome', 'Genomic Position'])
tissue_no_germline = tissue_bc_union[tissue_bc_union.iloc[:,[13]].isnull().any(axis=1)].iloc[:, :10]

blood_no_germline = blood_no_germline.rename(columns={'Reference Call_x': 'Reference Call', 'Alternative Call_x': 'Alternative Call', 'Allele Frequency_x': 'Allele Frequency',
    'Depth_x': 'Depth', 'P-Dot Notation_x': 'P-Dot Notation', 'C-Dot Notation_x': 'C-Dot Notation', 'Consequence(s)_x':'Consequence(s)'})
tissue_no_germline = tissue_no_germline.rename(columns={'Reference Call_x': 'Reference Call', 'Alternative Call_x': 'Alternative Call', 'Allele Frequency_x': 'Allele Frequency',
    'Depth_x': 'Depth', 'P-Dot Notation_x': 'P-Dot Notation', 'C-Dot Notation_x': 'C-Dot Notation', 'Consequence(s)_x':'Consequence(s)'})

patient_num=486

inter = pd.merge(blood_no_germline, tissue_no_germline, how='inner', on=['Gene', 'Chromosome', 'Genomic Position'])
union = pd.merge(blood_no_germline, tissue_no_germline, how='outer', on=['Gene', 'Chromosome', 'Genomic Position'])
blood_specific = union[union.iloc[:,[13]].isnull().any(axis=1)]
tissue_specific = union[union.iloc[:,[5]].isnull().any(axis=1)]
common_no_germline = inter
plt.figure()
plt.scatter(common_no_germline.iloc[:, 5], common_no_germline.iloc[:, 12],s=9, label='common')
plt.scatter(np.zeros((tissue_specific.shape[0])),tissue_specific.iloc[:,12],s=9, label='tissue specific')
plt.scatter(blood_specific.iloc[:,5], np.zeros((blood_specific.shape[0])),s=9, label='blood specific')
#plt.ylim(0,1)
#plt.xlim(0,1)
plt.xlabel("Blood VAF")
plt.ylabel("Tissue VAF")
plt.title(patient_num)
plt.legend()
#plt.savefig(directory / str(patient_num) + '_vaf.png')

import seaborn as sns

blood_no_germline['type'] = 'blood'
tissue_no_germline['type'] = 'tissue'
df_temp_all = None
df_temp = pd.concat([blood_no_germline, tissue_no_germline], axis=0)
if df_temp_all is None:
    df_temp_all = df_temp
else:
    df_temp_all = pd.concat([df_temp_all, df_temp], axis=0)

ax = sns.violinplot(data = df_temp, x='Allele Frequency', y='type')
plt.title(patient_num)
plt.savefig(str(patient_num) + '_vaf_violin.png')


with pd.ExcelWriter('patient_486.xlsx') as writer:
    blood.to_excel(writer, sheet_name='blood')
    tissue.to_excel(writer, sheet_name='tissue')
    blood_bc.to_excel(writer, sheet_name='blood_buffy_coat')
    tissue_bc.to_excel(writer, sheet_name='tissue_buffy_coat')
    blood_no_germline.to_excel(writer, sheet_name='blood_no_germline')
    tissue_no_germline.to_excel(writer, sheet_name='tissue_no_germline')
    common_no_germline.to_excel(writer, sheet_name='common_blood_tissue_no_germline')
