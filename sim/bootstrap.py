import pandas as pd 
import matplotlib.pyplot as plt
import numpy as np
import os
from pathlib import Path
import time
import sys
import argparse


'''
Complete:
1)Bootstrap method;
2)Input a parameter that ask for the directory contains *.xlsx files and bootstrap them;
3)Input a parameter that ask for how many times of bootstrapping you want to perform;
4)Output the Excel file for each patient in the format plyWGS can take in;
5)Perform the excel output for plyWGS input preparation for the original directory first;
6)Bootstrap each patient excel data and save each bootstrap in different directory;
7)Perform the excel output for plyWGS input preparation and save them in each bootstrap directory;
'''

############################
######### METHODS ##########
############################

#Take the AF column and the depth of reads to randomly choose a number between 0 to 1 
def bootstrap_va(AF_list, Depth_list, bootstrap_num):
      """
           This method takes the depth list and allele frequency list to bootstrap each sample 
           and create a new VAF;
           Then nest each iteration of the bootstrap into AF_lists;
      
      """
      
      AF_list_update = [] # Create an updated list of AF after bootstrapping

      for i in range(0,len(AF_list)):
            #(Bernouli Distribution)
            # Take its corresponding depth of reads to bootstrap the same number of reads;
            samples = np.random.binomial(n=Depth_list[i], p = AF_list[i], size = bootstrap_num)

            AF_list_update.append(samples/Depth_list[i]) 
      

      return AF_list_update


def bootstrap_va_dt(AF_list, Depth_list, bootstrap_num):
      """
           This method takes the depth list and allele frequency list to bootstrap each sample
           and create a new VAF;
           Then nest each iteration of the bootstrap into AF_lists;
      """

      AF_list_update = []  # Create an updated list of AF after bootstrapping
      total_depth = sum(Depth_list)
      count = 0
      while True:
          count += 1
          new_Depth_list = np.random.multinomial(n=total_depth, pvals=np.array(Depth_list)/total_depth, size=bootstrap_num)
          #print(new_Depth_list)
          if not np.any(new_Depth_list == 0): ### set the sampled depth greater than 0, otherwise causing error later
                break
          print(np.where(new_Depth_list == 0))
          if count >= 10:
              new_Depth_list[np.where(new_Depth_list == 0)] = 1
              break
      for i in range(len(AF_list)):
            AF_list_update_temp = []
            for j in range(bootstrap_num):
            # (Bernouli Distribution)
            # Take its corresponding depth of reads to bootstrap the same number of reads;
                  sample = np.random.binomial(n=new_Depth_list[j, i], p=AF_list[i], size=1)[0]
                  AF_list_update_temp.append(sample / new_Depth_list[j, i])
            AF_list_update.append(AF_list_update_temp)
      AF_list_update = np.array(AF_list_update)
      print(AF_list_update.shape, new_Depth_list.T.shape)
      return AF_list_update, new_Depth_list.T


def phyloWGS_output(input_path, output_path, num_bootstraps=0, surgery=False, sheet_type="union", num_tissue=2, num_blood=1, bootstrap_inds=0):
      """
            Access the blood and tissue sheets in every xlsx file in a directory 
            and output a dataframe that can be used for phyloWGS program;
      """


      #try:
      if surgery:
            sample_type = '_post'
      else:
            sample_type = ''

      if sheet_type == "union":
            sheet_name_list = ["union_blood_tissue_no_germline"]
            sheet_type_list = ["union"]
      elif sheet_type == "common":
            sheet_name_list = ["common_blood_tissue_no_germline"]
            sheet_type_list = ["common"]
      elif sheet_type == "all":
            sheet_name_list = ["union_blood_tissue_no_germline", "common_blood_tissue_no_germline"]
            sheet_type_list = ["union", "common"]
      for sheet_idx in range(len(sheet_name_list)):
            sheet_name = sheet_name_list[sheet_idx]
            sheet_type = sheet_type_list[sheet_idx]
            calls = pd.read_excel(input_path, sheet_name= sheet_name, index_col=0)
            print('calls:' , calls.shape)
            #pd.set_option('display.max_columns', None)


            count = 0
            if num_bootstraps == 0:

                  phylowgs_format = []
                  for i in range(calls.shape[0]):
                        flg = False
                        gene = calls["Gene"][i]
                        if not isinstance(gene, str):
                              gene = str(calls["Chromosome"][i]) + '_' + str(calls["Genomic Position"][i])

                        d_blood_list = []
                        a_blood_list = []
                        d_tissue_list = []
                        a_tissue_list = []
                        for b_idx in range(num_blood):
                              blood_no_germline = pd.read_excel(input_path, sheet_name='blood_no_germline_' + str(b_idx), index_col=0)
                              median_depth_blood = int(np.median(blood_no_germline["Depth_b{}".format(b_idx)]))

                              if np.isnan(calls['Depth_b' + str(b_idx)][i]):
                                    d_blood = median_depth_blood
                                    a_blood = d_blood
                              else:
                                    d_blood = int(calls['Depth_b'+ str(b_idx)][i])
                                    vaf_blood = calls['Allele Frequency_b'+ str(b_idx)][i]
                                    a_blood = int(np.round(d_blood * (1 - vaf_blood)))
                                    if vaf_blood > 0.9:
                                        flg = True
                              d_blood_list.append(str(d_blood))
                              a_blood_list.append(str(a_blood))

                        if flg:
                              continue
                        for t_idx in range(num_tissue):
                              tissue_no_germline = pd.read_excel(input_path, sheet_name='tissue_no_germline_' + str(t_idx),
                                                                 index_col=0)
                              median_depth_tissue = int(np.median(tissue_no_germline["Depth_t{}".format(t_idx)]))
                              if np.isnan(calls['Depth_t' + str(t_idx)][i]):
                                    d_tissue = median_depth_tissue
                                    a_tissue = d_tissue

                              else:
                                    d_tissue = int(calls['Depth_t' + str(t_idx)][i])
                                    vaf_tissue = calls['Allele Frequency_t'+ str(t_idx)][i]
                                    a_tissue = int(np.round(d_tissue * (1 - vaf_tissue)))
                                    if vaf_tissue > 0.9:
                                          flg=True
                              d_tissue_list.append(str(d_tissue))
                              a_tissue_list.append(str(a_tissue))
                        if flg:
                              continue
                        ### create phylowgs for unions
                        if num_blood != 0:
                              phylowgs_format.append({'id': 's' + str(count), 'gene': gene, 'a': ','.join(a_blood_list) + ',' +  ','.join(a_tissue_list),
                                                'd': ','.join(d_blood_list) + ',' + ','.join(d_tissue_list), 'mu_r': 0.999, 'mu_v': 0.499})
                        else:
                              phylowgs_format.append({'id': 's' + str(count), 'gene': gene,
                                                      'a': ','.join(a_tissue_list),
                                                      'd': ','.join(d_tissue_list),
                                                      'mu_r': 0.999, 'mu_v': 0.499})
                        count += 1
                  df_phylowgs = pd.DataFrame(phylowgs_format)

                  if sample_type == '_post':
                        if not os.path.exists(output_path / f'{sheet_type}/{sample_type}'):
                              os.mkdir(output_path / f'{sheet_type}/{sample_type}')
                        df_phylowgs.to_csv(output_path / f'{sheet_type}/{sample_type}/ssm_data{sample_type}.txt', sep='\t')
                  else:
                        if not os.path.exists(output_path / f'{num_blood}_{num_tissue}/{sheet_type}'):
                              os.makedirs(output_path / f'{num_blood}_{num_tissue}/{sheet_type}')
                        df_phylowgs.to_csv(output_path / f'{num_blood}_{num_tissue}/{sheet_type}/ssm_data{sample_type}.txt', sep='\t')

            else:
                  ###TODO: need to be compatible with multi-regional tumor samples
                  if bootstrap_inds == 0:
                        interpolate = 1
                        interval = num_bootstraps
                  else:
                        interpolate = 0
                        interval = 1
                  for bootstrap_idx in range(bootstrap_inds+interpolate, bootstrap_inds +interpolate + interval):
                        count = 0
                        phylowgs_format = []
                        for i in range(calls.shape[0]):
                              flg = False
                              gene = calls["Gene"][i]
                              if not isinstance(gene, str):
                                    gene = str(calls["Chromosome"][i]) + '_' + str(calls["Genomic Position"][i])
                              print(gene)
                              d_blood_list = []
                              a_blood_list = []
                              d_tissue_list = []
                              a_tissue_list = []
                              for b_idx in range(num_blood):
                                    blood_no_germline = pd.read_excel(input_path,
                                                                      sheet_name='blood_no_germline_' + str(b_idx),
                                                                      index_col=0)
                                    median_depth_blood = int(np.median(blood_no_germline["Depth_b{}".format(b_idx)]))

                                    if np.isnan(calls['Depth_b{}.Bootstrap.B{}'.format(b_idx, bootstrap_idx)][i]):
                                          d_blood = median_depth_blood
                                          a_blood = d_blood
                                    else:
                                          d_blood = int(calls['Depth_b{}.Bootstrap.B{}'.format(b_idx, bootstrap_idx)][i])
                                          vaf_blood = calls['Allele.Frequency_b{}.Bootstrap.B{}'.format(b_idx, bootstrap_idx)][i]
                                          a_blood = int(np.round(d_blood * (1 - vaf_blood)))
                                          if vaf_blood > 0.9:
                                                flg = True
                                    d_blood_list.append(str(d_blood))
                                    a_blood_list.append(str(a_blood))

                              if flg:
                                    continue
                              for t_idx in range(num_tissue):
                                    tissue_no_germline = pd.read_excel(input_path,
                                                                       sheet_name='tissue_no_germline_' + str(t_idx),
                                                                       index_col=0)
                                    median_depth_tissue = int(np.median(tissue_no_germline["Depth_t{}".format(t_idx)]))
                                    if np.isnan(calls['Depth_t{}.Bootstrap.T{}'.format(t_idx, bootstrap_idx)][i]):
                                          d_tissue = median_depth_tissue
                                          a_tissue = d_tissue

                                    else:
                                          d_tissue = int(calls['Depth_t{}.Bootstrap.T{}'.format(t_idx, bootstrap_idx)][i])
                                          vaf_tissue = calls['Allele.Frequency_t{}.Bootstrap.T{}'.format(t_idx, bootstrap_idx)][i]
                                          a_tissue = int(np.round(d_tissue * (1 - vaf_tissue)))
                                          if vaf_tissue > 0.9:
                                                flg = True
                                    d_tissue_list.append(str(d_tissue))
                                    a_tissue_list.append(str(a_tissue))
                              if flg:
                                    continue
                              ### create phylowgs for unions
                              if num_blood != 0:
                                    phylowgs_format.append({'id': 's' + str(count), 'gene': gene,
                                                            'a': ','.join(a_blood_list) + ',' + ','.join(a_tissue_list),
                                                            'd': ','.join(d_blood_list) + ',' + ','.join(d_tissue_list),
                                                            'mu_r': 0.999, 'mu_v': 0.499})
                              else:
                                    phylowgs_format.append({'id': 's' + str(count), 'gene': gene,
                                                            'a': ','.join(a_tissue_list),
                                                            'd': ','.join(d_tissue_list),
                                                            'mu_r': 0.999, 'mu_v': 0.499})

                              '''
                        median_depth_blood = int(np.median(blood_no_germline['Depth.Bootstrap.B{}'.format(bootstrap_idx)]))
                        median_depth_tissue = int(np.median(tissue_no_germline['Depth.Bootstrap.T{}'.format(bootstrap_idx)]))
                        for i in range(calls.shape[0]):
                              gene = calls["Gene"][i]
                              print(gene)
                              if not isinstance(gene, str):
                                    gene = calls["Chromosome"][i] + '_' + str(calls["Genomic Position"][i])
                              if np.isnan(calls['Depth_b{}.Bootstrap.B{}'.format(b_idx, bootstrap_idx)][i]):
                                    d_blood = median_depth_blood
                                    a_blood = d_blood
                                    vaf_blood = 0
                                    print(d_blood, a_blood)
                              else:
                                    d_blood = int(calls['Depth.Bootstrap.B{}'.format(bootstrap_idx)][i])
                                    vaf_blood = calls['Allele.Frequency.Bootstrap.B{}'.format(bootstrap_idx)][i]
                                    print(d_blood, vaf_blood)
                                    a_blood = int(np.round(d_blood * (1 - vaf_blood)))
                                    print(d_blood, vaf_blood, a_blood)
                              if np.isnan(calls['Depth.Bootstrap.T{}'.format(bootstrap_idx)][i]):
                                    d_tissue = median_depth_tissue
                                    a_tissue = d_tissue
                                    vaf_tissue = 0
                                    print(d_tissue, a_tissue)
                              else:
                                    d_tissue = int(calls['Depth.Bootstrap.T{}'.format(bootstrap_idx)][i])
                                    vaf_tissue = calls['Allele.Frequency.Bootstrap.T{}'.format(bootstrap_idx)][i]
                                    print(d_tissue, vaf_tissue)
                                    a_tissue = int(np.round(d_tissue * (1 - vaf_tissue)))
                                    print(d_tissue, vaf_tissue, a_tissue)
                              if vaf_blood > 0.9 or vaf_tissue > 0.9:
                                    continue
                              ### create phylowgs for unions
                              phylowgs_format.append(
                                    {'id': 's' + str(count), 'gene': gene, 'a': str(a_blood) + ',' + str(a_tissue),
                                     'd': str(d_blood) + ',' + str(d_tissue), 'mu_r': 0.999, 'mu_v': 0.499})
                                     '''
                              count += 1
                        df_phylowgs = pd.DataFrame(phylowgs_format)
                        os.makedirs(output_path / f'{num_blood}_{num_tissue}/{sheet_type}/bootstrap{bootstrap_idx}')
                        df_phylowgs.to_csv(output_path / f'{num_blood}_{num_tissue}/{sheet_type}/bootstrap{bootstrap_idx}/ssm_data_bootstrap{bootstrap_idx}.txt', index=False, sep='\t', header=True, columns=['id', 'gene', 'a', 'd', 'mu_r', 'mu_v'])


def merge_phyloWGS_input_bootstrap(dir1, dir2, num_bootstrap):
      for i in range(1, num_bootstrap + 1):
          file1 = f"{dir1}/bootstrap{i}/ssm_data_bootstrap{i}.txt"
          file2 = f"{dir2}/bootstrap{i}/ssm_data_bootstrap{i}.txt"
          merge_phyloWGS_input(file1, file2)
def merge_phyloWGS_input(path_file_1, path_file_2):
      with open(path_file_1, "r") as file1, open(path_file_2, "r") as file2:
            data1 = [line.strip().split('\t') for line in file1]
            data2 = [line.strip().split('\t') for line in file2]

      for line in range(1, len(data1)):
            print(line)
            data2[line][2] = data1[line][2] + ',' + data2[line][2]
            data2[line][3] = data1[line][3] + ',' + data2[line][3]

      with open(path_file_2, 'w') as f:
            for record in data2:
                  f.write('\t'.join(record) + '\n')

def search(item_list, column_1, column_2, column_3, df):
    """
        This function uses item_list to search for rows containing those items in the list 
        and export a dataframe composed of all rows;
    """
    all_list = df[[column_1, column_2, column_3]].values.tolist() #Store all samples under this column
    common_index = []
    for i in all_list:
        if i in item_list:
            common_index.append(True)
        else:
            common_index.append(False)

    df_common = df[common_index]

    df_common = df_common.reset_index(drop=True)
    return df_common


def output_bootstrap(dfs, num_bootstraps, used_num_blood=2, used_num_tissue=2,):
      Common_sheet_name = "common_blood_tissue_no_germline"
      Union_sheet_name = 'union_blood_tissue_no_germline'

      for blood_idx in range(used_num_blood):
            AF_name = "Allele Frequency_b{}".format(blood_idx)
            Depth_name = "Depth_b{}".format(blood_idx)
            df_blood = dfs['blood_no_germline_' + str(blood_idx)]
            # Access the blood_no_germline sheet with AF and Depth column;
            AF_list = df_blood[AF_name].tolist()
            Depth_list = df_blood[Depth_name].tolist()
            # Update the AF column using bootstrap;
            AF_lists_update, Depth_lists_update = bootstrap_va_dt(AF_list, Depth_list, num_bootstraps)
            col_names_af = ["Allele.Frequency_b{}.Bootstrap.B{}".format(blood_idx, i + 1) for i in range(num_bootstraps)]
            col_names_dt = ["Depth_b{}.Bootstrap.B{}".format(blood_idx, i + 1) for i in range(num_bootstraps)]
            # Make a list of column name:
            df_insert_af = pd.DataFrame(AF_lists_update, columns=col_names_af)
            df_insert_af = df_insert_af.reset_index(drop=True)
            df_insert_dt = pd.DataFrame(Depth_lists_update, columns=col_names_dt)
            df_insert_dt = df_insert_dt.reset_index(drop=True)
            # The column index number of Allele.Frequency, therefore, can insert at the right position
            column_index = df_blood.columns.get_loc(Depth_name)

            num = 0  # tracking the index number of column

            for idx in range(len(df_insert_af.columns)):
                  column_af = df_insert_af.columns[idx]
                  column_dt = df_insert_dt.columns[idx]
                  num += 1
                  df_blood.insert(num + column_index, column_af, df_insert_af[column_af].values, True)
                  num += 1
                  df_blood.insert(num + column_index, column_dt, df_insert_dt[column_dt].values, True)
            print("Pass blood")
      if used_num_blood != 0:
            df_blood_common = dfs['blood_no_germline_0'].copy(deep=True)
            for idx in range(1, used_num_blood):
                  df_blood_common = pd.merge(df_blood_common, dfs['blood_no_germline_{}'.format(idx)], how='inner',
                              on=['Gene', 'Chromosome', 'Genomic Position','Alternative Call','Reference Call', 'Consequence(s)'])

            df_blood_union = dfs['blood_no_germline_0'].copy(deep=True)
            for idx in range(1, used_num_blood):
                  df_blood_union = pd.merge(df_blood_union, dfs['blood_no_germline_{}'.format(idx)], how='outer',
                              on=['Gene', 'Chromosome', 'Genomic Position','Alternative Call','Reference Call', 'Consequence(s)'])
      # except:
      #       raise ("No blood sheet found!")

      for tissue_idx in range(used_num_tissue):
            AF_name = "Allele Frequency_t{}".format(tissue_idx)
            Depth_name = "Depth_t{}".format(tissue_idx)
            df_tissue = dfs['tissue_no_germline_' + str(tissue_idx)]
            # Access the tissue_no_germline sheet with AF and Depth column;
            AF_list = df_tissue[AF_name].tolist()
            Depth_list = df_tissue[Depth_name].tolist()
            # Update the AF column using bootstrap;
            AF_lists_update, Depth_lists_update = bootstrap_va_dt(AF_list, Depth_list, num_bootstraps)

            col_names_af = ["Allele.Frequency_t{}.Bootstrap.T{}".format(tissue_idx, i + 1) for i in range(num_bootstraps)]
            col_names_dt = ["Depth_t{}.Bootstrap.T{}".format(tissue_idx, i + 1) for i in range(num_bootstraps)]

            # Make a list of column name:
            df_insert_af = pd.DataFrame(AF_lists_update, columns=col_names_af)
            print("df_insert_af.shape", len(df_insert_af.columns))
            df_insert_af = df_insert_af.reset_index(drop=True)
            df_insert_dt = pd.DataFrame(Depth_lists_update, columns=col_names_dt)
            df_insert_dt = df_insert_dt.reset_index(drop=True)

            # The column index number of Allele.Frequency, therefore, can insert at the right position
            column_index = df_tissue.columns.get_loc(Depth_name)
            num = 0  # tracking the index number of column

            for idx in range(len(df_insert_af.columns)):
                  column_af = df_insert_af.columns[idx]
                  column_dt = df_insert_dt.columns[idx]
                  num += 1
                  df_tissue.insert(num + column_index, column_af, df_insert_af[column_af].values, True)
                  num += 1
                  df_tissue.insert(num + column_index, column_dt, df_insert_dt[column_dt].values, True)
                  print(idx, df_tissue.shape)
            print("Pass Tissue")
      df_tissue_common = dfs['tissue_no_germline_0'].copy(deep=True)
      for idx in range(1, used_num_tissue):
            df_tissue_common = pd.merge(df_tissue_common, dfs['tissue_no_germline_{}'.format(idx)], how='inner',
                                                    on=['Gene', 'Chromosome', 'Genomic Position','Alternative Call','Reference Call','Consequence(s)'])

      df_tissue_union = dfs['tissue_no_germline_0'].copy(deep=True)
      for idx in range(1, used_num_tissue):
            df_tissue_union = pd.merge(df_tissue_union, dfs['tissue_no_germline_{}'.format(idx)], how='outer',
                                       on=['Gene', 'Chromosome', 'Genomic Position','Alternative Call','Reference Call','Consequence(s)'])
      # except:
      if used_num_blood != 0:
            df_common = pd.merge(df_blood_common, df_tissue_common, how='inner', on=['Gene', 'Chromosome', 'Genomic Position','Alternative Call','Reference Call','Consequence(s)'])
            dfs[Common_sheet_name] = df_common

            df_union = pd.merge(df_blood_union, df_tissue_union, how='outer', on=['Gene', 'Chromosome', 'Genomic Position','Alternative Call','Reference Call','Consequence(s)'])
            dfs[Union_sheet_name] = df_union
      else:
            df_common = df_tissue_common
            dfs[Common_sheet_name] = df_common
            df_union = df_tissue_union
            dfs[Union_sheet_name] = df_union
      print("COMPLETE")
      print("-----------")
      print(dfs)
      return dfs

########################################################
def main_bootstrap(directory_path_input, num_bootstraps, used_blood, used_tissue='full', surgery=False):
      """
            This function is made for performing bootstrapping on the given real data from patients
            in the same directory.

            ** The parameter refers to the column & sheet name for the patient .xlsx file
               it already has default names for parameter but can be changed if needed as
               the given data name is not always consistent;
            directory_path_input: The path of the directory that you want to bootstrap
            num_iteration: Number of the bootstrapping iteration;
      """


      directory = Path(directory_path_input)  # Make the input in Path format;

      # Perform the bootstrapping for the input number of iteration;
      directory_bootstrap = Path(os.path.join(directory, "bootstrap"))
      #os.makedirs(directory_bootstrap)
      # Make a new directory for each iteration of bootstrapping;

      # print("[---BOOTSTRAP ",num+1," HAS STARTED!---]")
      start_time = time.time()
      for file in directory.glob("*.xlsx"):

            # The number of the patient
            try:
                  patient_num = int(file.stem.split('_')[1])
                  print("-----------")
                  print("Patient: ", patient_num)
            except IndexError:
                  patient_num = 0
                  patient_name = file.stem
                  patient_num = patient_name
                  print(patient_num)
            # Read the blood_no_germline and tissue_no germline sheet;
            # Store the whole file in a dictionary of dataframe;
            dfs = pd.read_excel(file, sheet_name=None, index_col=0)
            print(dfs)
            if used_tissue == 'full':
                  used_tissue_num = len(dfs)
            elif isinstance(used_tissue, int):
                  used_tissue_num = used_tissue
                  print(used_tissue_num)
            dfs = output_bootstrap(dfs, num_bootstraps, used_num_blood=used_blood, used_num_tissue=used_tissue_num)

            # Create a writer object to write the dataframes to a single Excel file
            if isinstance(patient_num, int):
                  path_new_bootstrap = Path(directory_bootstrap / f'Patient_{patient_num}_bootstrap.xlsx')
            else:
                  path_new_bootstrap = Path(directory_bootstrap / f'{patient_num}_bootstrap.xlsx')
            writer = pd.ExcelWriter(path_new_bootstrap, engine='openpyxl')

            # Write the bootstrap AF of each patient in this bootstrap directory
            for sheet_name, data in dfs.items():
                  data.to_excel(writer, index=True, sheet_name=sheet_name, engine='openpyxl')
            writer.close()
            if not os.path.exists(directory_bootstrap / f'{patient_num}'):
                  os.mkdir(directory_bootstrap / f'{patient_num}')
            #phyloWGS_output(path_new_bootstrap, directory_bootstrap / f'{patient_num}', num_bootstraps, surgery, "union", num_tissue=used_tissue, num_blood=used_blood)
      end_time = time.time()
      print(end_time - start_time)


def write_phylowgs_output(directory_path_input, num_bootstraps, used_blood, used_tissue='full', surgery=False, bootstrap_inds=0):
      """
            This function is made for performing bootstrapping on the given real data from patients
            in the same directory.

            ** The parameter refers to the column & sheet name for the patient .xlsx file
               it already has default names for parameter but can be changed if needed as
               the given data name is not always consistent;
            directory_path_input: The path of the directory that you want to bootstrap
            num_iteration: Number of the bootstrapping iteration;
      """

      directory = Path(directory_path_input)  # Make the input in Path format;

      # Perform the bootstrapping for the input number of iteration;
      directory_bootstrap = Path(os.path.join(directory, "bootstrap"))
      # os.makedirs(directory_bootstrap)
      # Make a new directory for each iteration of bootstrapping;

      # print("[---BOOTSTRAP ",num+1," HAS STARTED!---]")
      start_time = time.time()
      for file in directory.glob("*.xlsx"):

            # The number of the patient
            try:
                  patient_num = int(file.stem.split('_')[1])
                  print("-----------")
                  print("Patient: ", patient_num)
            except IndexError:
                  patient_num = 0
                  patient_name = file.stem
                  patient_num = patient_name
                  print(patient_num)
            # Read the blood_no_germline and tissue_no germline sheet;
            # Store the whole file in a dictionary of dataframe;
            dfs = pd.read_excel(file, sheet_name=None, index_col=0)
            # print(dfs)
            if used_tissue == 'full':
                  used_tissue_num = len(dfs)
            elif isinstance(used_tissue, int):
                  used_tissue_num = used_tissue
                  print(used_tissue_num)
            # dfs = output_bootstrap(dfs, num_bootstraps, used_blood, used_tissue)

            # Create a writer object to write the dataframes to a single Excel file
            if isinstance(patient_num, int):
                  path_new_bootstrap = Path(directory_bootstrap / f'Patient_{patient_num}_bootstrap.xlsx')
            else:
                  path_new_bootstrap = Path(directory_bootstrap / f'{patient_num}_bootstrap.xlsx')
            # writer = pd.ExcelWriter(path_new_bootstrap, engine='openpyxl')
            #
            # # Write the bootstrap AF of each patient in this bootstrap directory
            # for sheet_name, data in dfs.items():
            #       data.to_excel(writer, index=True, sheet_name=sheet_name, engine='openpyxl')
            # writer.close()
            # if not os.path.exists(directory_bootstrap / f'{patient_num}'):
            #       os.mkdir(directory_bootstrap / f'{patient_num}')
            phyloWGS_output(path_new_bootstrap, directory_bootstrap / f'{patient_num}', num_bootstraps, surgery,
                            "common", num_tissue=used_tissue_num, num_blood=used_blood, bootstrap_inds=bootstrap_inds)
      end_time = time.time()
      print(end_time - start_time)


def write_phylowgs_bootstrap(directory_path_input, num_bootstraps, surgery=False):
      """
            This function is made for performing bootstrapping on the given real data from patients
            in the same directory.

            ** The parameter refers to the column & sheet name for the patient .xlsx file
               it already has default names for parameter but can be changed if needed as
               the given data name is not always consistent;
            directory_path_input: The path of the directory that you want to bootstrap
            num_iteration: Number of the bootstrapping iteration;
      """


      directory = Path(directory_path_input)  # Make the input in Path format;

      # Perform the bootstrapping for the input number of iteration;
      directory_bootstrap = Path(os.path.join(directory, "bootstrap"))
      os.makedirs(directory_bootstrap, exist_ok=True)
      # Make a new directory for each iteration of bootstrapping;

      # print("[---BOOTSTRAP ",num+1," HAS STARTED!---]")
      start_time = time.time()
      for file in directory.glob("*.xlsx"):

            # The number of the patient
            patient_num = int(file.stem.split('_')[1])
            print("-----------")
            print("Patient: ", patient_num)
            # Read the blood_no_germline and tissue_no germline sheet;
            # Store the whole file in a dictionary of dataframe;


            # Create a writer object to write the dataframes to a single Excel file
            path_new_bootstrap = Path(directory_bootstrap / f'Patient_{patient_num}_bootstrap.xlsx')

            if not os.path.exists(directory_bootstrap / f'{patient_num}'):
                  os.mkdir(directory_bootstrap / f'{patient_num}')
            phyloWGS_output(path_new_bootstrap, directory_bootstrap / f'{patient_num}', num_bootstraps, surgery, "all")
      end_time = time.time()
      print(end_time - start_time)

def get_args(argv):
    parser = argparse.ArgumentParser(prog='bootstrap.py',)
    parser.add_argument('-d', '--directory', type=str, dest='directory')
    parser.add_argument('-b', '--num_bootstrap', type=int, dest='num_bootstrap', default=100)
    parser.add_argument('-i', '--bootstrap_inds', type=int, dest='bootstrap_inds', default=0)
    return vars(parser.parse_args(argv))

##################################################################
################   Main Function Implementation  #################
##################################################################
if __name__ == "__main__":
      argv = get_args(sys.argv[1:])
      bootstrap_num = argv["num_bootstrap"]
      bootstrap_inds = argv["bootstrap_inds"]
      directory = argv["directory"]
      write_phylowgs_output(directory, num_bootstraps=bootstrap_num, used_blood=0, used_tissue="full", bootstrap_inds=int(bootstrap_inds))
      #main_bootstrap(directory, bootstrap_num, used_blood=0, used_tissue="full" )
      #write_phylowgs_bootstrap("MOONSHOT2/filtered_xlsx", bootstrap_num)
