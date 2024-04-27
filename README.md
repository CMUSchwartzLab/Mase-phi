# Mase-phi-commented
Marker selection strategies for ctDNA guided by phylogenetic inference.

Here is a quick overview of how to use this github repository.

Steps:
1. First, liquid biopsy is conducted, which produces an EXCEL file with a many markers. The purpose of this project is to select the best markers that will tell us how the tumor is developing.
2. Run ([phyloWGS](https://github.com/morrislab/phylowgs)) to get organized trees of the tumor genes and how they mutated.
3. Then feed the output from phyloWGS to process_tracerx_bootstrap.py to get a summary and aggregate files of the trees.
4. Then input the EXCEL file and aggregate files to run_data.py to find the most common/important markers.
5. Use run_simu_select_adjust.py or run_simu_select_track.py for more detailed outputs to run_data.py.
6. Run iterative_pipeline_real.py to update the trees with the results of running ddPCR on the selected markers from steps 4 or 5.