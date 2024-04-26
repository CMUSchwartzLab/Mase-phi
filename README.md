# Mase-phi
Marker selection strategies for ctDNA guided by phylogenetic inference

Steps:
1. Liquid biopsy produces an EXCEL file with a mess of markers.
2. Run phyloWGS to get organized trees of the genes and how they mutated.
3. Then feed the output from phyloWGS to process_tracarx_bootstrap.py to get a summary and aggregate files of the trees.
4. Then input the EXCEL file and the aggregate file to run_data.py to find the most common markers.
5. Use run_simu_select_adjust.py or run_simu_select_track.py for more detailed outputs to run_data.py.
6. Run iterative_pipeline_real.py to update the trees with the results of running ddPCR on the selected markers from steps 4 or 5.