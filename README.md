# Mase-phi
Marker selection strategies for ctDNA guided by phylogenetic inference.

***********
Main Files:
***********
process_tracerx_bootstrap.py
    Generates summary and aggregate files based on trees generated with phyloWGS on the patient data (bootstrapped).
run_data.py
    Takes the summary and aggregate files and outputs an initial set of significant biomarkers. These biomarkers can be used in a biological assay for further inspection.
run_simu_select_adjust.py, run_simu_select_track.py
    These files take the output of a biological assay using the initial selected biomarkers, and futher update the empirical trees, track mutations over time, and output a new set of optimized biomarkers. These new biomarkers can be run in another round of biological assay.
iterative_pipeline_real.py
    With the output of the biological assay, this file produces phylogentic trees of the tumor DNA over specific points in time, which helps view the progression of the cancer over time.

***********************************
Other support files for the project
***********************************
adjust_tree_distribution.py
    Uses Bayesian inference to update tree distribution.
analyze.py
    Contains helper functions to help analyze trees.
evaluate_tree_fraction.py
    Generates a figure that shows the trackable fractions of simulated tumor clonal population in noisy liquid biopsy samples. The lines show the change of trackable fractions as the selected marker number increases, comparing marker selection optimized for the task against random marker selection.
evaluate.py
    Computes distances between true and estimated trees using CASet and DISC metrics. Collects evaluation results and stores them in a pandas DataFrame. Visualizes the evaluation results using seaborn's catplot.
optimize_fraction.py
    Uses Gurobi to solve integer linear programs meant to find the optimal fractions or weights to pick the best biomarkers.
optimize.py
    Contains helper functions to optimize phylogenetic trees based on the clonal frequencies of genetic mutations. 
visualize.py
    Generates the trees into graphs.

simulation (sim) folder
bootstrap.py
    Contains helper functions that take the patient data and "bootstraps" it, meaning that we generate similar data for testing. 
pipeline.py 
    Uses bootstrap.py functions to integrate clonal tree & variant mutation simulation with bootstrap in a pre- & post-surgery scenario.
simulate.py
    Simulates the ddPCR process.

******************
External Resources
******************
PhyloWGS [https://github.com/morrislab/phylowgs/blob/master/README.md]
CASET and DISC from this ([repository](https://bitbucket.org/oesperlab/stereodist/src/master/)). 
A ([Gurobi](https://www.gurobi.com/)) license to solve optimization problems.
Python3

********************************
Building and running the project
********************************

Steps:
1. First, liquid biopsy is conducted, which produces an EXCEL file with a many markers. The purpose of this project is to select the best markers that will tell us how the tumor is developing.
2. Run phyloWGS to get organized trees of the tumor genes and how they mutated.
3. Then feed the output from phyloWGS to process_tracerx_bootstrap.py to get a summary and aggregate files of the trees.
4. Then input the EXCEL file and aggregate files to run_data.py to create an estimated empirical tree distribution and find an initial optimal set of markers.
5. Then apply these markers in biological assays (eg ddPCR).
6. Take the output of the ddPCR and use run_simu_select_adjust.py to update the tree distributions and use run_simu_select_track.py to track the frequencies of the mutations. run_simu_select_adjust.py also outputs a file with more optimal markers that capture a range of markers and reduces redundancies.
7. Apply the new markers in biological assays (like ddPCR) again.
8. With the output of ddPCR, run iterative_pipeline_real.py to to estimate clonal fractions at each sampled timepoint. This should give the user an overview of how the tumor is potentially forming over time.