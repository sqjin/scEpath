# scEpath
Package of scEpath (a novel tool for analyzing single cell RNA-seq data)
===============

Overview
--------

This is a MATLAB Package of scEpath ("single-cell Energy path"). scEpath is a novel computational method for quantitatively measuring developmental potency and plasticity of single cells and transition probabilities between cell states, and inferring lineage relationships and pseudotemporal ordering from single-cell gene expression data. In addition, scEpath performs many downstream analyses including identification of the most important marker genes or transcription factors for given cell clusters or over pseudotime.

The rational of scEpath for inferring cellular trajectories is the [Waddington landscape](https://dev.biologists.org/content/142/19/3274)

<img src="https://github.com/sqjin/scEpath/blob/master/example_data/Waddington's%20landscape.jpg" width="500">


Check out [our paper on Bioinformatics](https://academic.oup.com/bioinformatics/article/34/12/2077/4838235) for the detailed methods and applications. 

![Overview of scEpath](https://github.com/sqjin/scEpath/blob/master/example_data/overview_scEpath.png)

Systems Requirements
--------------------

scEpath is independent of operating systems because it is written in Matlab. Basic requirement for running scEpath includes MATLAB and the Statistics toolbox. The pseudotime estimation step requires the R package "princurve" for principal curve analysis. In this case, both R and Matlab are required for running scEpath. 

This Package has been tested using MATLAB 2016a/b/2017a on Mac OS/64-bit Windows. 


Usage
-----

Unzip the package. Change the current directory in Matlab to the folder containing the scripts.

This directory includes the following main scripts:
1) scEpath_demo.m -- an example run of scEpath on a specific dataset
2) preprocessing.m -- do preprocessing of the input data (if applicable) 
3) constructingNetwork.m -- construct a gene-gene co-expression network
4) estimatingscEnergy.m -- estimate the single cell energy (scEnergy) for each cell
5) ECA.m -- prinpipal component analysis of energy matrix
6) clusteringCells.m -- perform unsupervised clustering of single cell data
7) addClusterInfo.m -- integrate clustering information
8) inferingLineage.m -- infer the cell lineage hierarchy
9) FindMDST.m -- find the minimal directed spanning tree in a directed graph
10) inferingPseudotime.m -- reconstruct pseudotime
11) smootheningExpr.m -- calculating the smooth version of expression level based on pseudotime
12) identify_pseudotime_dependent_genes.m -- identify pseudotime dependent marker genes
13) identify_keyTF.m -- identify key transcription factors responsible for cell fate decision
---------------------------
14) cluster_visualization.m -- visualize cells on two-dimensional space
15) lineage_visualization.m -- display cell lineage hierarchy with transition probability
16) scEnergy_comparison_visualization.m -- comparison of scEnergy among different clusters
17) landscape_visualization -- display energy landscape in 2-D contour plot and 3-D surface
18) plot_genes_in_pseudotime.m -- plot the temporal dynamics of individual gene along pseudotime
19) plot_rolling_wave.m -- create "rolling wave" showing the temporal pattern of pseudotime-dependent genes and display gene clusters showing similar patterns
20) plot_rolling_wave_TF.m -- create "rolling wave" showing the temporal pattern of key transcription factors


For each run, the final results of the analysis are deposited in the "results" directory:
1) results/figures, containing PDF figures of the analysis.
2) results/PDG_in_each_cluster, containing the identified pseudotime-dependent marker genes in each cluster
3) results/temporalfiles, containing intermediate results from the analysis.


Please refer to scEpath_demo.m for instructions on how to use this code.
Input Data are gene expression data matrix (rows are genes and columns are cells). 

If you have any problem or question using the package please contact suoqin.jin@uci.edu

