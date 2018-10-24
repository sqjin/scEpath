function scEpath_demo
% This is a demo showing how to running scEpath using mouse lung epithelial specification single cell RNA-seq data (LES data)
% This demo can reproduce all the figures when analyzing LES data in our paper
clc;clear;

%% add the path where the codes are located
path_codes = '/Users/XXX/Downloads/scEpath-master/';
addpath(genpath(path_codes)) % or simply use addpath(genpath('./')) if current directory is ./scEpath-master/

%%%% running the first several steps of scEpath to calculate the single cell energy and transition probabilities, and infer cell lineages and pseudotime
%% step 1: load data and (if applicable) cell attributes (such as cell type, culture condition, day captured)
iniData = importdata(fullfile('example_data','GSE52583_LESdata.mat'));
% iniData is a struct including three fields: data, genes, cells
% iniData.data: n x m data matrix (rows are genes and columns are cells)
% iniData.genes: n x 1 cell array contains the gene names
% iniData.cells: a table, rownames are cell barcodes and each column contains cell attributes (such as cell type, culture condition, day captured). The order shoule be the same with the columns of data matrix

%% step 2: do preprocessing
minCells = 3; minGenes = 100; logNormalize = 1; filterRibo = 1; % default parameters for QC (see functions for details):
proData = preprocessing(iniData,minCells, minGenes,logNormalize,filterRibo)
%   proData: a struct variable (store the data after QC) contains the same fields with iniData

%% step 3: construct a gene-gene co-expression network
% by default: the network will be constructed by choosing the highest threshold without a significant reduction in the total number of genes;
% if one would like to quickly construct a network using a given threshold,then set quick_construct = 1 and give a tau (e.g.0.4);
% quick_construct = 0; tau = [];
quick_construct = 1; tau = 0.4;
networkIfo = constructingNetwork(proData.data',quick_construct,tau); % a struct variable
% networkIfo.R: adjacency matrix (upper matrix) of the constructed network
% networkIfo.IDselect: the index of selected genes in the constructed network

%% step 4: estimate the single cell energy (scEnergy) for each cell
[scE,scEcell] = estimatingscEnergy(proData.data,networkIfo);
% scE: energy matrix, the same dimension with data
% scEcell; scEnergy of each cell
% perform principal component analysis of the energy matrix scE
ydata = ECA(scE);% ydata : m x nPC, nPC-D coordindates from dimension reduction
ydata = ydata(:,1:2); % In this demo, only the first two significant components are used. 

%% step 5: perform unsupervised clustering of single cell data
C = []; % set C to be empty if one would like to determine the number of clusters by eigengap; otherwise please provide the number of desired clusters
y = clusteringCells(proData.data,networkIfo,C); % the cluster label of each cell

clusterIfo = addClusterInfo(y);% user can also add external clustering results here. y can be either a numerical or cell array.
% clusterIfo.identity: the updated cluster label of each cell
% clusterIfo.idxCluster: a cell array, each cell contains the index of cells belong to the each cluster
rootNode = inferStart(clusterIfo.identity,scEcell)
numCluster = length(unique(clusterIfo.identity));

%% step 6: infer the cell lineage hierarchy
alpha = 0.01; theta1 = 0.8; % default parameters (see functions for details):
lineageIfo = inferingLineage(scEcell,ydata,clusterIfo,rootNode,alpha,theta1);
% lineageIfo.TP: transition probability TP
% lineageIfo.MDST: minimal directed spanning tree, i.e.,the inferred lineage
% lineageIfo.path: the node in each path
% lineageIfo.PDG: the constructed probabilistic directed graph
% lineageIfo.centroid: the centroid of each core cell

%% step 7: reconstruct pseudotime
% Replace the following line by the appropriate path for Rscript
% Rscript = '"C:\Program Files\R\R-3.4.0\bin\Rscript"'; % for 64-bit windows
Rscript = '"/usr/local/bin/Rscript"'; % for Mac OS
% User may also need to change the path of R script to execute inside inferingPseudotime.m if current directory is not ./scEpath-master/
pseudotimeIfo = inferingPseudotime(Rscript,ydata,lineageIfo,clusterIfo);
% pseudotimeIfo.pseudotime: a cell array, each cell gives pseudotime value for each path
% pseudotimeIfo.cellOrder: a cell array, each cell gives cell order for each path
% pseudotimeIfo.cellIndex: a cell array, each cell gives cell index in each path
% pseudotimeIfo.Pcurve: a cell array, each cell gives projection vectors of fitted principal curve in each path

%%%% visualization of results including clustering, cell lineage hierarchy, scEnergy comparison and energy landscape
%% plotting
colorCell = distinguishable_colors(numCluster);% colors for each cluster
group = clusterIfo.identity; % m x 1 numerical vector, the cluster assignment for each cell
class_labels = strcat('C',cellstr(num2str([1:numCluster]'))); % text annotations of each cluster

% visualize cells on two-dimensional space
% true_labs = []; % it can be empty
true_labs = proData.cells.Time;
true_labs = categorical(true_labs);
true_labs = reordercats(true_labs,{'E14.5','E16.5','E18.5','Adult'});
marker_size = scEcell*50;  % the size of individual cell (dots), propotational to the scEnergy of each cell (by default)
fig_width = 600;fig_height = 250;
cluster_visualization(ydata, group,class_labels, true_labs, marker_size,colorCell,fig_width,fig_height)

% display cell lineage hierarchy with transition probability
node_size = grpstats(full(scEcell),group,'median')*30; % the size of tree node, propotational to the median scEnergy of each cluster
showLoops = 1; fig_width = 200;
lineage_visualization(lineageIfo,class_labels,node_size,colorCell,showLoops,fig_width)

% comparison of scEnergy among different clusters using boxplot
scEnergy_comparison_visualization(scEcell,clusterIfo,class_labels,colorCell)

% display energy landscape in 2-D contour plot and 3-D surface
nlevels = 8; % contour levels in the contour plot
landscape_visualization(scEcell,ydata,clusterIfo,colorCell,nlevels)

%%%% running other steps of scEpath to perform downstream analyses,
%%%% including gene temporal dynamics along pseudotime, identification of pseudotime dependent genes,
%%%% creation of "rolling wave" and identification of key transcription factors responsible for cell fate decision
%% downstream analysis
% (1) calculating the smooth version of expression level based on pseudotime
nbins = 10; % the number of bins for dividing the pseudotime
smoothExprIfo = cell(1,length(pseudotimeIfo.pseudotime));
for j = 1:length(pseudotimeIfo.pseudotime)
    pseudotime = pseudotimeIfo.pseudotime{j};
    cellOrder = pseudotimeIfo.cellOrder{j};
    cellw = scEcell;
    smoothExprIfo{j} = smootheningExpr(proData.data,pseudotime,cellOrder,nbins,cellw);
end

% (2) plot individual gene temporal dynamics along pseudotime
my_genes = {'Cdk1','Sox11','Pdpn','Sftpc', 'Sftpb', 'Lyz2'};
color_by = clusterIfo.identity; % color cells according to the cell attributes such as inferred cluster
plot_genes_in_pseudotime(my_genes,proData,smoothExprIfo,pseudotimeIfo,lineageIfo,color_by,colorCell)

% (3) identify pseudotime dependent genes
% Note: results of the following analyses may be slightly different with the results presented in our paper because of the use of "random number generator" in identifying pseudotime dependent genes
sd_thresh = 0.5; sig_thresh = 0.01;nboot = 1000; % default parameters (see functions for details). nboot can be changed to 500 or less
PDGIfo = identify_pseudotime_dependent_genes(proData,smoothExprIfo,pseudotimeIfo,sd_thresh,sig_thresh,nboot);
%   PDGIfo.PDG : a cell array, each cell contains the identified pseudotime-dependent genes for each path
%   PDGIfo.allGenes: a cell array, each cell contains the genes with their calculated adjusted P-values and standard deviation

% (4) create "rolling wave" showing the temporal pattern of pseudotime dependent genes and identify gene clusters showing similar pattern
optimalK = 8; % the number of desired gene clusters
pathUsed = 1; % order the pseudotime-dependent genes based on their peak expression in the "pathUsed" branch. e.g. pathUsed = 1
PDGIfo = plot_rolling_wave(PDGIfo,smoothExprIfo,proData,optimalK,pathUsed);
% PDGIfo: update the PDG information (see functions for details).

% (5) identify key transcription factors responsible for cell fate decision
% load TF information downloaded from AnimalTFDB 2.0 (http://bioinfo.life.hust.edu.cn/AnimalTFDB/)
load TF_Ifo.mat
TF = TF_Ifo.mouse.Symbol; % TF for mouse
% TF = TF_Ifo.human.Symbol; % TF for mouse
thresh_SD_TF = 0.5; thresh_FC_TF = 1; % default parameters (see functions for details):
PDG_TFIfo = identify_keyTF(proData,PDGIfo.orderedPDG,smoothExprIfo,lineageIfo,clusterIfo,TF,thresh_SD_TF,thresh_FC_TF);
% PDG_TFIfo.PDG_TFSig : key TFs that may be responsible for cell fate decisions
% PDG_TFIfo.PDG_TF : all the TFs found in the pseudotime-dependent genes

% (6) create "rolling wave" showing the temporal pattern of key transcription factors
show_TF_names = 1; % to show the TF names as Yticklabels or not
plot_rolling_wave_TF(PDG_TFIfo,PDGIfo,optimalK,show_TF_names)


%% save all the variables in the workspace
save(fullfile('results','resultsLESdata.mat'));


