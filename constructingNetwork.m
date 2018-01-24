function networkIfo = constructingNetwork(data,quick_construct,thresh,thresh_percent,showFigure,fig_width,fig_height)
% construct a gene-gene co-expression network
% Inputs:
%    data: single cell data (rows are cells and columns are genes)
%    quick_construct: default=0,the network will be constructed by choosing the highest threshold without a significant reduction in the number of genes;
%    thresh: if quick_construct=1, thresh is the threshold for constructing a network. default = 0.1
%    thresh_percent: the percentage threshold for indicating a significant reduction in the number of genes
%    showFigure: boolean, to show the network metrics or not, default= 1 (true)
%    fig_width : the figure width, default=800
%    fig_height : the figure height, default=250
% Outputs:
%    networkIfo: network information for the constructed gene-gene network
%    networkIfo.R: adjacency matrix (upper matrix) of the constructed network
%    networkIfo.IDselect: the index of selected genes in the constructed network
%    networkIfo.deg: the connectivity (degree) of each node
%    networkIfo.reduction_percent; the percentage of reduction in the number of nodes with different threshold tau
%    networkIfo.tau; the range of threshold tau

if ~exist('thresh_percent','var') || isempty(thresh_percent)
    thresh_percent = 0.2; % 20% reduction in the number of nodes
end
if ~exist('showFigure','var') || isempty(showFigure)
    showFigure = 1;
end
if ~exist('fig_width', 'var') || isempty(fig_width)
    fig_width = 600;
end
if ~exist('fig_height', 'var') || isempty(fig_height)
    fig_height = 180;
end
R00 = corr(data,'Type','Spearman');
R00 = triu(R00,1);
R00 = sparse(R00);
R00 = abs(R00);
total_nodes = size(R00,1);
if quick_construct || total_nodes < 50
    if ~exist('thresh','var') || isempty(thresh),thresh = 0.1; end
    threshSig = thresh;
    [R,IDselect,deg] = calculatingAdjacencyMatrix(R00,threshSig);
    num_nodes = length(IDselect);num_edges = nnz(R);
    reduction_percent = (total_nodes-num_nodes)/total_nodes;
else
    [thresh,num_nodes,num_edges] = calculatingNetworkMetrics(R00);
    %     reduction_percent = [0,-diff(num_nodes)./num_nodes(1:end-1)];
    reduction_percent = (total_nodes-num_nodes)/total_nodes;
    idx = find(reduction_percent > thresh_percent,1);
    threshSig = thresh(idx-1);
    [R,IDselect,deg] = calculatingAdjacencyMatrix(R00,threshSig);
end
networkIfo.R = R; networkIfo.IDselect = IDselect; networkIfo.deg = deg; networkIfo.reduction_percent = reduction_percent;networkIfo.tau = thresh;

if showFigure
    hFig = figure('position', [600, 200, fig_width, fig_height]);
    subplot(1,3,1)
    plot(thresh,num_nodes,'k-o')
    ylim([0 total_nodes*1.05])
    xlim([0.05 0.85])
    xlabel('\tau','FontName','Arial','FontSize',10);
    ylabel('Number of nodes','FontName','Arial','FontSize',10);
    box on
    grid on
    
    subplot(1,3,2)
    plot(thresh,reduction_percent*100,'k-o')
    ylim([0 100])
    xlim([0.05 0.85])
    xlabel('\tau','FontName','Arial','FontSize',10);
    ylabel({'Percentage of reduction', 'in the number of nodes'},'FontName','Arial','FontSize',10);
    ytickformat('percentage')
    grid on
    box on
    
    subplot(1,3,3)
    plot(thresh,num_edges,'k-o')
    xlim([0.05 0.85])
    xlabel('\tau','FontName','Arial','FontSize',10);
    ylabel('Number of edges','FontName','Arial','FontSize',10);
    box on
    grid on
        
    folderName = fullfile('results','figures');
    if ~exist(folderName, 'dir')
        mkdir(folderName);
    end
    saveas(hFig,fullfile(folderName,'topological_metrics_gene_network.pdf'))
    
end


function [thresh,num_nodes,num_edges] = calculatingNetworkMetrics(R00)
thresh = 0.1:0.1:0.8;
num_nodes = zeros(1,length(thresh)); num_edges = zeros(1,length(thresh));
for i = 1:length(thresh)
    [R,IDselect] = calculatingAdjacencyMatrix(R00,thresh(i));
    num_edges(i) = nnz(R);
    num_nodes(i) = length(IDselect);
    sprintf('When tau is %.2f, the number of nodes is %d.',thresh(i),length(IDselect))
end


function [R,IDselect,deg] = calculatingAdjacencyMatrix(R,thresh)
R = R>thresh;
G = graph(R,'upper');
comp = conncomp(G,'OutputForm','cell');
[~,idx] = max(cellfun('length',comp));
IDselect = comp{idx};
R = R(IDselect,IDselect);
G = graph(R,'upper');
deg = degree(G);




