function [y, S] = clusteringCells(data,networkIfo,C,clusterRange,showFigure)
% perform unsupervised clustering of single cell data
% Inputs:
%    data: single cell data (rows are genes and columns are cells)
%    networkIfo: network information for the constructed gene-gene network
%    C: number of clusters, by default automatically choosing based on eigengap
%    clusterRange: a vector,the range of potential clusters when automatically choosing the number of clusters, default=2:10
%    showFigure: boolean, to show the similarity matrix or not, default= 1 (true)
% Outputs:
%    y is the cluster label of each cell
%    S is the similarity matrix
if ~exist('showFigure','var') || isempty(showFigure)
    showFigure = 1;
end
if ~exist('C','var')
    C = [];
end
if ~exist('clusterRange','var') || isempty(clusterRange)
    clusterRange = 2:10;
end
data = data(networkIfo.IDselect,:);
rng('default'); %%% for reproducibility
%% check if SIMLR has been successfully installed (https://github.com/BatzoglouLabSU/SIMLR)
try
    SIMLR(data(1:min(50,size(data,1)),1:20)',2);
catch
    INSTALL_SIMLR % compile external C program for the computations of SIMLR
end
%%
if isempty(C)
    for c = clusterRange
        [y, S] = SIMLR(data',c,10);
        K1 = EstimateNumberClusters(abs(S), clusterRange,0);
        if isequal(K1,c)
            K1 = EstimateNumberClusters(abs(S), clusterRange,1);
            break;
        end
    end
else
    [y, S] = SIMLR(data',C,10);
end

%% visualization of similarity matrix
if showFigure
    [~,labelOrdered] = sort(y);
    Ssym = (abs(S)+abs(S'))/2;
    SsymOrdered = Ssym(labelOrdered,labelOrdered);
    SsymOrdered(logical(eye(size(SsymOrdered)))) = 0;
    hFig2 = figure;
    imagesc(SsymOrdered); axis square;
    colormap hot;
    c = colorbar;
    c.Location = 'eastoutside';
    % c.Label.String = 'Similarity';
    % c.Label.FontSize = 8;%c.Label.FontWeight = 'bold';
    c.FontSize = 8;
    set(gca,'FontSize',8)
    xlim([0.5 length(y)+0.5]);
    xlabel('Cells','FontName','Arial','FontSize',10)
    ylabel('Cells','FontName','Arial','FontSize',10)
    folderName = fullfile('results','figures');
    if ~exist(folderName, 'dir')
        mkdir(folderName);
    end
    saveas(hFig2,fullfile(folderName,'cell_similarity_matrix.pdf'))
end


function [K1, K12,eigengap] = EstimateNumberClusters(W, NUMC,showFigure)
%%%This function estimates the number of clusters using eigengap

%W is the similarity graph
%NUMC is a vector which contains the possible choices of number of clusters.
%%K1 is the estimated best number of clusters according to eigen-gaps
%%K12 is the estimated SECOND best number of clusters according to eigen-gaps
% an example would be [K1, K12,eigengap] = Estimate_Number_of_Clusters_given_graph(W, [2:5]);

%%Note that this function can only give an estimate of the number of
%%clusters. How to determine the "OPTIMAL" number of clusters, is still an
%%open question so far.

if nargin < 2
    NUMC = 2:10;
end

if min(NUMC)==1
    warning('Note that we always assume there are more than one cluster.');
    NUMC(NUMC<=1) = [];
end
W = (W + W')/2;
W = W - diag(diag(W));

if ~isempty(NUMC)
    degs = sum(W, 2);
    D    = sparse(1:size(W, 1), 1:size(W, 2), degs);
    % compute unnormalized Laplacian
    L = D - W;
    degs(degs == 0) = eps;
    % calculate D^(-1/2)
    D = spdiags(1./(degs.^0.5), 0, size(D, 1), size(D, 2));
    % calculate normalized Laplacian
    L = D * L * D;
    % compute the eigenvectors corresponding to the k smallest
    % eigenvalues
    [U, eigenvalue] = eigs(L, max(NUMC)+1, eps);
    eigenvalue  = diag(eigenvalue);
    [a,b] = sort((eigenvalue),'ascend');
    eigenvalue = (eigenvalue(b));
    eigengap = abs(diff(eigenvalue));
    [tt1, t1] = sort(eigengap(NUMC),'descend');K1 = NUMC(t1(1));K12 = NUMC(t1(2));
    
end
if ~exist('showFigure','var')
    showFigure = 1;
end
if showFigure
    hFig1 = figure;
    plot(NUMC,eigengap(NUMC),'k-')
    hold on
    scatter(K1,tt1(1),20,'r','filled')
    line([K1,K1],[0 tt1(1)],'Color','red','LineStyle','--')
    xlabel('Number of clusters','FontName','Arial','FontSize',10)
    ylabel('Eigengap','FontName','Arial','FontSize',10)
    set(gca,'Xtick',NUMC)
    folderName = fullfile('results','figures');
    if ~exist(folderName, 'dir')
        mkdir(folderName);
    end
    saveas(hFig1,fullfile(folderName,'eigenGap.pdf'))
end




