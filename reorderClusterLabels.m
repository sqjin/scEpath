function clusterIfo = reorderClusterLabels(y,scEcell,clusterOrder,proData,markers)
% reorder the clusters such that the scEnergies decrease
% Inputs:
%   y : the cluster label of each cell
%   scEcell : m x 1 vector, single cell energy
%   clusterOrder : default= []; the clusters are reordered such that the first one has biggest average scEnergy (i.e., the fist one will be the initial state in cell lineages)
%   proData : the single cell data including genes and cell attributes
%   markers : one marker gene name. scEpath needs the prior knowledge to determine an initial state if there is no significant difference (p-value > 0.01) between the two metacells with the highest entropies.
% Outputs:
%   clusterIfo: a struct giving the cluster information
%   clusterIfo.identity: the updated cluster label of each cell
%   clusterIfo.idxCluster: a cell array, each cell contains the index of cells belong to the each cluster
if ~exist('clusterOrder','var')
    clusterOrder = [];
else
    clusterOrderAuto = clusterOrder;
end
if ~exist('markers','var')
    markers = [];
end

numCluster = length(unique(y));
if isempty(clusterOrder)
    EcellTypeMedian = zeros(1,numCluster);
    for i = 1:numCluster
        Q = quantile(scEcell(y == i),[0.25, 0.5, 0.75]);
        EcellTypeMedian(i) = 0.25*Q(1)+0.5*Q(2)+0.25*Q(3); % Tukey's trimean
    end
    [~,clusterOrderAuto] = sort(EcellTypeMedian,'descend');
end
idx = cell(1,numCluster);
for i = 1:numCluster
    idx{i} = find(y == clusterOrderAuto(i));
end
for i = 1:numCluster
    y(idx{i}) = i;
end

if isempty(clusterOrder)
    % judge whether the first two clusters have significant difference of energy
    p = ranksum(scEcell(y == 1),scEcell(y == 2));
    if p > 0.01
        if ~isempty(markers)
            [~,~,idxmarkers] = intersect(markers,proData.genes,'stable');
            if mean(proData.data(idxmarkers,y == 1)) < mean(proData.data(idxmarkers,y == 2))
                y1 = find(y == 1); y2 = find(y == 2);
                y(y1) = 2*ones(length(y1),1);y(y2) = 1*ones(length(y2),1);
            end
        else
            warning('Please provide a marker gene to determine the initial state!')
        end
    end
end

group = y;% assign new clustering labels to cells

idxCluster = cell(1,numCluster);
for i = 1:numCluster
    idxCluster{i} = find(group == i);
end

if ~isempty(markers)
    [~,~,idxmarkers] = intersect(markers,proData.genes,'stable');
    median_levels = grpstats(proData.data(idxmarkers,:),group,'median');
    fprintf('%s\n',['The median expression of ', markers, ' in each cluster is '])
    median_levels
end

clusterIfo.identity = group;
clusterIfo.idxCluster = idxCluster;

