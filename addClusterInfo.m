function clusterIfo = addClusterInfo(y,groupUniNames)
% Inputs:
%   y : the cluster label of each cell
%   groupUniNames : m x 1 cell array. If y is a cell array, then groupUniNames is needed.

% Outputs:
%   clusterIfo: a struct giving the cluster information
%   clusterIfo.identity: the updated cluster label of each cell
%   clusterIfo.idxCluster: a cell array, each cell contains the index of cells belong to the each cluster

if iscell(y) && exist('groupUniNames','var')
    y = categorical(y);y = reordercats(y,groupUniNames);
    y = grp2idx(y);
elseif iscell(y) && ~exist('groupUniNames','var')
    warning('Please input the cluster labels as a second input!')
end

numCluster = length(unique(y));
idxCluster = cell(1,numCluster);
for i = 1:numCluster
    idxCluster{i} = find(y == i);
end

clusterIfo.identity = y;
clusterIfo.idxCluster = idxCluster;

