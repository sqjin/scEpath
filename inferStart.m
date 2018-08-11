function rootNode = inferStart(y,scEcell,proData,markers)
% infer the start state for lineage inference 
% Inputs:
%   y : the cluster label of each cell
%   scEcell : m x 1 vector, single cell energy
%   proData : the single cell data object
%   markers : one marker gene name. scEpath needs the prior knowledge to determine an initial state if there is no significant difference (p-value > 0.01) between the two metacells with the highest entropies.
% Outputs:
%   rootNode: the inferred start state for lineage inference 

if ~exist('markers','var')
    markers = [];
end

numCluster = length(unique(y));

EcellTypeMedian = zeros(1,numCluster);
for i = 1:numCluster
    Q = quantile(scEcell(y == i),[0.25, 0.5, 0.75]);
    EcellTypeMedian(i) = 0.25*Q(1)+0.5*Q(2)+0.25*Q(3); % Tukey's trimean
end
[~,clusterOrder] = sort(EcellTypeMedian,'descend');



% judge whether the first two clusters have significant difference of energy
p = ranksum(scEcell(y == clusterOrder(1)),scEcell(y == clusterOrder(2)));
if p > 0.01
    if ~isempty(markers)
        [~,~,idxmarkers] = intersect(markers,proData.genes,'stable');
        if mean(proData.data(idxmarkers,y == clusterOrder(1))) < mean(proData.data(idxmarkers,y == clusterOrder(2)))
            rootNode = clusterOrder(2);
        else
            rootNode = clusterOrder(1);
        end
    else
        disp('Please provide a marker gene to determine the start state!')
        rootNode = NaN;
    end
else
    rootNode = clusterOrder(1);
end

if ~isempty(markers)
    [~,~,idxmarkers] = intersect(markers,proData.genes,'stable');
    median_levels = grpstats(full(proData.data(idxmarkers,:)),y,'median');
    fprintf('%s\n',['The median expression of ', markers, ' in each cluster is '])
    median_levels
end


