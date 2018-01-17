function lineageIfo = inferingLineage(sEcell,ydata,clusterIfo,alpha,thresh)
% reconstruct the pseudotime 
% Inputs: 
%   scEcell : m x 1 vector, single cell energy
%   ydata : m x nPC, nPC-D coordindates from dimension reduction
%   clusterIfo : a struct giving the cell cluster information
%   alpha : significance level for a two-sided Wilcoxon rank-sum test of the distribution of scEnergy between two metacell,default= 0.01
%   thresh : proportion of cells included in the metacell, default= 0.8
% Outputs:
%   lineageIfo: a struct giving the cell lineage information
%   lineageIfo.TP: transition probability TP
%   lineageIfo.MDST: minimal directed spanning tree, i.e.,the inferred lineage
%   lineageIfo.path: the node in each path
%   lineageIfo.PDG: the constructed probabilistic directed graph
%   lineageIfo.centroid: the centroid of each core cell
if ~exist('alpha', 'var') || isempty(alpha)
    alpha = 0.01;
end
if ~exist('thresh', 'var') || isempty(thresh)
    thresh = 0.8;
end

idxCluster = clusterIfo.idxCluster;
numCluster = length(unique(clusterIfo.identity));
%% (1) find the core cell in each cluster
% thresh = 0.8; % select the cells which include 80% energy in each cluster
centroidCoreCell = zeros(numCluster,size(ydata,2));FCoreCell = zeros(numCluster,1);
idxClusterCorecell = cell(1,numCluster);  ydataCoreCell = cell(1,numCluster);
if thresh < 1
    for i = 1:numCluster
        FcellCi = sEcell(idxCluster{i});
        ydataCi = ydata(idxCluster{i},:);
        [FcellCiOrdered,idx] = sort(FcellCi);
        idxCoreCellCi = find(cumsum(FcellCiOrdered)/sum(FcellCiOrdered) >= thresh,1);
        idxCoreCellCi = 1:idxCoreCellCi;
        index = idxCluster{i}(idx);
        idxClusterCorecell{i} = index(idxCoreCellCi);
        ydataCiOrdered = ydataCi(idx,:);
        centroidCoreCell(i,:) = mean(ydataCiOrdered(idxCoreCellCi,:));
        ydataCoreCell{i} = ydataCiOrdered(idxCoreCellCi,:);
        %     FCoreCell(i) = median(FcellCiOrdered(idx));
        Q = quantile(FcellCiOrdered(idxCoreCellCi),[0.25, 0.5, 0.75]);
        FCoreCell(i) = 0.25*Q(1)+0.5*Q(2)+0.25*Q(3); % Tukey's trimean
    end
else
    for i = 1:numCluster
        centroidCoreCell(i,:) = mean(ydata(idxCluster{i},:));
        ydataCoreCell{i} = ydata(idxCluster{i},:);
        Q = quantile(sEcell(idxCluster{i}),[0.25, 0.5, 0.75]);
        FCoreCell(i) = 0.25*Q(1)+0.5*Q(2)+0.25*Q(3); % Tukey's trimean
    end
end

%% (2) calculate the Transition Probability (TP)
% (2.1) based on distance in low-dimension space
D1 = pdist(centroidCoreCell);
epsilon = 1/2*max(D1);
P1 = exp(-D1.^2/epsilon^2);
P1 = squareform(P1);
wx = sum(P1,2); % sum of weights of each vertice
wg = sum(wx); % sum of all weights of the vertices
sd = wx/wg; % stationary distribution
% transition matrix
P1 = P1./(repmat(wx,1,size(P1,2)));
% symmetric transition matirx
P1 = sqrt(repmat(sd,1,size(P1,2))).*P1./(sqrt(repmat(sd',size(P1,1),1)));
% (2.2) based on the energy
% calculate the probability of each meatacell
PCorecell = exp(-FCoreCell)/sum(exp(-FCoreCell));
P2 = PCorecell;
% (2.3) define the likehood (probability) of transition by combing these two
TP = diag(1-P2)*P1+diag(P2);

%% (3) construct a probabilistic directed graph using the energy information
C = nchoosek(1:numCluster,2);
A = zeros(numCluster);
for i = 1:size(C,1)
    p = ranksum(sEcell(idxCluster{C(i,1)}),sEcell(idxCluster{C(i,2)}));    
    if FCoreCell(C(i,1)) > FCoreCell(C(i,2)) && p < alpha
        A(C(i,1),C(i,2)) = 1-TP(C(i,1),C(i,2))+eps;
    else
        A(C(i,2),C(i,1)) = 1-TP(C(i,2),C(i,1))+eps;
        A(C(i,1),C(i,2)) = 1-TP(C(i,1),C(i,2))+eps;   
    end
end

%% (4) find the minimal directed spanning tree
rootNode = 1;
MDST = FindMDST(A,rootNode,0);
MDST = digraph(MDST);

% find the node in each path
idxS = rootNode;
idxT = find(outdegree(MDST) == 0);
path = cell(1,length(idxT));
for i = 1:length(idxT)
    path{i} = shortestpath(MDST,idxS,idxT(i),'Method','unweighted');
end
% restore the edge weight using transition probability TP
MDST.Edges.Weight = 1-MDST.Edges.Weight;
MDST = addedge(MDST,1:numCluster,1:numCluster,P2); % add the propability of clusters

lineageIfo.TP = TP; % transition probability TP
lineageIfo.MDST = MDST; % minimal directed spanning tree, i.e.,the inferred lineage
lineageIfo.path = path; % the node in each path
lineageIfo.PDG = A; % the constructed probabilistic directed graph
lineageIfo.centroid = centroidCoreCell; % the centroid of each core cell

