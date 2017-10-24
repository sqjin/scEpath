function MDST = FindMDST(G,rootNode,showTree)
%% Find minimal directed spanning tree in directed graph
% Input G: adjacent matrix with edge weights for directed graph
%       rootNode: the index of root node.
% Output MDST: the weighted adjacent matrix of MDST (in sparse form)
% e.g.,G = [0 17 14 0;0 0 16 0;0 0 0 0;23 31 15 0]; 
if nargin < 3
    showTree = 1;
end

V = 1:size(G,1); % set of vertices (vertices vector)
E = incidence_to_3n(G);% convert weighted adjacent matrix (incidence  matrix)into 3 column format
GT = edmonds(V,E,rootNode); 
% GT is a struct, where GT.V is the vertices vector,GT.E is the edge list with weights
%                       GT.BV is the order of vertices in the MDST,GT.BE is the order of edges in the MDST

%% show the minimal weighted directed spanning tree
% output MDSTedgeOrder is the order of edges of MDST in the edges of G
MDSTedgeOrder = reconstruct_2(GT,showTree); 
%% extract the weighted adjacent matrix of MDST (in sparse form)
MDST = GT(1).E(MDSTedgeOrder,:);
MDST = sparse(MDST(:,1),MDST(:,2),MDST(:,3));
MDST(max(size(MDST)),max(size(MDST))) = 0;