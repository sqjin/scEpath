function ydata = ECA(E)
% perform principal component analysis of the energy matrix
% Inputs: 
%   E : energy matrix
% Outputs:
%   ydata : m x 2, 2-D coordindates from dimension reduction
Ezscored = zscore(full(E));
[~,score] = pca(Ezscored','NumComponents',3);
ydata = score(:,1:2); 