function ydata = ECA(E)
% perform principal component analysis of the energy matrix
% Inputs: 
%   E : energy matrix
% Outputs:
%   ydata : m x nPC, nPC-D coordindates from dimension reduction
Ezscored = zscore(full(E));
s = svd(Ezscored');
y = s;
beta = min(size(E,2)/size(E,1),size(E,1)/size(E,2));
coef = optimal_SVHT_coef(beta, 0)*median(s);
y( y < coef ) = 0;
nPC = nnz(y);
[~,score,~,~,explained] = pca(Ezscored','NumComponents',nPC);
ydata = score; 
