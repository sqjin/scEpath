function [E,Ecell] = estimatingscEnergy(data,networkIfo,normalization)
% estimate the single cell energy
% Inputs:
%    data: single cell data (rows are genes and columns are cells)
%    networkIfo: network information for the constructed gene-gene network
%    normalization: boolean, to normalize the single cell energy, default= 1 (true)
% Outputs:
%    E: energy matrix, the same dimension with data
%    Ecell; scEnergy of each cell
if ~exist('normalization','var') || isempty(normalization)
    normalization = 1;
end

% select the genes that form a giant connected component in the constructed network
data = data(networkIfo.IDselect,:);
% rescaling the data into [0 1]
emin = min(data);emax = max(data);
data = (data-repmat(emin,size(data,1),1))./repmat((emax-emin),size(data,1),1);

%% load network information
A = networkIfo.R+networkIfo.R';
A(logical(eye(size(A)))) = 1;
dataSumNeibor = zeros(size(data));
for i = 1:size(A,1)
    a = A(i,:) ~= 0;
    dataSumNeibor(i,:) = sum(data(a,:));
end
%% calculate the  scEnergy for each cell
E = -data.*log(data./dataSumNeibor);
E(isnan(E)) = 0;
E(isinf(E)) = 0;
Ecell = sum(E);
E = sparse(E);
if normalization
    Ecell = (Ecell/mean(Ecell)).^2./(1+(Ecell/mean(Ecell)).^2);
end
