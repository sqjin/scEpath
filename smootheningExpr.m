function smoothExprIfo = smootheningExpr(data,pseudotime,cellOrder,nbins,cellw)
% calculate the smooth gene expression based on the pseudotime
% Inputs:
%   data : the single cell data matrix
%   pseudotime : the ordered pseudotime
%   cellOrder : the cell order
%   nbins : the number of bins for dividing the pseudotime, default = 10
%   cellw : 1 x m vector, weights for cubic spline fitting, default = [1,...,1].
% Outputs: 
%   smoothExprIfo: a struct variable
%   smoothExprIfo.smoothExpr: smooth gene expression in each cell
%   smoothExprIfo.aveExpr: average gene expression in each bin
%   smoothExprIfo.avecellw: average weights in each bin
%   smoothExprIfo.nbins: the number of used bins
if ~exist('nbins','var') || isempty(nbins)
    nbins = 10;
end
if ~exist('cellw','var') || isempty(cellw)
    cellw = ones(1,size(data,2));
end
t0 = pseudotime(1);
t1 = pseudotime(end);
binHalfSize = (t1-t0)/2/nbins;
boundariesAndCenters = t0:binHalfSize:t1;
binCenters = boundariesAndCenters([1 2:2:end-1 end]);
binWidths = 0.08*(t1-t0);
expr = data(:,cellOrder); cellw = cellw(cellOrder);

aveExpr = zeros(size(expr,1),length(binCenters));
avecellw = zeros(1,length(binCenters));
index0 = [];
for index = 1:length(binCenters)
    cellLocations = pseudotime > binCenters(index) - binWidths & pseudotime < binCenters(index) + binWidths;
    if nnz(cellLocations) < 5
        index0 = [index0,index];
        continue;
    end
    Q = quantile(expr(:,cellLocations),[0.25, 0.5, 0.75],2);
    aveExpr(:,index) = 0.25*Q(:,1)+0.5*Q(:,2)+0.25*Q(:,3);
    avecellw(index) = mean(cellw(cellLocations));
end
if ~isempty(index0)
    for i = 1:length(index0)
        try
        aveExpr(:,index0(i)) = mean([aveExpr(:,index0(i)-1),aveExpr(:,index0(i)+1),aveExpr(:,index0(i)+2)],2);
        catch
        continue;
        end
    end
end
% fitting cubic curve
pp = csaps(binCenters,aveExpr,[],[],avecellw);
xii = linspace(min(binCenters),max(binCenters),length(pseudotime));
smoothExpr = fnval(pp,xii);
smoothExpr(smoothExpr < 0) = 0;
smoothExprIfo.smoothExpr = smoothExpr;
smoothExprIfo.aveExpr = aveExpr;
smoothExprIfo.avecellw = avecellw;
smoothExprIfo.nbins = nbins;
