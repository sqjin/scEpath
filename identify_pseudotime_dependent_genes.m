function PDGIfo = identify_pseudotime_dependent_genes(proData,smoothExprIfo,pseudotimeIfo,sd_thresh,sig_thresh,nboot)
% identify pseudotime dependent genes
% Inputs:
%   proData : the single cell data including genes and cell attributes
%   smoothExprIfo : the smooth expression information
%   pseudotimeIfo : the pseudotime information
%   sd_thresh : threshold of standard deviation, default = 0.5
%   sig_thresh: threshold of adjusted P-values, default = 0.01
%   nboot : the number of bootstraps, default = 1000
% Outputs:
%   PDGIfo: a struct giving the pseudotime dependent genes
%   PDGIfo.PDG : a cell array, each cell contains the identified pseudotime-dependent genes for each path
%   PDGIfo.allGenes: a cell array, each cell contains the genes with their calculated adjusted P-values and standard deviation
if ~exist('sd_thresh','var') || isempty(sd_thresh)
    sd_thresh = 0.05;
end
if ~exist('sig_thresh','var') || isempty(sig_thresh)
    sig_thresh = 0.01;
end
if ~exist('nboot','var') || isempty(nboot)
    nboot = 1000;
end

PDGPath = cell(1,length(smoothExprIfo));Tpath = PDGPath;
for j = 1:length(smoothExprIfo)
    pseudotime = pseudotimeIfo.pseudotime{j};
    cellOrder = pseudotimeIfo.cellOrder{j};
    smoothExprnull = smoothExprIfo{j}.aveExpr;
    nbins = smoothExprIfo{j}.nbins;
    SDnull = nanstd(smoothExprnull,[],2);
    SDboot = zeros(size(proData.data,1),nboot);
    for nE = 1:nboot
        index = randperm(length(cellOrder));
        bootcellOrder = cellOrder(index);
        smoothExprIfoboot = smootheningExpr(proData.data,pseudotime,bootcellOrder,nbins);
        smoothExprboot = smoothExprIfoboot.aveExpr;
        SDboot(:,nE) = nanstd(smoothExprboot,[],2);
    end
    nReject = sum(bsxfun(@ge,SDboot,SDnull),2);
    p = nReject/nboot;
    % compute the adjusted P-values (BH correction)
    padj = mafdr(p,'BHFDR',true);
    padj(SDnull < sd_thresh) = NaN; % can use a smaller variance threshold, e.g. 0.05
    Tpath{j} = table(proData.genes,padj,SDnull,'VariableNames',{'genes','padj','SD'});
    Tpath{j} = sortrows(Tpath{j},{'padj','SD'},{'ascend','descend'});
    PDG = Tpath{j}.genes(Tpath{j}.padj < sig_thresh);
    PDGPath{j} = PDG;
end

PDGIfo.PDG = PDGPath;
PDGIfo.allGenes = Tpath;