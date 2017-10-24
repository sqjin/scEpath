function PDG_TFIfo = identify_keyTF(proData,PDGordered,smoothExprIfo,lineageIfo,clusterIfo,TFnames,thresh_SD_TF,thresh_FC_TF,thresh_sig_TF)
% identify pseudotime dependent genes
% Inputs:
%   proData : the single cell data including genes and cell attributes
%   PDGordered: ordered pseudotime dependent genes
%   smoothExprIfo : the smooth expression information
%   lineageIfo : the lineage information
%   clusterIfo : the cluster information
%   TFnames: collected TF names
%   thresh_SD_TF : threshold of standard deviation, default = 0.5
%   thresh_FC_TF : threshold of log2fold change, default = 1
%   thresh_sig_TF: threshold of adjusted P-values, default = 0.01
% Outputs:
%   PDGIfo: a struct variable
%   PDG_TFIfo.PDG_TFSig : a cell array, each cell contains key TFs that may be responsible for cell fate decisions
%   PDG_TFIfo.PDG_TF : a cell array, each cell contains all the TFs found in the pseudotime-dependent genes

if ~exist('thresh_SD_TF','var') || isempty(thresh_SD_TF)
    thresh_SD_TF = 0.05;
end
if ~exist('thresh_FC_TF','var') || isempty(thresh_FC_TF)
    thresh_FC_TF = 1;
end
if ~exist('thresh_sig_TF','var') || isempty(thresh_sig_TF)
    thresh_sig_TF = 0.01;
end

%% find the TF in each cluster
PDG_TFname = intersect(PDGordered,TFnames,'stable');
PDG_TFname0 = PDG_TFname;
[PDG_TFname,~,idx] = intersect(PDG_TFname,proData.genes,'stable');
PDG_TFnameSigIndi = cell(1,length(smoothExprIfo));
for j = 1:length(smoothExprIfo)
    SD = nanstd(smoothExprIfo{j}.aveExpr,[],2);
    SD_TF = SD(idx);
    idxSD = SD_TF > thresh_SD_TF;
    PDG_TFname = PDG_TFname(idxSD);
    
    [PDG_TFname,~,idx] = intersect(PDG_TFname,proData.genes,'stable');
    PDG_TFdata = proData.data(idx,:);
    
    % select the TF with significant difference among these subpopulstions
    group = clusterIfo.identity;
    node = lineageIfo.path{j};
    C = zeros(length(node)-1,2);
    for i = 1:length(node)-1
        C(i,:) = [node(i),node(i+1)];
    end
    pFCSigIdx = [];
    for i = 1:size(C,1)
        pvalues = mattest(PDG_TFdata(:,group == C(i,1)), PDG_TFdata(:,group == C(i,2)),'VarType','unequal');
        padj = mafdr(pvalues,'BHFDR',true);
        pSigIdx = padj < thresh_sig_TF;
        log2FC0 = log2(mean(PDG_TFdata(:,group == C(i,1)),2)./mean(PDG_TFdata(:,group == C(i,2)),2));
        log2FC0(isinf(log2FC0)) = 0;log2FC0(isnan(log2FC0)) = 0;
        FCSigIdx = abs(log2FC0) > thresh_FC_TF;
        pFCSigIdx = [pFCSigIdx;find(pSigIdx & FCSigIdx)];
    end
    PDG_TFnameSigIndi{j} = intersect(PDG_TFname,PDG_TFname(pFCSigIdx),'stable');
end

PDG_TFnameSig = [];
for j = 1:length(PDG_TFnameSigIndi)
    PDG_TFnameSig = [PDG_TFnameSig;PDG_TFnameSigIndi{j}];
end
PDG_TFnameSig = unique(PDG_TFnameSig);
[PDG_TFnameSig] = intersect(PDG_TFname,PDG_TFnameSig,'stable');
PDG_TFSig = PDG_TFnameSig;

PDG_TFIfo.PDG_TFSig = PDG_TFSig; % key TFs that may be responsible for cell fate decisions
PDG_TFIfo.PDG_TF = PDG_TFname0; % all the TFs found in the pseudotime-dependent genes