function proData = preprocessing(iniData, minCells, minGenes, logNormalize,filterRibo)
% Preprocess scRNAseq data using the following steps
% (1) Filter out low-quality cells
% (2) Filter out low-expressed genes
% (3) log2 tranformation
% Inputs:
%   iniData: a struct variable, store the loaded data
%   minCells: a value, filter genes that are expressed in less than minCells cells; default = 3
%   minGenes: a value, filter low-quality cells in which the number of expressed genes is less than minGenes; default = 200
%   logNormalize: boolean, to do log2 normalization or not, default= 1 (true)
%   filterRibo: boolean, to filter not ERCC and Ribosomal genes or not, default= 1 (true)
% Outputs:
%   proData: a struct variable, store the data after QC
%   proData.data: a matrix giving the single cell data (rows are cells and columns are genes)
%   proData.genes: a cell array giving the gene names
%   proData.cells: a cell array, each cell giving cell attributes (such as cell type, culture condition, day captured)
if ~exist('minGenes','var') || isempty(minGenes)
    minGenes = 100;
end
if ~exist('minCells','var') || isempty(minCells)
    minCells = 3;
end
if ~exist('logNormalize','var') || isempty(logNormalize)
    logNormalize = 1;
end
if ~exist('filterRibo','var') || isempty(filterRibo)
    filterRibo = 1;
end

data0 = iniData.data; gene0 = iniData.genes;
%% filter ERCC genes and Ribosomal genes
if filterRibo == 1
    idxERCC = find(cellfun(@isempty,strfind(gene0,'ERCC-')) == 0);
    idxRPL = find(cellfun(@isempty,strfind(gene0,'Rpl')) == 0 | cellfun(@isempty,strfind(gene0,'RPL')) == 0);
    idxRPS = find(cellfun(@isempty,strfind(gene0,'Rps')) == 0 | cellfun(@isempty,strfind(gene0,'RPS')) == 0);
    idxMRPL = find(cellfun(@isempty,strfind(gene0,'Mrpl')) == 0 | cellfun(@isempty,strfind(gene0,'MRPL')) == 0);
    idxMRPS = find(cellfun(@isempty,strfind(gene0,'Mrps')) == 0 | cellfun(@isempty,strfind(gene0,'MRPS')) == 0);
    idxGeneFiltered = [idxERCC;idxRPL;idxRPS;idxMRPL;idxMRPS];
    gene0(idxGeneFiltered) = []; data0(idxGeneFiltered,:) = [];
end
%% filter cells that have expressed genes less than #100
dataTemp = data0;
dataTemp(find(data0 > 0)) = 1;
msum = sum(dataTemp,1);
data0(:,msum < minGenes) = [];
iniData.cells(msum < minGenes,:) = [];
%% filter genes that only express less than #3 cells
dataTemp = data0;
dataTemp(find(data0 > 0)) = 1;
nsum = sum(dataTemp,2);
data0(nsum < minCells,:) = [];gene0(nsum < minCells) = [];%%

%% normalization
if logNormalize
    data0 = log2(data0+1);
end

proData.data = data0;
proData.genes = gene0;
proData.cells = iniData.cells;
