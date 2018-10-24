function pseudotimeIfo = inferingPseudotime(Rscript,ydata,lineageIfo,clusterIfo,theta)
% reconstruct the pseudotime
% Inputs:
%   Rscript : path for Rscript
%   ydata : m x 2, 2-D coordindates from dimension reduction
%   lineageIfo : a struct giving the cell lineage information (returned from the function inferingLineage.m)
%   clusterIfo : a struct giving the cell cluster information
%   theta : proportion of cells used to fitting the principal curve,default= 0.75
% Outputs:
%   pseudotimeIfo: a struct giving the pseudotime information
%   pseudotimeIfo.pseudotime: a cell array, each cell gives pseudotime value for each path
%   pseudotimeIfo.cellOrder: a cell array, each cell gives cell order for each path
%   pseudotimeIfo.cellIndex: a cell array, each cell gives cell index in each path
%   pseudotimeIfo.Pcurve: a cell array, each cell gives projection vectors of fitted principal curve in each path
if ~exist('theta', 'var') || isempty(theta)
    theta = 0.75;
end
ydata = ydata(:,1:2);
centroidCoreCell = lineageIfo.centroid(:,1:2); path = lineageIfo.path;
idxCluster = clusterIfo.idxCluster; group = clusterIfo.identity;

% identify the main cells in each group to infer pesudotime
numCluster = size(centroidCoreCell,1);
mainCell = cell(1,numCluster);outCell = mainCell;
for i = 1:numCluster
    d = pdist([centroidCoreCell(i,:);ydata(idxCluster{i},:)]);
    d = d(1:length(idxCluster{i}));
    thresh = quantile(d,theta);
    mainCell{i} = idxCluster{i}(d <= thresh);
    outCell{i} = idxCluster{i}(d > thresh);
end

filefolder = fullfile('results','temporalfiles');
if ~exist(fullfile(pwd,filefolder),'dir')
    mkdir(filefolder)
end
cellIndexpath = cell(1,length(path));
for j = 1:length(path)
    ydataFpathj = [];cellwpathj = [];centroidCoreCellpathj = [];ydataFOutCellpathj = [];
    maincellIndexpathj = [];outcellIndexpathj = [];
    for i = 1:length(path{j})
        ydataFpathj = [ydataFpathj;ydata(mainCell{path{j}(i)},:)];
        centroidCoreCellpathj = [centroidCoreCellpathj;centroidCoreCell(path{j}(i),:)];
        ydataFOutCellpathj = [ydataFOutCellpathj;ydata(outCell{path{j}(i)},:)];
        maincellIndexpathj = [maincellIndexpathj;mainCell{path{j}(i)}];
        outcellIndexpathj = [outcellIndexpathj;outCell{path{j}(i)}];
    end
    cellIndexpath{j} = [maincellIndexpathj;outcellIndexpathj];
    dlmwrite(fullfile(filefolder,'ydata.txt'),ydataFpathj,'delimiter','\t','precision','%.4f');
    dlmwrite(fullfile(filefolder,'ydataOutCell.txt'),ydataFOutCellpathj,'delimiter','\t','precision','%.4f');
    dlmwrite(fullfile(filefolder,'pathLength.txt'),length(path{j}),'delimiter','\t');
    % Calling R's principal curve
    RscriptFileName = ' ./principal_curve_fitting.R '; % the full path of the R script to execute. Default: current directory is ./scEpath-master/
    % User can change the default path to the path where codes are located, e.g., RscriptFileName = ' /Users/XXX/Downloads/scEpath-master/principal_curve_fitting.R '; 
    eval([' system([', '''', Rscript, RscriptFileName, '''', ' filefolder]);']);
    try
        pseudotimeMainCell = importdata(fullfile(filefolder,'PcurveLambdaMainCell.txt'));
        projectionsPcurveMaincell = importdata(fullfile(filefolder,'PcurveProjectionValueMainCell.txt'));
        pseudotimeOutCell = importdata(fullfile(filefolder,'PcurveLambdaOutCell.txt'));
        projectionsPcurveOutCell = importdata(fullfile(filefolder,'PcurveProjectionValueOutCell.txt'));
    catch
        error(sprintf('Error!!! Please check if the R package "princurve" has been successfully installed!'))
    end
    dlmwrite(fullfile(filefolder,['PcurveLambdaMainCellPath' num2str(j) '.txt']),pseudotimeMainCell,'delimiter','\t','precision','%.4f');
    dlmwrite(fullfile(filefolder,['PcurveProjectionValueMainCellPath' num2str(j) '.txt']),projectionsPcurveMaincell,'delimiter','\t','precision','%.4f');
    dlmwrite(fullfile(filefolder,['PcurveLambdaOutCellPath' num2str(j) '.txt']),pseudotimeOutCell,'delimiter','\t','precision','%.4f');
    dlmwrite(fullfile(filefolder,['PcurveProjectionValueOutCellPath' num2str(j) '.txt']),projectionsPcurveOutCell,'delimiter','\t','precision','%.4f');
end

% rescontruct the pseudotime
projectionsPcurvePath = cell(1,length(path));cellOrderPath = cell(1,length(path));pseudotimePath = cell(1,length(path));
for j = 1:length(path)
    pseudotimeMainCell = importdata(fullfile(filefolder,['PcurveLambdaMainCellPath' num2str(j) '.txt']));
    projectionsPcurveMaincell = importdata(fullfile(filefolder,['PcurveProjectionValueMainCellPath' num2str(j) '.txt']));
    pseudotimeOutCell = importdata(fullfile(filefolder,['PcurveLambdaOutCellPath' num2str(j) '.txt']));
    projectionsPcurveOutCell = importdata(fullfile(filefolder,['PcurveProjectionValueOutCellPath' num2str(j) '.txt']));
    pseudotime = [pseudotimeMainCell;pseudotimeOutCell];
    projectionsPcurve = [projectionsPcurveMaincell;projectionsPcurveOutCell];
    [~, cellOrder] = sort(pseudotime);% ordering the projection index
    % judge which is the start state
    s = cellIndexpath{j}(cellOrder(1:floor(0.2*length(cellOrder))));
    t = cellIndexpath{j}(cellOrder(floor((1-0.2)*length(cellOrder)):end));
    sNum = length(intersect(s,find(group == lineageIfo.rootNode))); tNum = length(intersect(t,find(group == lineageIfo.rootNode)));
    if sNum < tNum 
        cellOrder = flipud(cellOrder);
        pseudotime = max(pseudotime(cellOrder))-pseudotime(cellOrder);
    else
        pseudotime = pseudotime(cellOrder);
    end
    projectionsPcurve = projectionsPcurve(cellOrder,:);
    cellOrder = cellIndexpath{j}(cellOrder);
    pseudotime = (pseudotime- min(pseudotime))/(max(pseudotime)-min(pseudotime));
    projectionsPcurvePath{j} = projectionsPcurve;
    cellOrderPath{j} = cellOrder;
    pseudotimePath{j} = pseudotime;
end

pseudotimeIfo.pseudotime = pseudotimePath; % a cell array, each cell gives pseudotime value for each path
pseudotimeIfo.cellOrder = cellOrderPath; % a cell array, each cell gives cell order for each path
pseudotimeIfo.cellIndex = cellIndexpath; % a cell array, each cell gives cell index in each path
pseudotimeIfo.Pcurve = projectionsPcurvePath; % a cell array, each cell gives projection vectors of fitted principal curve in each path
