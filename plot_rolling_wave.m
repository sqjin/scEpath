function PDGIfo = plot_rolling_wave(PDGIfo,smoothExprIfo,proData,optimalK,pathUsed,fig_width,fig_height)
% Inputs:
%    PDGIfo.PDG:a cell array, each cell gives the identified pseudotime-dependent genes in each branch
%    smoothExprIfo:a cell array, each cell gives the smoothed expression levels of all genes in each branch
%    optimalK: the number of desired genes clusters
%    pathUsed:order the pseudotime-dependent genes based on their peak expression in the "pathUsed" branch. e.g. pathUsed = 1
%    fig_width : the figure width
%    fig_height : the figure height
% Outputs:
%    PDGIfo: updated information of the pseudotime dependent genes
%    PDGIfo.smoothExpr: a cell array, each cell contains the smooth expression of pseudotime dependent genes in each branch
%    PDGIfo.PDGall: a cell array contains the combined pseudotime dependent genes in all the branches
%    PDGIfo.orderedPDG: a cell array contains the ordered pseudotime dependent genes
%    PDGIfo.PDGOrder: 1 x optimalK cell array, each cell contains the index of ordered pseudotime dependent genes
%    gene list in each identified gene cluster. .txt files are saved in the folder called PDGclusters
% Figures:
%      (i) Rolling wave plot showing the normalized-smoothed expression pattern of pseudotime-dependent genes
%      (ii)  Average expression pattern of the identified gene clusters along pseudotime in each branch
if ~exist('optimalK','var') || isempty(optimalK)
    optimalK = 8;
end
if ~exist('pathUsed','var') || isempty(pathUsed)
    pathUsed = 1;
end

PDGPath = PDGIfo.PDG;
% combine the pseudotime dependent genes (PDG) in each lineage(path)
PDGall = [];
for j = 1:length(PDGPath)
    PDGall = [PDGall; PDGPath{j}];
end
PDGall = unique(PDGall,'stable');
PDG = PDGall;
[PDG,~,idxPDG] = intersect(PDG,proData.genes,'stable');
smoothCurve = [];smoothCurveIndi = cell(1,length(PDGPath));
for j = 1:length(PDGPath)
    smoothCurveIndi{j} = smoothExprIfo{j}.smoothExpr(idxPDG,:);
    smoothCurve = [smoothCurve, smoothCurveIndi{j}];
end
PDGIfo.smoothExpr = smoothCurveIndi;
PDGIfo.PDGall = PDGall;

folderName = fullfile('results','figures');
if ~exist(folderName, 'dir')
    mkdir(folderName);
end

%% performing hierarchical clustering of temporal dynamics using smoothed gene expression
smoothCurveZ = zscore(smoothCurve,[],2);
Z = linkage(smoothCurveZ,'ward','correlation');
groupDEG = cluster(Z,'maxclust',optimalK);
thresh = Z(end-optimalK+2,3)-eps;
cgo = clustergram(smoothCurveZ,'Standardize',3,'Linkage','ward','RowPDist','correlation','Cluster','column','ColumnLabels',{},'RowLabels',{},'Colormap',redbluecmap,'Dendrogram',thresh);
hFig1 = plot(cgo);
saveas(hFig1,fullfile(folderName,'hierarchical_clustering_pseudotime_depedent_genes.pdf'))

clear smoothCurve
% order the genes based on the peak order in pseudotime
smoothCurve = smoothCurveIndi{pathUsed};
xii = linspace(0,1,size(smoothCurve,2));
locMax = zeros(1,optimalK);locDecay = locMax;
groupCenter = zeros(optimalK,size(smoothCurve,2));
for i = 1:optimalK
    Q = quantile(smoothCurve(groupDEG == i,:),[0.25, 0.5, 0.75]);
    groupCenter(i,:) = 0.25*Q(1,:)+0.5*Q(2,:)+0.25*Q(3,:);
    [groupMax,idxMax] = max(groupCenter(i,:));
    [~,idxDecay] = find(abs(groupCenter(i,:)-(min(groupCenter(i,:))+range(groupCenter(i,:))/2))<0.05*max(groupCenter(i,:)),1,'last');
    locMax(i) = xii(idxMax);
    if ~isempty(idxDecay)
        locDecay(i) = xii(idxDecay);
    else
        locDecay(i) = locMax(i);
    end
end

[locMaxS,idx1] = sort(locMax);idx2 = find(diff(locMaxS) <= 0.1);
for ii = 1:length(idx2)
    locMax(idx1(idx2(ii))) = mean([locMax(idx1(idx2(ii))),locMax(idx1(idx2(ii)+1))]);
    locMax(idx1(idx2(ii)+1)) = locMax(idx1(idx2(ii)));
end

[count,locMaxUni] = hist(locMax,unique(locMax));
idx = find(count > 1);
if ~isempty(idx)
    for i = 1:length(idx)
        idxEqual = find(ismember(locMax,locMaxUni(idx(i))));
        
        [locSort,idxDecay] = sort(locDecay(idxEqual));
        if i == 1
            locMaxCorrect = locMax;
            locMaxCorrect(idxEqual(idxDecay)) = locMax(idxEqual(idxDecay))+[1:length(idxEqual)]*(10^-4);
        else
            locMaxCorrect(idxEqual(idxDecay)) = locMaxCorrect(idxEqual(idxDecay))+[1:length(idxEqual)]*(10^-4);
        end
    end
else
    locMaxCorrect = locMax;
end

[~,idxMax] = sort(locMaxCorrect);

idxGroupgene = cell(1,optimalK);
for i = 1:optimalK
    idxGroupgene{i} = find(groupDEG == idxMax(i));
end
for i = 1:optimalK
    groupDEG(idxGroupgene{i}) = i;
end

idxGeneOrder = cell(1,optimalK);
for i = 1:optimalK
    locMax = zeros(1,length(idxGroupgene{i}));locDecay = locMax;
    for j = 1:length(idxGroupgene{i})
        [~,idxMax] = max(smoothCurve(idxGroupgene{i}(j),:));
        [~,idxDecay] = find(abs(smoothCurve(idxGroupgene{i}(j),:)-(min(smoothCurve(idxGroupgene{i}(j),:))+range(smoothCurve(idxGroupgene{i}(j),:))/2))<0.05*max(smoothCurve(idxGroupgene{i}(j),:)),1,'last');
        locMax(j) = xii(idxMax);%locDecay(j) = xii(idxDecay);
    end
    [locSort,idxMax] = sort(locMax);
    %     idxEqual = find(diff(locSort) == 0)+1; idxEqual = [idxEqual(1)-1,idxEqual];
    %     idxEqual = idxMax(idxEqual);
    %     [locSort,idxDecay] = sort(locDecay(idxEqual));
    locMaxCorrect = locMax;
    %     locMaxCorrect(idxEqual) = locMax(idxEqual(idxDecay))+[1:length(idxEqual)]*(10^-4);
    [~,idx] = sort(locMaxCorrect);
    idxGeneOrder{i} = idxGroupgene{i}(idx);
end

%¡¡save the gene in each cluster
folderName = fullfile('results','PDG_in_each_cluster');
if ~exist(folderName, 'dir')
    mkdir(folderName);
end
PDGordered = [];PDGorderIdx = [];
for i = 1:optimalK
    PDGordered = [PDGordered; PDG(idxGeneOrder{i})];
    PDGorderIdx = [PDGorderIdx; idxGeneOrder{i}];
    fileName = fullfile(folderName, ['PDGcluster' num2str(i) '.txt']);
    T = cell2table(PDG(idxGeneOrder{i}));
    writetable(T,fileName,'Delimiter','\t','WriteVariableNames',false)
end
PDGIfo.orderedPDG = PDGordered;
PDGIfo.PDGOrder = idxGeneOrder;


if ~exist('fig_width', 'var') || isempty(fig_width)
    fig_width = 150*length(smoothExprIfo);
end
if ~exist('fig_height', 'var') || isempty(fig_height)
    fig_height = 80*optimalK;
end

% plot the heatmap
colorClusters = jet(optimalK);
hFig1 = figure('position', [600, 0, fig_width, fig_height]);
for j = 1:length(smoothCurveIndi)
    smoothCurve = smoothCurveIndi{j};
    dataHeatmapDEG = [];
    for i = 1:optimalK
        dataHeatmapDEG = [dataHeatmapDEG; smoothCurve(idxGeneOrder{i},:)];
    end
    subplot(1,10*length(smoothCurveIndi),[1+8*(j-1) 8+8*(j-1)])
    dataHeatmapDEGZscore = zscore(dataHeatmapDEG,[],2);
    imagesc(dataHeatmapDEGZscore);
    xlabel('Pseudotime','FontName','Arial','FontSize',10)
    set(gca,'FontName','Arial','FontSize',8)
    set(gca,'Xtick',[1 size(smoothCurve,2)/2 size(smoothCurve,2)]);
    set(gca,'XtickLabel',0:0.5:1);
    set(gca,'Ytick',[]);
    if length(PDGordered) < 100
     set(gca,'Ytick',1:1:length(PDGordered));
     set(gca,'YtickLabel',PDGordered,'FontName','Arial','FontSize',8)
    end

    colormap(jet)
    title(['Path ' num2str(j)],'FontName','Arial','FontSize',10);
    if j == length(smoothCurveIndi)
        c = colorbar;
        c.Location = 'east';
        c.FontSize = 8;
        %         c.Position = [0.85 .1 .02 .2];
        c.Position = [0.8 .11 .02 .2];
    end
    caxis([-3 3])
    if j == 1 & length(smoothCurveIndi) > 1
        set(gca, 'XDir', 'reverse')
    end
    if j == length(smoothCurveIndi)
        subplot(1,10*length(smoothCurveIndi),[1+8*j 8*j+2])
        ybar = [];
        for i = optimalK:-1:1
            ybar = [ybar,nnz(groupDEG == i)];
        end
        H = bar([ybar;nan(1,length(ybar))],'stacked','BarWidth',0.0004,'ShowBaseLine','off','LineWidth',0.01);
        % colormap(h,flipud(colormap(jet)))
        k = 0;
        for i = optimalK:-1:1
            k = k+1;
            set(H(i),'facecolor',colorClusters(k,:))
        end
        axis tight
        xlim([1 1.001])
        set(gca,'xtick',1)
        axis off
    end
end
folderName = fullfile('results','figures');
if ~exist(folderName, 'dir')
    mkdir(folderName);
end
saveas(hFig1,fullfile(folderName,'rolling_wave_pseudotime_depedent_genes.pdf'))


hFig3 = figure('position', [600, 0, fig_width/length(smoothExprIfo), fig_height]);
value_max = zeros(length(smoothCurveIndi),optimalK); value_min = value_max;
for j = 1:length(smoothCurveIndi)
    smoothCurve = smoothCurveIndi{j};
    xii = linspace(0,1,size(smoothCurve,2));
    groupCenter = zeros(optimalK,size(smoothCurve,2));
    for i = 1:optimalK
        H(i) = subplot(optimalK,1,i);
        hold on
        Q = quantile(smoothCurve(groupDEG == i,:),[0.25, 0.5, 0.75]);
        groupCenter(i,:) = 0.25*Q(1,:)+0.5*Q(2,:)+0.25*Q(3,:);
        value_max(j,i) = max(groupCenter(i,:));value_min(j,i) = min(groupCenter(i,:));
        
        if j == 1
            if i == 1
                h1 = plot(xii,groupCenter(i,:),'LineWidth',1.5,'color',colorClusters(i,:),'LineStyle','-');
            else
                plot(xii,groupCenter(i,:),'LineWidth',1.5,'color',colorClusters(i,:),'LineStyle','-');
            end
        else
            if i == 1
                h2 = plot(xii,groupCenter(i,:),'LineWidth',1.5,'color',colorClusters(i,:),'LineStyle','-.');
            else
                plot(xii,groupCenter(i,:),'LineWidth',1.5,'color',colorClusters(i,:),'LineStyle','-.');
            end
        end
      %  set(gca,'Xtick',[]);
        box on
        set(gca,'FontName','Arial','FontSize',8)
        ylim([0.95*min(value_min(:,i)) 1.05*max(value_max(:,i))])
        yticks('auto')
        yticklabels('auto')
        box off
        set(gca,'Xtick',[]);

        if i == optimalK
            set(gca,'Xtick',[min(xii) range(xii)/2 max(xii)]);
            set(gca,'XtickLabel',0:0.5:1);
            xlabel('Pseudotime','FontName','Arial','FontSize',10)
        end
        % set(gca, 'XDir', 'reverse')
    end
    legendText{j} = ['Path ', num2str(j)];
end
if length(smoothCurveIndi) > 1
    legend([h1,h2],legendText)
end
saveas(hFig3,fullfile(folderName,'temporal_patterns_pseudotime_depedent_genes.pdf'))

