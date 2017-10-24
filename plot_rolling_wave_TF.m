function plot_rolling_wave_TF(PDG_TFIfo,PDGIfo,optimalK,show_TF_names,fig_width,fig_height)
% Inputs:
%      PDGIfo.PDG:a cell array, each cell gives the identified pseudotime-dependent genes in each branch
%      PDG_TFIfo:a struct variable contains identified TF
%      optimalK: the number of desired genes clusters
%      show_TF_names: boolean, to show the TF names as Yticklabels or not, default= 1 (true)
%    fig_width : the figure width
%    fig_height : the figure height

% Outputs: one figure
%      Rolling wave plot showing the normalized-smoothed expression pattern of key TF
if ~exist('optimalK','var') || isempty(optimalK)
    optimalK = 8;
end
if ~exist('show_TF_names','var') || isempty(show_TF_names)
    show_TF_names = 1;
end
if ~exist('fig_width', 'var') || isempty(fig_width)
    fig_width = 200*length(PDGIfo.smoothExpr);
end
if ~exist('fig_height', 'var') || isempty(fig_height)
    fig_height = 80*length(PDGIfo.PDGOrder);
end

PDG_TFname = PDG_TFIfo.PDG_TFSig;
smoothCurveIndi = PDGIfo.smoothExpr;
idxGeneOrder = PDGIfo.PDGOrder;

% plot the heatmap
colorClusters = jet(optimalK);
hFig = figure('position', [600, 0, fig_width, fig_height]);

for j = 1:length(smoothCurveIndi)
    smoothCurve = smoothCurveIndi{j};
    dataHeatmapPDG = [];
    for i = 1:optimalK
        dataHeatmapPDG = [dataHeatmapPDG; smoothCurve(idxGeneOrder{i},:)];
    end
    [PDG_TFname,~,idx] = intersect(PDG_TFname,PDGIfo.orderedPDG,'stable');
    DEG_TFsmoothdata = dataHeatmapPDG(idx,:);
    if j == 1
        subplot(1,10*length(smoothCurveIndi),[3+8*(j-1) 8+8*(j-1)])
    else
        subplot(1,10*length(smoothCurveIndi),[1+8*(j-1) 8+8*(j-1)-2])
    end
    PDG_TFsmoothdataZscore = zscore(DEG_TFsmoothdata,[],2);
    imagesc(PDG_TFsmoothdataZscore);
    xlabel('Pseudotime','FontName','Arial','FontSize',10)
    set(gca,'FontName','Arial','FontSize',8)
    set(gca,'Xtick',[1 size(smoothCurve,2)/2 size(smoothCurve,2)]);
    set(gca,'XtickLabel',0:0.5:1);
    set(gca,'Ytick',[]);
    if show_TF_names & j == 1
        set(gca,'Ytick',1:1:length(PDG_TFname));
        set(gca,'YtickLabel',PDG_TFname(1:1:end),'FontName','Arial','FontSize',8)
    else
        set(gca,'Ytick',[]);
    end
    colormap(jet)
    title(['Path ' num2str(j)],'FontName','Arial','FontSize',10);
    if j == length(smoothCurveIndi)
        c = colorbar;
        c.Location = 'east';
        %         c.Label.String = 'TF expression';
        %         c.Label.FontSize = 8;%c.Label.FontWeight = 'bold';
        c.FontSize = 8;
        if length(smoothCurveIndi) > 1
            c.Position = [0.7 .11 .02 .2];
            %
        else
            c.Position = [0.81 .11 .02 .2];
        end
    end
    caxis([-3 3])
    if j == 1 & length(smoothCurveIndi) > 1
        set(gca, 'XDir', 'reverse')
    end
    if j == length(smoothCurveIndi)
        if length(smoothCurveIndi) > 1
            subplot(1,10*length(smoothCurveIndi),[1+8*j-2 8*j+2-2])
        else
            subplot(1,10*length(smoothCurveIndi),[1+8*j 8*j+2])
        end
        
        ybar = [];
        for i = optimalK:-1:1
            ybar = [ybar,length(intersect(PDG_TFname,PDGIfo.PDGall(idxGeneOrder{i})))];
        end
        H = bar([ybar;nan(1,length(ybar))],'stacked','BarWidth',0.2,'ShowBaseLine','off','LineWidth',0.01);
        xlim([1 1.001])
        set(gca,'xtick',1)
        axis tight
        axis off
        k = 0;
        for i = optimalK:-1:1
            k = k+1;
            if ybar(i) ~= 0
                set(H(i),'facecolor',colorClusters(k,:))
            end
        end
    end
end
folderName = fullfile('results','figures');
if ~exist(folderName, 'dir')
    mkdir(folderName);
end
saveas(hFig,fullfile(folderName,'rolling_wave_critical_TF.pdf'))
