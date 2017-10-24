function cluster_visualization(ydata, group,class_labels, true_labs, marker_size,colorCell,fig_width,fig_height)
% display cell lineage hierarchy
% Inputs:
%   ydata : m x 2, 2-D coordindates from dimension reduction
%   group : m x 1 numerical vector, the cluster assignment for each cell
%   class_labels : N x 1 cell strings, text annotations of each cluster in the figure legend
%   true_labs : N x 1 numerical vector or cell strings or categorical variables, text annotations (prior knowledge of cells, e.g. time information) of each cluster in the figure legend
%   marker_size : the size of individual cell, default= 8
%   colorCell : the color of each cluster, default is randomly generated
%   fig_width : the figure width, default=600
%   fig_height : the figure height, default=250
% Outputs:
%   one figure
if ~exist('marker_size', 'var') || isempty(marker_size)
    marker_size = 8;
end
if ~exist('colorCell', 'var') || isempty(colorCell)
    colorCell = distinguishable_colors(length(unique(group))+1);% colors for each cluster
    if length(unique(group)) >= 3
        colorCell(4,:) = []; % the fourth color is black
    end
end
if exist('true_labs', 'var') & ~isempty(true_labs)
    colorCell2 = distinguishable_colors(length(unique(true_labs))+1);% colors for each cluster
    if length(unique(true_labs)) >= 3
        colorCell2(4,:) = []; % the fourth color is black
    end
end
if ~exist('fig_width', 'var') || isempty(fig_width)
    fig_width = 600;
end
if ~exist('fig_height', 'var') || isempty(fig_height)
    fig_height = 250;
end

figure('position', [600, 200, fig_width, fig_height])

if ~exist('true_labs', 'var') || isempty(true_labs)
    numCluster = length(unique(group));
    h = zeros(1,numCluster);
    for i = 1:length(unique(group))
        hold on
        h(i) = scatter(ydata(group==i,1),ydata(group==i,2),marker_size(group==i),colorCell(i,:),'filled');
    end
    legend(h,class_labels,'Location','best')
    set(gca,'Xtick',[]);set(gca,'Ytick',[])
    axis tight;box off
    set(gcf,'color','w');
    xlabel('Component 1','FontName','Arial','FontSize',10);
    ylabel('Component 2','FontName','Arial','FontSize',10);
else
    subplot(1,2,1); gscatter(ydata(:,1),ydata(:,2),true_labs,colorCell2);
    alpha(0.5)
    set(gca,'Xtick',[]);set(gca,'Ytick',[])
    axis tight;box off
    set(gcf,'color','w');
    xlabel('Component 1','FontName','Arial','FontSize',10);
    ylabel('Component 2','FontName','Arial','FontSize',10);
    subplot(1,2,2);
    numCluster = length(unique(group));
    class_labels = cell(1,numCluster);h = zeros(1,numCluster);
    for i = 1:length(unique(group))
        hold on
        h(i) = scatter(ydata(group==i,1),ydata(group==i,2),marker_size(group==i),colorCell(i,:),'filled');
        class_labels{i} = ['C' num2str(i)];
    end
    legend(h,class_labels,'Location','best')
    set(gca,'Xtick',[]);set(gca,'Ytick',[])
    axis tight;box off
    set(gcf,'color','w');
    xlabel('Component 1','FontName','Arial','FontSize',10);
    ylabel('Component 2','FontName','Arial','FontSize',10);
end
folderName = fullfile('results','figures');
if ~exist(folderName, 'dir')
    mkdir(folderName);
end
saveas(gcf,fullfile(folderName,'cluster_2Dvisualization.pdf'))