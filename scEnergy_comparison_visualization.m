function scEnergy_comparison_visualization(scEcell,clusterIfo,class_labels,colorCell,show_pvalue)
% display cell lineage hierarchy
% Inputs:
%   scEcell : m x 1 vector, single cell energy
%   clusterIfo : a struct giving the cell cluster information
%   class_labels : N x 1 cell strings, text annotations of each cluster in the figure legend
%   colorCell : the color of each box, default= auto
%   show_pvalue : boolean, to show the legend or not, default= 1 (true)
% Outputs:
%   one figure
if ~exist('show_pvalue', 'var') || isempty(show_pvalue)
    show_pvalue = 1;
end
group = clusterIfo.identity;
numCluster = length(unique(group));
idxCluster = clusterIfo.idxCluster;

figure
if ~exist('colorCell', 'var') || isempty(colorCell)
    boxplot(scEcell,group);
else
    boxplot(scEcell,group, 'colors',colorCell);
end
set(gca,'xtickLabel',class_labels)
set(gca,'FontName','Arial','FontSize',10)
ylabel('scEnergy','FontSize',10,'FontName','Arial');
xlabel('Clusters','FontSize',10,'FontName','Arial');
ylim([min(scEcell)-0.05,max(scEcell)+0.05])
if show_pvalue
    for i = 1:numCluster-1
        p = ranksum(scEcell(idxCluster{i}),scEcell(idxCluster{i+1}));
        if p < 0.01
            text(i+0.35,max(scEcell(idxCluster{i+1}))+0.02,['P=' num2str(p,'%3.0e')],'FontSize',10,'FontWeight','bold')
        else
            text(i+0.35,max(scEcell(idxCluster{i+1}))+0.02,['P=' num2str(p,'%.2f')],'FontSize',10,'FontWeight','bold')
        end
    end
end
folderName = fullfile('results','figures');
if ~exist(folderName, 'dir')
    mkdir(folderName);
end
saveas(gcf,fullfile(folderName,'scEnergy_comparison.pdf'))