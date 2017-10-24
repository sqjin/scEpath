function plot_genes_in_pseudotime(my_genes,proData,smoothExprIfo,pseudotimeIfo,lineageIfo,color_by,colorCell,marker_size,fig_width, fig_height)
% plot individual gene temporal dynamics along pseudotime
% Inputs: 
%   my_genes : a cell array, the gene names
%   proData : the single cell data including genes and cell attributes
%   smoothExprIfo : the smooth expression information
%   pseudotimeIfo : the pseudotime information
%   lineageIfo : the cell lineage information
%   color_by: 1 x m array, cell attributes such as inferred cluster, culture condition, day captured 
%   colorCell : the color of each cluster, default=[0 0 0]
%   marker_size: the size of dots, default = 6
%   fig_width : the figure width
%   fig_height : the figure height
% Outputs:
%   one figure
if ~exist('marker_size', 'var') || isempty(marker_size)
    marker_size = 6;
end
if ~exist('fig_width', 'var') || isempty(fig_width)
    fig_width = 150*length(smoothExprIfo);
end
if ~exist('fig_height', 'var') || isempty(fig_height)
    fig_height = 100*length(my_genes);
end

[my_genes,~,idxmy_genes] = intersect(my_genes,proData.genes,'stable');

figure('position', [600, 0, fig_width, fig_height])

for j = 1:length(smoothExprIfo)
    pseudotime = pseudotimeIfo.pseudotime{j};
    cellOrder = pseudotimeIfo.cellOrder{j};
    smoothCurveTemp = smoothExprIfo{j}.smoothExpr(idxmy_genes,:);
    dataMarkers = proData.data(idxmy_genes,cellOrder);
    group = color_by;
    xii = linspace(min(pseudotime),max(pseudotime),length(pseudotime));
    
    for i = 1:length(my_genes)
        subplot(length(my_genes),length(smoothExprIfo),length(smoothExprIfo)*(i-1)+j)
        gscatter(pseudotime, dataMarkers(i,:),group(cellOrder)',colorCell(lineageIfo.path{j},:),[],marker_size);
        hold on
        plot(xii,smoothCurveTemp(i,:),'k-','Linewidth',1)
        set(gca,'Xtick',[]);
        set(gca,'xlabel',[]);
        set(gca,'FontName','Arial','FontSize',8)
        if i == length(my_genes)
            xlabel('Pseudotime','FontName','Arial','FontSize',10)
            set(gca,'Xtick',0:0.5:1);
        end
        if i == 1
            title(['Path ',num2str(j)],'FontName','Arial','FontSize',10)
        end
        if j == 1
            ylabel(my_genes{i},'FontName','Arial','FontSize',10)
            % title(markers{i},'FontName','Arial','FontSize',10,'FontWeight','bold')
        end
        axis tight
        legend off
        yl = ylim(gca);
        yl(1) = min([0,yl(1)]);
        ylim(gca, yl);
    end
end
folderName = fullfile('results','figures');
if ~exist(folderName, 'dir')
    mkdir(folderName);
end
saveas(gcf,fullfile(folderName,'terporal_dynamics_individual_genes.pdf'))