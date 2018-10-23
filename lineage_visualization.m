function lineage_visualization(lineageIfo,class_labels,node_size,colorCell,showLoops,fig_width)
% display cell lineage hierarchy
% Inputs: 
%   lineageIfo : a struct giving the cell lineage information (returned from the function inferingLineage.m)
%   class_labels : N x 1 cell strings, text annotations of each cluster in the figure legend
%   node_size : the size of tree node, default= 12
%   colorCell : the color of each tree node, default=[0 0 0]
%   showLoops : boolean, to show the self-loops of lineage tree, default= 1 (true)
%   fig_width : the figure width, default=150
% Outputs:
%   one figure
if ~exist('node_size', 'var') || isempty(node_size)
    node_size = 12;
end
if ~exist('colorCell', 'var') || isempty(colorCell)
    colorCell = [0 0 0];
end
if ~exist('fig_width', 'var') || isempty(fig_width)
    fig_width = 150;
end
if ~exist('showLoops', 'var') || isempty(showLoops)
    showLoops = 1;
end

figure('position', [600, 0, fig_width, fig_width*2.5])

MDST = lineageIfo.MDST;
if ~showLoops
    MDST = rmedge(MDST, 1:numnodes(MDST), 1:numnodes(MDST));   
end

h = plot(MDST,'EdgeAlpha',0.8,'EdgeColor','k','LineWidth',1,'NodeColor',colorCell);
layout(h,'layered','Sources',lineageIfo.rootNode)
set(h,'MarkerSize',node_size)
set(h,'EdgeLabel',round(MDST.Edges.Weight*10^2)/10^2)
set(h,'NodeLabel',[])
set(h,'ArrowSize',8)
xd = get(h, 'XData');
yd = get(h, 'YData');
text(xd-0.4, yd, class_labels, 'FontSize',10, 'FontWeight','bold','Interpreter','none', 'HorizontalAlignment','left', 'VerticalAlignment','middle')
axis off
set(gcf,'color','w');
folderName = fullfile('results','figures');
if ~exist(folderName, 'dir')
    mkdir(folderName);
end
saveas(gcf,fullfile(folderName,'inferred_lineage.pdf'))
