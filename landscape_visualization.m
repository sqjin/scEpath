function landscape_visualization(scEcell,ydata,clusterIfo,colorCell,nlevels, viewpoint1,viewpoint2,fig_width, fig_height)
% display energy landscape in 2-D contour plot and 3-D surface
% Inputs:
%   scEcell : m x 1 vector, single cell energy
%   ydata : m x 2, 2-D coordindates from dimension reduction
%   clusterIfo : a struct giving the cell cluster information
%   colorCell : the color of each cluster, default is randomly generated
%   nLevels : contour levels in the contour plot, default=8
%   viewpoint: Viewpoint specification of contour plot
%   fig_width : the figure width
%   fig_height : the figure height
% Outputs:
%   two figures
if ~exist('colorCell', 'var') || isempty(colorCell)
    colorCell = distinguishable_colors(length(unique(clusterIfo.identity))+1);% colors for each cluster
    if length(unique(clusterIfo.identity)) >= 3
        colorCell(4,:) = []; % the fourth color is black
    end
end
if ~exist('nlevels', 'var') || isempty(nlevels)
    nlevels = 8;
end
if ~exist('viewpoint1', 'var') || isempty(viewpoint1)
    viewpoint1 = [90 90];
end
if ~exist('viewpoint2', 'var') || isempty(viewpoint2)
    viewpoint2 = [137 34];
end
if ~exist('fig_width', 'var') || isempty(fig_width)
    fig_width = 300;
end
if ~exist('fig_height', 'var') || isempty(fig_height)
    fig_height = 300;
end

group = clusterIfo.identity;idxCluster = clusterIfo.idxCluster;
numCluster = length(unique(group));
outliers = zeros(numCluster,2);
idxCenter = cell(numCluster,1);
idxDataCenter = [];
for i = 1:numCluster
    Q = quantile(scEcell(group == i),[0.25 0.75]);
    r = iqr(scEcell(group == i));
    k = 1.5;
    outliers(i,:) = [Q(1) - k*r, Q(2) + k*r];
    idx1 = find(scEcell(group == i) > outliers(i,1));
    idx2 = find(scEcell(group == i) < outliers(i,2));
    idxCenter{i} = intersect(idx1,idx2);idxCenter{i} = idxCluster{i}(idxCenter{i});
    idxDataCenter = [idxDataCenter;idxCenter{i}];
end
warning('off','all')
ydataCenter = ydata(idxDataCenter,:); scEcellCenter = scEcell(idxDataCenter);
bpoints = convhull(ydataCenter(:,1),ydataCenter(:,2));
ydatab = ydataCenter(bpoints,:);  scEcellb = scEcellCenter(bpoints)+0.1;
for i = 1:size(ydatab,1)
    for j = 1:size(ydatab,2)
        if ydatab(i,j) < 0
            ydatab(i,j) = ydatab(i,j)-range(ydataCenter(:))/5;
        else
            ydatab(i,j) = ydatab(i,j)+range(ydataCenter(:))/5;
        end
    end
end
ydatac = [ydataCenter;ydatab];scEcellc = [scEcellCenter,scEcellb];
sf = fit([ydatac(:,1), ydatac(:,2)],scEcellc','linearinterp');
x1range = [min(ydata(:,1)),max(ydata(:,1))]; x2range = [min(ydata(:,2)),max(ydata(:,2))];
x = linspace(x1range(1),x1range(2),500);
y = linspace(x2range(1),x2range(2),500);
[X,Y] = meshgrid(x,y);
Z = sf(X,Y);

hFig = figure('position', [600, 200, fig_width*1.3, fig_height]);
contourf(X,Y,Z,nlevels);
colormap(flipud(colormap(pink)))
c = colorbar;
c.Location = 'eastoutside';
c.Label.String = 'scEnergy';
c.Label.FontSize = 10;%c.Label.FontWeight = 'bold';
c.FontSize = 8;
hold on
% for i = 1:numCluster
%     scatter(ydata(group==i,1),ydata(group==i,2),8,colorCell(i,:),'filled');
% end
for i = 1:numCluster
    x1 = ydata(group == i,1);x2 = ydata(group == i,2);
    U = scEcell(group == i);
    for k = 1:length(x1)
        if U(k) < outliers(i,1) || U(k) > outliers(i,2)
            continue;
        end
        generate_three_dimension_ball(x1(k),x2(k),U(k),colorCell(i,:),x1range,x2range) %use another function to draw the 3D-balls
    end
    %     camlight right
    %     lighting phong
end
camlight right
lighting phong
set(gca,'Xtick',[]);set(gca,'Ytick',[])
xlabel('Component 1','FontSize',10,'FontName','Arial')
ylabel('Component 2','FontSize',10,'FontName','Arial')
title('Contour plot','FontSize',10,'FontName','Arial')
view(viewpoint1)

folderName = fullfile('results','figures');
if ~exist(folderName, 'dir')
    mkdir(folderName);
end
saveas(hFig,fullfile(folderName,'contour_plot_energy_landscape.pdf'))


hFig2 = figure('position', [600, 200, fig_width, fig_height]);
plot(sf)
colormap(flipud(colormap(gray)))
shading interp;
alpha(0.5) % this parameter is to make the landscape transparent
view(viewpoint2)
xlabel('Component 1','FontSize',10,'FontName','Arial')
ylabel('Component 2','FontSize',10,'FontName','Arial')
zlabel('scEnergy','FontSize',10,'FontName','Arial')
set(gca,'FontSize',8)

hold on
for i = 1:numCluster
    x1 = ydata(group == i,1);x2 = ydata(group == i,2);
    U = scEcell(group == i);
    for k = 1:length(x1)
        if U(k) < outliers(i,1) || U(k) > outliers(i,2)
            continue;
        end
        generate_three_dimension_ball(x1(k),x2(k),U(k),colorCell(i,:),x1range,x2range) %use another function to draw the 3D-balls
    end
    % camlight right
    % lighting phong
end
camlight right
lighting phong
axis([min(ydataCenter(:,1))-5,max(ydataCenter(:,1))+5, min(ydataCenter(:,2))-5, max(ydataCenter(:,2))+5, -0.05, 1])

folderName = fullfile('results','figures');
if ~exist(folderName, 'dir')
    mkdir(folderName);
end
% saveas(hFig2,fullfile(folderName,'3D_energy_landscape.pdf'))
print(fullfile(folderName,'3D_energy_landscape'),'-depsc','-tiff','-r300')


