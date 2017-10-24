function networkIfo = constructingSFNetwork(data,quick_construct,thresh,showFigure,fig_width,fig_height)
% construct a gene-gene co-expression network
% Inputs:
%    data: single cell data (rows are cells and columns are genes)
%    quick_construct: default=0,the network will be constructed based on scale-free creteria of biological networks;
%    thresh: if quick_construct=1, thresh is the threshold for constructing a network. default = 0.1
%    showFigure: boolean, to show the network metric (e.g.signed R^2, average connectivity and degree discontribution) or not, default= 1 (true)
%    fig_width : the figure width, default=800
%    fig_height : the figure height, default=250
% Outputs:
%    networkIfo: network information for the constructed gene-gene network
%    networkIfo.R: adjacency matrix (upper matrix) of the constructed network
%    networkIfo.IDselect: the index of selected genes in the constructed network
%    networkIfo.deg: the connectivity (degree) of each node
%    networkIfo.signedR2; the signed R^2 of the constructed network
if ~exist('showFigure','var') || isempty(showFigure)
    showFigure = 1;
end
if ~exist('fig_width', 'var') || isempty(fig_width)
    fig_width = 600;
end
if ~exist('fig_height', 'var') || isempty(fig_height)
    fig_height = 180;
end

R00 = corr(data,'Type','Spearman');
R00 = triu(R00,1);
R00 = sparse(R00);
R00 = abs(R00);

if quick_construct || size(R00,1) < 50
    if ~exist('thresh','var') || isempty(thresh),thresh = 0.1; end
    threshSig = thresh;
    [R,IDselect,deg,k,pk,signedR2] = calculatingDegreeDistribution(R00,threshSig);
    signedR2T = signedR2;avDeg = mean(deg);
else
    [thresh,signedR2T,avDeg] = determineThresh(R00);
    idx = find(signedR2T > 0.8,1);
    threshSig = thresh(idx);
    [R,IDselect,deg,k,pk,signedR2] = calculatingDegreeDistribution(R00,threshSig);
end
networkIfo.R = R; networkIfo.IDselect = IDselect; networkIfo.deg = deg; networkIfo.signedR2 = signedR2;

if showFigure
    hFig = figure('position', [600, 200, fig_width, fig_height]);
    subplot(1,3,1)
    scatter(thresh,signedR2T,'k')
    hold on
    if length(signedR2T) > 1
        line([0.1 0.9],[signedR2T(idx) signedR2T(idx)],'Color','r')
    end
    ylim([-1 1])
    xlim([0.1 0.9])
    xlabel('\tau','FontName','Arial','FontSize',10);
    ylabel('Signed R^2','FontName','Arial','FontSize',10);
    box on
    
    subplot(1,3,2)
    scatter(thresh,avDeg,'k')
    xlabel('\tau','FontName','Arial','FontSize',10);
    ylabel('Average connectivity','FontName','Arial','FontSize',10);
    xlim([0.1 0.9])
    box on
    
    mdl = fitlm(log10(k(:)),log10(pk(:)));
    beta = mdl.Coefficients.Estimate;
    subplot(1,3,3)
    scatter(log10(k(:)),log10(pk(:)),'k')
    hold on;
    x = [min(log10(k)):0.01:max(log10(k))];
    y = beta(1) + beta(2)*x;
    plot(x,y,'r');
    hold off;
    title(['\tau = ' num2str(threshSig) ', R^2 = ' num2str(signedR2,2)])
    xlabel('log10(k)','FontName','Arial','FontSize',10);
    ylabel('log10(p(k))','FontName','Arial','FontSize',10);
    xlim([min(x)-0.02,max(x)+0.02])
    box on
    
    folderName = fullfile('results','figures');
    if ~exist(folderName, 'dir')
        mkdir(folderName);
    end
    saveas(hFig,fullfile(folderName,'topological_metrics_gene_network.pdf'))
    
end



function [thresh,signedR2,avDeg] = determineThresh(R00)
thresh = [0.1:0.1:0.5 0.55:0.05:0.85];
signedR2 = zeros(1,length(thresh)); avDeg = signedR2;
for i = 1:length(thresh)
    try
        [~,~,deg,~,~,signedR2(i)] = calculatingDegreeDistribution(R00,thresh(i));
        avDeg(i) = mean(deg);
    catch
        break;
    end
    
end


