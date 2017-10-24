function [R,IDselect,deg,k,pk,signedR2] = calculatingDegreeDistribution(R,thresh)
R = R>thresh;
G = graph(R,'upper');
comp = conncomp(G,'OutputForm','cell');
[~,idx] = max(cellfun('length',comp));
IDselect = comp{idx};
R = R(IDselect,IDselect);
G = graph(R,'upper');
deg = degree(G);
if length(deg) > 100
    nbins = 10;
else
    nbins = 5;
end
[pk,k] = hist(deg,nbins);
idx = pk==0; pk(idx) = [];k(idx) = [];
pk = pk/length(deg);
r = corrcoef(log10(k(:)),log10(pk(:)));r = r(1,2);
signedR2 = -sign(r)*r^2;