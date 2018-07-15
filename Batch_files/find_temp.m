function [clust_num temp auto_sort] = find_temp(tree,clu,par) %,spikes,ipermut)


min_clus = par.min_clus;

c_ov = par.c_ov;
elbow_min = par.elbow_min;
clu = clu(1:end-1,3:end)+1; %first dim temp

y = tree(1:end,5);
maxdiff = max(diff(tree(1:end,6:end)),[],2);
maxdiff(maxdiff<0)=0;
prop = (y(2:end)+maxdiff.*(maxdiff>0))./y(1:end-1);
aux = find(prop<elbow_min,1,'first')+1; %percentaje of the rest

tree = tree(1:end-1,5:end);
clus = zeros(size(tree));
clus(tree(:,:) >= min_clus)=1; %only check the ones that cross the thr

dt = diff(tree);
clus = clus & [ones(size(clus(1,:)));dt(1:end,:)>min_clus];

for ii = 1:size(clus,1)
    detect = find(clus(ii,:),1,'last');
    if ~isempty(detect)
        clus(ii,1:detect)=1;
    end
end

auto_sort.elbow = size(tree,1);
if ~isempty(aux)
    clus(aux:end,1:end)=0;
    auto_sort.elbow = aux;
end
auto_sort.peaks = clus;

for ti = size(clus,1):-1:1
    detect = find(clus(ti,:));
    for ci = 1:length(detect) %the clusters removed aren't detected
        cl = (clu(ti,:) == detect(ci));
        for tj = ti-1:-1:1
            toremove = find(clus(tj,:));
            for j = 1:length(toremove)
                totest = (clu(tj,:) == toremove(j));
                if nnz(cl & totest)/min(nnz(totest),nnz(cl)) >= c_ov
                    clus(tj,toremove(j))=0;
                end
            end
        end
    end
end

[temp clust_num]=find(clus);
