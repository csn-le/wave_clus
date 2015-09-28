function [templates, maxdist, pointdist] = build_templates(classes,features)
% function [templates maxdist] = build_templates(classes,inspk,pointdist)
max_class = max(classes);
feature_dim = size(features,2);
templates = zeros(max_class, feature_dim);
maxdist   = zeros(1,max_class);
pointdist = zeros(max_class,feature_dim);
for i=1:max_class,
    fi = features(classes==i,:);
    templates(i,:) = mean(fi);
    maxdist(i)     = sqrt(sum(var(fi,1)));   % the 1 means that we want sum(x-m)^2/N, not N-1
                                             % maxdist is the std dev of
                                             % the euclidean distance from
                                             % mean.
    pointdist(i,:)   = sqrt(var(fi,1));      % the std dev of the variation along each dimension.
    
end

