function index = nearest_neighbor(x,vectors,maxdist,varargin)
% function index = nearest_neigbor(x,vectors,maxdist,pointdist*,pointlimit*,k*)
% x is a row vector
% pointdist (optional) - vector of standard deviations
% pointlimit (optional) - upper bound on number of points outside pointdist
% k (optional) - number of points used for nearest neighbor

% Find the distance to all neighbors. Consider only those neighbors where
% the point falls in the radius of possibility for that point. Find the
% nearest possible neighbor.
% Return 0 if there is no possible nearest neighbor.

distances = sqrt(sum((ones(size(vectors,1),1)*x - vectors).^2,2)');
conforming = find(distances < maxdist);
if( length(varargin) > 0 )
    pointdist = varargin{1};
    if( length(varargin) > 1 )
        pointlimit = varargin{2};
    else
        pointlimit = Inf;
    end
    pointwise_conforming = [];
    for i=1:size(vectors,1),
        if( sum( abs(x-vectors(i,:)) > pointdist(i,:) ) < pointlimit )  % number of deviations from pointdist allowed.
            pointwise_conforming = [pointwise_conforming i];
        end
    end
    conforming = intersect(conforming, pointwise_conforming);
end
if( length( conforming ) == 0 )
    index = 0;
else
    if( length(varargin) > 2 )
        k = varargin{3};
        [y i] = sort(distances(conforming)); % k-nearest neighbors
        i = i(1:min(length(i),k));
    else
        [y i] = min(distances(conforming));   
    end
    index = conforming(i);
end
