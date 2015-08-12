function index = nearest_mahal(x,mu,sigma)
% function index = nearest_mahal(x,mu,sigma)
% x is a vector
% mu(i,:) is the mean of the ith Gaussian
% sigma(:,:,i) is the covariance of the ith Gaussian
% 
% Returns the index of the Gaussian closest (by the Mahalanobis distance)
% to x.

N = size(mu,1);  % number of Gaussians
d = [];
if( N == 0 )
    index = 0;
else
    for i=1:N,
        d(i) = (x-mu(i,:))*inv(sigma(:,:,i))*(x-mu(i,:))';
    end
    [m index] = min(d);
end
