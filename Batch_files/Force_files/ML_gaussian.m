function index = ML_gaussian(x,mu,inv_sigma)
% function index = ML_gaussian(x,mu,sigma)
% x is a vector drawn from some multivariate gaussian
% mu(i,:) is the mean of the ith Gaussian
% sigma(:,:,i) is the covariance of the ith Gaussian
% 
% Returns the index of the Gaussian with the highest value of p(x).

N = size(mu,1);  % number of Gaussians

if( N == 0 )
    index = 0;
else
    for i=1:N,
        % leave out factor of 1/(2*pi)^(N/2) since it doesn't affect argmax
        p(i) = sqrt(det(inv_sigma(:,:,i)))*exp(-0.5*(x-mu(i,:))*inv_sigma(:,:,i)*(x-mu(i,:))');
    end
    [m index] = max(p);
end
