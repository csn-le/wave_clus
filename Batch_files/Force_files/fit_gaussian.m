function [mu, sigma] = fit_gaussian(x,class)

N = max(class);
mu = zeros(N,size(x,2));
sigma = zeros(size(x,2),size(x,2),N);

for i=1:N,
    mu(i,:) = mean(x(class==i,:));
    sigma(:,:,i) = cov(x(class==i,:));
end


    

