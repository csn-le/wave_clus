function [mu, inv_sigma] = fit_gaussian(x,class)

N = max(class);
mu = zeros(N,size(x,2));
inv_sigma = zeros(size(x,2),size(x,2),N);

for i=1:N,
    mu(i,:) = mean(x(class==i,:));
    inv_sigma(:,:,i) = inv(cov(x(class==i,:)));
end


    

