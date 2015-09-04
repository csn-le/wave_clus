function [KSmax] = test_ks(x)
% 
% Calculates the CDF (expcdf)
%[y_expcdf,x_expcdf]=cdfcalc(x);

yCDF = [];
xCDF = [];
x = x(~isnan(x));
n = length(x);
x = sort(x(:));
% Get cumulative sums
yCDF = (1:n)' / n;
% Remove duplicates; only need final one with total count
notdup = ([diff(x(:)); 1] > 0);
x_expcdf = x(notdup);
y_expcdf = [0; yCDF(notdup)];

%
% The theoretical CDF (theocdf) is assumed to be normal  
% with unknown mean and sigma

zScores  =  (x_expcdf - mean(x))./std(x);

%theocdf  =  normcdf(zScores , 0 , 1);
mu = 0; 
sigma = 1; 
theocdf = 0.5 * erfc(-(zScores-mu)./(sqrt(2)*sigma));


%
% Compute the Maximum distance: max|S(x) - theocdf(x)|.
%

delta1    =  y_expcdf(1:end-1) - theocdf;   % Vertical difference at jumps approaching from the LEFT.
delta2    =  y_expcdf(2:end)   - theocdf;   % Vertical difference at jumps approaching from the RIGHT.
deltacdf  =  abs([delta1 ; delta2]);

KSmax =  max(deltacdf);
