function y = fix_filter(x)
% fix_filter (in case signal processing toolbox is not available).
% filters data x with an eliptic passband between [300 3000] Hz.

a = [1.0000 -2.3930  2.0859 -0.9413 0.2502];
b = [0.1966 -0.0167 -0.3598 -0.0167 0.1966];

x = x(:);  
len = size(x,1);  
b = b(:).';
a = a(:).';
nb = length(b);
na = length(a);
nfilt = max(nb,na);

nfact = 3*(nfilt-1);  % length of edge transients

rows = [1:nfilt-1  2:nfilt-1  1:nfilt-2];
cols = [ones(1,nfilt-1) 2:nfilt-1  2:nfilt-1];
data = [1+a(2) a(3:nfilt) ones(1,nfilt-2)  -ones(1,nfilt-2)];
sp = sparse(rows,cols,data);
zi = sp \ ( b(2:nfilt).' - a(2:nfilt).'*b(1) );

y = [2*x(1)-x((nfact+1):-1:2);x;2*x(len)-x((len-1):-1:len-nfact)];

if exist('FilterM','file')
    y = FilterM(b,a,y,[zi*y(1)]);
else
    y = filter(b,a,y,[zi*y(1)]);
end
y = y(length(y):-1:1);

%second filter, in the other way
if exist('FilterM','file')
    y = FilterM(b,a,y,[zi*y(1)]);
else
    y = filter(b,a,y,[zi*y(1)]);
end
y = y(length(y):-1:1);

y([1:nfact len+nfact+(1:nfact)]) = [];

y = y.';   
