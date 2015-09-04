function [c,l] = fix_wavedec(x,n)
% Does a haar wavelet decomposition with n scales.
% Avoids using wavelet toolbox


Lo_D = [ 0.7071 0.7071];
Hi_D = [-0.7071 0.7071];

s = size(x); x = x(:)'; 
c = []; 
l = [length(x)];

dwtEXTM = 'sym';
shift = 0;

for k = 1:n
    lf = length(Lo_D);
    lx = length(x);
    lenEXT = lf-1; lenKEPT = lx+lf-1;
   
    I = getSymIndices(lx,1);
    y  = x(I);
    
    x = convdown(y,Lo_D,lenKEPT,shift);
    d = convdown(y,Hi_D,lenKEPT,shift);
    
    c     = [d c];            % store detail
    l     = [length(d) l];    % store length
end

% Last approximation.
c = [x c];
l = [length(x) l];

if s(1)>1, c = c'; l = l'; end


%-----------------------------------------------------%
% Internal Function(s)
%-----------------------------------------------------%
function y = convdown(x,f,lenKEPT,shift)

y = conv2(x(:)',f(:)'); if size(x,1)>1 , y = y'; end

sx = length(y);
begInd = 1;
[first,last,ok] = GetFirstLast(sx,begInd,lenKEPT);
if ok , y = y(first(1):last(1)); end

y = y(2-rem(shift,2):2:end);

%-----------------------------------------------------%
%----------------------------------------------------------------------------%
function I = getSymIndices(lx,lf)

I = [lf:-1:1 , 1:lx , lx:-1:lx-lf+1];
if lx<lf
    K = (I<1);
    I(K) = 1-I(K);
    J = (I>lx);
    while any(J)
        I(J) = 2*lx+1-I(J);
        K = (I<1);
        I(K) = 1-I(K);
        J = (I>lx);
    end
end
%----------------------------------------------------------------------------%
%----------------------------------------------------------------------------%
function [first,last,ok] = GetFirstLast(sx,begInd,varargin)

oneDIM = isequal(begInd,1);
s = varargin{1}(:)';
if ~oneDIM
    K  = find(s>sx);
    s(K) = sx(K);
    m = find((s < 0) | (s ~= fix(s)));
    ok = isempty(m);
else
    ok = (s>=0) & (s<sx) & (s == fix(s));
end
if ok==0 , first = begInd; last = s; return; end

nbarg = length(varargin);
if nbarg<2, o = 'c'; else , o = lower(varargin{2}); end

err = 0;
if ischar(o(1))
    switch o(1)
        case 'c'
            d = (sx-s)/2;
            if nbarg<3
                if length(o)>1 , side = o(2:end); else , side = 'l'; end
            else
                side = varargin{3};
            end
            if oneDIM
                [first,last] = GetFirst1D(side,sx,d);
            else
                if length(side)<2 , side(2) = 0; end
                for k = 1:2
                    [first(k),last(k)] = GetFirst1D(side(k),sx(k),d(k));
                end
            end

        case {'l','u'} , first = begInd; last = s;
        case {'r','d'} , first = sx-s+1; last = sx;
        otherwise      , err = 1;
    end
else
    first = o; last = first+s-1;
    if ~isequal(first,fix(first)) | any(first<1) | any(last>sx)
        err = 1;
    end
end
if err
    errargt(mfilename,'invalid argument','msg');
    error('*');
end
%----------------------------------------------------------------------------%
function [first,last] = GetFirst1D(side,s,d)

switch side
  case {'u','l','0',0} , first = 1+floor(d); last = s-ceil(d);
  case {'d','r','1',1} , first = 1+ceil(d);  last = s-floor(d);
  otherwise    , first = 1+floor(d); last = s-ceil(d);  % Default is left side
end
%----------------------------------------------------------------------------%

