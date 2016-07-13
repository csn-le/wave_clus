function [spikes1] = int_spikes(spikes,par)
%Interpolates with cubic splines to improve alignment.

w_pre = par.w_pre;
w_post = par.w_post;
ls = w_pre + w_post;
detect = par.detection;
int_factor = par.int_factor;
nspk = size(spikes,1);
extra = (size(spikes,2)-ls)/2;

s = 1:size(spikes,2);
ints = 1/int_factor:1/int_factor:size(spikes,2);

intspikes = zeros(1,length(ints));
spikes1 = zeros(nspk,ls);

switch detect
    case 'pos'
        for i=1:nspk
            intspikes(:) = spline(s,spikes(i,:),ints);
%             [maxi iaux] = max((intspikes(w_pre*int_factor:w_pre*int_factor+2*int_factor))); 
            [maxi iaux] = max(intspikes((w_pre+extra-1)*int_factor:(w_pre+extra+1)*int_factor)); 
%             iaux = iaux + w_pre*int_factor -1;
            iaux = iaux + (w_pre+extra-1)*int_factor -1;
            spikes1(i,1:ls)= intspikes(iaux-w_pre*int_factor+int_factor:int_factor:iaux+w_post*int_factor);
        end
    case 'neg'
        for i=1:nspk
            intspikes(:) = spline(s,spikes(i,:),ints);
            [maxi iaux] = min(intspikes((w_pre+extra-1)*int_factor:(w_pre+extra+1)*int_factor)); 
            iaux = iaux + (w_pre+extra-1)*int_factor -1;
            spikes1(i,1:ls)= intspikes(iaux-w_pre*int_factor+int_factor:int_factor:iaux+w_post*int_factor);
        end
    case 'both'
        for i=1:nspk
            intspikes(:) = spline(s,spikes(i,:),ints);
            [maxi iaux] = max(abs(intspikes((w_pre+extra-1)*int_factor:(w_pre+extra+1)*int_factor))); 
            iaux = iaux + (w_pre+extra-1)*int_factor -1;
            spikes1(i,1:ls)= intspikes(iaux-w_pre*int_factor+int_factor:int_factor:iaux+w_post*int_factor);
        end
end
