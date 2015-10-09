function [spikes1] = int_spikes_pol(spikes,nch,thr,par); 
%Interpolates with cubic splines to improve alignment.

w_pre = par.w_pre;
w_post = par.w_post;
ref = ceil(par.ref_ms/1000 * par.sr);
detect = par.detection;
int_factor = par.int_factor;
nspk = size(spikes,1);
s = 1:size(spikes,2)/nch; % length of the single spike
ints = 1/int_factor:1/int_factor:size(spikes,2)/nch; 
intspikes = zeros(1,length(ints));
lspk = w_pre + w_post; % length of a spike
spikes1 = zeros(nspk,lspk*nch);
awin = par.alignment_window;  
lsa = w_pre + w_post + 2*awin; % spike length with extra samples for alignment

switch detect
    case 'pos'
        for i=1:nspk
            for j=1:nch
                intspikes(:) = spline(s,spikes(i,1+(j-1)*lsa:j*lsa),ints); %makes the interpolation and intspikes have double the points in this case
                [maxi iaux] = max( intspikes( w_pre*int_factor:(w_pre+2*awin)*int_factor ) ); 
                iaux = iaux + w_pre*int_factor - 1;
                spikes1(i,(j-1)*lspk+w_pre:-1:(j-1)*lspk+1) = intspikes(iaux:-int_factor:iaux-w_pre*int_factor+int_factor); %20 samples 1 of 2 (diezmado)
                spikes1(i,(j-1)*lspk+w_pre+1:j*lspk) = intspikes(iaux+int_factor:int_factor:iaux+w_post*int_factor); %44 samples 1 of 2 (diezmado)
            end
        end

    case 'neg'
        for i=1:nspk
            for j=1:nch
                intspikes(:) = spline(s,spikes(i,1+(j-1)*lsa:j*lsa),ints); %makes the interpolation and intspikes have double the points in this case
                [maxi iaux] = min( intspikes( w_pre*int_factor:(w_pre+2*awin)*int_factor ) ); 
                iaux = iaux + w_pre*int_factor - 1;
                spikes1(i,(j-1)*lspk+w_pre:-1:(j-1)*lspk+1) = intspikes(iaux:-int_factor:iaux-w_pre*int_factor+int_factor); %20 samples 1 of 2 (diezmado)
                spikes1(i,(j-1)*lspk+w_pre+1:j*lspk) = intspikes(iaux+int_factor:int_factor:iaux+w_post*int_factor); %44 samples 1 of 2 (diezmado)
            end
        end

    case 'both'
        for i=1:nspk
            for j=1:nch
                intspikes(:) = spline(s,spikes(i,1+(j-1)*lsa:j*lsa),ints); %makes the interpolation and intspikes have double the points in this case
                [maxi iaux] = max( abs(intspikes( w_pre*int_factor:(w_pre+2*awin)*int_factor )) ); 
                iaux = iaux + w_pre*int_factor - 1;
                spikes1(i,(j-1)*lspk+w_pre:-1:(j-1)*lspk+1) = intspikes(iaux:-int_factor:iaux-w_pre*int_factor+int_factor); %20 samples 1 of 2 (diezmado)
                spikes1(i,(j-1)*lspk+w_pre+1:j*lspk) = intspikes(iaux+int_factor:int_factor:iaux+w_post*int_factor); %44 samples 1 of 2 (diezmado)
            end
        end
end
