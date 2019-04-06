function [spikes1] = int_spikes_pol(spikes,nch,thr,par); 
%Interpolates with cubic splines to improve alignment.

w_pre = par.w_pre;
w_post = par.w_post;
detect = par.detection;
int_factor = par.int_factor;
nspk = size(spikes,1);
s = 1:size(spikes,2)/nch; % length of the single spike
ints = 1/int_factor:1/int_factor:size(spikes,2)/nch; 
lspk = w_pre + w_post; % length of a spike
awin = par.alignment_window;  
lsa = w_pre + w_post + 2*awin; % spike length with extra samples for alignment

spikes1 = zeros(nspk,lspk*nch);
for j=1:nch
    intspikes = spline(s,spikes(:,1+(j-1)*lsa:j*lsa),ints);
    switch detect
        case 'pos'
            [maxi iaux] = max(intspikes(:, w_pre*int_factor:(w_pre+2*awin)*int_factor ),[],2); 
        case 'neg'
            [maxi iaux] = min(intspikes(:, w_pre*int_factor:(w_pre+2*awin)*int_factor ),[],2); 
        case 'both'
            [maxi iaux] = max(abs(intspikes(:, w_pre*int_factor:(w_pre+2*awin)*int_factor )),[],2); 
    end
    
    iaux = iaux + w_pre*int_factor - 1;

    for si=1:nspk
        spikes1(si,(j-1)*lspk+1:j*lspk) = intspikes(si,iaux(si)-w_pre*int_factor+int_factor:int_factor:iaux(si)+w_post*int_factor);
    end
end
