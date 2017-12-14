function [index,xf,thr] = index_detect_pol(x,par)
% return the spike times within a channel

%PARAMETERS
sr=par.sr;
w_pre=par.w_pre;
w_post=par.w_post;
ref = ceil(par.ref_ms/1000 * par.sr);
detect =  par.detection;
stdmin =  par.stdmin;
stdmax =  par.stdmax;
fmin_detect =  par.detect_fmin;
fmax_detect =  par.detect_fmax;
fmin_sort =  par.sort_fmin;
fmax_sort =  par.sort_fmax;
awin =  par.alignment_window;  

% HIGH-PASS FILTER OF THE DATA
if exist('ellip')                         %Checks for the signal processing toolbox
    [b_detect,a_detect] = ellip(par.detect_order,0.1,40,[fmin_detect fmax_detect]*2/sr);
    [b,a] = ellip(par.sort_order,0.1,40,[fmin_sort fmax_sort]*2/sr);
    if exist('FiltFiltM','file')
        xf_detect = FiltFiltM(b_detect,a_detect,x);
        xf = FiltFiltM(b,a,x);
    else
        xf_detect = filtfilt(b_detect,a_detect,x);
        xf = filtfilt(b,a,x);
    end
else
    xf=fix_filter(x);                   %Does a bandpass filtering between [300 3000] without the toolbox.
    xf_detect = xf;
end
clear x;

noise_std_detect = median(abs(xf_detect))/0.6745;
noise_std_sorted = median(abs(xf))/0.6745;
thr = stdmin * noise_std_detect;        %thr for detection is based on detected settings.
thrmax = stdmax * noise_std_sorted;     %thrmax for artifact removal is based on sorted settings.
index = [];

% LOCATE SPIKE TIMES
switch detect
    case 'pos'
        nspk = 0;
        xaux = find(xf_detect(w_pre+awin+2:end-w_post-awin-2) > thr) + w_pre+awin+1;
        xaux0 = 0;
        for i=1:length(xaux)
            if xaux(i) >= xaux0 + ref
                [maxi, iaux]=max((xf(xaux(i):xaux(i)+floor(ref/2)-1)));    %introduces alignment
                nspk = nspk + 1;
                index(nspk) = iaux + xaux(i) - 1; % we selec the max value among the values over the threshold
                % minimum value for index is "w_pre+2"
                xaux0 = index(nspk);
            end
        end
    case 'neg'
        nspk = 0;
        xaux = find(xf_detect(w_pre+awin+2:end-w_post-awin-2) < -thr) + w_pre+awin+1;
        xaux0 = 0;
        for i=1:length(xaux)
            if xaux(i) >= xaux0 + ref
                [maxi, iaux]=min((xf(xaux(i):xaux(i)+floor(ref/2)-1)));    %introduces alignment
                nspk = nspk + 1;
                index(nspk) = iaux + xaux(i) -1;
                xaux0 = index(nspk);
            end
        end
    case 'both'
        nspk = 0;
        xaux = find(abs(xf_detect(w_pre+awin+2:end-w_post-awin-2)) > thr) + w_pre+awin+1;
        xaux0 = 0;
        for i=1:length(xaux)
            if xaux(i) >= xaux0 + ref
                [maxi, iaux]=max(abs(xf(xaux(i):xaux(i)+floor(ref/2)-1)));    %introduces alignment
                nspk = nspk + 1;
                index(nspk) = iaux + xaux(i) -1;
                xaux0 = index(nspk);
            end
        end

end

xf(length(xf)+1:length(xf)+w_post+awin)=0; %add 44 more points just in case the last point is a spike and we have to align it
aux=[];
for i=1:nspk                          %Eliminates artifacts
    if max(abs( xf(index(i)-w_pre:index(i)+w_post) )) > thrmax
        aux(i)=1;
    end
end
% aux = find(spikes(:,w_pre)==0);       %erases indexes that were artifacts, every spike over thrmax is erased
index(aux==1)=[];