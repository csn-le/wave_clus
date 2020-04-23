function [spikes,thr,index] = amp_detect(x, par)
% Detect spikes with amplitude thresholding. Uses median estimation.
% Detection is done with filters set by fmin_detect and fmax_detect. Spikes
% are stored for sorting using fmin_sort and fmax_sort. This trick can
% eliminate noise in the detection but keeps the spikes shapes for sorting.


sr = par.sr;
w_pre = par.w_pre;
w_post = par.w_post;


if isfield(par,'ref_ms')
    ref = floor(par.ref_ms * par.sr/1000);
else
    ref = par.ref; %for retrocompatibility
end

detect = par.detection;
stdmin = par.stdmin;
stdmax = par.stdmax;

if par.sort_order > 0
    xf = filt_signal(x,par.sort_order,par.sort_fmin,par.sort_fmax,par.sr);
else
    xf = x;
end
if par.detect_order > 0
    xf_detect = filt_signal(x,par.detect_order,par.detect_fmin,par.detect_fmax,par.sr);
else
    xf_detect = x;
end

noise_std_detect = median(abs(xf_detect))/0.6745;
noise_std_sorted = median(abs(xf))/0.6745;
thr = stdmin * noise_std_detect;        %thr for detection is based on detect settings.
thrmax = stdmax * noise_std_sorted;     %thrmax for artifact removal is based on sorted settings.

index = [];
sample_ref = floor(ref/2);
% LOCATE SPIKE TIMES
switch detect
    case 'pos'
        nspk = 0;
        xaux = find(xf_detect(w_pre+2:end-w_post-2-sample_ref) > thr) +w_pre+1;
        xaux0 = 0;
        for i=1:length(xaux)
            if xaux(i) >= xaux0 + ref
                [aux_unused, iaux] = max((xf(xaux(i):xaux(i)+sample_ref-1)));    %introduces alignment
                nspk = nspk + 1;
                index(nspk) = iaux + xaux(i) -1;
                xaux0 = index(nspk);
            end
        end
    case 'neg'
        nspk = 0;
        xaux = find(xf_detect(w_pre+2:end-w_post-2-sample_ref) < -thr) +w_pre+1;
        xaux0 = 0;
        for i=1:length(xaux)
            if xaux(i) >= xaux0 + ref
                [aux_unused, iaux] = min((xf(xaux(i):xaux(i)+sample_ref-1)));    %introduces alignment
                nspk = nspk + 1;
                index(nspk) = iaux + xaux(i) -1;
                xaux0 = index(nspk);
            end
        end
    case 'both'
        nspk = 0;
        xaux = find(abs(xf_detect(w_pre+2:end-w_post-2-sample_ref)) > thr) +w_pre+1;
        xaux0 = 0;
        for i=1:length(xaux)
            if xaux(i) >= xaux0 + ref
                [aux_unused, iaux] = max(abs(xf(xaux(i):xaux(i)+sample_ref-1)));    %introduces alignment
                nspk = nspk + 1;
                index(nspk) = iaux + xaux(i) -1;
                xaux0 = index(nspk);
            end
        end
end

% SPIKE STORING (with or without interpolation)
ls = w_pre+w_post;
spikes = zeros(nspk,ls+4);

xf(length(xf)+1:length(xf)+w_post)=0;

for i=1:nspk                          %Eliminates artifacts
    if max(abs( xf(index(i)-w_pre:index(i)+w_post) )) < thrmax
        spikes(i,:)=xf(index(i)-w_pre-1:index(i)+w_post+2);
    end
end
aux = find(spikes(:,w_pre)==0);       %erases indexes that were artifacts
spikes(aux,:)=[];
index(aux)=[];

switch par.interpolation
    case 'n'
        spikes(:,end-1:end)=[];       %eliminates borders that were introduced for interpolation
        spikes(:,1:2)=[];
    case 'y'
        %Does interpolation
        spikes = int_spikes(spikes,par);
end
end

function filtered = filt_signal(x,order,fmin,fmax,sr)
    %HIGH-PASS FILTER OF THE DATA
    if exist('ellip','file')                         %Checks for the signal processing toolbox
        [b,a] = ellip(order,0.1,40,[fmin fmax]*2/sr);
        if exist('FiltFiltM','file')
            filtered = FiltFiltM(b, a, x); 
        else
            filtered = filtfilt(b, a, x);      
        end
    else
        filtered = fix_filter(x);    %Does a bandpass filtering between [300 3000] without the toolbox.
    end
end
