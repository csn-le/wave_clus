function [spikes,thr,index] = amp_detect(x, par)
% Detect spikes with amplitude thresholding. Uses median estimation.
% Detection is done with filters set by fmin_detect and fmax_detect. Spikes
% are stored for sorting using fmin_sort and fmax_sort. This trick can
% eliminate noise in the detection but keeps the spikes shapes for sorting.


sr = par.sr;
w_pre = par.w_pre;
w_post = par.w_post;
ref = par.ref;
detect = par.detection;
stdmin = par.stdmin;
stdmax = par.stdmax;
fmin_detect = par.detect_fmin;
fmax_detect = par.detect_fmax;
fmin_sort = par.sort_fmin;
fmax_sort = par.sort_fmax;


%HIGH-PASS FILTER OF THE DATA
if exist('ellip','file')                         %Checks for the signal processing toolbox
    [b_detect,a_detect] = ellip(par.detect_order,0.1,40,[fmin_detect fmax_detect]*2/sr);
    [b,a] = ellip(par.sort_order,0.1,40,[fmin_sort fmax_sort]*2/sr);
%     [z_det,p_det,k_det] = ellip(par.detect_order,0.1,40,[fmin_detect fmax_detect]*2/sr);
%     [z,p,k] = ellip(par.sort_order,0.1,40,[fmin_sort fmax_sort]*2/sr);
%     
%     [SOS,G] = zp2sos(z,p,k);
%     [SOS_det,G_det] = zp2sos(z_det,p_det,k_det);
    if exist('FiltFiltM','file')
    	xf_detect = FiltFiltM(b_detect, a_detect, x);
        xf = FiltFiltM(b, a, x); 
    else
        xf_detect = filtfilt(b_detect, a_detect, x);
        xf = filtfilt(b, a, x);
%         xf_detect = filtfilt(SOS_det, G_det, x);
%         xf = filtfilt(SOS,G, x);
        
    end
    
else
    xf = fix_filter(x);                   %Does a bandpass filtering between [300 3000] without the toolbox.
    xf_detect = xf;
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
                if iaux == 1 && ~any((xf(xaux(i)+1:xaux(i)+sample_ref))>thr)
                    continue
                end
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
                if iaux == 1 && ~any((xf(xaux(i)+1:xaux(i)+sample_ref))<thr)
                    continue
                end
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
                if iaux == 1 && ~any((abs(xf(xaux(i)+1:xaux(i)+sample_ref)))>thr)
                    continue
                end
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


