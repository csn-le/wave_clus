function xf_detect = spike_detection_filter(x, par)
%this function filter the signal, using the detection filter. Is used in the
%readInData class. 

sr = par.sr;
fmin_detect = par.detect_fmin;
fmax_detect = par.detect_fmax;


% HIGH-PASS FILTER OF THE DATA
if exist('ellip','file')                         %Checks for the signal processing toolbox
    [b,a] = ellip(par.detect_order,0.1,40,[fmin_detect fmax_detect]*2/sr);
    [b_notch,a_notch] = butter(par.detect_order,[55 65]*2/sr, 'stop');
    if exist('FiltFiltM','file')
        xf_detect = FiltFiltM(b_notch, a_notch, x);
        xf_detect = FiltFiltM(b, a, xf_detect);
    else
        xf_detect = filtfilt(b_notch, a_notch, x);
        xf_detect = filtfilt(b, a, xf_detect);
    end
else
    xf_detect = fix_filter(x);                   %Does a bandpass filtering between [300 3000] without the toolbox.
end
