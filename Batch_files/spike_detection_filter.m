function xf_detect = spike_detection_filter(x, par)
%this function filter the signal, using the detection filter. Is used in the
%readInData class. 

sr = par.sr;
fmin_detect = par.detect_fmin;
fmax_detect = par.detect_fmax;


% HIGH-PASS FILTER OF THE DATA
if exist('ellip','file')                         %Checks for the signal processing toolbox
    [b,a] = ellip(par.detect_order,0.1,40,[fmin_detect fmax_detect]*2/sr);
    if exist('FiltFiltM','file')
        xf_detect = FiltFiltM(b, a, x);
    else
        xf_detect = filtfilt(b, a, x);
    end
else
    xf_detect = fix_filter(x);                   %Does a bandpass filtering between [300 3000] without the toolbox.
end
