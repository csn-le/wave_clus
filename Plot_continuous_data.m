function Plot_continuous_data(xf_detect, sr_sub,handles)
% PLOT CONTINUOUS DATA

detect = handles.par.detection;
stdmin = handles.par.stdmin;
stdmax = handles.par.stdmax;

lx = length(xf_detect);

noise_std_detect = median(abs(xf_detect))/0.6745;

thr = stdmin * noise_std_detect;        %thr for detection is based on detect settings.
%thrmax = stdmax * noise_std_sorted;    %thr
thrmax = stdmax * noise_std_detect;     %aprox thrmax for plotting

cla(handles.cont_data);
hold(handles.cont_data, 'on');

plot(handles.cont_data, (1:lx)/sr_sub, xf_detect)

switch detect
    case 'pos'
        plot(handles.cont_data, [0 lx/sr_sub], [thr thr],'-r')
    case 'neg'
        plot(handles.cont_data, [0 lx/sr_sub], [-thr -thr],'-r')
    case 'both'
        plot(handles.cont_data, [0 lx/sr_sub], [thr thr],'-r')
        plot(handles.cont_data, [0 lx/sr_sub], [-thr -thr],'-r')
end

ypmax = min(thrmax, max(xf_detect));
ypmin = max(-thrmax/2, min(xf_detect));

%axis(handles.cont_data, [0 lx/sr_sub ypmin ypmax])
ylim(handles.cont_data, [ypmin ypmax])
xlim(handles.cont_data, [0 lx/sr_sub])

%xlabel('Time (sec)')
