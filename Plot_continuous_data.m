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

axes(handles.cont_data)
cla
hold on

plot((1:lx)/sr_sub, xf_detect)

switch detect
    case 'pos'
        line([0 lx/sr_sub], [thr thr],'color','r')
    case 'neg'
        line([0 lx/sr_sub], [-thr -thr],'color','r')
    case 'both'
        line([0 lx/sr_sub], [thr thr],'color','r')
        line([0 lx/sr_sub], [-thr -thr],'color','r')
end

axis([0 lx/sr_sub -thrmax/2 thrmax])
%xlabel('Time (sec)')
