function Plot_continuous_data(xf_detect, sr_sub,handles)
% Function that plot the segment of continuous data. 
% It will try to plot a minute, maximun.

detect = handles.par.detection;
stdmin = handles.par.stdmin;
stdmax = handles.par.stdmax;

lx = length(xf_detect);

noise_std_detect = median(abs(xf_detect))/0.6745;

thr = stdmin * noise_std_detect;        %thr for detection is based on detect settings.
throut  = stdmax * noise_std_detect;    %thrmax
max_ylim = 15 * noise_std_detect;     %ylim for plotting. aprox thrmax

cla(handles.cont_data);
hold(handles.cont_data, 'on');

plot(handles.cont_data, (1:lx)/sr_sub, xf_detect)

switch detect
    case 'pos'
        plot(handles.cont_data, [0 lx/sr_sub], [thr thr],'-r')
        ylim(handles.cont_data, [-max_ylim/2 max_ylim])
    case 'neg'
        plot(handles.cont_data, [0 lx/sr_sub], [-thr -thr],'-r')
        ylim(handles.cont_data, [-max_ylim max_ylim/2])
    case 'both'
        plot(handles.cont_data, [0 lx/sr_sub], [thr thr],'-r')
        plot(handles.cont_data, [0 lx/sr_sub], [-thr -thr],'-r')
        ylim(handles.cont_data, [-max_ylim max_ylim])
end

xlim(handles.cont_data, [0 lx/sr_sub])
yl = ylim(handles.cont_data);
plot(handles.cont_data, [0 lx/sr_sub], [throut throut],':r')
plot(handles.cont_data, [0 lx/sr_sub], [-throut -throut],':r')
ylim(handles.cont_data,yl);

%xlabel('Time (sec)')
