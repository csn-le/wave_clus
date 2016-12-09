function Plot_continuous_data(xf_detect, sr_sub,handles, spikesIdx, spikesClass, colors)
% Function that plot the segment of continuous data. 
% It will try to plot a minute, maximun.
% Optional arguments:
%   - spikesIdx - indices of samples of detected spikes
%   - spikesClass - classes of detected spikes
%   - colors - cluster colors

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

if nargin==6
    % mark positions of detected spikes

    hold on

    % update ylim to leave more space for markers and 
    % compute y-position where to plot spike markers
    ylo=yl;
    switch detect
        case 'pos'
            yl=[yl(1) yl(2)*1.2];
            y1=[yl(2) ylo(2)];
            y2=[];
        case 'neg'
            yl=[yl(1)*1.2 yl(2)];
            y1=[yl(1) ylo(1)];
            y2=[];
        case 'both'
            yl=yl*1.2;
            y1=[yl(2) ylo(2)];
            y2=[yl(1) ylo(1)];
    end
    ylim(handles.cont_data,yl);

    % plot
    for i=0:max(spikesClass)
        idx=find(spikesClass==i);
        n=length(idx);
        if ~isempty(idx)
            if i>0
                clr=colors(1+mod(i-1,length(colors)),:);
                % make the color a little bit lighter to ease
                % its visual perception on the white background
                clr=min([1 1 1],clr+[.2 .2 .2]);
            else
                % cluster 0 plotted in black
                clr=[0 0 0];
            end
            % samples of detected spikes
            xx=spikesIdx(idx);
            % force them to form a row vector
            if size(xx,1)>1 xx=xx'; end
            plot(handles.cont_data, repmat(xx,2,1)/1000,repmat(y1',1,n),...
                'color',clr,'LineWidth',3);
            if ~isempty(y2)
                plot(handles.cont_data, repmat(xx,2,1)/1000,repmat(y2',1,n),...
                    'color',clr,'LineWidth',3);
            end
        end
    end
    hold off
end

%xlabel('Time (sec)')
