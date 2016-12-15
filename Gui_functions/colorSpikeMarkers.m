function colorSpikeMarkers(spikeMarkers, spikesClass, colors)
% Function that colors spike markers (tiny lines in the trace plot
% created by Plot_continuous_data) according to the current clusters.
% Arguments:
%   spikesMarkers - 1x2 cell array of arrays of handles to spike markers
%   spikesClass - classes of spikes
%   colors - cluster colors
%

    for i=0:max(spikesClass)
        idx=find(spikesClass==i);
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
            set(spikeMarkers{1}(idx),'color',clr);
            if ~isempty(spikeMarkers{2})
                set(spikeMarkers{2}(idx),'color',clr);
            end
        end
    end
end
