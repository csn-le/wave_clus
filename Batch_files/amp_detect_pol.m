function [spikes,index] = amp_detect_pol(x,handles)
%detect spikes in a tetrode

%PARAMETERS
sr = handles.par.sr; % ~28KHz
w_pre = handles.par.w_pre;
w_post = handles.par.w_post;
ref = handles.par.ref; % refractory period (value in counts)
detect = handles.par.detection;
stdmin = handles.par.stdmin;
stdmax = handles.par.stdmax;
fmin_detect = handles.par.detect_fmin;
fmax_detect = handles.par.detect_fmax;
fmin_sort = handles.par.sort_fmin;
fmax_sort = handles.par.sort_fmax;
awin = handles.awin;  

% SPIKE DETECTION FOR EACH CHANNEL
% Detects all spike times including all channels
% there will be spike times very close
% probably belonged to the same spike at 
% different channels. 
index = [];
for i=1:size(x,1)
   [indexch,xf(i,:),thr(i)] = index_detect_pol(x(i,:),handles); % spike times and filtered data
   index = [index indexch];
   clear indexch;
end

% ALIGNMENT OF THE POLYTRODE CHANNELS
% Only spike times separated more than the refractory period 
% will be taken. The minimum distance between two 
% different spikes allowed by the detection algorithm 
% is the refractory period 
index = sort(index);
index( find(abs(diff(index))<ref/2) ) = []; 



% SPIKE CONCATENATION: POLY-SPIKE
for i=1:length(index)
    for j=1:size(x,1)
        spikes(i,1+(j-1)*(w_pre+w_post+2*awin):j*(w_pre+w_post+2*awin)) = xf(j,index(i)-w_pre-awin+1:index(i)+w_post+awin);
    end
end

% POLY-SPIKE INTERPOLATION
% interpolation and downsampling is done separately in every channel within
% each poly-spike which does the interchannel alignment
switch handles.par.interpolation
    case 'n'
        spikes(:,end-1:end)=[];       %eliminates borders that were introduced for interpolation 
        spikes(:,1:2)=[];
    case 'y'
        %Does interpolation
        spikes = int_spikes_pol(spikes,size(x,1),thr,handles);   
end









    