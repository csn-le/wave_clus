function [spikes] = spike_alignment(spikes, handles)


% spike_alignment_batch
% function for alligning spikes detected outside wave_clus.
% INPUT: Matrix of spikes in an ascii file. Spikes must be grouped in a
% matrix called 'spikes' with a spike per row.
% OUTPUT: matrix of spikes alligned with the maximum in sample set by the
% parameters. This is necessary for a correct clustering of the spikes.

% handles.par.w_pre=20;                       %number of pre-event data points stored
% handles.par.w_post=44;                      %number of post-event data points stored
% handles.par.detection = 'pos';              %type of threshold
% handles.par.interpolation = 'y';            %interpolation for alignment
% handles.par.int_factor = 2;                 %interpolation factor 
% handles.par.sr = 24000;                     %sampling frequency, in Hz.
% handles.par.alignment_window = 10;          %number of sample points around the maximum 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

w_pre = handles.par.w_pre;                      %number of pre-event data points stored
w_post = handles.par.w_post;                    %number of post-event data points stored
sr = handles.par.sr;
align_window = handles.par.alignment_window;    %number of sample points for re-alignment
detec_threshold = handles.par.detection;        %type of threshold

ls = w_pre + w_post;


% Introduces first alignment
spikes1 = zeros(size(spikes,1),ls+(2*align_window)+4);
if (size(spikes,2)< (ls + align_window)+2)
    diff_size = ls+align_window+2-size(spikes,2);
    spikes1(:,1:align_window+2) = -spikes(:,align_window+2:-1:1);
    spikes1(:,1+align_window:align_window+size(spikes,2)) = spikes;
    spikes1(:,1+align_window+size(spikes,2)+2:end) = -spikes(:,end:-1:end-diff_size+1);
else
    spikes1(:,1:align_window+2) = -spikes(:,align_window:-1:1);
    spikes1(:,1+align_window:(2*align_window+ls)) = spikes(1:align_window+ls);
end
correct_times = zeros(size(spikes,1),1);
spikes2 = zeros(size(spikes,1),ls+4);
for i=1:size(spikes1,1)

    if strcmp(handles.par.detection, 'pos')
        [maxi iaux] = max(spikes1(i,w_pre+2:w_pre+2*align_window+1));    %introduces alignment
    else 
        [mini iaux] = min(spikes1(i,w_pre+2:w_pre+2*align_window+1));    %introduces alignment
    end;
    if iaux > 1
        spikes2(i,:) = spikes1(i,iaux-1:iaux+ls+2);
    else 
        spikes2(i,:) = spikes1(i,iaux:iaux+ls+3);
    end
    correct_times(i) = (iaux+w_pre-align_window+1)*1/sr;  %corrects spike-times.
end

spikes = spikes2;

switch handles.par.interpolation
case 'n'
    spikes(:,end-1:end)=[];       %eliminates borders that were introduced for interpolation 
    spikes(:,1:2)=[];
case 'y'
    %Does interpolation
    spikes = int_spikes(spikes,handles);   
end;

