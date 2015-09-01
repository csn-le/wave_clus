function Get_spikes_NSX_new(channels,notchfilter,detection_type)
% Gets spikes from continuous data. This batch is to be used with Neuroport data.
% Saves spikes and spike times corresponding to the timestamps.
% index_ts %absolute times (in msec starting from the start of the recording)
% MJI, 27/7/2010
% detection_type 'both', 'neg' or 'pos'
% notchfilter set to 0 to avoid using notches before detection. If set to 1, it needs "check_lfp_power_NSX" to be
% run first (so that the notch frequencies can be computed)
% HGR FEB 2013
% Now saving handles.par
% Includes the option of using notch filters for line cancellation before doing spike detection
% If the channels are arranged in a 2 by n matrix, it uses the signals from
% the difference between the channels on each column. If it is a row vector it
% uses the signals from the individual channels

if nargin <1
    load channels
end
if ~exist('detection_type','var'), detection_type='both'; end
if ~exist('notchfilter','var'), notchfilter = 0; end

direc_resus_base = pwd;
% if ~exist(fullfile(direc_resus_base,'test'),'dir')
%         mkdir(direc_resus_base,'test');
% end
% direc_resus = fullfile(direc_resus_base,'test');
direc_resus = direc_resus_base;
if notchfilter
    mkdir(direc_resus,'notch');
    direc_resus = fullfile(direc_resus,'notch');
    if exist('stimulus.mat','file')
        copyfile('stimulus.mat',fullfile(direc_resus,'stimulus.mat'));
    else
        fprintf('stimulus has not been created yet\n')
    end
    if exist('sites.mat','file')
        copyfile('sites.mat',fullfile(direc_resus,'sites.mat'));
    else
        fprintf('sites has not been created yet\n')
    end
    if exist('finalevents.mat','file')
        copyfile('finalevents.mat',fullfile(direc_resus,'finalevents.mat'));
    else
        fprintf('finalevents has not been created yet\n')
    end
end

handles.par.w_pre=20;                       %number of pre-event data points stored (default 20)
handles.par.w_post=44;                      %number of post-event data points stored (default 44)
handles.par.detection = detection_type;     %type of threshold ('pos', 'neg', 'both', default 'neg')
handles.par.stdmin = 5;                     %minimum threshold (default 5)
handles.par.stdmax = 50;                    %maximum threshold (default 50)
handles.par.interpolation = 'y';            %interpolation for alignment (default 'y')
handles.par.int_factor = 2;                 %interpolation factor (default 2)
handles.par.detect_fmin = 300;              %high pass filter for detection (default 300)
handles.par.detect_fmax = 3000;             %low pass filter for detection (default 1000)
handles.par.sort_fmin = 300;                %high pass filter for sorting (default 300)
handles.par.sort_fmax = 3000;               %low pass filter for sorting (default 3000)
handles.par.segments_length = 5;            %length of segments in which the data is cutted (default 5min) 
handles.par.notchfilter = notchfilter;               

if size(channels,1)==2
    diff_chs = 1;
else
    diff_chs = 0;
end
    
for k= 1:size(channels,2)
    tic
    index_all=[];
    spikes_all=[];
    threshold_all=[];
    channel1=channels(1,k)
    filename1=sprintf('NSX%d.NC5',channel1);    
    if diff_chs
        channel2=channels(2,k)
        filename2=sprintf('NSX%d.NC5',channel2);
    end
    
    %Read TimeStamps
	load('NSX_TimeStamps');
   
    %sr and lts read from NSX_TimeStamps
    min_ref_per=1.5;                                    %detector dead time (in ms)
    ref = floor(min_ref_per*sr/1000);                  %number of counts corresponding the dead time
    handles.par.sr = sr; handles.par.ref = ref;
       
    %That's for cutting the data into pieces
    handles.par.segments=ceil(lts/(handles.par.segments_length  * sr * 60)); %number of segments in which data is cut
    segmentLength = floor (lts/handles.par.segments);
    tsmin = 1 : segmentLength :lts;
    tsmin = tsmin(1: handles.par.segments);
    tsmax = tsmin - 1;
    tsmax = tsmax (2:end);
    tsmax = [tsmax, lts];
    recmin = tsmin; %in samples
    recmax = tsmax; %in samples
    tsmin = TimeStamps(tsmin); %in microsec
    tsmax = TimeStamps(tsmax); %in microsec
        
    for j=1:length(tsmin)
        % LOAD NSX DATA
        Samples = read_NC5(filename1,recmin(j),recmax(j));
        if(handles.par.notchfilter~=0)
            [Samples,~,~,handles.par.que_armon1,handles.par.qsnew1] = mete_notches(Samples,sr,channel1,handles.par.detect_fmin,direc_resus_base); 
        end
        
        if diff_chs
            Samples2 = read_NC5(filename2,recmin(j),recmax(j));
            if(handles.par.notchfilter~=0)
                [Samples2,~,~,handles.par.que_armon2,handles.par.qsnew2] = mete_notches(Samples2,sr,channel2,handles.par.detect_fmin,direc_resus_base);
            end
            Samples = Samples - Samples2;
        end
        x=Samples';clear Samples;        
                
        % SPIKE DETECTION WITH AMPLITUDE THRESHOLDING
        [spikes,thr,index]  = amp_detect(x,handles);       %detection with amp. thresh.
        index_ts=TimeStamps(recmin(j)+index);
        index_all = [index_all index_ts];
        spikes_all = [spikes_all; spikes];
        threshold_all=[threshold_all thr];
    end
                    
    index_ts = index_all/1000; %absolute times (spikes in msec starting from the start of the recording)
    spikes = spikes_all;
    par = handles.par;
    
    if diff_chs
        save(fullfile(direc_resus,['NSX' num2str(channel1) '-' num2str(channel2) '_spikes']),'spikes','index_ts','threshold_all','par');
    else
        save(fullfile(direc_resus,['NSX' num2str(channel1) '_spikes']),'spikes','index_ts','threshold_all','par');
    end
    toc
end   
