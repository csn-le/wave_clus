function Get_spikes_CSC_new(channels,TimeStamps)
% function Get_spikes_CSC_new(channel,TimeStamps)
% Gets spikes from all channels in the .mat file channels. This batch is
% to be used with Neuralynx data.
% Saves spikes and spike times corresponding to the Neuralynx timestamps.
% Modified 27/3/08: Notch filter option has been added to remove 2000KHz
% noise that appears in some sessions.
% 
if nargin <1
    load channels
end

handles.par.w_pre=20;                       %number of pre-event data points stored (default 20)
handles.par.w_post=44;                      %number of post-event data points stored (default 44)
handles.par.detection = 'pos';              %type of threshold (default 'pos')
handles.par.stdmin = 5;                     %minimum threshold (default 5)
handles.par.stdmax = 50;                    %maximum threshold (default 50)
handles.par.interpolation = 'y';            %interpolation for alignment (default 'y')
handles.par.int_factor = 2;                 %interpolation factor (default 2)
handles.par.detect_fmin = 300;              %high pass filter for detection (default 300)
handles.par.detect_fmax = 1000;             %low pass filter for detection (default 1000)
handles.par.sort_fmin = 300;                %high pass filter for sorting (default 300)
handles.par.sort_fmax = 3000;               %low pass filter for sorting (default 3000)
handles.par.segments_length = 5;            %length of segments in which the data is cutted (default 5min)
handles.par.notchfilter = 0;                %if set to 1, applies a notch filter to remove 2000KHz noise (default 0) 

for k= 1:length(channels)

    tic
    index_all=[];
    spikes_all=[];
    threshold_all=[];
    channel=channels(k)
    filename=sprintf('CSC%d.Ncs',channel);
    f=fopen(filename,'r','l');

	if nargin <2
        fseek(f,16384,'bof'); % Skip Header, put pointer to the first record
        TimeStamps=fread(f,inf,'int64',(4+4+4+2*512)); %Read all TimeStamps
	end

	%Find the starting of the recording and gets sampling frequency
 	time0 = TimeStamps(1); 
 	timeend = TimeStamps(end);
%	dt = TimeStamps(2)-TimeStamps(1);
	lts = length(TimeStamps);  
    dt=min(diff(TimeStamps));
    if dt<=0, %check for corrupt TimeStamps that would violate the monotonic increase
        warning('corrupt TimeStamps - attempting correction')
        tsdiff=diff(TimeStamps);
        dt=median(tsdiff);
        ind=find(tsdiff<=0); %corrupt TimeStamps are usually too low
        TimeStamps(ind+1)=TimeStamps(ind)+dt;
        clear tsdiff;
        dt=min(diff(TimeStamps));
        if dt<=0, error('automated correction of corrupt TimeStamps failed - try manually!'); end
    end
	sr = 512*1e6/dt;
	min_ref_per=1.5;                                    %detector dead time (in ms)
	ref = floor(min_ref_per *sr/1000);                  %number of counts corresponding the dead time
	handles.par.sr = sr; handles.par.ref = ref;
	handles.par.segments = ceil((lts*dt) / (handles.par.segments_length * 1e6 * 60)); %number of segments in which data is cut
%	handles.par.segments = ceil((timeend - time0) /(handles.par.segments_length * 1e6 * 60)); %number of segments in which data is cut
	
	%That's for cutting the data into pieces
	segmentLength = floor (lts/handles.par.segments);
	tsmin = 1 : segmentLength :lts;
	tsmin = tsmin(1: handles.par.segments);
	tsmax = tsmin - 1;
	tsmax = tsmax (2:end);
	tsmax = [tsmax, lts];
	recmax=tsmax;
	recmin=tsmin;
	tsmin = TimeStamps(int64(tsmin));
	tsmax = TimeStamps(int64(tsmax));
	
    fseek(f,16384+8+4+4+4,'bof'); % put pointer to the beginning of data
    
    %GETS THE GAIN AND CONVERTS THE DATA TO MICRO V.
    eval(['scale_factor=textread(''CSC' num2str(channel) '.Ncs'',''%s'',43);']);
    
    for j=1:length(tsmin)
        % LOAD CSC DATA
        Samples=fread(f,512*(recmax(j)-recmin(j)+1),'512*int16=>int16',8+4+4+4);
        x=double(Samples(:))'; clear Samples;
        x=x*str2num(scale_factor{43})*1e6;
        if(handles.par.notchfilter==1)
            [b,a]=ellip(2,0.5,20,[1999 2001]*2/sr,'stop');
            x=filtfilt(b,a,x);
        end
        
        % SPIKE DETECTION WITH AMPLITUDE THRESHOLDING
        [spikes,thr,index]  = amp_detect(x,handles);       %detection with amp. thresh.
%        index=index*1e6/sr+tsmin(j);
        index_ts=TimeStamps(recmin(j)+floor(index/512))'+mod(index,512)/512*dt;
        index_all = [index_all index_ts];
        spikes_all = [spikes_all; spikes];
        threshold_all=[threshold_all thr];
    end
    fclose(f);
    %index = (index_all-time0)/1000; %this was our old convention, which corrects breaks in the recording a posteriori
    index_ts = index_all/1000; %this is the new convention, saving absolute times
    spikes = spikes_all;
    eval(['save CSC' num2str(channel) '_spikes spikes index_ts threshold_all']); %saves Sc files
    toc
end   
