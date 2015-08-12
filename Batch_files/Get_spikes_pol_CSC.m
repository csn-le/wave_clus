function Get_spikes_pol_CSC(polytrodes)


% PARAMETERS
handles.par.w_pre = 20;                       %number of pre-event data points stored (default 20)
handles.par.w_post = 44;                      %number of post-event data points stored (default 44)
% handles.par.detection = 'both';              %type of threshold (default 'pos')
handles.par.detection = 'pos';              
% handles.par.detection = 'neg';
handles.par.stdmin = 5;                     %minimum threshold (default 5)
handles.par.stdmax = 50;                    %maximum threshold (default 50)
handles.par.interpolation = 'y';            %interpolation for alignment (default 'y')
handles.par.int_factor = 2;                 %interpolation factor (default 2)
handles.par.detect_fmin = 300;              %high pass filter for detection (default 300)
handles.par.detect_fmax = 1000;             %low pass filter for detection (default 1000)
handles.par.sort_fmin = 300;                %high pass filter for sorting (default 300)
handles.par.sort_fmax = 3000;               %low pass filter for sorting (default 3000)
handles.par.segments_length = 5;            %length of segments in which the data is cutted (default 5min).
handles.par.features = 'wav';               %choice of spike features wav,pca (default 'wav')
handles.par.inputs = 40;                    %number of inputs to the clustering (default 10)
handles.par.scales = 6;                     %scales for wavelet decomposition (default 4)
handles.awin = 4; % alignment window

ls = handles.par.w_pre + handles.par.w_post; % length of the spike

for k=1:length(polytrodes)

    polytrode = polytrodes(k);

    % LOAD POLYTRODE CHANNELS
    pol = strcat('polytrode',num2str(polytrode),'.txt');
    channels = textread(pol,'%s');
    filename = channels(1);
    filename = cell2mat(filename);
    
    % Get TimeStamps
    f = fopen(filename,'r','l');
    fseek(f,16384,'bof'); % Skip Header, put pointer to the first record
    TimeStamps = fread(f,inf,'int64',(4+4+4+2*512)); %Read all TimeStamps in microseconds
    fclose(f);
    
    % Find the starting of the recording and gets sampling frequency
    time0 = TimeStamps(1);                              %time in microseconds
    timeend = TimeStamps(end);
    sr = 512*1e6/(TimeStamps(2)-TimeStamps(1));         %seconds^(-1), ~28KHz sampling rate (~36microseconds every sample)
    min_ref_per = 1.5;                                    %minimum refractory period (in ms)
    ref = floor(min_ref_per *sr/1000);                  %number of counts corresponding the minimum refractory period
    handles.par.sr = sr; handles.par.ref = ref;
    handles.par.segments = ceil((timeend - time0) / (handles.par.segments_length * 1e6 * 60)); %number of segments in which data is cutted
    
    %That's for cutting the data into pieces
    lts = length(TimeStamps);
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
    clear TimeStamps;
        
    index_all = [];
    spikes_all = [];
    for j=1:length(tsmin)        %that's for cutting the data into pieces
        poly_spikes = [];
        for i=1:length(channels)
            if j==1
                filename = char(channels(i));
                f(i) = fopen(filename,'r','l');
                fseek(f(i),16384+8+4+4+4,'bof'); % put pointer to the beginning of data

                %GETS THE GAIN AND CONVERTS THE DATA TO MICRO V.
                eval(['scale_factor=textread(''CSC' filename(4) '.Ncs'',''%s'',41);']);
            end
            
            % LOAD CSC DATA
            Samples = fread(f(i),512*(recmax(j)-recmin(j)+1),'512*int16=>int16',8+4+4+4);
            x(i,:) = double(Samples(:))'; clear Samples;
            x(i,:) = x(i,:)*str2num(scale_factor{41})*1e6;
        end
        % SPIKE DETECTION WITH AMPLITUDE THRESHOLDING
        [spikes, index] = amp_detect_pol(x,handles); clear x;
        size(spikes,1)
        % JOIN SEGMENTS
        index=(index+tsmin(j))*1e6/sr;
        index_all = [index_all index]; clear index;
        if length(spikes)~=0
            for i=1:length(channels)
                poly_spikes = [poly_spikes spikes(:,1+(i-1)*ls:i*ls)];
            end
            spikes_all = [spikes_all; poly_spikes];
        end
    end
    clear spikes;
    clear poly_spikes;
    index = index_all/1000; % time in msecs

    spikes = zeros(size(spikes_all,1),size(spikes_all,2));
    spikes = spikes_all;
    clear spikes_all;

    outfile = strcat(pol(1:end-4),'_spikes');
    save(outfile, 'spikes', 'index')
    clear spikes;

end





