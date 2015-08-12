function Get_spikes_pol(polytrodes)


% PARAMETERS
handles.par.w_pre = 20;                       %number of pre-event data points stored (default 20)
handles.par.w_post = 44;                      %number of post-event data points stored (default 44)
% handles.par.detection = 'both';              %type of threshold (default 'pos')
handles.par.detection = 'pos';              
% handles.par.detection = 'neg';
handles.par.stdmin = 5;                     %minimum threshold (default 5)
handles.par.stdmax = 100;                    %maximum threshold (default 50)
handles.par.interpolation = 'y';            %interpolation for alignment (default 'y')
handles.par.int_factor = 2;                 %interpolation factor (default 2)
handles.par.detect_fmin = 300;              %high pass filter for detection (default 300)
handles.par.detect_fmax = 3000;             %low pass filter for detection (default 3000)
handles.par.sort_fmin = 300;                %high pass filter for sorting (default 300)
handles.par.sort_fmax = 3000;               %low pass filter for sorting (default 3000)
handles.par.segments = 10;                  %length of segments in which the data is cutted 
handles.par.features = 'wav';               %choice of spike features wav,pca (default 'wav')
handles.par.inputs = 20;                    %number of inputs to the clustering (default 10)
handles.par.scales = 6;                     %scales for wavelet decomposition (default 4)
handles.par.sr = 24414.1;                     % sampling frequency in Hz
min_ref_per = 2;   %minimum refractory period (in ms)
handles.par.ref = floor(min_ref_per *handles.par.sr/1000); %number of counts corresponding the minimum refractory period
handles.awin = 8; % alignment window

ref = handles.par.ref;
sr = handles.par.sr;
ls = handles.par.w_pre + handles.par.w_post; % length of the spike


for k=1:length(polytrodes)

    polytrode = polytrodes(k);

    % LOAD POLYTRODE CHANNELS
    pol = strcat('polytrode',num2str(polytrode),'.txt');
    channels = textread(pol,'%s');
    filename = channels(1);
    filename = cell2mat(filename);
    filename = strcat(pol(1:end-4));
    
    index_all = [];
    spikes_all = [];
    for j=1:handles.par.segments        %that's for cutting the data into pieces
        poly_spikes = [];
        for i=1:length(channels)
            % LOAD CONTINUOUS DATA
            file_to_cluster = channels(i);
            eval(['load ' char(file_to_cluster) ';']);
            tsmin = (j-1)*floor(length(data)/handles.par.segments)+1;
            tsmax = j*floor(length(data)/handles.par.segments);
            x(i,:)=data(tsmin:tsmax); clear data;
            x = double(x);
        end
        % SPIKE DETECTION WITH AMPLITUDE THRESHOLDING
        [spikes, index] = amp_detect_pol(x,handles); clear x;
       
        % JOIN SEGMENTS
        index=(index+tsmin)*1e6/sr;
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

    outfile = strcat(filename,'_spikes');
    save(outfile, 'spikes', 'index')
    clear spikes;

end





