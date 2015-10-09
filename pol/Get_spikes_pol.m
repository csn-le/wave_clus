function Get_spikes_pol(input, par)

sr = par.sr;
ls = par.w_pre + par.w_post; % length of the spike

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
    for j=1:par.segments        %that's for cutting the data into pieces
        poly_spikes = [];
        for i=1:length(channels)
            % LOAD CONTINUOUS DATA
            file_to_cluster = channels(i);
            eval(['load ' char(file_to_cluster) ';']);
            tsmin = (j-1)*floor(length(data)/par.segments)+1;
            tsmax = j*floor(length(data)/par.segments);
            x(i,:)=data(tsmin:tsmax); clear data;
            x = double(x);
        end
        % SPIKE DETECTION WITH AMPLITUDE THRESHOLDING
        [spikes, index] = amp_detect_pol(x,par); clear x;
       
        % JOIN SEGMENTS
        index=index*1e6/sr+tsmin;
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





