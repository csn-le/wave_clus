function Get_spikes_pol(polytrodes, par_input)
% function Get_spikes_pol(polytrodes, par_input)
% Make polytrode spikes detecting and concatenating the spikes of all the
% channels in the polytrodesN.txt used. 
% Saves spikes, spike times (in ms) and used parameters in filename_spikes.mat.

% input must be a vector with the numbers(N) of polytrodesN.txt to use.

% par_input must be a cell with some of the detecction parameters. All the
% parameters included will overwrite the parameters load from set_parameters()



for k = 1:length(polytrodes)
    par = set_parameters();
    if exist('par_input','var')
        par = update_parameters(par,par_input,'relevant');
    end

    sr = par.sr;
    ls = par.w_pre + par.w_post; % length of the spike
    par.reset_results = true;
    polytrode = polytrodes(k);
    index_all = [];
    spikes_all = [];
    
    % LOAD POLYTRODE CHANNELS
    pol = strcat('polytrode',num2str(polytrode),'.txt');
    out_filename = strcat(pol(1:end-4));
    electrodes = textread(pol,'%s');
    n_channels = length(electrodes);
    
    data_handler_ch{n_channels} = [];
    par_ch{n_channels} = [];
    min_num_seg = inf;
    
    for i=1:n_channels
        par_ch{i} = par;
        par_ch{i}.filename = electrodes{i};
        data_handler_ch{i} = readInData(par_ch{i});
        if ~data_handler_ch{i}.with_raw
        	ME = MException('MyComponent:FileError', 'The file %s doesn''t have raw data',electrodes{i});
            throw(ME)
        end
        
        par_ch{i} = data_handler_ch{i}.par;
        min_num_seg =  min(data_handler_ch{i}.max_segments, min_num_seg);
    end
    
    

    
    index = [];
    spikes_all = [];
    
    for n = 1:min_num_seg   %that's for cutting the data into pieces
        x(1,:) = data_handler_ch{1}.get_segment();
        poly_spikes = [];
        x_len = length(x(1,:));
        x(n_channels,:) = zeros(1,x_len); %pre-allocating memory
        for i=2:n_channels
            xaux = data_handler_ch{i}.get_segment();
            if length(xaux) == x_len
                x(i,:) = xaux;
            elseif length(xaux)>x_len
                x(i,:) = xaux(1:x_len);
            elseif length(xaux)<x_len
                x_len = length(xaux);
                x = x(1:end,1:x_len);
                x(i,:) = xaux;
            end
        end       
    
        [spikes, new_index] = amp_detect_pol(x,par); clear x;
        index = [index data_handler_ch{1}.index2ts(new_index)];
        if length(spikes)~=0
            for i=1:n_channels
                poly_spikes = [poly_spikes spikes(:,1+(i-1)*ls:i*ls)];
            end
            spikes_all = [spikes_all; poly_spikes];
        end

        
    end

    
    clear spikes;
    clear poly_spikes;


    spikes = spikes_all;
    clear spikes_all;

    current_par = par;
    par = struct;
    par.channels = n_channels;
    par = update_parameters(par, current_par, 'detect');
    save([out_filename '_spikes'], 'spikes', 'index','par')
    clear spikes;

end





