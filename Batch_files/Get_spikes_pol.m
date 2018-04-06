function Get_spikes_pol(polytrodes, varargin)
% function Get_spikes_pol(polytrodes, par_input)
% Make polytrode spikes detecting and concatenating the spikes of all the
% channels in the polytrodeN.txt used. 
% Saves spikes, spike times (in ms) and used parameters in filename_spikes.mat.

% input must be a vector with the numbers(N) of polytrodeN.txt to use.

% par_input must be a struct with some of the detecction parameters. All the
% parameters included will overwrite the parameters load from set_parameters()


%default config
par_input = struct;
parallel = false;

%optinal inputs
nvar = length(varargin);
for v = 1:2:nvar
    if strcmp(varargin{v},'par')
        if (nvar>=v+1) && isstruct(varargin{v+1})
            par_input = varargin{v+1};
        else
            error('Error in ''par'' optional input.')
        end
    elseif strcmp(varargin{v},'parallel')
        if (nvar>=v+1) && islogical(varargin{v+1})
            parallel = varargin{v+1};
        else
            error('Error in ''parallel'' optional input.')
        end
    end
end

run_parfor = parallel;

if parallel == true
    if exist('matlabpool','file')
        try
            matlabpool('open');
        catch
            parallel = false;
        end
    else
        poolobj = gcp('nocreate'); % If no pool, do not create new one.
        if isempty(poolobj)
            parpool
        else
            parallel = false;
        end
    end
end

init_date = now;
if run_parfor == true
    parfor j = 1:length(polytrodes)
        get_spikes_pol_single(polytrodes(j), par_input);
        disp(sprintf('%d of %d ''spikes'' files done.',count_new_sp_files(init_date, polytrodes),length(polytrodes)))
    end
else
    for j = 1:length(polytrodes)
        get_spikes_pol_single(polytrodes(j), par_input);
        disp(sprintf('%d of %d ''spikes'' files done.',count_new_sp_files(init_date, polytrodes),length(polytrodes)))
    end
    
end

if parallel == true
    if exist('matlabpool','file')
        matlabpool('close')
    else
        poolobj = gcp('nocreate');
        delete(poolobj);
    end
end


end

function get_spikes_pol_single(polytrode, par_input)
    par = set_parameters();
    par = update_parameters(par,par_input,'detect');


    sr = par.sr;
    ls = par.w_pre + par.w_post; % length of the spike
    par.reset_results = true;
    index_all = [];
    spikes_all = [];
    
    % LOAD POLYTRODE CHANNELS
    pol = strcat('polytrode',num2str(polytrode),'.txt');
    out_filename = pol(1:end-4);
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
    thr = cell(min_num_seg,1);
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
    
        [spikes, new_index,thr{n}] = amp_detect_pol(x,par); clear x;
        index = [index data_handler_ch{1}.index2ts(new_index)];
        if length(spikes)~=0
            for i=1:n_channels
                poly_spikes = [poly_spikes spikes(:,1+(i-1)*ls:i*ls)];
            end
            spikes_all = [spikes_all; poly_spikes];
        end

        
    end

    spikes = spikes_all;
    current_par = par;
    par = struct;
    par = update_parameters(par, current_par, 'detect');
    par.detection_date =  datestr(now);
	par.channels = n_channels;
    thr = cell2mat(thr);
	try
		save([out_filename '_spikes'], 'spikes', 'index','par','thr')
	catch
		save([out_filename '_spikes'], 'spikes', 'index','par','thr','-v7.3')
	end
end

function counter = count_new_sp_files(initial_date, polytrodes)
counter = 0;
for i = 1:length(polytrodes)
    polytrode = polytrodes(i);
    sp_file = strcat('polytrode',num2str(polytrode),'_spikes.mat');
    FileInfo = dir(sp_file);
    if length(FileInfo)==1 && (FileInfo.datenum > initial_date)
        counter = counter + 1;
    end
end
end
