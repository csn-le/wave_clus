function Get_spikes(input, par_input)
% function Get_spikes(input, par_input)


% Saves spikes and spike times.

if isnumeric(input) || strcmp(input,'all')
    filenames = {};
    se = supported_wc_extensions();
    dirnames = dir();
    dirnames = {dirnames.name};
    
    for i = 1:length(dirnames)
        fname = dirnames{i};
        [unused, f, ext] = fileparts(fname);
        ext = lower(ext(2:end));
        if any(strcmp(ext,se)) 
            if strcmp(ext,'mat')
                sprintf('Skipped file ''%s''. The ''.mat'' files should be added by name.\n',fname);
                continue
            end
            if strcmp(input,'all')
                filenames = [filenames {fname}];
            else
                aux = regexp(fname, '\d+', 'match');
                if ismember(str2num(aux{1}),input)
                    filenames = [filenames {fname}];   
                end
            end
        end
    end
    
elseif ischar(input) && length(input) > 4
    if  strcmp (input(end-3:end),'.txt')
        filenames =  textread(input,'%s');
    else
        filenames = {input};
    end
    
elseif iscellstr(input)
    filenames = input;
else
    ME = MException('MyComponent:noValidInput', 'Invalid input arguments');
    throw(ME)
end

for fnum = 1:length(filenames)
    filename = filenames{fnum};
    par = set_parameters();
    par.cont_segment = true;
    par.filename = filename;
    par.reset_results = true;

    par.cont_segment = true;  %false to don't save the sample in spikes
    
    data_handler = readInData(par);
    
    par = data_handler.par;
    if exist('par_input','var')
        par = update_parameters(par,par_input,'detect');
    end

    if data_handler.with_spikes            %data have some time of _spikes files
        [spikes, index] = data_handler.load_spikes(); 
        if ~data_handler.with_wc_spikes
            [spikes] = spike_alignment(spikes,par);
        end
    else    
        index = [];
        spikes = [];
        for n = 1:data_handler.max_segments
            x = data_handler.get_segment();
                %<----  Add here extra processing of the signal (x)
            [new_spikes, temp_aux_th, new_index]  = amp_detect(x, par);
            index = [index data_handler.index2ts(new_index)]; %new_index to ms
            spikes = [spikes; new_spikes];
        end
    end

    current_par = par;
    par = struct;
    par = update_parameters(par, current_par, 'detect');
    
    
    %<----  Add here auxiliar parameters

    if current_par.cont_segment
        [psegment, sr_psegment] = data_handler.get_signal_sample();
        save([data_handler.nick_name '_spikes'], 'spikes', 'index', 'par','psegment','sr_psegment')
        clear psegment
    else
        save([data_handler.nick_name '_spikes'], 'spikes', 'index', 'par')
    end

    clear spikes
    
end
end
