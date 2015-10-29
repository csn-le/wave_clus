function Get_spikes(input, parallel, par_input)
% function Get_spikes(input, par_input)
% Saves spikes, spike times (in ms), used parameters and a sample segment 
% of the continuous data (optional) in filename_spikes.mat.
%input must be: 
%               A .txt file with the names of the files to use.
%               A matlab cell with the names of the files to use.
%               A vector, in this case the function will proccess all the
%                   supported files with that numbers in the folder
%                   (except .mat files). 
%                   (ipunt=2 don't implies 20 or viceversa)
%               'all', in this case the functions will process all the
%                   supported files in the folder (except .mat files).
%par_input must be a struct with some of the detecction parameters. All the
%parameters included will overwrite the parameters load from set_parameters()



if isnumeric(input) || any(strcmp(input,'all'))  %cases for numeric or 'all' input
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
                aux = regexp(f, '\d+', 'match');
                if ismember(str2num(aux{1}),input)
                    filenames = [filenames {fname}];   
                end
            end
        end
    end
    
elseif ischar(input) && length(input) > 4  %case for .txt input
    if  strcmp (input(end-3:end),'.txt')
        filenames =  textread(input,'%s');
    else
        filenames = {input};
    end
    
elseif iscellstr(input)   %case for cell input
    filenames = input;
else
    ME = MException('MyComponent:noValidInput', 'Invalid input arguments');
    throw(ME)
end

if exist('parallel','var') && parallel == true
    if exist('matlabpool','file')
        if matlabpool('size') > 0
            parallel = false;
        else
            matlabpool('open');
        end
    else
        poolobj = gcp('nocreate'); % If no pool, do not create new one.
        if isempty(poolobj)
            parallel = false;
        else
            parpool
        end
    end
end

if ~exist('par_input','var')
	par_input = struct;
end

parfor fnum = 1:length(filenames)
    filename = filenames{fnum};
    try
        get_spikes_single(filename, par_input);
    catch MExc
        warning(MExc.message);
        continue
    end
    
end


if exist('parallel','var') && parallel == true
    if exist('matlabpool','file')
        matlabpool('close')
    else
        poolobj = gcp('nocreate');
        delete(poolobj);
    end
 
end

end


function get_spikes_single(filename, par_input)
    
    par = set_parameters();
    par.filename = filename;
    par.reset_results = true;  %if true,  don't load times_ or _spikes files
    par.cont_segment = true;  %false doesn't save the segment of the continuous data in the spikes file
    data_handler = readInData(par);
    par = data_handler.par;
    par = update_parameters(par,par_input,'detect')
    return
    if data_handler.with_spikes            %data have some type of _spikes files
        [spikes, index] = data_handler.load_spikes(); 
        if ~data_handler.with_wc_spikes
            [spikes] = spike_alignment(spikes,par);
        end
    else    
        index = [];
        spikes = [];
        threshold = [];
        for n = 1:data_handler.max_segments
            x = data_handler.get_segment();
                %<----  Add here extra processing of the signal (x)
            [new_spikes, aux_th, new_index]  = amp_detect(x, par);
            index = [index data_handler.index2ts(new_index)]; %new_index to ms
            spikes = [spikes; new_spikes];
            threshold = [threshold, aux_th];
        end
    end

    current_par = par;
    par = struct;
    par = update_parameters(par, current_par, 'detect');
    
    
    %<----  Add here auxiliar parameters

    if current_par.cont_segment
        [psegment, sr_psegment] = data_handler.get_signal_sample();
        
        save([data_handler.nick_name '_spikes'], 'spikes', 'index', 'par','psegment','sr_psegment')
    else
        save([data_handler.nick_name '_spikes'], 'spikes', 'index', 'par')
    end

    if ~data_handler.with_spikes
        save([data_handler.nick_name '_spikes'],'threshold','-append')
    end
end
