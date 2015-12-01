function Get_spikes(input, varargin)
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
% optional argument 'par' and the next input must be a struct with some of
%       the detecction parameters. All the parameters included will 
%       overwrite the parameters load from set_parameters()
% optional argument 'parallel' : true for use parallel computing



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

if parallel
    if exist('matlabpool','file')
        try
            matlabpool('open');
        catch
            parallel = false;
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


init_date = now;
parfor fnum = 1:length(filenames)
    filename = filenames{fnum};
    get_spikes_single(filename, par_input);
    disp(sprintf('%d of %d ''spikes'' files done.',count_new_sp_files(init_date, filenames),length(filenames)))

    
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


function get_spikes_single(filename, par_input)
    
    par = set_parameters();
    par.filename = filename;
    par.reset_results = true;  %if true,  don't load times_ or _spikes files
    par.cont_segment = true;  %false doesn't save the segment of the continuous data in the spikes file
    try 
        data_handler = readInData(par);
    catch MExc
        warning(MExc.message);
        return
    end
    
    par = data_handler.par;
    par = update_parameters(par,par_input,'detect');
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

    if current_par.cont_segment && data_handler.with_raw
        [psegment, sr_psegment] = data_handler.get_signal_sample();
        
        save([data_handler.nick_name '_spikes'], 'spikes', 'index', 'par','psegment','sr_psegment')
    else
        save([data_handler.nick_name '_spikes'], 'spikes', 'index', 'par')
    end

    if ~data_handler.with_spikes
        save([data_handler.nick_name '_spikes'],'threshold','-append')
    end
end

function counter = count_new_sp_files(initial_date, filenames)
counter = 0;
for i = 1:length(filenames)
    fname = filenames{i};
    [unu, fname] = fileparts(fname);
    FileInfo = dir([fname '_spikes.mat']);
    if length(FileInfo)==1 && (FileInfo.datenum > initial_date)
        counter = counter + 1;
    end
end
end

