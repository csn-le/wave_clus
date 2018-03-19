function Get_spikes(input, varargin)

% PROGRAM Get_spikes.
% Detect spikes and save them in a file.
%
% Saves spikes, spike times (in ms), used parameters and a sample segment 
% of the continuous data (optional) in filename_spikes.mat.
% input must be: 
%               A .txt file with the names of the files to use.
%               A matlab cell with the names of the files to use.
%               A vector with channel numbers. In this case the function will proccess all the
%                   supported files (except .mat files) located in the folder with those
%                   channel numbers (e.g., CSC1.Ncs or NSX4.NC5)                    
%               'all', in this case the functions will process all the
%                   supported files in the folder (except .mat files).
% optional argument 'par' and the next input must be a struct with some of
%       the detecction parameters. All the parameters included in the structure 
%       par will overwrite the parameters loaded from set_parameters()
% optional argument 'parallel' : true for use parallel computing

% Example
% param.stdmin = 4;
% param.sr = 24000;
% param.detection = 'pos';
% Get_spikes([1:8 33:56],'parallel',true,'par',param);


%default config
par_input = struct;
parallel = false;

%search for optional inputs
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



% get a cell of filenames from the input
if isnumeric(input) || any(strcmp(input,'all'))  % cases for numeric or 'all' input
    filenames = {};
    se = supported_wc_extensions();
    dirnames = dir();
    dirnames = {dirnames.name};
    
    for i = 1:length(dirnames)
        fname = dirnames{i};
        [unused, f, ext] = fileparts(fname);
        ext = lower(ext(2:end));
        if any(strcmp(ext,se)) 
            if strcmp(input,'all')
                if strcmp(ext,'mat')
                    warning('Skipped file ''%s''. The ''.mat'' files should be added by name.\n',fname);
                    continue
                end
                filenames = [filenames {fname}];
            else
                aux = regexp(f, '\d+', 'match');
                if ~isempty(aux) && ismember(str2num(aux{1}),input)
                    if strcmp(ext,'mat')
                        warning('Skipped file ''%s''. The ''.mat'' files should be added by name.\n',fname);
                        continue
                    end
                    filenames = [filenames {fname}];   
                end
            end
        end
    end
    
elseif ischar(input) && length(input) > 4  
    if  strcmp (input(end-3:end),'.txt')   % case for .txt input
        filenames =  textread(input,'%s');
    else
        filenames = {input};               % case for cell input
    end
    
elseif iscellstr(input)   %case for cell input
    filenames = input;
else
    ME = MException('MyComponent:noValidInput', 'Invalid input arguments');
    throw(ME)
end

run_parfor = parallel;

% open parallel pool, if parallel input is true
if parallel
    if exist('matlabpool','file')   % old versions of Matlab have the matlabpool function
        try 
            matlabpool('open');     % If a matlabpool is already open 
        catch                       % will throw a error and catch here.
            parallel = false;
        end
    else
        poolobj = gcp('nocreate');
        if isempty(poolobj) % If already a pool, do not create new one.
            parpool
        else
            parallel = false;
        end
    end
end


init_date = now;


if run_parfor == true
    parfor fnum = 1:length(filenames)
        filename = filenames{fnum};
        get_spikes_single(filename, par_input);
        disp(sprintf('%d of %d ''spikes'' files done.',count_new_sp_files(init_date, filenames),length(filenames)))
    end
else
    for fnum = 1:length(filenames)
        filename = filenames{fnum};
        get_spikes_single(filename, par_input);
        disp(sprintf('%d of %d ''spikes'' files done.',count_new_sp_files(init_date, filenames),length(filenames)))
    end
end

% if a pool was open, close it
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
    par = update_parameters(par,par_input,'relevant');
    try 
        data_handler = readInData(par);
    catch MExc
        warning(MExc.message);
        return
    end
    par = data_handler.par;
    par = update_parameters(par,par_input,'detect');
    data_handler.par = par;
    threshold = [];
    if data_handler.with_spikes            %data have some type of _spikes files
        [spikes, index] = data_handler.load_spikes(); 
        if ~data_handler.with_wc_spikes
            spikes = spike_alignment(spikes,par);
            % spikes = int_spikes(spikes,par);
        else
            disp([filename ': Using detected spikes'])
        end
    else    
        index = [];
        spikes = [];
        
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
    par.detection_date =  datestr(now);
    
    %<----  Add here auxiliar parameters

    if current_par.cont_segment && data_handler.with_raw
        [psegment, sr_psegment] = data_handler.get_signal_sample();
        try
			save([data_handler.nick_name '_spikes'], 'spikes', 'index', 'par','psegment','sr_psegment','threshold')
		catch
			save([data_handler.nick_name '_spikes'], 'spikes', 'index', 'par','psegment','sr_psegment','threshold','-v7.3')
		end
    else
		try
			save([data_handler.nick_name '_spikes'], 'spikes', 'index', 'par')
		catch
			save([data_handler.nick_name '_spikes'], 'spikes', 'index', 'par','-v7.3')
		end
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

