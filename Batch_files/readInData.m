classdef readInData < handle
    %This calss will handle different types of data, and present a easy and
    %common interface.
	properties
        par             % parameters used, this will be updated with the data readed.
        nick_name       % name of the file without extension.
        with_raw        % boolean flag. True if the raw data was found and supported.
        with_spikes     % boolean flag. True if a file with spikes was readed.
        with_results    % boolean flag. True if 'times_' file was found.
        with_gui_status % boolean flag. True if the 'times_' includes GUI data.
        n_to_read       % number of segments of raw data not readed.
        sample_signal   % segment of continuous downsampled signal that wave_clus plot
        max_segments    % total number of segments of raw data to read.
        file_reader     % object member of the class that handle the file extension.
        with_wc_spikes  % boolean flag. True if '_spikes' file was found.
        with_psegment   % boolean flag. True if the continuous segment was found in the '_spikes' file.
        with_spc        % boolean flag. True if the '.dg_01' files created for the SPC algorithm were found.
    end 
	methods 
        function obj = readInData(par_ui)
            [unused, fnam, ext] = fileparts(par_ui.filename);
            ext = lower(ext);
            obj.par = par_ui;
            obj.nick_name = fnam;
            obj.n_to_read = 1;
            obj.with_results = false;
            obj.with_raw = false;
            obj.with_spikes = false;
            obj.with_wc_spikes = false;
            obj.with_psegment = false;
            obj.with_gui_status = false;
            obj.with_spc = false;
            with_par = false;
            results_selected = false;

            if length(fnam)>6 && strcmp(fnam(1:6),'times_') && strcmp(ext,'.mat') %if a 'times' file was selected.
                obj.with_results = true;
                results_selected = true;
                obj.nick_name = fnam(7:end);
                obj.with_wc_spikes = true;
            end
            if length(fnam)>7 && strcmp(fnam(end-6:end),'_spikes') && (strcmp(ext,'.mat') || isempty(ext))%if a 'spikes' file was selected.
                obj.with_wc_spikes = true;
                results_selected =true;
                obj.nick_name = fnam(1:end-7);
            end
            
            
            if (~isfield(par_ui,'reset_results')) || (~ par_ui.reset_results) ||  obj.with_results
                %Search for previous results
                if exist(['times_' obj.nick_name '.mat'],'file')
                    finfo = whos('-file',['times_' obj.nick_name '.mat']);
                    if ~ismember('spikes',{finfo.name})
                        ME = MException('MyComponent:FileError', 'Coultn''t find spikes variable in ''_times'' file');
                        throw(ME)
                    end
                   
                    if exist(['data_' obj.nick_name '.dg_01.lab'],'file') && exist(['data_' obj.nick_name '.dg_01'],'file')
                        obj.with_spc = true;
                        obj.with_results = true;
                    end
                    if  obj.with_results
                        if ismember('par',{finfo.name})
                            load(['times_' obj.nick_name '.mat'],'par');
                            obj.par = update_parameters(obj.par, par, 'relevant',true);
                            with_par = true;
                        end
                        if ismember('gui_status',{finfo.name})  
                            obj.with_gui_status = true;
                        end
                    end
                end
            end
            if (~isfield(par_ui,'reset_results')) || (~ par_ui.reset_results) || obj.with_wc_spikes

                %Search for previously detected spikes
                if exist([obj.nick_name '_spikes.mat'],'file')
                    obj.with_wc_spikes = true;
                    obj.with_spikes = true;
                    finfo = whos('-file', [obj.nick_name '_spikes.mat']);
                    if ~ismember('spikes',{finfo.name})
                        ME = MException('MyComponent:FileError', 'Coultn''t find spikes variable in ''_spikes'' file');
                        throw(ME)
                    end
                    if ismember('par',{finfo.name}) && ~ with_par 
                        load([obj.nick_name '_spikes.mat'],'par'); 
                        obj.par = update_parameters(obj.par,par,'detect',true);
                        with_par = true;
                    end
                    if ismember('psegment',{finfo.name})
                    	obj.with_psegment = true;
                    end
                    
                end
            end
            % Search raw data
            if exist([ext(2:end) '_wc_reader'],'file')
                obj.file_reader = eval([ext(2:end) '_wc_reader(obj.par,obj.par.filename)']);
                [sr, obj.max_segments, obj.with_raw, with_spikes] = obj.file_reader.get_info();
                obj.with_spikes = obj.with_spikes || with_spikes;
                if ~with_par                                                                                            %if didn't load sr from previous results 
                    if isempty(sr)
                        disp('Wave_clus didn''t find a sampling rate in file. It will use the input parameter or set_parameters.m')  %use default sr (from set_parameters) 
                    else
                        obj.par.sr = sr;                                                                                %load sr from raw data
                    end
                end
            else
                if ~(obj.with_results || obj.with_wc_spikes)
                    ME = MException('MyComponent:noSuchExt', 'File type ''%s'' isn''t supported',ext);
                    throw(ME)
                elseif results_selected
                    disp ('Wave_clus data selected. Raw data wasn''t loaded.')
                else
                    disp ('File type isn''t supported.')
                    disp ('Using Wave_clus data found.')
                end
            end 
            
            if ~isfield(obj.par,'channels')
                 obj.par.channels = 1;
            end

            if isfield(obj.par,'sr') && isfield(obj.par,'ref_ms')
                obj.par.ref = floor(obj.par.ref_ms *obj.par.sr/1000);
            end
            
        end
        
        % method for change the given par replacing with the parameters read from the data
        function par = update_par(obj,par) 
            load_par_names = fieldnames(obj.par);
            for i= 1:length(load_par_names)
                par.(load_par_names{i}) = obj.par.(load_par_names{i});
            end
        end
            
            
        function [spikes, index_ts] = load_spikes(obj)
            if ~ (obj.with_spikes || obj.with_wc_spikes)
                ME = MException('MyComponent:noSpikesFound', 'Wave_Clus couldn''t find a file with spikes');
                throw(ME)
            end
            
            if obj.with_wc_spikes                               %wc data have priority
                load([obj.nick_name '_spikes.mat']);
                if ~ exist('index_ts','var')                    %for retrocompatibility
                    index_ts = index;
                end
            else
                [spikes, index_ts] = obj.file_reader.load_spikes();
            end
            
        end
        
        
        function [clu, tree, spikes, index, inspk, ipermut, classes, forced, Temp] = load_results(obj)
        	
            if ~ obj.with_results
            	ME = MException('MyComponent:noClusFound', 'This file don''t have a associated ''times_%s.mat'' file',obj.nick_name);
            	throw(ME)
            end
            load(['times_' obj.nick_name '.mat']);
            if ~exist('ipermut','var')
            	ipermut = [];
            end
            if ~exist('forced','var')
            	forced = false(size(spikes,1), 1);
            end
           
            
            % cluster_class(:,1);
            index = cluster_class(1:end,2);
            classes = cluster_class(1:end,1);
            
            if ~exist('Temp','var')
            	Temp = ones(max(classes));
            end
            if obj.with_spc
                clu = load(['data_' obj.nick_name '.dg_01.lab']);
                tree = load(['data_' obj.nick_name '.dg_01']);
            else
                clu = [];
                tree = [];
            end
        end
        
        % method for load boolean vector for rejected spikes. Vector only used in develop mode.
        function [rejected] = load_rejected(obj)
        	
            if ~ obj.with_results
            	ME = MException('MyComponent:noClusFound', 'This file don''t have a associated ''times_%s.mat'' file',obj.nick_name);
            	throw(ME)
            end
            load(['times_' obj.nick_name '.mat']);
            if ~exist('rejected','var')
            	rejected = false(1,size(spikes,1));
            end
        end
        

        function [original_classes, current_temp,auto_sort_info] = get_gui_status(obj)
            load(['times_' obj.nick_name '.mat'],'gui_status','cluster_class');
            current_temp = gui_status.current_temp;
            if isempty(current_temp) || (current_temp == -1 && obj.with_spc)
                current_temp=1;
            end
            original_classes = gui_status.original_classes;
            auto_sort_info = [];
            if isfield(gui_status,'auto_sort_info')
                auto_sort_info = gui_status.auto_sort_info;
            end
        end
            
            
        function x = get_segment(obj)
            
            if ~ obj.with_raw
                ME = MException('MyComponent:noRawFound', 'Wave_Clus couldn''t find the raw data');
                throw(ME)
            end
            if obj.n_to_read > obj.max_segments
                ME = MException('MyComponent:allFileReaded', 'The raw data is already fully loaded');
                throw(ME)
            end
            
            x = obj.file_reader.get_segment(obj.n_to_read);
            
            if ~ isa(x,'double')
                x = double(x);
            end

            if ~obj.with_psegment && obj.n_to_read == 1 && obj.par.cont_segment
                lplot = min(floor(60*obj.par.sr), length(x));
                xf_detect = spike_detection_filter(x(1:lplot), obj.par);
                
                max_samples = obj.par.cont_plot_samples;
                sub = max(1,floor(lplot/max_samples));
                obj.sample_signal.xd_sub = xf_detect(1:sub:end) ;
                obj.sample_signal.sr_sub = obj.par.sr/sub ;
            end
            
            obj.n_to_read = obj.n_to_read + 1;
        end
        
        
        function index_ts = index2ts(obj,index)
            index_ts = obj.file_reader.index2ts(index,obj.n_to_read-1);
        end
        
        
        function [xd_sub, sr_sub] = get_signal_sample(obj)
            if obj.with_psegment
                load([obj.nick_name '_spikes.mat'],'psegment','sr_psegment');
                xd_sub = psegment;
                sr_sub = sr_psegment;
            else
                if obj.n_to_read == 1
                    disp('Segment read only for plotting');
                    obj.get_segment();
                end
                xd_sub = obj.sample_signal.xd_sub;
                sr_sub = obj.sample_signal.sr_sub;
                obj.sample_signal.xd_sub = [];
            end
        end

  
	end

end
