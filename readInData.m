classdef readInData < handle
	properties
        par
        nick_name
        with_raw
        with_spikes
        with_results
        n_to_read
        sample_signal
        max_segments
        file_reader
        with_wc_spikes
        
    end 
	methods 
        function obj = readInData(par)
            [~, fnam, ext] = fileparts(par.filename);
            obj.par = par;
            obj.nick_name = fnam;
            obj.n_to_read = 1;
            obj.with_results = false;
            obj.with_spikes = false;
            obj.with_wc_spikes = false;
            
            if exist([ext(2:end) '_reader'],'file')
                obj.file_reader = eval([ext(2:end) '_reader(par,par.filename)']);
                [sr, obj.max_segments, obj.with_raw, obj.with_spikes] = obj.file_reader.get_info();
                if isempty(sr)
                    disp('Wave_clus din''t find a sampling rate in file. It will use the set in set_parameters.m')
                else
                    obj.par.sr = sr;
                end
            else
                ME = MException('MyComponent:noSuchExt', 'File type ''%s'' is not supported',ext);
                throw(ME)
            end
            
            %Search for previously detected spikes
            if exist([obj.nick_name '_spikes.mat'],'file')
%                 spf_names = whos('-file','[obj.nick_name '_times.mat'].mat')
%                 ismember('index_ts', {spf_names.name})
                obj.with_wc_spikes = true;
            end
            

            %Search for previous results
            if exist(['data_' obj.nick_name '.dg_01.lab'],'file') && exist(['data_' obj.nick_name '.dg_01'],'file') && exist(['times_' obj.nick_name '.mat'],'file')
                obj.with_results = true;
            end
            
            obj.par.ref = floor(obj.par.ref_ms *obj.par.sr/1000);
        end
        
        
        function [spikes, index_ts] = load_spikes(obj)
            if ~ (obj.with_spikes || obj.with_wc_spikes)
                ME = MException('MyComponent:noSpikesFound', 'Wave_Clus couldn''t find a file with spikes');
                throw(ME)
            end
            
            if obj.with_wc_spikes  %wc data have priority
                load([obj.nick_name '_times.mat']);
                if ~ exist('index_ts','var')  %for retrocompatibility
                    index_ts = index;
                end
            else
                [spikes, index_ts] = obj.file_reader.load_spikes();
            end
            
        end
        
        
        function [clu, tree, spikes, index, inspk, ipermut] = load_results(obj)
        	
            if ~ obj.with_results
            	ME = MException('MyComponent:noClusFound', 'This file don''t have a associated ''times_%s.mat'' file',obj.nick_name);
            	throw(ME)
            end
            load(['times_' obj.nick_name '.mat']);
            
            % cluster_class(:,1);
            index = cluster_class(:,2);
         	clu = load(['data_' obj.nick_name '.dg_01.lab']);
         	tree = load(['data_' obj.nick_name '.dg_01']);
        end

        
        function x = get_segment(obj)
            
            if ~ obj.with_raw
                ME = MException('MyComponent:noRawFound', 'Wave_Clus couldn''t find the raw data',ext);
                throw(ME)
            end
            if obj.n_to_read > obj.max_segments
                ME = MException('MyComponent:allFileReaded', 'The raw data is already fully loaded',ext);
                throw(ME)
            end
            
            x = obj.file_reader.get_segment(obj.n_to_read);
            
            if ~ isa(x,'double')
                x = double(x);
            end

            if obj.n_to_read == 1 && obj.par.show_signal
                lplot = min(floor(60*obj.par.sr), length(x));
                xf_detect = spike_detection_filter(x(1:lplot), obj.par);
                
                max_samples = 100000;
                sub = floor(lplot/max_samples);
                obj.sample_signal.xd_sub = xf_detect(1:sub:end) ;
                obj.sample_signal.sr_sub = obj.par.sr/sub ;
            end
            
            obj.n_to_read = obj.n_to_read + 1;
        end
        
        
        function index_ts = index2ts(obj,index)
            index_ts = obj.file_reader.index2ts(index,obj.n_to_read-1);
        end
        
        
        function [xd_sub, sr_sub] = get_signal_sample(obj)
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




%    case 'nev data (pre-clustered)'                                   %nev files matlab files
%         if length(filename) == 15
%             channel = filename(4);
%         else
%             channel = filename(4:5);
%         end
%         f=fopen(filename,'r','l');
% 
%         sr = 30000
%         handles.par.sr = sr;                        % sampling rate (in Hz).
%         handles.par.ref = floor(handles.par.ref_ms *sr/1000);     % conversion to datapoints
% 
%         %Load spikes and parameters
%         %Load spikes and parameters
%         eval(['load times_' num2str(channel) ';']);
%         index=cluster_class(:,2)';
% 
%         %Load clustering results
%         %fname = [handles.par.fname '_' filename(1:end-4)];         %filename for interaction with SPC
%         fname = [handles.par.fname '_ch' channel];         %filename for interaction with SPC
%         clu = load([fname '.dg_01.lab']);
%         tree = load([fname '.dg_01']);
%         handles.par.fnamespc = fname;
%         handles.par.fnamesave = fname;
%         
%         USER_DATA{3} = index;
% 
%         %Load continuous data (for ploting)
% %         if ~strcmp(filename(1:4),'poly')
% %             load(filename);
% %             if length(data)> 60*handles.par.sr
% %                 x=data(1:60*handles.par.sr)'; 
% %             else
% %                 x=data(1:length(data))'; 
% %             end
% %             [spikes,thr,index] = amp_detect_wc(x,handles);                   %Detection with amp. thresh.
% %         end     
%         
%     case 'NSX data (pre-clustered)'                                   %nev files matlab files
%         channel = filename(4:4+length(filename)-8);
%         f=fopen(filename,'r','l');
%         
%         %Load spikes and parameters
%         eval(['load times_NSX' num2str(channel) ';']);
%         index=cluster_class(:,2)';
%         
%         handles.par = par;      %Load parameters
% 
%         %Load clustering results
%         %fname = [handles.par.fname '_' filename(1:end-4)];         %filename for interaction with SPC
% %         fname = [handles.par.fname '_ch' channel];         %filename for interaction with SPC                                                                
% %         fname = [handles.par.fname channel];         %filename for interaction with SPC
%         fname = handles.par.fname;         %filename for interaction with SPC
%         
%         clu = load([fname '.dg_01.lab']);
%         tree = load([fname '.dg_01']);
%         handles.par.fnamespc = fname;
%         handles.par.fnamesave = fname;
%         USER_DATA{3} = index;
% 
%     case '.nse'
%         if length(filename) == 7
%             channel = filename(3);
%         else
%             channel = filename(3:4);
%         end
%         [index, Samples] = Nlx2MatSE(['Sc' num2str(channel) '.Nse'],1,0,0,0,1,0);
%         spikes(:,:)= Samples(:,1,:); clear Samples; spikes = spikes';
%         sr = 24000
%         handles.par.sr = sr;                        % sampling rate (in Hz).
%         handles.par.ref = floor(handles.par.ref_ms *sr/1000);     % conversion to datapoints
%         handles.nsegment = 1;
%         axes(handles.cont_data); cla
%         
%         [spikes] = spike_alignment(spikes,handles);
%         set(handles.file_name,'string','Calculating spike features ...');
%         [inspk] = wave_features_wc(spikes,handles);                 %Extract spike features.
% 
%         if handles.par.permut == 'y'
%             if handles.par.match == 'y';
%                 naux = min(handles.par.max_spk,size(inspk,1));
%                 ipermut = randperm(length(inspk));
%                 ipermut(naux+1:end) = [];
%                 inspk_aux = inspk(ipermut,:);
%             else
%                 ipermut = randperm(length(inspk));
%                 inspk_aux = inspk(ipermut,:);
%             end
%         else
%             if handles.par.match == 'y';
%                 naux = min(handles.par.max_spk,size(inspk,1));
%                 inspk_aux = inspk(1:naux,:);
%             else
%                 inspk_aux = inspk;
%             end
%         end
%             
%         %Interaction with SPC
%         set(handles.file_name,'string','Running SPC ...');
%         handles.par.fname_in = 'tmp_data';
%         fname_in = handles.par.fname_in;
%         save([fname_in],'inspk_aux','-ascii');                         %Input file for SPC
%         handles.par.fname = [handles.par.fname '_wc'];             %Output filename of SPC
%         handles.par.fnamesave = [handles.par.fname '_ch' ...
%                 num2str(channel)];                                 %filename if "save clusters" button is pressed
%         handles.par.fnamespc = handles.par.fname;
%         USER_DATA{3} = index/1000;
%         
%     case 'Sc data (pre-clustered)'
%         channel = filename(3:end-4);
%         [index, Samples] = Nlx2MatSE(['Sc' channel '.Nse'],1,0,0,0,1,0);
%         index = index/1000;
%         spikes(:,:)= Samples(:,1,:); clear Samples; spikes = spikes';
%         sr = 24000
%         handles.par.sr = sr;                        % sampling rate (in Hz).
%         handles.par.ref = floor(handles.par.ref_ms *sr/1000);     % conversion to datapoints
%         
%         axes(handles.cont_data); cla
% 
%         [spikes] = spike_alignment(spikes,handles);
%         
%         %Load clustering results
%         fname = [handles.par.fname '_ch' channel];         %filename for interaction with SPC
%         clu = load([fname '.dg_01.lab']);
%         tree = load([fname '.dg_01']);
%         handles.par.fnamespc = fname;
%         handles.par.fnamesave = handles.par.fnamespc; 
% 
%         USER_DATA{3} = index;  
%     
%     case '.mat'            % ASCII matlab files
%         sr = 20000
%         handles.par.sr = sr;                        % sampling rate (in Hz).
%         handles.nsegment = 1;
%         handles.par.ref = floor(handles.par.ref_ms *sr/1000);     % conversion to datapoints
%         index_all=[];
%         spikes_all=[];
%         for j=1:handles.par.segments                                %that's for cutting the data into pieces
%             % LOAD CONTINUOUS DATA
%             load(filename);
%             x=data(:)';
%             tsmin = (j-1)*floor(length(data)/handles.par.segments)+1;
%             tsmax = j*floor(length(data)/handles.par.segments);
%             x=data(tsmin:tsmax); clear data; 
%             handles.nsegment = j;                                      %flag for ploting only in the 1st loop
%             
%             % SPIKE DETECTION WITH AMPLITUDE THRESHOLDING
%             set(handles.file_name,'string','Detecting spikes ...');
%             [spikes,thr,index]  = amp_detect_wc(x,handles,true);        %detection with amp. thresh.
%             index=index+tsmin-1;
%             
%             index_all = [index_all index];
%             spikes_all = [spikes_all; spikes];
%         end
%         index = index_all *1e3/handles.par.sr;                     %spike times in ms.
%         spikes = spikes_all;
%         
%         USER_DATA{2}=spikes;
%         USER_DATA{3}=index;
%         set(handles.file_name,'string','Calculating spike features ...');
%         [inspk] = wave_features_wc(spikes,handles);                %Extract spike features.
% 
%         if handles.par.permut == 'y'
%             if handles.par.match == 'y';
%                 naux = min(handles.par.max_spk,size(inspk,1));
%                 ipermut = randperm(length(inspk));
%                 ipermut(naux+1:end) = [];
%                 inspk_aux = inspk(ipermut,:);
%             else
%                 ipermut = randperm(length(inspk));
%                 inspk_aux = inspk(ipermut,:);
%             end
%         else
%             if handles.par.match == 'y';
%                 naux = min(handles.par.max_spk,size(inspk,1));
%                 inspk_aux = inspk(1:naux,:);
%             else
%                 inspk_aux = inspk;
%             end
%         end
%             
%         %Interaction with SPC
%         set(handles.file_name,'string','Running SPC ...');
%         handles.par.fname_in = 'tmp_data';
%         fname_in = handles.par.fname_in;
%         save([fname_in],'inspk_aux','-ascii');                         %Input file for SPC
%         handles.par.fnamesave = [handles.par.fname '_' ...
%                 filename(1:end-4)];                                %filename if "save clusters" button is pressed
% %         handles.par.fname = [handles.par.fname '_wc'];             %Output filename of SPC
% %         handles.par.fnamespc = handles.par.fname;
%         handles.par.fnamespc = handles.par.fname;
%         handles.par.fname = [handles.par.fname '_wc'];             %Output filename of SPC
%         
%         
%         
%     case 'ASCII (pre-clustered)'                                   %ASCII matlab files
%         %In case of polytrode data 
%         if strcmp(filename(1:5),'times')
%             filename = filename(7:end);
%         end
%         sr = 24000
%         handles.par.sr = sr;                                       % sampling rate (in Hz).
%         handles.par.ref = floor(handles.par.ref_ms *sr/1000);     % conversion to datapoints
%                 
%         %Load spikes and parameters
%         eval(['load times_' filename ';']);
%         index=cluster_class(:,2)';
% 
%         %Load clustering results
%         fname = [handles.par.fname '_' filename(1:end-4)];         %filename for interaction with SPC
%         clu = load([fname '.dg_01.lab']);
%         tree = load([fname '.dg_01']);
%         handles.par.fnamespc = fname;
%         handles.par.fnamesave = fname;
%         
%         USER_DATA{3} = index;
% 
%         %Load continuous data (for ploting)
%         if ~strcmp(filename(1:4),'poly')
%             load(filename);
%             lplot = min(length(data),floor(60*handles.par.sr));
%             [spikes,thr,index] = amp_detect_wc(data(1:lplot), handles,false);                   %Detection with amp. thresh.
%         end
%         
%     case 'ASCII spikes'
%         sr = 30000
%         handles.par.sr = sr;                        % sampling rate (in Hz).
%         handles.par.ref = floor(handles.par.ref_ms *sr/1000);     % conversion to datapoints
%         axes(handles.cont_data); cla
%         
%         %Load spikes
%         load(filename);
%                 
%         [spikes] = spike_alignment(spikes,handles);
%         set(handles.file_name,'string','Calculating spike features ...');
%         [inspk] = wave_features_wc(spikes,handles);                      %Extract spike features.
%         
%         if handles.par.permut == 'y'
%             if handles.par.match == 'y';
%                 naux = min(handles.par.max_spk,size(inspk,1));
%                 ipermut = randperm(length(inspk));
%                 ipermut(naux+1:end) = [];
%                 inspk_aux = inspk(ipermut,:);
%             else
%                 ipermut = randperm(length(inspk));
%                 inspk_aux = inspk(ipermut,:);
%             end
%         else
%             if handles.par.match == 'y';
%                 naux = min(handles.par.max_spk,size(inspk,1));
%                 inspk_aux = inspk(1:naux,:);
%             else
%                 inspk_aux = inspk;
%             end
%         end
%             
%         %Interaction with SPC
%         set(handles.file_name,'string','Running SPC ...');
%         handles.par.fname_in = 'tmp_data';
%         fname_in = handles.par.fname_in;
%         save([fname_in],'inspk_aux','-ascii');                      %Input file for SPC
%         handles.par.fnamesave = [handles.par.fname '_' ...
%                 filename(1:end-4)];                             %filename if "save clusters" button is pressed
%         handles.par.fname = [handles.par.fname '_wc'];          %Output filename of SPC
%         handles.par.fnamespc = handles.par.fname;
% 
%         USER_DATA{3} = index(:)';        
%         
%     case 'ASCII spikes (pre-clustered)' 
%         %         with_ascci_wc_old = 1;
%         with_ascci_wc_old = 0;
%         
%         [filename, pathname] = uigetfile('*.mat','Select file');
%         if ~with_ascci_wc_old
%             % BEGIN NEW %
%             filename = filename(1:end-11);  % removes the "_spikes.mat" part of the filename
%             % END NEW %
%         end
%         set(handles.file_name,'string',['Loading:    ' pathname filename]);
%         cd(pathname);
%         if with_ascci_wc_old
%             % BEGIN OLD %
% 
%             sr = 24000
%             handles.par.sr = sr;                        % sampling rate (in Hz).
%             handles.par.ref = floor(handles.par.ref_ms *sr/1000);     % conversion to datapoints
%             axes(handles.cont_data); cla
%             % END OLD %
%         end
%         %Load spikes and parameters
%         eval(['load times_' filename ';']);
%         index=cluster_class(:,2)';
%         
%         if ~with_ascci_wc_old
%             % BEGIN NEW %
%             par.filename = [filename '.BLA'];   
%             handles.par = par;      %Load parameters
%             % END NEW %
%         end 
%         
%         %Load clustering results
%         if with_ascci_wc_old
%             % fname REPLACED %
%             fname = [handles.par.fname '_' filename(1:end-4)];               %filename for interaction with SPC
%         else
%             fname = handles.par.fname;         %filename for interaction with SPC
%         end
%         clu = load([fname '.dg_01.lab']);
%         tree = load([fname '.dg_01']);
%         handles.par.fnamespc = fname;
%         handles.par.fnamesave = fname;
% 
%         USER_DATA{3} = index(:)';
%               
% end