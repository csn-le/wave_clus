classdef readInData < handle
	properties
        par
        file_type
        nick_name
        with_raw
        with_spikes
        with_clust
        n_to_read
        sample_signal
        max_segments
        %channel
        raw_reader
    end 
	methods 
        function obj = readInData(par)
            [~, fnam, ext] = fileparts(par.filename);
            obj.par = par;
            obj.nick_name = fnam;
            obj.n_to_read = 1;
            obj.par.sr =[];
            
            if strcmpi(ext,'.mat')
                %intento encontrar el raw
                %o si es el raw lo cargo como esta
                disp('searching raw data')
                obj.with_raw = true;            %maybe
                raw_filename = par.filename;    %perphaps
            else
               raw_filename = par.filename;
               obj.with_raw = true;
            end

            if obj.with_raw
                if exist([ext(2:end) '_reader'],'file')
                    obj.raw_reader = eval([ext(2:end) '_reader(par,raw_filename)']);
                    [sr, obj.max_segments] = obj.raw_reader.get_info();
                    if isempty(obj.par.sr)
                        obj.par.sr = sr;
                    end
                else
                    ME = MException('MyComponent:noSuchExt', 'File type ''%s'' is not supported',ext);
                    throw(ME)
                end
                
            end
            
            
            obj.par.ref = floor(obj.par.ref_ms *obj.par.sr/1000);
            
        end
        
        function [x t0]= get_segment(obj)
            
            if ~ obj.with_raw
                ME = MException('MyComponent:noRawFound', 'Wave_Clus couldn''t find the raw data',ext);
                throw(ME)
            end
            if obj.n_to_read > obj.max_segments
                ME = MException('MyComponent:allFileReaded', 'The raw data is already fully loaded',ext);
                throw(ME)
            end
            
            [x t0] = obj.raw_reader.get_segment(obj.n_to_read);
            
            if ~ isa(x,'double')
                x = double(x);
            end

            if obj.n_to_read == 1 && obj.par.show_signal
                lplot = min(floor(60*obj.par.sr), length(x));
                xf_detect = spike_detection_filter(x(1:lplot), obj.par);
                
                max_samples = 20000;
                sub = floor(lplot/max_samples);
                obj.sample_signal.xd_sub = xf_detect(1:sub:end) ;
                obj.sample_signal.sr_sub = obj.par.sr/sub ;
            end
            
            obj.n_to_read = obj.n_to_read + 1;
            
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




% 
% switch lower(ext)        
%     case 'CSC data (pre-clustered)'                                 %Neuralynx (CSC files)
%         channel = filename(4:end-4);
% 
%         f=fopen(filename,'r','l');
%         fseek(f,16384,'bof');                                       %Skip Header, put pointer to the first record
%         TimeStamps=fread(f,inf,'int64',(4+4+4+2*512));              %Read all TimeStamps
%         time0 = TimeStamps(1); 
%         timeend = TimeStamps(end);
%         sr = 512*1e6/(TimeStamps(2)-TimeStamps(1));
%         clear TimeStamps;
% 
%         handles.par.sr = sr;                                        % sampling rate (in Hz).
%         handles.par.ref = floor(handles.par.ref_ms *sr/1000);       % conversion to datapoints
% 
%                
%         %Load spikes and parameters
%         eval(['load times_CSC' channel ';']);
%         index=cluster_class(:,2)';
% 
%         %Load clustering results
%         fname = [handles.par.fname '_ch' channel];         %filename for interaction with SPC
%         clu=load([fname '.dg_01.lab']);
%         tree=load([fname '.dg_01']);
%         handles.par.fnamespc = fname;
%         handles.par.fnamesave = handles.par.fnamespc;
%                         
%         USER_DATA{3} = index;
%  
%         % LOAD CSC DATA (for plotting)
%         fseek(f,16384+8+4+4+4,'bof');                               %put pointer to the beginning of data
%         Samples=fread(f,ceil(sr*60),'512*int16=>int16',8+4+4+4);  
%         x=double(Samples(:))';
%         clear Samples; 
%         fclose(f);
% 
%         %GETS THE GAIN AND CONVERTS THE DATA TO MICRO V.
%         scale_factor = textread(['CSC' channel '.ncs'],'%s',43);
%         if(str2num(scale_factor{41})*1e6 > 0.5)
%             num_scale_factor=str2num(scale_factor{43});
%         else
%             num_scale_factor=str2num(scale_factor{41});
%         end
%         x=x*num_scale_factor*1e6;
%     
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