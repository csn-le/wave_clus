function varargout = wave_clus(varargin)
% WAVE_CLUS M-file for wave_clus.fig
%      WAVE_CLUS, by itself, creates a new WAVE_CLUS or raises the existing
%      singleton*.
%
%      H = WAVE_CLUS returns the handle to a new WAVE_CLUS or the handle to
%      the existing singleton*.
%
%      WAVE_CLUS('Property','Value',...) creates a new WAVE_CLUS using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to wave_clus_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      WAVE_CLUS('CALLBACK') and WAVE_CLUS('CALLBACK',hObject,...) call the
%      local function named CALLBACK in WAVE_CLUS.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help wave_clus

% Last Modified by GUIDE v2.5 14-Jan-2009 10:01:11

% JMG101208
% USER_DATA DEFINITIONS
% USER_DATA{1} = par;
% USER_DATA{2} = spikes;
% USER_DATA{3} = index;
% USER_DATA{4} = clu;
% USER_DATA{5} = tree;
% USER_DATA{7} = inspk;
% USER_DATA{6} = classes(:)'
% USER_DATA{8} = temp
% USER_DATA{9} = classes(:)', backup for non-forced classes
% USER_DATA{10} = clustering_results
% USER_DATA{11} = clustering_results_bk
% USER_DATA{12} = ipermut, indexes of the previously permuted spikes for clustering taking random number of points 
% USER_DATA{13} - USER_DATA{19}, for future changes
% USER_DATA{20} - USER_DATA{42}, fix clusters



% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @wave_clus_OpeningFcn, ...
                   'gui_OutputFcn',  @wave_clus_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
                   'gui_Callback',   []);
if nargin & isstr(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before wave_clus is made visible.
function wave_clus_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for wave_clus
handles.output = hObject;
handles.datatype ='CSC data (pre-clustered)';
set(handles.isi1_accept_button,'value',1);
set(handles.isi2_accept_button,'value',1);
set(handles.isi3_accept_button,'value',1);
set(handles.spike_shapes_button,'value',1);
set(handles.force_button,'value',0);
set(handles.plot_all_button,'value',1);
set(handles.plot_average_button,'value',0);
set(handles.fix1_button,'value',0);
set(handles.fix2_button,'value',0);
set(handles.fix3_button,'value',0);


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes wave_clus wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = wave_clus_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

clus_colors = [0 0 1; 1 0 0; 0 0.5 0; 0 0.75 0.75; 0.75 0 0.75; 0.75 0.75 0; 0.25 0.25 0.25];
set(0,'DefaultAxesColorOrder',clus_colors)


% --- Executes on button press in load_data_button.
function load_data_button_Callback(hObject, eventdata, handles)
set(handles.isi1_accept_button,'value',1);
set(handles.isi2_accept_button,'value',1);
set(handles.isi3_accept_button,'value',1);
set(handles.isi1_reject_button,'value',0);
set(handles.isi2_reject_button,'value',0);
set(handles.isi3_reject_button,'value',0);
set(handles.isi1_nbins,'string','Auto');
set(handles.isi1_bin_step,'string','Auto');
set(handles.isi2_nbins,'string','Auto');
set(handles.isi2_bin_step,'string','Auto');
set(handles.isi3_nbins,'string','Auto');
set(handles.isi3_bin_step,'string','Auto');
set(handles.isi0_nbins,'string','Auto');
set(handles.isi0_bin_step,'string','Auto');
set(handles.force_button,'value',0);
set(handles.force_button,'string','Force');
set(handles.fix1_button,'value',0);
set(handles.fix2_button,'value',0);
set(handles.fix3_button,'value',0);

switch char(handles.datatype)
    case 'Simulator'  
        [filename, pathname] = uigetfile('C_*.mat','Select file');
        set(handles.file_name,'string',['Loading:    ' pathname filename]);
        cd(pathname);
        %load(filename);                                 %Load data
        load([pathname filename]);                      %Load data
        x.data=data;
        x.sr=1000./samplingInterval;
        
        handles.par = set_parameters_simulation(x.sr,filename,handles);     % Load parameters
        set(handles.min_clus_edit,'string',num2str(handles.par.min_clus));

        [spikes,thr,index] = amp_detect_wc(x.data,handles);     % Detection with amp. thresh.
        [inspk] = wave_features_wc(spikes,handles);             % Extract spike features.
        
        %Interaction with SPC
        set(handles.file_name,'string','Running SPC ...');
        handles.par.fname_in = 'tmp_data';
        fname_in = handles.par.fname_in;
                         
        if handles.par.permut == 'y'
            if handles.par.match == 'y';
                naux = min(handles.par.max_spk,size(inspk,1));
                ipermut = randperm(length(inspk));
                ipermut(naux+1:end) = [];
                inspk_aux = inspk(ipermut,:);
            else
                ipermut = randperm(length(inspk));
                inspk_aux = inspk(ipermut,:);
            end
        else
            if handles.par.match == 'y';
                naux = min(handles.par.max_spk,size(inspk,1));
                inspk_aux = inspk(1:naux,:);
            else
                inspk_aux = inspk;
            end
        end
        
        save([fname_in],'inspk_aux','-ascii');                      %Input file for SPC
        handles.par.fname = [handles.par.fname '_wc'];          %Output filename of SPC
        handles.par.fnamespc = handles.par.fname;
        handles.par.fnamesave = handles.par.fnamespc;
        [clu,tree] = run_cluster(handles);
        USER_DATA = get(handles.wave_clus_figure,'userdata');
        
        if exist('ipermut')
            clu_aux = zeros(size(clu,1),length(index)) + 1000;
            for i=1:length(ipermut)
                clu_aux(:,ipermut(i)+2) = clu(:,i+2);
            end
            clu_aux(:,1:2) = clu(:,1:2);
            clu = clu_aux; clear clu_aux
            USER_DATA{12} = ipermut;
        end
        
        USER_DATA{2}=spikes;
        USER_DATA{3}=index;
        USER_DATA{4} = clu;
        USER_DATA{5} = tree;
        USER_DATA{7} = inspk;
        set(handles.wave_clus_figure,'userdata',USER_DATA)
        
    case 'CSC data'                                              %Neuralynx (CSC files)
        [filename, pathname] = uigetfile('*.Ncs','Select file');
        set(handles.file_name,'string',['Loading:    ' pathname filename]);
        cd(pathname);
        if length(filename) == 8
            channel = filename(4);
        else
            channel = filename(4:5);
        end
        f=fopen(filename,'r','l');
        fseek(f,16384,'bof');                                     %Skip Header, put pointer to the first record
        TimeStamps=fread(f,inf,'int64',(4+4+4+2*512));            %Read all TimeStamps
        fseek(f,16384+8+4+4+4,'bof');                             %put pointer to the beginning of data
        time0 = TimeStamps(1); 
        timeend = TimeStamps(end);
        delta_time=(TimeStamps(2)-TimeStamps(1));
        sr = 512*1e6/delta_time;
        handles.par = set_parameters_CSC(sr,filename,handles);     % Load parameters
        set(handles.min_clus_edit,'string',num2str(handles.par.min_clus));
        
        %Load continuous data 
        if strcmp(handles.par.tmax,'all')                          %Loads all data
            index_all=[];
            spikes_all=[];
            lts = length(TimeStamps);
            %Segments the data in par.segments pieces
            handles.par.segments = ceil((timeend - time0) / ...
                (handles.par.segments_length * 1e6 * 60));         %number of segments in which data is cutted
            segmentLength = floor (lts/handles.par.segments);
            tsmin = 1 : segmentLength :lts;
            tsmin = tsmin(1:handles.par.segments);
            tsmax = tsmin - 1;
            tsmax = tsmax (2:end);
            tsmax = [tsmax, lts];
            recmax=tsmax;
	        recmin=tsmin;
            tsmin = TimeStamps(int64(tsmin));
            tsmax = TimeStamps(int64(tsmax));

            for j=1:length(tsmin)
                
                Samples=fread(f,512*(recmax(j)-recmin(j)+1),'512*int16=>int16',8+4+4+4);
                x=double(Samples(:))';
                clear Samples;
               
                %GETS THE GAIN AND CONVERTS THE DATA TO MICRO V.
                eval(['scale_factor=textread(''CSC' num2str(channel) '.Ncs'',''%s'',41);']);
                x=x*str2num(scale_factor{41})*1e6;
                
                handles.flag = j;                                   %flag for plotting only in the 1st loop
                [spikes,thr,index]  = amp_detect_wc(x,handles);     %detection with amp. thresh.
                index = index*1e6/sr+tsmin(j);
                index_all = [index_all index];
                spikes_all = [spikes_all; spikes];
            end
            index = (index_all-time0)/1000;
            spikes = spikes_all;
            USER_DATA = get(handles.wave_clus_figure,'userdata');
            USER_DATA{2}=spikes;
            USER_DATA{3}=index;
            set(handles.wave_clus_figure,'userdata',USER_DATA);
        else                                                        %Loads a data segment
            tsmin = time0 + handles.par.tmin*1e6;                   %min time to read (in micro-sec)
            tsmax = time0 + handles.par.tmax*1e6;                   %max time to read (in micro-sec)
            index_tinitial = find(tsmin > TimeStamps);
            if isempty(index_tinitial) ==1;
                index_tinitial = 0;
            else
                index_tinitial = index_tinitial(end);
            end    
            index_tfinal = find(tsmax < TimeStamps);
            if isempty(index_tfinal) ==1;
                index_tfinal = timeend;
            else
                index_tfinal = index_tfinal(1);
            end    
            fseek(f,16384+8+4+4+4+index_tinitial,'bof');            %put pointer to the correct time
            Samples=fread(f,512*(index_tfinal-index_tinitial+1),'512*int16=>int16',8+4+4+4);
            x=double(Samples(:))';
            clear Samples;
            [spikes,thr,index] = amp_detect_wc(x,handles);          %Detection with amp. thresh.
        end
        fclose(f);
        
        [inspk] = wave_features_wc(spikes,handles);                 %Extract spike features.
        
        if handles.par.permut == 'y'
            if handles.par.match == 'y';
                naux = min(handles.par.max_spk,size(inspk,1));
                ipermut = randperm(length(inspk));
                ipermut(naux+1:end) = [];
                inspk_aux = inspk(ipermut,:);
            else
                ipermut = randperm(length(inspk));
                inspk_aux = inspk(ipermut,:);
            end
        else
            if handles.par.match == 'y';
                naux = min(handles.par.max_spk,size(inspk,1));
                inspk_aux = inspk(1:naux,:);
            else
                inspk_aux = inspk;
            end
        end
        
        %Interaction with SPC
        set(handles.file_name,'string','Running SPC ...');
        fname_in = handles.par.fname_in;
        save([fname_in],'inspk_aux','-ascii');                      %Input file for SPC
        handles.par.fnamesave = [handles.par.fname '_ch' ...
                num2str(channel)];                                  %filename if "save clusters" button is pressed
        handles.par.fnamespc = handles.par.fname;
        handles.par.fname = [handles.par.fname '_wc'];              %Output filename of SPC 
        [clu,tree] = run_cluster(handles);
        USER_DATA = get(handles.wave_clus_figure,'userdata');
        
        if exist('ipermut')
            clu_aux = zeros(size(clu,1),length(index)) + 1000;
            for i=1:length(ipermut)
                clu_aux(:,ipermut(i)+2) = clu(:,i+2);
            end
            clu_aux(:,1:2) = clu(:,1:2);
            clu = clu_aux; clear clu_aux
            USER_DATA{12} = ipermut;
        end
        
        USER_DATA{4} = clu;
        USER_DATA{5} = tree;
        USER_DATA{7} = inspk;
        set(handles.wave_clus_figure,'userdata',USER_DATA)
        
        
    case 'CSC data (pre-clustered)'                                 %Neuralynx (CSC files)
        [filename, pathname] = uigetfile('*.Ncs','Select file');
        set(handles.file_name,'string',['Loading:    ' pathname filename]);
        cd(pathname);
        if length(filename) == 8
            channel = filename(4);
        else
            channel = filename(4:5);
        end
        f=fopen(filename,'r','l');
        fseek(f,16384,'bof');                                       %Skip Header, put pointer to the first record
        TimeStamps=fread(f,inf,'int64',(4+4+4+2*512));              %Read all TimeStamps
        time0 = TimeStamps(1); 
        timeend = TimeStamps(end);
        sr = 512*1e6/(TimeStamps(2)-TimeStamps(1));
        clear TimeStamps;
        handles.par = set_parameters_CSC(sr,filename,handles);      %Load parameters
        set(handles.min_clus_edit,'string',num2str(handles.par.min_clus));
        
        %Load spikes and parameters
        eval(['load times_CSC' num2str(channel) ';']);
        index=cluster_class(:,2)';

        %Load clustering results
        fname = [handles.par.fname '_ch' num2str(channel)];         %filename for interaction with SPC
        clu=load([fname '.dg_01.lab']);
        tree=load([fname '.dg_01']);
        handles.par.fnamespc = fname;
        handles.par.fnamesave = handles.par.fnamespc;
        
        USER_DATA = get(handles.wave_clus_figure,'userdata');
        
        if exist('ipermut')
            clu_aux = zeros(size(clu,1),length(index)) + 1000;
            for i=1:length(ipermut)
                clu_aux(:,ipermut(i)+2) = clu(:,i+2);
            end
            clu_aux(:,1:2) = clu(:,1:2);
            clu = clu_aux; clear clu_aux
            USER_DATA{12} = ipermut;
        end
        
        USER_DATA{2} = spikes;
        USER_DATA{3} = index;
        USER_DATA{4} = clu;
        USER_DATA{5} = tree;
        if exist('inspk');
            USER_DATA{7} = inspk;
        end
        set(handles.wave_clus_figure,'userdata',USER_DATA)
 
        % LOAD CSC DATA (for plotting)
        fseek(f,16384+8+4+4+4,'bof');                               %put pointer to the beginning of data
        Samples=fread(f,ceil(sr*60),'512*int16=>int16',8+4+4+4);  
        x=double(Samples(:))';
        clear Samples; 
        fclose(f);

        %GETS THE GAIN AND CONVERTS THE DATA TO MICRO V.
        eval(['scale_factor=textread(''CSC' num2str(channel) '.Ncs'',''%s'',41);']);
        x=x*str2num(scale_factor{41})*1e6;

        [spikes,thr,index] = amp_detect_wc(x,handles);              %Detection with amp. thresh.
       
    case 'Sc data'
        [filename, pathname] = uigetfile('*.Nse','Select file');
        set(handles.file_name,'string',['Loading:    ' pathname filename]);
        cd(pathname);
        if length(filename) == 7
            channel = filename(3);
        else
            channel = filename(3:4);
        end
        eval(['[index, Samples] = Nlx2MatSE(''Sc' num2str(channel) '.Nse'',1,0,0,0,1,0);']);
        spikes(:,:)= Samples(:,1,:); clear Samples; spikes = spikes';
        handles.par = set_parameters_Sc(filename,handles);          %Load parameters
        set(handles.min_clus_edit,'string',num2str(handles.par.min_clus));
        axes(handles.cont_data); cla
        
        [spikes] = spike_alignment(spikes,handles);
        
        [inspk] = wave_features_wc(spikes,handles);                 %Extract spike features.

        if handles.par.permut == 'y'
            if handles.par.match == 'y';
                naux = min(handles.par.max_spk,size(inspk,1));
                ipermut = randperm(length(inspk));
                ipermut(naux+1:end) = [];
                inspk_aux = inspk(ipermut,:);
            else
                ipermut = randperm(length(inspk));
                inspk_aux = inspk(ipermut,:);
            end
        else
            if handles.par.match == 'y';
                naux = min(handles.par.max_spk,size(inspk,1));
                inspk_aux = inspk(1:naux,:);
            else
                inspk_aux = inspk;
            end
        end
            
        %Interaction with SPC
        set(handles.file_name,'string','Running SPC ...');
        handles.par.fname_in = 'tmp_data';
        fname_in = handles.par.fname_in;
        save([fname_in],'inspk_aux','-ascii');                         %Input file for SPC
        handles.par.fname = [handles.par.fname '_wc'];             %Output filename of SPC
        handles.par.fnamesave = [handles.par.fname '_ch' ...
                num2str(channel)];                                 %filename if "save clusters" button is pressed
        handles.par.fnamespc = handles.par.fname;
        [clu,tree] = run_cluster(handles);
        USER_DATA = get(handles.wave_clus_figure,'userdata');
        
        if exist('ipermut')
            clu_aux = zeros(size(clu,1),length(index)) + 1000;
            for i=1:length(ipermut)
                clu_aux(:,ipermut(i)+2) = clu(:,i+2);
            end
            clu_aux(:,1:2) = clu(:,1:2);
            clu = clu_aux; clear clu_aux
            USER_DATA{12} = ipermut;
        end
        
        USER_DATA{2} = spikes;
        USER_DATA{3} = index/1000;
        USER_DATA{4} = clu;
        USER_DATA{5} = tree;
        USER_DATA{7} = inspk;
        set(handles.wave_clus_figure,'userdata',USER_DATA)
        
    case 'Sc data (pre-clustered)'
        [filename, pathname] = uigetfile('*.Nse','Select file');
        set(handles.file_name,'string',['Loading:    ' pathname filename]);
        cd(pathname);
        if length(filename) == 7
            channel = filename(3);
        else
            channel = filename(3:4);
        end
        eval(['[index, Samples] = Nlx2MatSE(''Sc' num2str(channel) '.Nse'',1,0,0,0,1,0);']);
        spikes(:,:)= Samples(:,1,:); clear Samples; spikes = spikes';
        handles.par = set_parameters_Sc(filename,handles);          %Load parameters
        set(handles.min_clus_edit,'string',num2str(handles.par.min_clus));
        axes(handles.cont_data); cla

        [spikes] = spike_alignment(spikes,handles);
        
        %Load clustering results
        fname = [handles.par.fname '_ch' num2str(channel)];         %filename for interaction with SPC
        clu=load([fname '.dg_01.lab']);
        tree=load([fname '.dg_01']);
        handles.par.fnamespc = fname;
        handles.par.fnamesave = handles.par.fnamespc; 
        
        USER_DATA = get(handles.wave_clus_figure,'userdata');
        
        if exist('ipermut')
            clu_aux = zeros(size(clu,1),length(index)) + 1000;
            for i=1:length(ipermut)
                clu_aux(:,ipermut(i)+2) = clu(:,i+2);
            end
            clu_aux(:,1:2) = clu(:,1:2);
            clu = clu_aux; clear clu_aux
            USER_DATA{12} = ipermut;
        end
        
        USER_DATA{2} = spikes;
        USER_DATA{3} = index/1000;
        USER_DATA{4} = clu;
        USER_DATA{5} = tree;
        
        set(handles.wave_clus_figure,'userdata',USER_DATA)
           
    
    case 'ASCII'            % ASCII matlab files
        [filename, pathname] = uigetfile('*.mat','Select file');
        set(handles.file_name,'string',['Loading:    ' pathname filename]);
        cd(pathname);
        handles.par = set_parameters_ascii(filename,handles);       %Load parameters
        set(handles.min_clus_edit,'string',num2str(handles.par.min_clus));
        
        index_all=[];
        spikes_all=[];
        for j=1:handles.par.segments                                %that's for cutting the data into pieces
            % LOAD CONTINUOUS DATA
            load(filename);
            x=data(:)';
            tsmin = (j-1)*floor(length(data)/handles.par.segments)+1;
            tsmax = j*floor(length(data)/handles.par.segments);
            x=data(tsmin:tsmax); clear data; 
            handles.flag = 1;                                      %flag for ploting only in the 1st loop
            
            % SPIKE DETECTION WITH AMPLITUDE THRESHOLDING
            [spikes,thr,index]  = amp_detect_wc(x,handles);        %detection with amp. thresh.
            index=index+tsmin-1;
            
            index_all = [index_all index];
            spikes_all = [spikes_all; spikes];
        end
        index = index_all *1e3/handles.par.sr;                     %spike times in ms.
        spikes = spikes_all;
        
        USER_DATA = get(handles.wave_clus_figure,'userdata');
        USER_DATA{2}=spikes;
        USER_DATA{3}=index;
        set(handles.wave_clus_figure,'userdata',USER_DATA);

        [inspk] = wave_features_wc(spikes,handles);                %Extract spike features.

        if handles.par.permut == 'y'
            if handles.par.match == 'y';
                naux = min(handles.par.max_spk,size(inspk,1));
                ipermut = randperm(length(inspk));
                ipermut(naux+1:end) = [];
                inspk_aux = inspk(ipermut,:);
            else
                ipermut = randperm(length(inspk));
                inspk_aux = inspk(ipermut,:);
            end
        else
            if handles.par.match == 'y';
                naux = min(handles.par.max_spk,size(inspk,1));
                inspk_aux = inspk(1:naux,:);
            else
                inspk_aux = inspk;
            end
        end
            
        %Interaction with SPC
        set(handles.file_name,'string','Running SPC ...');
        handles.par.fname_in = 'tmp_data';
        fname_in = handles.par.fname_in;
        save([fname_in],'inspk_aux','-ascii');                         %Input file for SPC
        handles.par.fnamesave = [handles.par.fname '_' ...
                filename(1:end-4)];                                %filename if "save clusters" button is pressed
        handles.par.fname = [handles.par.fname '_wc'];             %Output filename of SPC
        handles.par.fnamespc = handles.par.fname;
        
        [clu,tree] = run_cluster(handles);
        USER_DATA = get(handles.wave_clus_figure,'userdata');
        
        if exist('ipermut')
            clu_aux = zeros(size(clu,1),length(index)) + 1000;
            for i=1:length(ipermut)
                clu_aux(:,ipermut(i)+2) = clu(:,i+2);
            end
            clu_aux(:,1:2) = clu(:,1:2);
            clu = clu_aux; clear clu_aux
            USER_DATA{12} = ipermut;
        end
        
        USER_DATA{4} = clu;
        USER_DATA{5} = tree;
        USER_DATA{7} = inspk;
        set(handles.wave_clus_figure,'userdata',USER_DATA)
        
        
    case 'ASCII (pre-clustered)'                                   %ASCII matlab files
        [filename, pathname] = uigetfile('*.mat','Select file');
        set(handles.file_name,'string',['Loading:    ' pathname filename]);
        cd(pathname);
        
        %In case of polytrode data 
        if strcmp(filename(1:5),'times')
            filename = filename(7:end);
            handles.par = set_parameters_pol(filename,handles);      %Load parameters
        else
            handles.par = set_parameters_ascii(filename,handles);      %Load parameters
        end
        
        set(handles.min_clus_edit,'string',num2str(handles.par.min_clus));
        
        %Load spikes and parameters
        eval(['load times_' filename ';']);
        index=cluster_class(:,2)';

        %Load clustering results
        fname = [handles.par.fname '_' filename(1:end-4)];         %filename for interaction with SPC
        clu=load([fname '.dg_01.lab']);
        tree=load([fname '.dg_01']);
        handles.par.fnamespc = fname;
        handles.par.fnamesave = fname;

        USER_DATA = get(handles.wave_clus_figure,'userdata');
        
        if exist('ipermut')
            clu_aux = zeros(size(clu,1),length(index)) + 1000;
            for i=1:length(ipermut)
                clu_aux(:,ipermut(i)+2) = clu(:,i+2);
            end
            clu_aux(:,1:2) = clu(:,1:2);
            clu = clu_aux; clear clu_aux
            USER_DATA{12} = ipermut;
        end
        
        USER_DATA{2} = spikes;
        USER_DATA{3} = index;
        USER_DATA{4} = clu;
        USER_DATA{5} = tree;
        if exist('inspk');
            USER_DATA{7} = inspk;
        end
        set(handles.wave_clus_figure,'userdata',USER_DATA)

        %Load continuous data (for ploting)
        if ~strcmp(filename(1:4),'poly')
            load(filename);
            if length(data)> 60*handles.par.sr
                x=data(1:60*handles.par.sr)'; 
            else
                x=data(1:length(data))'; 
            end
            [spikes,thr,index] = amp_detect_wc(x,handles);                   %Detection with amp. thresh.
        end
        
    case 'ASCII spikes'
        [filename, pathname] = uigetfile('*.mat','Select file');
        set(handles.file_name,'string',['Loading:    ' pathname filename]);
        cd(pathname);
        handles.par = set_parameters_ascii_spikes(filename,handles);     %Load parameters
        set(handles.min_clus_edit,'string',num2str(handles.par.min_clus));
        axes(handles.cont_data); cla
        
        %Load spikes
        load(filename);
                
        [spikes] = spike_alignment(spikes,handles);
        
        [inspk] = wave_features_wc(spikes,handles);                      %Extract spike features.
        
        if handles.par.permut == 'y'
            if handles.par.match == 'y';
                naux = min(handles.par.max_spk,size(inspk,1));
                ipermut = randperm(length(inspk));
                ipermut(naux+1:end) = [];
                inspk_aux = inspk(ipermut,:);
            else
                ipermut = randperm(length(inspk));
                inspk_aux = inspk(ipermut,:);
            end
        else
            if handles.par.match == 'y';
                naux = min(handles.par.max_spk,size(inspk,1));
                inspk_aux = inspk(1:naux,:);
            else
                inspk_aux = inspk;
            end
        end
            
        %Interaction with SPC
        set(handles.file_name,'string','Running SPC ...');
        handles.par.fname_in = 'tmp_data';
        fname_in = handles.par.fname_in;
        save([fname_in],'inspk_aux','-ascii');                      %Input file for SPC
        handles.par.fnamesave = [handles.par.fname '_' ...
                filename(1:end-4)];                             %filename if "save clusters" button is pressed
        handles.par.fname = [handles.par.fname '_wc'];          %Output filename of SPC
        handles.par.fnamespc = handles.par.fname;
        [clu,tree] = run_cluster(handles);
        USER_DATA = get(handles.wave_clus_figure,'userdata');
        
        if exist('ipermut')
            clu_aux = zeros(size(clu,1),length(index)) + 1000;
            for i=1:length(ipermut)
                clu_aux(:,ipermut(i)+2) = clu(:,i+2);
            end
            clu_aux(:,1:2) = clu(:,1:2);
            clu = clu_aux; clear clu_aux
            USER_DATA{12} = ipermut;
        end
         
        USER_DATA{2} = spikes;
        USER_DATA{3} = index(:)';
        USER_DATA{4} = clu;
        USER_DATA{5} = tree;
        USER_DATA{7} = inspk;
        set(handles.wave_clus_figure,'userdata',USER_DATA)
        
        
    case 'ASCII spikes (pre-clustered)'
        [filename, pathname] = uigetfile('*.mat','Select file');
        set(handles.file_name,'string',['Loading:    ' pathname filename]);
        cd(pathname);
        handles.par = set_parameters_ascii_spikes(filename,handles);     %Load parameters
        set(handles.min_clus_edit,'string',num2str(handles.par.min_clus));
        axes(handles.cont_data); cla

        %Load spikes and parameters
        eval(['load times_' filename ';']);
        index=cluster_class(:,2)';

        %Load clustering results
        fname = [handles.par.fname '_' filename(1:end-4)];               %filename for interaction with SPC
        clu=load([fname '.dg_01.lab']);
        tree=load([fname '.dg_01']);
        handles.par.fnamespc = fname;
        handles.par.fnamesave = fname;
      
        USER_DATA = get(handles.wave_clus_figure,'userdata');
        
        if exist('ipermut')
            clu_aux = zeros(size(clu,1),length(index)) + 1000;
            for i=1:length(ipermut)
                clu_aux(:,ipermut(i)+2) = clu(:,i+2);
            end
            clu_aux(:,1:2) = clu(:,1:2);
            clu = clu_aux; clear clu_aux
            USER_DATA{12} = ipermut;
        end
        
        USER_DATA{2} = spikes;
        USER_DATA{3} = index(:)';
        USER_DATA{4} = clu;
        USER_DATA{5} = tree;
        USER_DATA{7} = inspk;
        set(handles.wave_clus_figure,'userdata',USER_DATA)
        
end    

temp=find_temp(tree,handles);                                   %Selects temperature.
set(handles.file_name,'string',[pathname filename]);

if size(clu,2)-2 < size(spikes,1);
    classes = clu(temp(end),3:end)+1;
    if ~exist('ipermut')
        classes = [classes(:)' zeros(1,size(spikes,1)-handles.par.max_spk)];
    end
else
    classes = clu(temp(end),3:end)+1;
end

guidata(hObject, handles);
USER_DATA = get(handles.wave_clus_figure,'userdata');
USER_DATA{6} = classes(:)';
USER_DATA{8} = temp(end);
USER_DATA{9} = classes(:)';                                     %backup for non-forced classes.

% definition of clustering_results
clustering_results = []; 
clustering_results(:,1) = repmat(temp,length(classes),1); % GUI temperatures
clustering_results(:,2) = classes'; % GUI classes 
clustering_results(:,3) = repmat(temp,length(classes),1); % original temperatures 
clustering_results(:,4) = classes'; % original classes 
clustering_results(:,5) = repmat(handles.par.min_clus,length(classes),1); % minimum number of clusters
clustering_results_bk = clustering_results; % old clusters for undo actions
USER_DATA{10} = clustering_results;
USER_DATA{11} = clustering_results_bk;
handles.force = 0;
handles.merge = 0;
handles.reject = 0;
handles.undo = 0;
handles.minclus = handles.par.min_clus;
set(handles.wave_clus_figure,'userdata',USER_DATA);
handles.setclus = 0;

% mark clusters when new data is loaded
plot_spikes(handles);

USER_DATA = get(handles.wave_clus_figure,'userdata');
clustering_results = USER_DATA{10};
mark_clusters_temperature_diagram(handles,tree,clustering_results);

% --- Executes on button press in change_temperature_button.
function change_temperature_button_Callback(hObject, eventdata, handles)
axes(handles.temperature_plot)
hold off
[temp aux]= ginput(1);                                          %gets the mouse input
temp = round((temp-handles.par.mintemp)/handles.par.tempstep);
if temp < 1; temp=1;end                                         %temp should be within the limits
if temp > handles.par.num_temp; temp=handles.par.num_temp; end
min_clus = round(aux);
set(handles.min_clus_edit,'string',num2str(min_clus));

USER_DATA = get(handles.wave_clus_figure,'userdata');
par = USER_DATA{1};
par.min_clus = min_clus;
clu = USER_DATA{4};
classes = clu(temp,3:end)+1;
tree = USER_DATA{5};
USER_DATA{1} = par;
USER_DATA{6} = classes(:)';
USER_DATA{8} = temp;
USER_DATA{9} = classes(:)';                                     %backup for non-forced classes.


handles.minclus = min_clus;
clustering_results = USER_DATA{10};
clustering_results_bk = USER_DATA{11};
set(handles.wave_clus_figure,'userdata',USER_DATA);
temperature=handles.par.mintemp+temp*handles.par.tempstep;

switch par.temp_plot
    case 'lin'
        plot([handles.par.mintemp handles.par.maxtemp-handles.par.tempstep],[par.min_clus par.min_clus],'k:',...
            handles.par.mintemp+(1:handles.par.num_temp)*handles.par.tempstep, ...
            tree(1:handles.par.num_temp,5:size(tree,2)),[temperature temperature],[1 tree(1,5)],'k:')
    case 'log'
         semilogy([handles.par.mintemp handles.par.maxtemp-handles.par.tempstep], ...
            [par.min_clus par.min_clus],'k:',...
            handles.par.mintemp+(1:handles.par.num_temp)*handles.par.tempstep, ...
            tree(1:handles.par.num_temp,5:size(tree,2)),[temperature temperature],[1 tree(1,5)],'k:')
end
xlim([0 handles.par.maxtemp])
xlabel('Temperature'); 
if par.temp_plot == 'log' 
    set(get(gca,'ylabel'),'vertical','Cap');
else
    set(get(gca,'ylabel'),'vertical','Baseline');
end
ylabel('Clusters size');

handles.setclus = 0;
handles.force = 0;
handles.merge = 0;
handles.reject = 0;
handles.undo = 0;
plot_spikes(handles);
USER_DATA = get(handles.wave_clus_figure,'userdata');
clustering_results = USER_DATA{10};
clustering_results_bk = USER_DATA{11};
mark_clusters_temperature_diagram(handles,tree,clustering_results)
set(handles.wave_clus_figure,'userdata',USER_DATA);

set(handles.fix1_button,'value',0);
set(handles.fix2_button,'value',0);
set(handles.fix3_button,'value',0);
for i=4:par.max_clus
    eval(['par.fix' num2str(i) '=0;']);
end    
USER_DATA{1} = par;
set(handles.wave_clus_figure,'userdata',USER_DATA);


% --- Change min_clus_edit     
function min_clus_edit_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_figure,'userdata');
par = USER_DATA{1};
par.min_clus = str2num(get(hObject, 'String'));
clu = USER_DATA{4};
temp = USER_DATA{8};
classes = clu(temp,3:end)+1;
tree = USER_DATA{5};
USER_DATA{1} = par;
USER_DATA{6} = classes(:)';
USER_DATA{9} = classes(:)';                                     %backup for non-forced classes.
clustering_results = USER_DATA{10};
clustering_results(:,5) = par.min_clus;
set(handles.wave_clus_figure,'userdata',USER_DATA);

mark_clusters_temperature_diagram(handles,tree,clustering_results)
handles.setclus = 0;
handles.force = 0;
handles.merge = 0;
handles.undo = 0;
handles.reject = 0;
handles.minclus = par.min_clus;
plot_spikes(handles);
USER_DATA = get(handles.wave_clus_figure,'userdata');
clustering_results = USER_DATA{10};
set(handles.wave_clus_figure,'userdata',USER_DATA);
mark_clusters_temperature_diagram(handles,tree,clustering_results)

set(handles.force_button,'value',0);
set(handles.force_button,'string','Force');
set(handles.fix1_button,'value',0);
set(handles.fix2_button,'value',0);
set(handles.fix3_button,'value',0);
for i=4:par.max_clus
    eval(['par.fix' num2str(i) '=0;']);
end    


% --- Executes on button press in save_clusters_button.
function save_clusters_button_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_figure,'userdata');
spikes = USER_DATA{2};
par = USER_DATA{1};
classes = USER_DATA{6};

% Classes should be consecutive numbers
i=1;
while i<=max(classes)
    if isempty(classes(find(classes==i)))
        for k=i+1:max(classes)
            classes(find(classes==k))=k-1;
        end
    else
        i=i+1;
    end
end

%Saves clusters
cluster_class=zeros(size(spikes,1),2);
cluster_class(:,1) = classes(:);
cluster_class(:,2) = USER_DATA{3}';

outfile=['times_' par.filename(1:end-4)];

currentver = version;
currentver = currentver(1);
switch currentver
    case {'7'}
        if isempty(USER_DATA{7}) && isempty(USER_DATA{12})
            exec_line = strcat('save',' ''',outfile,'''',' cluster_class',' par',' spikes',' -v6;');
        elseif isempty(USER_DATA{7}) && ~isempty(USER_DATA{12})
            ipermut = USER_DATA{12};
            exec_line = strcat('save',' ''',outfile,'''',' cluster_class',' par',' spikes',' ipermut',' -v6;');
        elseif ~isempty(USER_DATA{7}) && isempty(USER_DATA{12})
            inspk = USER_DATA{7};
            exec_line = strcat('save',' ''',outfile,'''',' cluster_class',' par',' spikes',' inspk',' -v6;');
        else
            inspk = USER_DATA{7};
            ipermut = USER_DATA{12};
            exec_line = strcat('save',' ''',outfile,'''',' cluster_class',' par',' spikes',' inspk',' ipermut',' -v6;');
        end
        
    otherwise
        if isempty(USER_DATA{7}) && isempty(USER_DATA{12})
            exec_line = strcat('save',' ''',outfile,'''',' cluster_class',' par',' spikes');
        elseif isempty(USER_DATA{7}) && ~isempty(USER_DATA{12})
            ipermut = USER_DATA{12};
            exec_line = strcat('save',' ''',outfile,'''',' cluster_class',' par',' spikes',' ipermut');
        elseif ~isempty(USER_DATA{7}) && isempty(USER_DATA{12})
            inspk = USER_DATA{7};
            exec_line = strcat('save',' ''',outfile,'''',' cluster_class',' par',' spikes',' inspk');
        else
            inspk = USER_DATA{7};
            ipermut = USER_DATA{12};
            exec_line = strcat('save',' ''',outfile,'''',' cluster_class',' par',' spikes',' inspk',' ipermut');
        end
end

eval(exec_line);

if(~strcmp(handles.par.fnamespc,handles.par.fnamesave))
    copyfile([handles.par.fnamespc '.dg_01.lab'], [handles.par.fnamesave '.dg_01.lab']);
    copyfile([handles.par.fnamespc '.dg_01'], [handles.par.fnamesave '.dg_01']);
end

%Save figures
h_figs=get(0,'children');
h_fig = findobj(h_figs,'tag','wave_clus_figure');
h_fig1 = findobj(h_figs,'tag','wave_clus_aux');
h_fig2= findobj(h_figs,'tag','wave_clus_aux1');
h_fig3= findobj(h_figs,'tag','wave_clus_aux2');
h_fig4= findobj(h_figs,'tag','wave_clus_aux3');
h_fig5= findobj(h_figs,'tag','wave_clus_aux4');
h_fig6= findobj(h_figs,'tag','wave_clus_aux5');
if strcmp(outfile(7:9),'CSC')
    if ~isempty(h_fig)
        figure(h_fig); set(gcf,'papertype','usletter','paperorientation','portrait','paperunits','inches')
        set(gcf,'paperposition',[.25 .25 10.5 7.8])
        eval(['print(h_fig,''-djpeg'',''fig2print_' outfile(10:end) ''')' ]);
    end
    if ~isempty(h_fig1)
        figure(h_fig1); set(gcf,'papertype','usletter','paperorientation','portrait','paperunits','inches')
        set(gcf,'paperposition',[.25 .25 10.5 7.8])
        eval(['print(h_fig1,''-djpeg'',''fig2print_' outfile(10:end) 'a' ''')' ]);
    end
    if ~isempty(h_fig2)
        figure(h_fig2); set(gcf,'papertype','usletter','paperorientation','portrait','paperunits','inches')
        set(gcf,'paperposition',[.25 .25 10.5 7.8])
        eval(['print(h_fig2,''-djpeg'',''fig2print_' outfile(10:end) 'b' ''')' ]);
    end
    if ~isempty(h_fig3)
        figure(h_fig3); set(gcf,'papertype','usletter','paperorientation','portrait','paperunits','inches')
        set(gcf,'paperposition',[.25 .25 10.5 7.8])
        eval(['print(h_fig3,''-djpeg'',''fig2print_' outfile(10:end) 'c' ''')' ]);
    end
    if ~isempty(h_fig4)
        figure(h_fig4); set(gcf,'papertype','usletter','paperorientation','portrait','paperunits','inches')
        set(gcf,'paperposition',[.25 .25 10.5 7.8])
        eval(['print(h_fig4,''-djpeg'',''fig2print_' outfile(10:end) 'd' ''')' ]);
    end
    if ~isempty(h_fig5)
        figure(h_fig5); set(gcf,'papertype','usletter','paperorientation','portrait','paperunits','inches')
        set(gcf,'paperposition',[.25 .25 10.5 7.8])
        eval(['print(h_fig5,''-djpeg'',''fig2print_' outfile(10:end) 'e' ''')' ]);
    end
    if ~isempty(h_fig6)
        figure(h_fig6); set(gcf,'papertype','usletter','paperorientation','portrait','paperunits','inches')
        set(gcf,'paperposition',[.25 .25 10.5 7.8])
        eval(['print(h_fig6,''-djpeg'',''fig2print_' outfile(10:end) 'f' ''')' ]);
    end
    
elseif strcmp(outfile(7:9),'pol')
    if ~isempty(h_fig)
        figure(h_fig); set(gcf,'papertype','usletter','paperorientation','portrait','paperunits','inches')
        set(gcf,'paperposition',[.25 .25 10.5 7.8])
        eval(['print(h_fig,''-djpeg'',''fig2print_' outfile(7:end) ''')' ]);
    end
    if ~isempty(h_fig1)
        figure(h_fig1); set(gcf,'papertype','usletter','paperorientation','portrait','paperunits','inches')
        set(gcf,'paperposition',[.25 .25 10.5 7.8])
        eval(['print(h_fig1,''-djpeg'',''fig2print_' outfile(7:end) 'a' ''')' ]);
    end
    if ~isempty(h_fig2)
        figure(h_fig2); set(gcf,'papertype','usletter','paperorientation','portrait','paperunits','inches')
        set(gcf,'paperposition',[.25 .25 10.5 7.8])
        eval(['print(h_fig2,''-djpeg'',''fig2print_' outfile(7:end) 'b' ''')' ]);
    end
    if ~isempty(h_fig3)
        figure(h_fig3); set(gcf,'papertype','usletter','paperorientation','portrait','paperunits','inches')
        set(gcf,'paperposition',[.25 .25 10.5 7.8])
        eval(['print(h_fig3,''-djpeg'',''fig2print_' outfile(7:end) 'c' ''')' ]);
    end
    if ~isempty(h_fig4)
        figure(h_fig4); set(gcf,'papertype','usletter','paperorientation','portrait','paperunits','inches')
        set(gcf,'paperposition',[.25 .25 10.5 7.8])
        eval(['print(h_fig4,''-djpeg'',''fig2print_' outfile(7:end) 'd' ''')' ]);
    end
    if ~isempty(h_fig5)
        figure(h_fig5); set(gcf,'papertype','usletter','paperorientation','portrait','paperunits','inches')
        set(gcf,'paperposition',[.25 .25 10.5 7.8])
        eval(['print(h_fig5,''-djpeg'',''fig2print_' outfile(7:end) 'e' ''')' ]);
    end
    if ~isempty(h_fig6)
        figure(h_fig6); set(gcf,'papertype','usletter','paperorientation','portrait','paperunits','inches')
        set(gcf,'paperposition',[.25 .25 10.5 7.8])
        eval(['print(h_fig6,''-djpeg'',''fig2print_' outfile(7:end) 'f' ''')' ]);
    end
    
else
    if ~isempty(h_fig)
        figure(h_fig); set(gcf,'papertype','usletter','paperorientation','portrait','paperunits','inches')
        set(gcf,'paperposition',[.25 .25 10.5 7.8])
        eval(['print(h_fig,''-djpeg'',''fig2print_' outfile(7:end)  ''')' ]);
    end
    if ~isempty(h_fig1)
        figure(h_fig1); set(gcf,'papertype','usletter','paperorientation','portrait','paperunits','inches')
        set(gcf,'paperposition',[.25 .25 10.5 7.8])
        eval(['print(h_fig1,''-djpeg'',''fig2print_' outfile(7:end) 'a' ''')' ]);
    end
    if ~isempty(h_fig2)
        figure(h_fig2); set(gcf,'papertype','usletter','paperorientation','portrait','paperunits','inches')
        set(gcf,'paperposition',[.25 .25 10.5 7.8])
        eval(['print(h_fig2,''-djpeg'',''fig2print_' outfile(7:end) 'b' ''')' ]);
    end
    if ~isempty(h_fig3)
        figure(h_fig3); set(gcf,'papertype','usletter','paperorientation','portrait','paperunits','inches')
        set(gcf,'paperposition',[.25 .25 10.5 7.8])
        eval(['print(h_fig3,''-djpeg'',''fig2print_' outfile(7:end) 'c' ''')' ]);
    end
    if ~isempty(h_fig4)
        figure(h_fig4); set(gcf,'papertype','usletter','paperorientation','portrait','paperunits','inches')
        set(gcf,'paperposition',[.25 .25 10.5 7.8])
        eval(['print(h_fig4,''-djpeg'',''fig2print_' outfile(7:end) 'd' ''')' ]);
    end
    if ~isempty(h_fig5)
        figure(h_fig5); set(gcf,'papertype','usletter','paperorientation','portrait','paperunits','inches')
        set(gcf,'paperposition',[.25 .25 10.5 7.8])
        eval(['print(h_fig5,''-djpeg'',''fig2print_' outfile(7:end) 'e' ''')' ]);
    end
    if ~isempty(h_fig6)
        figure(h_fig6); set(gcf,'papertype','usletter','paperorientation','portrait','paperunits','inches')
        set(gcf,'paperposition',[.25 .25 10.5 7.8])
        eval(['print(h_fig6,''-djpeg'',''fig2print_' outfile(7:end) 'f' ''')' ]);
    end
end

% --- Executes on selection change in data_type_popupmenu.
function data_type_popupmenu_Callback(hObject, eventdata, handles)
aux = get(hObject, 'String');
aux1 = get(hObject, 'Value');
handles.datatype = aux(aux1);
guidata(hObject, handles);


% --- Executes on button press in set_parameters_button.
function set_parameters_button_Callback(hObject, eventdata, handles)
helpdlg('Check the set_parameters files in the subdirectory Wave_clus\Parameters_files');


%SETTING OF FORCE MEMBERSHIP
% --------------------------------------------------------------------
function force_button_Callback(hObject, eventdata, handles)
%set(gcbo,'value',1);
USER_DATA = get(handles.wave_clus_figure,'userdata');
par = USER_DATA{1};
spikes = USER_DATA{2};
classes = USER_DATA{6};
inspk = USER_DATA{7};

% Fixed clusters are not considered for forcing
if get(handles.fix1_button,'value') ==1     
    fix_class = USER_DATA{20}';
    classes(fix_class)=-1;
end
if get(handles.fix2_button,'value') ==1     
    fix_class = USER_DATA{21}';
    classes(fix_class)=-1;
end
if get(handles.fix3_button,'value') ==1     
    fix_class = USER_DATA{22}';
    classes(fix_class)=-1;
end
% Get fixed clusters from aux figures
for i=4:par.max_clus
    eval(['fixx = par.fix' num2str(i) ';']);
    if fixx == 1
        fix_class = USER_DATA{22+i-3}';
        classes(fix_class)=-1;
    end
end

switch par.force_feature
    case 'spk'
        f_in  = spikes(find(classes~=0 & classes~=-1),:);
        f_out = spikes(find(classes==0),:);
    case 'wav'
        if isempty(inspk)
            [inspk] = wave_features_wc(spikes,handles);        % Extract spike features.
            USER_DATA{7} = inspk;
        end
        f_in  = inspk(find(classes~=0 & classes~=-1),:);
        f_out = inspk(find(classes==0),:);
end

class_in = classes(find(classes~=0 & classes~=-1));
class_out = force_membership_wc(f_in, class_in, f_out, handles);
classes(find(classes==0)) = class_out;
  
USER_DATA{6} = classes(:)';
set(handles.wave_clus_figure,'userdata',USER_DATA)

clustering_results = USER_DATA{10};
handles.minclus = clustering_results(1,5);
handles.setclus = 1;
handles.force = 1;
handles.merge = 0;
handles.reject = 0;
handles.undo = 0;

plot_spikes(handles);

USER_DATA = get(handles.wave_clus_figure,'userdata');
clustering_results = USER_DATA{10};
set(handles.wave_clus_figure,'userdata',USER_DATA);

set(handles.fix1_button,'value',0);
set(handles.fix2_button,'value',0);
set(handles.fix3_button,'value',0);
for i=4:par.max_clus
    eval(['par.fix' num2str(i) '=0;']);
end    



% PLOT ALL PROJECTIONS BUTTON
% --------------------------------------------------------------------
function Plot_all_projections_button_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_figure,'userdata');
par = USER_DATA{1};
if strcmp(par.filename(1:4),'poly')
    Plot_amplitudes(handles)
else
    Plot_all_features(handles)
end
% --------------------------------------------------------------------

% fix1 button --------------------------------------------------------------------
function fix1_button_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_figure,'userdata');
classes = USER_DATA{6};
fix_class = find(classes==1);
if get(handles.fix1_button,'value') ==1
    USER_DATA{20} = fix_class;
else
    USER_DATA{20} = [];
end
set(handles.wave_clus_figure,'userdata',USER_DATA)
h_figs=get(0,'children');
h_fig4 = findobj(h_figs,'tag','wave_clus_aux');
h_fig3 = findobj(h_figs,'tag','wave_clus_aux1');
h_fig2 = findobj(h_figs,'tag','wave_clus_aux2');
h_fig1 = findobj(h_figs,'tag','wave_clus_aux3');
set(h_fig4,'userdata',USER_DATA)
set(h_fig3,'userdata',USER_DATA)
set(h_fig2,'userdata',USER_DATA)
set(h_fig1,'userdata',USER_DATA)


% fix2 button --------------------------------------------------------------------
function fix2_button_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_figure,'userdata');
classes = USER_DATA{6};
fix_class = find(classes==2);
if get(handles.fix2_button,'value') ==1
    USER_DATA{21} = fix_class;
else
    USER_DATA{21} = [];
end
set(handles.wave_clus_figure,'userdata',USER_DATA)
h_figs=get(0,'children');
h_fig4 = findobj(h_figs,'tag','wave_clus_aux');
h_fig3 = findobj(h_figs,'tag','wave_clus_aux1');
h_fig2 = findobj(h_figs,'tag','wave_clus_aux2');
h_fig1 = findobj(h_figs,'tag','wave_clus_aux3');
set(h_fig4,'userdata',USER_DATA)
set(h_fig3,'userdata',USER_DATA)
set(h_fig2,'userdata',USER_DATA)
set(h_fig1,'userdata',USER_DATA)


% fix3 button --------------------------------------------------------------------
function fix3_button_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_figure,'userdata');
classes = USER_DATA{6};
fix_class = find(classes==3);
if get(handles.fix3_button,'value') ==1
    USER_DATA{22} = fix_class;
else
    USER_DATA{22} = [];
end
set(handles.wave_clus_figure,'userdata',USER_DATA)
h_figs=get(0,'children');
h_fig4 = findobj(h_figs,'tag','wave_clus_aux');
h_fig3 = findobj(h_figs,'tag','wave_clus_aux1');
h_fig2 = findobj(h_figs,'tag','wave_clus_aux2');
h_fig1 = findobj(h_figs,'tag','wave_clus_aux3');
set(h_fig4,'userdata',USER_DATA)
set(h_fig3,'userdata',USER_DATA)
set(h_fig2,'userdata',USER_DATA)
set(h_fig1,'userdata',USER_DATA)


%SETTING OF SPIKE FEATURES OR PROJECTIONS
% --------------------------------------------------------------------
function spike_shapes_button_Callback(hObject, eventdata, handles)
set(gcbo,'value',1);
set(handles.spike_features_button,'value',0);
USER_DATA = get(handles.wave_clus_figure,'userdata');
cluster_results = USER_DATA{10};
handles.setclus = 1;
handles.force = 0;
handles.merge = 0;
handles.reject = 0;
handles.undo = 0;
handles.minclus = cluster_results(1,5);
plot_spikes(handles);
% -------------------------------------------------------------------
function spike_features_button_Callback(hObject, eventdata, handles)
set(gcbo,'value',1);
set(handles.spike_shapes_button,'value',0);
USER_DATA = get(handles.wave_clus_figure,'userdata');
cluster_results = USER_DATA{10};
handles.setclus = 1;
handles.force = 0;
handles.merge = 0;
handles.reject = 0;
handles.undo = 0;
handles.minclus = cluster_results(1,5);
plot_spikes(handles);


%SETTING OF SPIKE PLOTS
% --------------------------------------------------------------------
function plot_all_button_Callback(hObject, eventdata, handles)
set(gcbo,'value',1);
set(handles.plot_average_button,'value',0);
USER_DATA = get(handles.wave_clus_figure,'userdata');
cluster_results = USER_DATA{10};
handles.setclus = 1;
handles.force = 0;
handles.merge = 0;
handles.reject = 0;
handles.undo = 0;
handles.minclus = cluster_results(1,5);
plot_spikes(handles);
% --------------------------------------------------------------------
function plot_average_button_Callback(hObject, eventdata, handles)
set(gcbo,'value',1);
set(handles.plot_all_button,'value',0);
USER_DATA = get(handles.wave_clus_figure,'userdata');
cluster_results = USER_DATA{10};
handles.setclus = 1;
handles.force = 0;
handles.merge = 0;
handles.reject = 0;
handles.undo = 0;
handles.minclus = cluster_results(1,5);
plot_spikes(handles);


%SETTING OF ISI HISTOGRAMS
% --------------------------------------------------------------------
function isi1_nbins_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_figure,'userdata');
par = USER_DATA{1};
cluster_results = USER_DATA{10};
par.nbins1 = str2num(get(hObject, 'String'));
USER_DATA{1} = par;
set(handles.wave_clus_figure,'userdata',USER_DATA);
handles.setclus = 1;
handles.force = 0;
handles.merge = 0;
handles.reject = 0;
handles.undo = 0;
handles.minclus = cluster_results(1,5);
plot_spikes(handles)
% --------------------------------------------------------------------
function isi1_bin_step_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_figure,'userdata');
par = USER_DATA{1};
cluster_results = USER_DATA{10};
par.bin_step1 = str2num(get(hObject, 'String'));
USER_DATA{1} = par;
set(handles.wave_clus_figure,'userdata',USER_DATA);
handles.setclus = 1;
handles.force = 0;
handles.merge = 0;
handles.reject = 0;
handles.undo = 0;
handles.minclus = cluster_results(1,5);
plot_spikes(handles)
% --------------------------------------------------------------------
function isi2_nbins_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_figure,'userdata');
par = USER_DATA{1};
cluster_results = USER_DATA{10};
par.nbins2 = str2num(get(hObject, 'String'));
USER_DATA{1} = par;
set(handles.wave_clus_figure,'userdata',USER_DATA);
handles.setclus = 1;
handles.force = 0;
handles.merge = 0;
handles.reject = 0;
handles.undo = 0;
handles.minclus = cluster_results(1,5);
plot_spikes(handles)
% --------------------------------------------------------------------
function isi2_bin_step_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_figure,'userdata');
par = USER_DATA{1};
cluster_results = USER_DATA{10};
par.bin_step2 = str2num(get(hObject, 'String'));
USER_DATA{1} = par;
set(handles.wave_clus_figure,'userdata',USER_DATA);
handles.setclus = 1;
handles.force = 0;
handles.merge = 0;
handles.reject = 0;
handles.undo = 0;
handles.minclus = cluster_results(1,5);
plot_spikes(handles)
% --------------------------------------------------------------------
function isi3_nbins_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_figure,'userdata');
par = USER_DATA{1};
cluster_results = USER_DATA{10};
par.nbins3 = str2num(get(hObject, 'String'));
USER_DATA{1} = par;
set(handles.wave_clus_figure,'userdata',USER_DATA);
handles.setclus = 1;
handles.force = 0;
handles.merge = 0;
handles.reject = 0;
handles.undo = 0;
handles.minclus = cluster_results(1,5);
plot_spikes(handles)
% --------------------------------------------------------------------
function isi3_bin_step_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_figure,'userdata');
par = USER_DATA{1};
cluster_results = USER_DATA{10};
par.bin_step3 = str2num(get(hObject, 'String'));
USER_DATA{1} = par;
set(handles.wave_clus_figure,'userdata',USER_DATA);
handles.setclus = 1;
handles.force = 0;
handles.merge = 0;
handles.reject = 0;
handles.undo = 0;
handles.minclus = cluster_results(1,5);
plot_spikes(handles)
% --------------------------------------------------------------------
function isi0_nbins_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_figure,'userdata');
par = USER_DATA{1};
cluster_results = USER_DATA{10};
par.nbins0 = str2num(get(hObject, 'String'));
USER_DATA{1} = par;
set(handles.wave_clus_figure,'userdata',USER_DATA);
handles.setclus = 1;
handles.force = 0;
handles.merge = 0;
handles.reject = 0;
handles.undo = 0;
handles.minclus = cluster_results(1,5);
plot_spikes(handles)
% --------------------------------------------------------------------
function isi0_bin_step_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_figure,'userdata');
par = USER_DATA{1};
cluster_results = USER_DATA{10};
par.bin_step0 = str2num(get(hObject, 'String'));
USER_DATA{1} = par;
set(handles.wave_clus_figure,'userdata',USER_DATA);
handles.setclus = 1;
handles.force = 0;
handles.merge = 0;
handles.reject = 0;
handles.undo = 0;
handles.minclus = cluster_results(1,5);
plot_spikes(handles)



%SETTING OF ISI BUTTONS

% --------------------------------------------------------------------
function isi1_accept_button_Callback(hObject, eventdata, handles)
set(gcbo,'value',1);
set(handles.isi1_reject_button,'value',0);

% --------------------------------------------------------------------
function isi1_reject_button_Callback(hObject, eventdata, handles)
set(gcbo,'value',1);
set(handles.isi1_accept_button,'value',0);
USER_DATA = get(handles.wave_clus_figure,'userdata');
classes = USER_DATA{6};
tree = USER_DATA{5};
classes(find(classes==1))=0;
USER_DATA{6} = classes;
USER_DATA{9} = classes;

clustering_results = USER_DATA{10};
handles.undo = 0;
handles.force = 0;
handles.merge = 0;
handles.reject = 1; 
handles.setclus = 1;
handles.minclus = clustering_results(1,5);
set(handles.wave_clus_figure,'userdata',USER_DATA);
plot_spikes(handles)

USER_DATA = get(handles.wave_clus_figure,'userdata');
clustering_results = USER_DATA{10};
mark_clusters_temperature_diagram(handles,tree,clustering_results)
set(handles.wave_clus_figure,'userdata',USER_DATA);

set(gcbo,'value',0);
set(handles.isi1_accept_button,'value',1);

% --------------------------------------------------------------------
function isi2_accept_button_Callback(hObject, eventdata, handles)
set(gcbo,'value',1);
set(handles.isi2_reject_button,'value',0);

% --------------------------------------------------------------------
function isi2_reject_button_Callback(hObject, eventdata, handles)
set(gcbo,'value',1);
set(handles.isi2_accept_button,'value',0);
USER_DATA = get(handles.wave_clus_figure,'userdata');
classes = USER_DATA{6};
tree = USER_DATA{5};
classes(find(classes==2))=0;
USER_DATA{6} = classes;
USER_DATA{9} = classes;

clustering_results = USER_DATA{10};
handles.undo = 0;
handles.force = 0;
handles.merge = 0;
handles.reject = 1; 
handles.setclus = 1;
handles.minclus = clustering_results(1,5);
set(handles.wave_clus_figure,'userdata',USER_DATA);
plot_spikes(handles)

USER_DATA = get(handles.wave_clus_figure,'userdata');
clustering_results = USER_DATA{10};
mark_clusters_temperature_diagram(handles,tree,clustering_results)
set(handles.wave_clus_figure,'userdata',USER_DATA);

set(gcbo,'value',0);
set(handles.isi2_accept_button,'value',1);

% --------------------------------------------------------------------
function isi3_accept_button_Callback(hObject, eventdata, handles)
set(gcbo,'value',1);
set(handles.isi3_reject_button,'value',0);

% --------------------------------------------------------------------
function isi3_reject_button_Callback(hObject, eventdata, handles)
set(gcbo,'value',1);
set(handles.isi3_accept_button,'value',0);
USER_DATA = get(handles.wave_clus_figure,'userdata');
classes = USER_DATA{6};
tree = USER_DATA{5};

ilab = find(classes==3);

classes(find(classes==3))=0;
USER_DATA{6} = classes;
USER_DATA{9} = classes;

clustering_results = USER_DATA{10};
handles.undo = 0;
handles.force = 0;
handles.merge = 0;
handles.reject = 1; 
handles.setclus = 1;
handles.minclus = clustering_results(1,5);
set(handles.wave_clus_figure,'userdata',USER_DATA);
plot_spikes(handles)

USER_DATA = get(handles.wave_clus_figure,'userdata');
clustering_results = USER_DATA{10};
mark_clusters_temperature_diagram(handles,tree,clustering_results)
set(handles.wave_clus_figure,'userdata',USER_DATA);

set(gcbo,'value',0);
set(handles.isi3_accept_button,'value',1);

if isempty(ilab)
      nlab = imread('filelist.xlj','jpg'); figure('color','k'); image(nlab); axis off; set(gcf,'NumberTitle','off');
end

% --- Executes on button press in undo_button.
function undo_button_Callback(hObject, eventdata, handles)
handles.force = 0;
handles.merge = 0;
handles.undo = 1;
handles.setclus = 0;
handles.reject = 0;
USER_DATA = get(handles.wave_clus_figure,'userdata');
par = USER_DATA{1};
clustering_results_bk = USER_DATA{11};
USER_DATA{6} = clustering_results_bk(:,2); % old gui classes
USER_DATA{10} = clustering_results_bk;
handles.minclus = clustering_results_bk(1,5);
USER_DATA{8} = clustering_results_bk(1,1); % old gui temperatures
set(handles.wave_clus_figure,'userdata',USER_DATA)
plot_spikes(handles) % plot_spikes updates USER_DATA{11}
USER_DATA = get(handles.wave_clus_figure,'userdata');
tree = USER_DATA{5};
clustering_results = USER_DATA{10};
clustering_results_bk = USER_DATA{11};

mark_clusters_temperature_diagram(handles,tree,clustering_results_bk)
min_clus = handles.minclus;
set(handles.min_clus_edit,'string',num2str(min_clus));
set(handles.wave_clus_figure,'userdata',USER_DATA)
set(handles.wave_clus_figure,'userdata',USER_DATA)
set(handles.fix1_button,'value',0);
set(handles.fix2_button,'value',0);
set(handles.fix3_button,'value',0);
for i=4:par.max_clus
    eval(['par.fix' num2str(i) '=0;']);
end  
USER_DATA{1} = par;
set(handles.wave_clus_figure,'userdata',USER_DATA);



% --- Executes on button press in merge_button.
function merge_button_Callback(hObject, eventdata, handles)
handles.force = 0;
handles.merge = 1;
handles.undo = 0;
handles.setclus = 0;
handles.reject = 0;
USER_DATA = get(handles.wave_clus_figure,'userdata');
par = USER_DATA{1};
clustering_results = USER_DATA{10};
handles.minclus = clustering_results(1,5);
plot_spikes(handles)
USER_DATA = get(handles.wave_clus_figure,'userdata');
tree = USER_DATA{5};
clustering_results = USER_DATA{10};
mark_clusters_temperature_diagram(handles,tree,clustering_results)
set(handles.wave_clus_figure,'userdata',USER_DATA)
set(handles.fix1_button,'value',0);
set(handles.fix2_button,'value',0);
set(handles.fix3_button,'value',0);
for i=4:par.max_clus
    eval(['par.fix' num2str(i) '=0;']);
end  
USER_DATA{1} = par;
set(handles.wave_clus_figure,'userdata',USER_DATA);


% --- Executes on button press in Plot_polytrode_channels_button.
function Plot_polytrode_button_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_figure,'userdata');
par = USER_DATA{1};
if strcmp(par.filename(1:9),'polytrode')
    Plot_polytrode(handles)
elseif strcmp(par.filename(1:13),'C_sim_script_')
    handles.simname = par.filename;
    Plot_simulations(handles)
end





% --- Executes during object creation, after setting all properties.
function isi1_nbins_CreateFcn(hObject, eventdata, handles)
% hObject    handle to isi1_nbins (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function isi1_bin_step_CreateFcn(hObject, eventdata, handles)
% hObject    handle to isi1_bin_step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function isi2_nbins_CreateFcn(hObject, eventdata, handles)
% hObject    handle to isi2_nbins (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function isi2_bin_step_CreateFcn(hObject, eventdata, handles)
% hObject    handle to isi2_bin_step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function isi3_nbins_CreateFcn(hObject, eventdata, handles)
% hObject    handle to isi3_nbins (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function isi3_bin_step_CreateFcn(hObject, eventdata, handles)
% hObject    handle to isi3_bin_step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function isi0_nbins_CreateFcn(hObject, eventdata, handles)
% hObject    handle to isi0_nbins (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function isi0_bin_step_CreateFcn(hObject, eventdata, handles)
% hObject    handle to isi0_bin_step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function min_clus_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to min_clus_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton13.
function pushbutton13_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


