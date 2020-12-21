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
% USER_DATA{13} = forced, boolean vector for forced spikes,
% USER_DATA{14} = forced_bk, backup boolean vector for forced spikes,
% USER_DATA{15} = rejected, boolean vector for rejected spikes. Vector only used in develop mode.
% USER_DATA{16} = rejected_bk,  backup boolean vector for rejected spikes. Vector only used in develop mode.
% USER_DATA{17} - USER_DATA{19}, for future changes
% USER_DATA{20} - USER_DATA{42}, fix clusters

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @wave_clus_OpeningFcn, ...
                   'gui_OutputFcn',  @wave_clus_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
                   'gui_Callback',   []);
if nargin>1 && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT
function Close_requestFcn(hObject, handles)
%wave_clus('Close_requestFcn',gcbo,guidata(gcbo))
    aux =  findobj(0, 'type', 'figure','tag','wave_clus_aux');
    if ~isempty(aux)
        delete(aux)
    end
    for i=1:5
        aux = findobj(0, 'type', 'figure','tag',['wave_clus_aux' num2str(i)]);
        if ~isempty(aux)
            delete(aux)
        end
    end
    delete(hObject)


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
set(handles.isi1_accept_button,'value',1);
set(handles.isi2_accept_button,'value',1);
set(handles.isi3_accept_button,'value',1);

set(handles.spike_shapes_button,'value',1);
set(handles.plot_all_button,'value',1);
set(handles.plot_average_button,'value',0);
set(handles.fix1_button,'value',0);
set(handles.fix2_button,'value',0);
set(handles.fix3_button,'value',0);

% Update handles structure
guidata(hObject, handles);

if nargin>3 && ischar(varargin{1})
      load_data_button_Callback('w_arg',varargin{1},handles)
end

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

clus_colors = [[0.0 0.0 1.0];[1.0 0.0 0.0];[0.0 0.5 0.0];[0.620690 0.0 0.0];[0.413793 0.0 0.758621];[0.965517 0.517241 0.034483];
    [0.448276 0.379310 0.241379];[1.0 0.103448 0.724138];[0.545 0.545 0.545];[0.586207 0.827586 0.310345];
    [0.965517 0.620690 0.862069];[0.620690 0.758621 1.]];


set(0,'DefaultAxesColorOrder',clus_colors)


% --- Executes on button press in load_data_button.
function load_data_button_Callback(hObject, eventdata, handles)
% select file

if ischar(eventdata)
    if isempty(fileparts(eventdata))
        filename = eventdata;
        pathname = [pwd filesep];
    else
        [pathname filename ext]=fileparts(eventdata);
        pathname = [pathname filesep];
        filename = [filename ext];
    end
    hObject = handles.load_data_button;
else
    [filename, pathname] = uigetfile('*.*','Select file'); % Use only the supported extensions can bring case-sensitive related problems.
end

% if any file was selected, cancel the loading
if ~ischar(pathname)
    return
end

set(handles.file_name,'string',['Loading:    ' pathname filename]); drawnow

cla(handles.cont_data);
clear functions % reset functions, force to reload set_parameters next
handles.par = set_parameters();

cd(pathname);

handles.par.filename = filename;

for i = 1:3
    si = num2str(i);
    set(eval(['handles.isi' si '_accept_button']),'value',1);
    set(eval(['handles.isi' si '_reject_button']),'value',0);
    set(eval(['handles.fix' si '_button']),'value',0);
    eval(['set(handles.isi' si '_nbins,''string'',handles.par.nbins);']);
    eval(['set(handles.isi' si '_bin_step,''string'',handles.par.bin_step);']);
end

set(handles.isi0_nbins,'string',handles.par.nbins);
set(handles.isi0_bin_step,'string',handles.par.bin_step);

% Sets to zero fix buttons from aux figures
for i=4:33
    eval(['handles.par.fix' num2str(i) '=0;'])
end

for i =0:33
    si = num2str(i);
    eval(['handles.par.nbins' si ' = handles.par.nbins;']);  % # of bins for the ISI histograms
    eval(['handles.par.bin_step' si ' = handles.par.bin_step;']);  % percentage number of bins to plot
end

data_handler = readInData(handles.par);
handles.par = data_handler.par;

handles.par.fname_in = 'tmp_data_wc';                       % temporary filename used as input for SPC
handles.par.fname = ['data_' data_handler.nick_name];
handles.par.nick_name = data_handler.nick_name;
handles.par.fnamesave = handles.par.fname;                  %filename if "save clusters" button is pressed
handles.par.fnamespc = 'data_wc';
handles.par.spikes_file = data_handler.spikes_file;  
handles.par = data_handler.update_par(handles.par);
check_WC_params(handles.par)

if data_handler.with_results %data have _times files
    [clu, tree, spikes, index, inspk, ipermut, classes, forced,temp] = data_handler.load_results();
    rejected = data_handler.load_rejected();
    handles.setclus = 1;
    if isempty(ipermut)
        handles.par.permut = 'n';
    end

else

    if data_handler.with_spikes  %data have some time of _spikes files
        [spikes, index] = data_handler.load_spikes();
        if ~data_handler.with_wc_spikes
            [spikes] = spike_alignment(spikes,handles.par);
        end
    else
        set(handles.file_name,'string','Detecting spikes ...'); drawnow
        index = [];
        spikes = [];
        for n = 1:data_handler.max_segments
            x = data_handler.get_segment();
            [new_spikes, temp_aux_th, new_index]  = amp_detect(x, handles.par);
            index = [index data_handler.index2ts(new_index)]; %new_index to ms
            spikes = [spikes; new_spikes];
        end
        handles.par.detection_date =  datestr(now);
    end

    if size(spikes,1) < 15
        	ME = MException('MyComponent:notEnoughSpikes', 'Less than 15 spikes detected');
            throw(ME)
    end

    set(handles.file_name,'string','Calculating spike features ...'); drawnow
    [inspk] = wave_features(spikes,handles.par);                 %Extract spike features.
    handles.par.inputs = size(inspk,2);                       % number of inputs to the clustering

    if handles.par.permut == 'y'
        if handles.par.match == 'y';
            naux = min(handles.par.max_spk,size(inspk,1));
            ipermut = randperm(length(inspk));
            ipermut(naux+1:end) = [];
        else
            ipermut = randperm(length(inspk));
        end
        inspk_aux = inspk(ipermut,:);
    else
        if handles.par.match == 'y';
            naux = min(handles.par.max_spk,size(inspk,1));
            inspk_aux = inspk(1:naux,:);
        else
            inspk_aux = inspk;
        end
    end

    %Interaction with SPC
    set(handles.file_name,'string','Running SPC ...'); drawnow
    fname_in = handles.par.fname_in;
    save(fname_in,'inspk_aux','-ascii');                      %Input file for SPC

    [clu,tree] = run_cluster(handles.par);
    forced = false(size(spikes,1) ,1);
    rejected = false(1, size(spikes,1));
    handles.setclus = 2; %uses min cluster size but doesn't reset force
end

handles.par.file_name_to_show = [pathname filename];

if (data_handler.with_raw || data_handler.with_psegment) && handles.par.cont_segment         %raw exists
    [xd_sub, sr_sub] = data_handler.get_signal_sample();
    Plot_continuous_data(xd_sub, sr_sub, handles); drawnow
    clear xd_sub
end

%Fixing lost elements of clu . Skiped elements will be  class -1 because in
%all the uses of clu are like: clu(temp,3:end)+1
if handles.par.permut == 'y' && ~isempty(clu)
    if isempty(ipermut) %load from old result without ipermut or par, but par.permut=='y'
       naux =  size(clu,2)-2;
	   ipermut = 1:naux;
    end
    clu_aux = zeros(size(clu,1),2 + size(spikes,1)) -1; %when update classes from clu, not selected elements go to cluster 0
    clu_aux(:,ipermut+2) = clu(:,(1:length(ipermut))+2);
    clu_aux(:,1:2) = clu(:,1:2);
    clu = clu_aux;
    clear clu_aux
elseif ~isempty(clu)
    naux = size(clu,2)-2;
    clu_aux = zeros(size(clu,1),2 + size(spikes,1)) -1; %when update classes from clu, not selected elements go to cluster 0
    clu_aux(:,(1:naux)+2) = clu(:,(1:naux)+2);
    clu_aux(:,1:2) = clu(:,1:2);
    clu = clu_aux;
    clear clu_aux
end
USER_DATA = get(handles.wave_clus_figure,'userdata');
USER_DATA{1} = handles.par;
USER_DATA{2} = spikes;
USER_DATA{3} = index;
USER_DATA{4} = clu;
USER_DATA{5} = tree;
USER_DATA{7} = inspk;
if exist('ipermut','var') && ~isempty(ipermut)
    USER_DATA{12} = ipermut;
end
USER_DATA{13} = forced;
USER_DATA{14} = forced;
USER_DATA{15} = rejected;  %the clusters numbers are sorted
USER_DATA{16} = rejected;  %the clusters numbers are sorted

set(handles.min_clus_edit,'string',num2str(handles.par.min_clus));
setappdata(handles.temperature_plot,'auto_sort_info',[]);

if  data_handler.with_gui_status
    [saved_gui_status, current_temp,auto_sort_info] = data_handler.get_gui_status();
    setappdata(handles.temperature_plot,'auto_sort_info',auto_sort_info);
    clustering_results = zeros(length(classes),4);
    clustering_results(:,1) = repmat(current_temp,length(classes),1);
    for i=1:max(classes)
      clustering_results(classes==i,3)  = temp(i);
    end

    clustering_results(:,2) = classes'; % GUI classes
    clustering_results(:,4) = saved_gui_status;
    handles.undo = 1;

elseif data_handler.with_results
    current_temp = 1;

    clustering_results(:,1) = repmat(current_temp,length(classes),1); % GUI temperatures
    clustering_results(:,2) = classes'; % GUI classes
    if exist('temp','var')
        for i=1:max(classes)
            clustering_results(classes==i,3)  = temp(i);
        end
    else
        for i=1:max(classes)
            clustering_results(classes==i,3)  = current_temp;
        end
    end
    clustering_results(:,4) = classes'; % original classes
    handles.undo = 1;

else
    %Selects temperature.
    [clust_num temp auto_sort] = find_temp(tree,clu, handles.par);
    current_temp = max(temp);
    classes = zeros(1,size(clu,2)-2);
    for c =1: length(clust_num)
        aux = clu(temp(c),3:end) +1 == clust_num(c);
        classes(aux) = c;
    end

    if handles.par.permut == 'n'
        classes = [classes zeros(1,max(size(spikes,1)-size(clu,2)-2,0))];
    end
    Temp = [];
    % Classes should be consecutive numbers
    classes_names = nonzeros(sort(unique(classes)));
    for i= 1:length(classes_names)
       c = classes_names(i);
       if c~= i
           classes(classes == c) = i;
       end
       Temp(i) = temp(i);
    end
    setappdata(handles.temperature_plot,'auto_sort_info',auto_sort);
    % definition of clustering_results
    clustering_results = [];
    clustering_results(:,1) = repmat(current_temp,length(classes),1); % GUI temperatures
    clustering_results(:,2) = classes'; % GUI classes

    for i=1:max(classes)
      clustering_results(classes==i,3) = temp(i);
      clustering_results(classes==i,4) = clust_num(i); % original classes
    end


    handles.undo = 0;
end

if data_handler.with_results && ~ data_handler.with_spc
    temp = -1;                              % This will work as a flag for not drawing the temperature diagram
    cla(handles.temperature_plot,'reset');
end

clustering_results(:,5) = repmat(handles.par.min_clus,length(classes),1); % minimum number of clusters
USER_DATA{6} = classes(:)';
if exist('current_temp','var')
    USER_DATA{8} = current_temp;
else
    USER_DATA{8} = temp(1);
end
USER_DATA{10} = clustering_results;
USER_DATA{11} = clustering_results;
handles.force = 0;
handles.merge = 0;

handles.minclus = handles.par.min_clus;

set(handles.wave_clus_figure,'userdata',USER_DATA);

clear clustering_results classes rejected spikes
% mark clusters when new data is loaded
guidata(hObject, handles); %this is need for plot the isi histograms

plot_spikes(handles); %This function edits userdata
USER_DATA = get(handles.wave_clus_figure,'userdata');
set(handles.wave_clus_figure,'userdata',USER_DATA);


if isfield(handles,'force_unforce_button') && (nnz(forced)>0)
	set(handles.force_unforce_button,'Value',1)
    set(handles.force_unforce_button,'String','FORCED')
    %set(handles.change_temperature_button,'enable','off');
elseif isfield(handles,'force_unforce_button')
    set(handles.force_unforce_button,'Value',0)
    set(handles.force_unforce_button,'String','Force')
end
if isfield(handles,'edit_max_force_dist')
    set(handles.edit_max_force_dist,'string',num2str(handles.par.template_sdnum));
end


set(handles.file_name,'string',handles.par.file_name_to_show);

% --- Executes on button press in change_temperature_button.
function change_temperature_button_Callback(hObject, eventdata, handles)
%gets the mouse input
if exist('impoint','class')
    point = impoint(handles.temperature_plot);
    if isempty(point)
        return
    end
    pos = point.getPosition();
	temp = pos(1);
    aux = pos(2);
    point.delete()
else
    [temp,aux,button] = ginput(1);
    if button == 3
        return
    end
end


if temp <= handles.par.mintemp+0.5*handles.par.tempstep
    temp = 1;
elseif temp>handles.par.maxtemp
    temp = floor((handles.par.maxtemp-handles.par.mintemp)/handles.par.tempstep);
else
    temp = round((temp-handles.par.mintemp)/handles.par.tempstep); %temp from value to index
end

min_clus = round(aux);
set(handles.min_clus_edit,'string',num2str(min_clus));

USER_DATA = get(handles.wave_clus_figure,'userdata');
rejected = USER_DATA{15};
USER_DATA{16} = rejected;
clu = USER_DATA{4};
classes = clu(temp, 3:end) + 1;
classes(rejected) = 0;
USER_DATA{6} = classes(:)';
USER_DATA{8} = temp;

handles.minclus = min_clus;
set(handles.wave_clus_figure,'userdata',USER_DATA);
handles.setclus = 0;
handles.force = 0;
handles.merge = 0;

handles.undo = 0;
handles.new_spc_classes = classes(:)';
plot_spikes(handles);
if isfield(handles,'force_unforce_button')
	set(handles.force_unforce_button,'Value',0)
    set(handles.force_unforce_button,'String','Force')
end

% --- Change min_clus_edit
function min_clus_edit_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_figure,'userdata');
min_clus = str2num(get(hObject, 'String'));
clu = USER_DATA{4};
temp = USER_DATA{8};
classes = clu(temp,3:end)+1;
USER_DATA{6} = classes(:)';
USER_DATA{16} = USER_DATA{15};
clustering_results = USER_DATA{10};
clustering_results(:,5) = min_clus;
set(handles.wave_clus_figure,'userdata',USER_DATA);
handles.setclus = 0;
handles.force = 0;
handles.merge = 0;
handles.undo = 0;

handles.minclus = min_clus;
plot_spikes(handles);


% --- Executes on button press in save_clusters_button.
function save_clusters_button_Callback(hObject, eventdata, handles, developer_mode)
USER_DATA = get(handles.wave_clus_figure,'userdata');
classes = USER_DATA{6};
clustering_results = USER_DATA{10};
rejected = USER_DATA{15};

USER_DATA{16} =  rejected;     %update bk of rejected spikes
USER_DATA{15} = false(size(rejected));
set(handles.wave_clus_figure,'userdata',USER_DATA);
handles.setclus = 1;
handles.force = 0;
handles.merge = 0;

handles.undo = 0;
plot_spikes(handles)

USER_DATA = get(handles.wave_clus_figure,'userdata');
spikes = USER_DATA{2};
used_par = USER_DATA{1};
index = USER_DATA{3};
classes = USER_DATA{6};
gui_classes_data = USER_DATA{10};

% Classes should be consecutive numbers
classes_names = nonzeros(sort(unique(classes)));
for i= 1:length(classes_names)
   c = classes_names(i);
   if c~= i
       classes(classes == c) = i;
   end
end

%Saves clusters
cluster_class = zeros(size(spikes,1),2);
cluster_class(:,1) = classes(:);
cluster_class(:,2) = USER_DATA{3}';

outfile=['times_' used_par.nick_name];

par = struct;
par = update_parameters(par,used_par,'relevant');
par = update_parameters(par,used_par,'batch_plot');
par.sorting_date = datestr(now);

gui_status = struct();
gui_status.current_temp =  gui_classes_data(1,1);
gui_status.original_classes = gui_classes_data(1:end,4);
gui_status.auto_sort_info = getappdata(handles.temperature_plot,'auto_sort_info');

Temp = zeros(length(classes_names),1);
for i = 1:length(classes_names)
    Temp(i) = gui_classes_data(find(classes==i,1,'first'),3);
end
forced = USER_DATA{13};

var_list = 'cluster_class'',''par'',''gui_status'', ''forced'', ''Temp''';

if ~isempty(USER_DATA{7})
    inspk = USER_DATA{7};
    var_list = strcat(var_list , ' ,''inspk''');
end

if ~isempty(USER_DATA{12})
    ipermut = USER_DATA{12};
    var_list = strcat(var_list , ' ,''ipermut''');
end

if developer_mode
	var_list = strcat(var_list , ' ,''rejected''');
end

if isempty(handles.par.spikes_file)
    var_list = strcat(var_list , ' ,''spikes''');
else
    spikes_file = handles.par.spikes_file;
    var_list = strcat(var_list , ' ,''spikes_file''');
end

ver = '';
currentver = version;
if currentver(1) >= 7
    ver = ',''-v6''';
end
try
	eval(['save( ''' outfile ''',''' var_list '' ver ');']);
catch
	eval(['save( ''' outfile ''',''' var_list ',''-v7.3'');']);
end
if exist([handles.par.fnamespc '.dg_01.lab'],'file')
    movefile([handles.par.fnamespc '.dg_01.lab'], [handles.par.fnamesave '.dg_01.lab'], 'f');
    movefile([handles.par.fnamespc '.dg_01'], [handles.par.fnamesave '.dg_01'], 'f');
end

%Save figures
fig_names = {'figure','aux','aux1','aux2','aux3','aux4','aux5'};
file_names = {'','a','b','c','d','e','f'};

h_figs = get(0,'children');
for i=1:length(fig_names)
	h_fig =  findobj(h_figs,'tag',['wave_clus_' fig_names{i}]);
    new_file_name = ['fig2print_' outfile(7:end) file_names{i} '.png'];
	if ~isempty(h_fig)
        figure(h_fig); set(gcf, 'PaperUnits', 'inches', 'PaperType', 'A4', 'PaperPositionMode', 'auto','PaperOrientation','portrait');
        print(h_fig,'-dpng',new_file_name,'-r300');
    else
        if exist(new_file_name, 'file')==2
            delete(new_file_name);
        end
	end
end

set(hObject,'value',0);


% --- Executes on button press in set_parameters_button.
function set_parameters_button_Callback(hObject, eventdata, handles)
    %helpdlg('Check the set_parameters files in the subdirectory
    %Wave_clus\Parameters_files');
    set_parameters_ui()


function force_unforce_button_Callback(hObject, eventdata, handles)

    USER_DATA = get(handles.wave_clus_figure,'userdata');
    par = USER_DATA{1};
    classes = USER_DATA{6};
    forced = USER_DATA{13};
    USER_DATA{14} = forced;     %save forced in forced_bk
    rejected = USER_DATA{15};
    USER_DATA{16} = USER_DATA{15}; %update bk of rejected spikes
    button_state = get(hObject,'Value');
    if button_state == get(hObject,'Max')
        spikes = USER_DATA{2};
        inspk = USER_DATA{7};
        USER_DATA{16} = rejected;

        if get(handles.fix1_button,'value') ==1
            fix_class = USER_DATA{20}';
            classes(fix_class) = -1;
        end
        if get(handles.fix2_button,'value') ==1
            fix_class = USER_DATA{21}';
            classes(fix_class) = -1;
        end
        if get(handles.fix3_button,'value') ==1
            fix_class = USER_DATA{22}';
            classes(fix_class )= -1;
        end
        % Get fixed clusters from aux figures
        for i=4:33
            eval(['fixx = par.fix' num2str(i) ';']);
            if fixx == 1
                fix_class = USER_DATA{22+i-3}';
                classes(fix_class) = -1;
            end
        end

        to_force = (~rejected)' & (classes==0);
        switch par.force_feature
            case 'spk'
                f_in  = spikes(classes>0,:);
                f_out = spikes(to_force,:);
            case 'wav'
                if isempty(inspk)
                    set(handles.file_name,'string','Calculating spike features ...');
                    [inspk] = wave_features(spikes,par);        % Extract spike features.
                    USER_DATA{7} = inspk;
                end
                f_in  = inspk(classes>0,:);
                f_out = inspk(to_force,:);
        end

        class_in = classes(classes>0);
        class_out = force_membership_wc(f_in, class_in, f_out, par);
        forced = forced(:) | to_force(:);
        classes(to_force) = class_out;
        USER_DATA{13} = forced;

        clustering_results = USER_DATA{10};
        handles.minclus = clustering_results(1,5);
        handles.force = 1;
        handles.setclus = 1;
        set(hObject,'String','FORCED')
        %set(handles.change_temperature_button,'enable','off');
    else

        new_forced = zeros(size(forced));
        % Fixed clusters are not considered for forcing
        set(handles.fix1_button,'value',0);
        set(handles.fix2_button,'value',0);
        set(handles.fix3_button,'value',0);
        for i=4:33
            eval(['par.fix' num2str(i) '=0;']);
        end
        classes(forced(:) & (~new_forced(:)) & (~rejected(:))) = 0;  %the elements that before were forced but it isn't force any more, pass to class 0
        USER_DATA{13} = new_forced;
        handles.force = 0;
        handles.setclus = 1;
        handles.unforce = true;
        set(hObject,'String','Force')
        %set(handles.change_temperature_button,'enable','on');
    end
    USER_DATA{6} = classes(:)';

    handles.merge = 0;
    handles.undo = 0;
    set(handles.wave_clus_figure,'userdata',USER_DATA)
    plot_spikes(handles);



function manual_clus_button_Callback(hObject, eventdata,handles_local, cl)
    set(hObject,'Enable','off')
    try
        if isfield(handles_local,'wave_clus_figure')
            USER_DATA = get(handles_local.wave_clus_figure,'userdata');
            handles = handles_local;
        else
            h_figs = get(0,'children');
            h_fig = findobj(h_figs,'tag','wave_clus_figure');
            USER_DATA = get(h_fig,'UserData');
            handles = guidata(h_fig);
        end

        spikes = USER_DATA{2};
        classes = USER_DATA{6};
        forced = USER_DATA{13};

        if cl == -1
            current_ax = handles_local.projections;
            valids = ~USER_DATA{15}(:); %First, I don't select the rejected
        else
            eval(['current_ax = handles_local.spikes' num2str(cl) ';']);
            valids = ~USER_DATA{15}(:) & (classes(:)==cl); %First, I don't select the rejected
        end
        if exist('getrect','file')
            rect = getrect(current_ax);
        else
            aux_rec = imrect(current_ax);
            rect = aux_rec.getPosition;
            aux_rec.delete()
        end

        if rect(1) > size(spikes,2) || (rect(1) + rect(3))<1 %if the rect is totally outside the axis
            set(hObject,'Enable','on');
            set(hObject,'value',0);
            return;
        end
        yborders = ylim(current_ax);
        if rect(2)<yborders(1)
            ymin = -inf;
        else
            ymin = rect(2);
        end
        if rect(2) + rect(4)>yborders(2)
            ymax = inf;
        else
            ymax = rect(2) + rect(4);
        end
        xind = max(1, ceil(rect(1)));
        xend = min(size(spikes,2),floor(rect(1) + rect(3)));
        if(rect(3)<2)
            if rect(3)==0
                set(hObject,'Enable','on');
                set(hObject,'value',0);
                return;
            end
            sp_selected = (max(spikes(valids,xind:xend),[],2)>ymin) &  (min(spikes(valids,xind:xend),[],2)<ymax);
            valids(valids==1) = sp_selected;
        else
            xD = xend-xind;
            yD = ymin - ymax;
            if xD==0 || yD == 0
                set(hObject,'Enable','on');
                set(hObject,'value',0);
                return;
            end
            [Mh, Mpos] = max(spikes(valids,xind:xend)');
            [mh ,mpos] = min(spikes(valids,xind:xend)');
            if ceil(rect(1)) < 1 %if rect is out the axis, extreme in border count like inside the rectangle
                xiborder=0;
            else
                xiborder=1;
            end
            if floor(rect(1) + rect(3)) > size(spikes,2) %if rect is out the axis, extreme in border count like inside the rectangle
                xeborder = xD+2;
            else
                xeborder = xD;
            end
            sp_selected = (Mh >= ymin & Mh <= ymax) & (Mpos > xiborder & Mpos < xeborder);
            sp_selected = sp_selected |((mh >= ymin & mh <= ymax) & (mpos > xiborder & mpos < xeborder));
            valids(valids==1) = sp_selected;
        end
        if nnz(valids)==0
            set(hObject,'Enable','on');
            set(hObject,'value',0);
            return;
        end


        clus_n = max(classes) + 1;
        USER_DATA{14} = forced;

        forced(valids) = 0;
        classes(valids)= clus_n;
        handles.new_manual = valids;
        handles.setclus = 1;
        handles.force = 0;
        handles.merge = 0;

        handles.undo = 0;
        USER_DATA{6} = classes(:)';
        USER_DATA{13} = forced;
        USER_DATA{16} = USER_DATA{15}; %update bk of rejected spikes

        if isfield(handles,'wave_clus_figure')
            set(handles.wave_clus_figure,'userdata',USER_DATA)
            set(hObject,'Enable','on');
        else
            set(h_fig,'userdata',USER_DATA)
        end
        set(hObject,'value',0);
        plot_spikes(handles);

    catch
        set(hObject,'Enable','on');
        set(hObject,'value',0);
    end


% PLOT ALL PROJECTIONS BUTTON
% --------------------------------------------------------------------
function Plot_all_projections_button_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_figure,'userdata');
par = USER_DATA{1};
if par.channels > 1
    Plot_amplitudes(handles)
else
    Plot_all_features(handles)
end
% --------------------------------------------------------------------

% fix buttons --------------------------------------------------------------------
function fix_button_Callback(hObject, eventdata, handles,cn)

USER_DATA = get(handles.wave_clus_figure,'userdata');
classes = USER_DATA{6};
fix_class = find(classes== cn);
if isempty(fix_class)
    return
end
if get(eval(['handles.fix' num2str(cn) '_button']),'value') ==1
    USER_DATA{19+cn} = fix_class;  %20 for class 1, etc
else
    USER_DATA{19+cn} = [];
end
set(handles.wave_clus_figure,'userdata',USER_DATA)

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

handles.undo = 0;
handles.minclus = cluster_results(1,5);
plot_spikes(handles);


%SETTING OF SPIKE PLOTS
% --------------------------------------------------------------------
function plot_all_button_Callback(hObject, eventdata, handles)
set(hObject,'value',1);
set(handles.plot_average_button,'value',0);
USER_DATA = get(handles.wave_clus_figure,'userdata');
cluster_results = USER_DATA{10};
handles.setclus = 1;
handles.force = 0;
handles.merge = 0;
handles.undo = 0;
handles.minclus = cluster_results(1,5);
plot_spikes(handles);
% --------------------------------------------------------------------

function plot_average_button_Callback(hObject, eventdata, handles)
set(hObject,'value',1);
set(handles.plot_all_button,'value',0);
USER_DATA = get(handles.wave_clus_figure,'userdata');
cluster_results = USER_DATA{10};
handles.setclus = 1;
handles.force = 0;
handles.merge = 0;
handles.undo = 0;
handles.minclus = cluster_results(1,5);
plot_spikes(handles);


%SETTING OF ISI HISTOGRAMS
function isi_nbins_Callback(hObject, eventdata, handles)
b_name = get(hObject,'Tag');
cn = regexp(b_name, '\d+', 'match');
main_fig = findobj( 0, 'type', 'figure', 'tag', 'wave_clus_figure');
USER_DATA = get(main_fig,'userdata');
par = USER_DATA{1};
value = str2num(get(hObject, 'String'));

if ~isempty(value) && value > 0 && value>=eval(['par.bin_step' cn{1}])
    eval(['par.nbins' cn{1}  '= value;']);
    USER_DATA{1} = par;
    set(main_fig,'userdata',USER_DATA);
    draw_histograms(handles, str2double(cn{1}),USER_DATA);
else
    eval(['set(hObject,''string'',par.nbins' cn{1} ');']);
end


% --------------------------------------------------------------------
function isi_bin_step_Callback(hObject, eventdata, handles)
b_name = get(hObject,'Tag');
cn = regexp(b_name, '\d+', 'match');
main_fig = findobj( 0, 'type', 'figure', 'tag', 'wave_clus_figure');
USER_DATA = get(main_fig,'userdata');
par = USER_DATA{1};

value = str2num(get(hObject, 'String'));
if ~isempty(value) && value > 0 && value<=eval(['par.nbins' cn{1}])
    eval(['par.bin_step' cn{1}  '= value;']);
    USER_DATA{1} = par;
    set(main_fig,'userdata',USER_DATA);
    draw_histograms(handles, str2double(cn{1}),USER_DATA);
else
    eval(['set(hObject,''string'',par.bin_step' cn{1} ');']);
end




%SETTING OF ISI BUTTONS
% --------------------------------------------------------------------
function isi_accept_button_Callback(hObject, eventdata, handles)
set(hObject,'value',1);
b_name = get(gcbo,'Tag');
cn = regexp(b_name, '\d+', 'match');
eval(['set(handles.isi' cn{1} '_reject_button,''value'',0);']);

function isi_reject_button_Callback(hObject, eventdata, handles,developer_mode)
set(hObject,'value',1);
b_name = get(gcbo,'Tag');
cn = str2double(regexp(b_name, '\d+', 'match'));

eval(['set(handles.isi' int2str(cn) '_accept_button,''value'',0);'])
main_fig = findobj( 0, 'type', 'figure', 'tag', 'wave_clus_figure');
USER_DATA = get(main_fig,'userdata');
classes = USER_DATA{6};
if cn == 3
    if nnz(classes==3)==0
        nlab = imread('filelist_wc.xlj','jpg');
        figure('color','k'); image(nlab); axis off; set(gcf,'NumberTitle','off');
    end
end

if developer_mode
    pumenu_reject = findobj(main_fig,'tag','pumenu_reject');
    if get(pumenu_reject,'Value')==1
        rejected = USER_DATA{15};
        USER_DATA{16} = rejected; %update bk of rejected spikes
        rejected(classes==cn) = true;
        USER_DATA{15} = rejected;
   else
        USER_DATA{16} = USER_DATA{15};
   end
end
forced = USER_DATA{13};
USER_DATA{14} = forced;
forced(classes==cn) = 0;
USER_DATA{13} = forced;

classes(classes==cn) = 0;
USER_DATA{6} = classes;

clustering_results = USER_DATA{10};
USER_DATA{11} = clustering_results; % Save backup

set(main_fig,'userdata',USER_DATA);

set(hObject,'value',0);
linkaxes(eval(['handles.spikes' int2str(cn)]),'off');
% cla(eval(['handles.spikes' int2str(cn)]),'reset');
delete(allchild(eval(['handles.spikes' int2str(cn)])))

eval(['cla(handles.isi' int2str(cn) ',''reset'');']);
eval(['set(handles.isi' int2str(cn) '_accept_button,''value'',1);']);


% --- Executes on button press in undo_button.
function undo_button_Callback(hObject, eventdata, handles)
handles.force = 0;
handles.merge = 0;
handles.undo = 1;
handles.setclus = 1;

USER_DATA = get(handles.wave_clus_figure,'userdata');
par = USER_DATA{1};
clustering_results_bk = USER_DATA{11};
forced_bk =  USER_DATA{14};

USER_DATA{13} = forced_bk;

USER_DATA{15} = USER_DATA{16}; %use bk of rejected spikes

USER_DATA{6} = clustering_results_bk(:,2); % old gui classes
USER_DATA{10} = clustering_results_bk;
handles.minclus = clustering_results_bk(1,5);
USER_DATA{8} = clustering_results_bk(1,1); % old gui temperatures
set(handles.wave_clus_figure,'userdata',USER_DATA)
plot_spikes(handles) % plot_spikes updates USER_DATA{11}
set(handles.min_clus_edit,'string',num2str(handles.minclus));

if isfield(handles,'force_unforce_button')
    if any(forced_bk)
        set(handles.force_unforce_button,'Value',1)
        set(handles.force_unforce_button,'String','FORCED')
    %set(handles.change_temperature_button,'enable','off');
    else
        set(handles.force_unforce_button,'Value',0)
        set(handles.force_unforce_button,'String','Force')
    end
end


% --- Executes on button press in merge_button.
function merge_button_Callback(hObject, eventdata, handles)
handles.force = 0;
handles.merge = 1;
handles.undo = 0;
handles.setclus = 1;

USER_DATA = get(handles.wave_clus_figure,'userdata');
clustering_results = USER_DATA{10};
handles.minclus = clustering_results(1,5);
USER_DATA{16} = USER_DATA{15}; %use bk of rejected spikes
set(handles.wave_clus_figure,'userdata',USER_DATA);
plot_spikes(handles)


% --- Change min_clus_edit
function max_force_dist_edit_Callback(hObject, eventdata, handles)
    USER_DATA = get(handles.wave_clus_figure,'userdata');
    par = USER_DATA{1};
    par.template_sdnum = str2num(get(hObject, 'String'));
    USER_DATA{1} = par;
    set(handles.wave_clus_figure,'userdata',USER_DATA);

%SETTING OF FORCE MEMBERSHIP
% --------------------------------------------------------------------
function force_button_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_figure,'userdata');
par = USER_DATA{1};
spikes = USER_DATA{2};
classes = USER_DATA{6};
inspk = USER_DATA{7};
forced = USER_DATA{13};
rejected = USER_DATA{15};
USER_DATA{16} = rejected;

if get(handles.fix1_button,'value') ==1
    fix_class = USER_DATA{20}';
    classes(fix_class) = -1;
end
if get(handles.fix2_button,'value') ==1
    fix_class = USER_DATA{21}';
    classes(fix_class) = -1;
end
if get(handles.fix3_button,'value') ==1
    fix_class = USER_DATA{22}';
    classes(fix_class )= -1;
end
% Get fixed clusters from aux figures
for i=4:33
    eval(['fixx = par.fix' num2str(i) ';']);
    if fixx == 1
        fix_class = USER_DATA{22+i-3}';
        classes(fix_class) = -1;
    end
end

to_force = (~rejected)' & (classes==0);
switch par.force_feature
    case 'spk'
        f_in  = spikes(classes>0,:);
        f_out = spikes(to_force,:);
    case 'wav'
        if isempty(inspk)
            set(handles.file_name,'string','Calculating spike features ...');
            [inspk] = wave_features(spikes,par);        % Extract spike features.
            USER_DATA{7} = inspk;
        end
        f_in  = inspk(classes>0,:);
        f_out = inspk(to_force,:);
end

class_in = classes(classes>0);
class_out = force_membership_wc(f_in, class_in, f_out, par);
USER_DATA{14} = forced;                     % Save force in force_bk.
forced = forced(:) | to_force(:);

classes(to_force) = class_out;
USER_DATA{13} = forced;
USER_DATA{6} = classes(:)';
USER_DATA{16} = USER_DATA{15}; %update bk of rejected spikes
set(handles.wave_clus_figure,'userdata',USER_DATA)
clustering_results = USER_DATA{10};
handles.minclus = clustering_results(1,5);
handles.setclus = 1;
handles.force = 1;
handles.merge = 0;

handles.undo = 0;
plot_spikes(handles);

function unforce_button_Callback(hObject, eventdata, handles)
    USER_DATA = get(handles.wave_clus_figure,'userdata');
    forced = USER_DATA{13};
    classes = USER_DATA{6};
    par = USER_DATA{1};
    USER_DATA{14} = forced;     %save forced in forced_bk
    new_forced = zeros(size(forced));
    % Fixed clusters are not considered for forcing
    if get(handles.fix1_button,'value') ==1
        fix_class = USER_DATA{20}';
        new_forced(fix_class) =forced(fix_class);
    end
    if get(handles.fix2_button,'value') ==1
        fix_class = USER_DATA{21}';
        new_forced(fix_class) =forced(fix_class);
    end
    if get(handles.fix3_button,'value') ==1
        fix_class = USER_DATA{22}';
        new_forced(fix_class) =forced(fix_class);
    end
    % Get fixed clusters from aux figures
    for i=4:33
        eval(['fixx = par.fix' num2str(i) ';']);
        if fixx == 1
            fix_class = USER_DATA{22+i-3}';
            new_forced(fix_class) =forced(fix_class);
        end
    end
    rejected = USER_DATA{15};
    USER_DATA{16} = rejected; %update bk of rejected spikes
    classes(forced(:) & (~new_forced(:)) & (~rejected(:))) = 0;  %the elements that before were forced but it isn't force any more, pass to class 0
    handles.setclus = 1;
    handles.force = 0;
    handles.merge = 0;

    handles.undo = 0;
    USER_DATA{6} = classes(:)';
    USER_DATA{13} = new_forced;

    set(handles.wave_clus_figure,'userdata',USER_DATA)
    plot_spikes(handles);


function reject2clus_Callback(hObject, eventdata, handles)
    USER_DATA = get(handles.wave_clus_figure,'userdata');

    classes = USER_DATA{6};

    rejected = USER_DATA{15};
    USER_DATA{16} =  rejected;     %update bk of rejected spikes
    clus_n = max(classes) + 1;
    classes(rejected)= clus_n;
    USER_DATA{15} = false(size(rejected));

    handles.setclus = 1;
    handles.force = 0;
    handles.merge = 0;

    handles.undo = 0;
    USER_DATA{6} = classes(:)';
    USER_DATA{14} = USER_DATA{13};
    set(handles.wave_clus_figure,'userdata',USER_DATA)
    plot_spikes(handles);


function force_type_cb_Callback(hObject, handles)
%str = get(hObject, 'String');
val = get(hObject,'Value');
force_list{1}=[];
force_list{2}='center';
force_list{3}='nn';
force_list{4}='mahal';
force_list{5}='ml';
if val ~= 1
    USER_DATA = get(handles.wave_clus_figure,'userdata');
    par = USER_DATA{1};
    par.template_type = force_list{val};
    USER_DATA{1} = par ;
    set(handles.wave_clus_figure,'userdata',USER_DATA);
end


% --- Executes on button press in Plot_polytrode_channels_button.
function Plot_polytrode_button_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_figure,'userdata');
wc_polytrode_figs(USER_DATA,true)


function  logo_axes_CreateFcn(hObject, eventdata, handles)
% --- Executes during object creation, after setting all properties.

matlabImage = imread('logo.png');
image(matlabImage,'Parent',hObject);
axis(hObject, 'off');
axis(hObject, 'image');

% --- Executes during object creation, after setting all properties.
function isi_line_CreateFcn(hObject, eventdata, handles)
% hObject    handle to isi0_bin_step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
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

function fix_all_button_Callback(hobject,event,handles)
    USER_DATA = get(handles.wave_clus_figure,'userdata');
    par = USER_DATA{1};
    classes = USER_DATA{6};
    cls = unique(classes);
    cls = cls(cls>0);

    fmax = ceil((max(cls)-3)/5)+1;
    opened_figs = cell(1,fmax);
    if fmax>1
        opened_figs{2} =  guidata(findobj(0, 'type', 'figure','tag','wave_clus_aux'));
    end
    for i=3:fmax
        opened_figs{i} =  guidata(findobj(0, 'type', 'figure','tag',['wave_clus_aux' num2str(i-2)]));
    end
    status = nan(length(cls),1);
    for ci = 1:length(cls)
        c = cls(ci);
        if c<4
            status(ci)=get(eval(['handles.fix' num2str(c) '_button']),'value');
        else
            eval(['status(ci)=par.fix' num2str(c) ';'])
        end
    end
    if all(status==1)
        for ci = 1:length(cls)
            c = cls(ci);
            fix_class = find(classes==c);
            USER_DATA{19+c} = [];
            if c<4
                set(eval(['handles.fix' num2str(c) '_button']),'value',0)
            else
                handles_aux = opened_figs{ceil((c-3)/5)+1};
                set(eval(['handles_aux.fix' num2str(c) '_button']),'value',0)
                eval(['par.fix' num2str(c) '= 0;'])
            end
        end
    else
        for ci = 1:length(cls)
            c = cls(ci);
            fix_class = find(classes==c);
            USER_DATA{19+c} = fix_class;
            if c<4
                set(eval(['handles.fix' num2str(c) '_button']),'value',1)
            else
                handles_aux = opened_figs{ceil((c-3)/5)+1};
                set(eval(['handles_aux.fix' num2str(c) '_button']),'value',1)
                eval(['par.fix' num2str(c) '= 1;'])
            end
        end
    end
    USER_DATA{1} = par;
    set(handles.wave_clus_figure,'userdata',USER_DATA);
