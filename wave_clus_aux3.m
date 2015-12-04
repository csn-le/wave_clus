function varargout = wave_clus_aux3(varargin)
% WAVE_CLUS_AUX1 M-file for wave_clus_aux.fig
%      WAVE_CLUS_AUX, by itself, creates a new WAVE_CLUS_AUX or raises the existing
%      singleton*.
%
%      H = WAVE_CLUS_AUX returns the handle to a new WAVE_CLUS_AUX or the handle to
%      the existing singleton*.
%
%      WAVE_CLUS_AUX('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in WAVE_CLUS_AUX.M with the given input arguments.
%
%      WAVE_CLUS_AUX('Property','Value',...) creates a new WAVE_CLUS_AUX or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before wave_clus_aux_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to wave_clus_aux_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help wave_clus_aux

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @wave_clus_aux3_OpeningFcn, ...
                   'gui_OutputFcn',  @wave_clus_aux3_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && isstr(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before wave_clus_aux is made visible.
function wave_clus_aux3_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to wave_clus_aux (see VARARGIN)

% Choose default command line output for wave_clus_aux
handles.output = hObject;
set(handles.isi19_accept_button,'value',1);
set(handles.isi20_accept_button,'value',1);
set(handles.isi21_accept_button,'value',1);
set(handles.isi22_accept_button,'value',1);
set(handles.isi23_accept_button,'value',1);
set(handles.fix19_button,'value',0);
set(handles.fix20_button,'value',0);
set(handles.fix21_button,'value',0);
set(handles.fix22_button,'value',0);
set(handles.fix23_button,'value',0);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes wave_clus_aux wait for user response (see UIRESUME)
% uiwait(handles.wave_clus_aux);


% --- Outputs from this function are returned to the command line.
function varargout = wave_clus_aux3_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

h_figs=get(0,'children');
h_fig = findobj(h_figs,'tag','wave_clus_figure');
USER_DATA = get(h_fig,'UserData');
par = USER_DATA{1};

for i = 19:23
	si = num2str(i);
	set(eval(['handles.isi' si '_accept_button']),'value',1);
	set(eval(['handles.isi' si '_reject_button']),'value',0);
	set(eval(['handles.fix' si '_button']),'value',0);
	
	eval(['set(handles.isi' si '_nbins,''string'',par.nbins' si ');']);
	eval(['set(handles.isi' si '_bin_step,''string'',par.bin_step' si ');']);
	
end

USER_DATA{1} = par;
set(handles.wave_clus_aux3,'userdata',USER_DATA)
plot_spikes_aux(handles,3)


% Change nbins
% -------------------------------------------------------------
function isi_nbins_Callback(hObject, eventdata, handles)
b_name = get(gcbo,'Tag');
cn = regexp(b_name, '\d+', 'match');
USER_DATA = get(handles.wave_clus_aux3,'userdata');
par = USER_DATA{1};
eval(['par.nbins' cn{1}  '= str2num(get(hObject, ''String''));']);
USER_DATA{1} = par;
set(handles.wave_clus_aux3,'userdata',USER_DATA);
draw_histograms(handles,  str2double(cn{1}),USER_DATA);

% Change bin steps
% --------------------------------------------------------------------
function isi_bin_step_Callback(hObject, eventdata, handles)
b_name = get(gcbo,'Tag');
cn = regexp(b_name, '\d+', 'match');
USER_DATA = get(handles.wave_clus_aux3,'userdata');
par = USER_DATA{1};
eval(['par.bin_step' cn{1}  '= str2num(get(hObject, ''String''));']);
USER_DATA{1} = par;
set(handles.wave_clus_aux3,'userdata',USER_DATA);
draw_histograms(handles, str2double(cn{1}),USER_DATA);

% Reject buttons

% --------------------------------------------------------------------
function isi_reject_button_Callback(hObject, eventdata, handles, developer_mode)
set(hObject,'value',1);
b_name = get(gcbo,'Tag');
cn = str2double(regexp(b_name, '\d+', 'match'));

eval(['set(handles.isi' int2str(cn) '_accept_button,''value'',0);'])
USER_DATA = get(handles.wave_clus_aux3,'userdata');
classes = USER_DATA{6};

h_figs = get(0,'children');

if developer_mode
    pumenu_reject = findobj(h_figs,'tag','pumenu_reject');
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
clustering_results(:,2) = classes;
USER_DATA{10} = clustering_results; 

h_fig{1} = findobj(h_figs,'tag','wave_clus_figure');
h_fig{2} = findobj(h_figs,'tag','wave_clus_aux');
h_fig{3} = findobj(h_figs,'tag','wave_clus_aux2');
h_fig{4} = findobj(h_figs,'tag','wave_clus_aux1');
h_fig{5} = findobj(h_figs,'tag','wave_clus_aux4');
h_fig{6} = findobj(h_figs,'tag','wave_clus_aux5');
set(handles.wave_clus_aux3,'userdata',USER_DATA);

for h = h_fig
    set(h{1},'userdata',USER_DATA)
end

set(hObject,'value',0);
eval(['cla(handles.spikes' int2str(cn) ',''reset'');']);
eval(['cla(handles.isi' int2str(cn) ',''reset'');']);
eval(['set(handles.isi' int2str(cn) '_accept_button,''value'',1);']);


% FIX buttons
% --------------------------------------------------------
function fix_button_Callback(hObject, eventdata, handles, cn)
USER_DATA = get(handles.wave_clus_aux3,'userdata');
par = USER_DATA{1};
classes = USER_DATA{6};
fix_class = find(classes==cn);

if get(eval(['handles.fix' num2str(cn) '_button']),'value') ==1
    USER_DATA{19+cn} = fix_class;
    eval(['par.fix' num2str(cn) '= 1;'])
else
    USER_DATA{19+cn} = [];
    eval(['par.fix' num2str(cn) '= 0;'])
end
USER_DATA{1} = par;
h_figs=get(0,'children');
h_fig = findobj(h_figs,'tag','wave_clus_figure');
h_fig1 = findobj(h_figs,'tag','wave_clus_aux1');
h_fig2 = findobj(h_figs,'tag','wave_clus_aux2');
h_fig3 = findobj(h_figs,'tag','wave_clus_aux');
h_fig4 = findobj(h_figs,'tag','wave_clus_aux4');
h_fig5 = findobj(h_figs,'tag','wave_clus_aux5');
set(handles.wave_clus_aux3,'userdata',USER_DATA);
set(h_fig,'userdata',USER_DATA)
set(h_fig1,'userdata',USER_DATA)
set(h_fig2,'userdata',USER_DATA)
set(h_fig3,'userdata',USER_DATA)
set(h_fig4,'userdata',USER_DATA)
set(h_fig5,'userdata',USER_DATA)

