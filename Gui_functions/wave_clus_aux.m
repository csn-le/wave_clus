function varargout = wave_clus_aux(varargin)
% WAVE_CLUS_AUX M-file for wave_clus_aux.fig
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

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @wave_clus_aux_OpeningFcn, ...
                   'gui_OutputFcn',  @wave_clus_aux_OutputFcn, ...
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
function wave_clus_aux_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to wave_clus_aux (see VARARGIN)

% Choose default command line output for wave_clus_aux
handles.output = hObject;
set(handles.isi4_accept_button,'value',1);
set(handles.isi5_accept_button,'value',1);
set(handles.isi6_accept_button,'value',1);
set(handles.isi7_accept_button,'value',1);
set(handles.isi8_accept_button,'value',1);
set(handles.fix4_button,'value',0);
set(handles.fix5_button,'value',0);
set(handles.fix6_button,'value',0);
set(handles.fix7_button,'value',0);
set(handles.fix8_button,'value',0);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes wave_clus_aux wait for user response (see UIRESUME)
% uiwait(handles.wave_clus_aux);


% --- Outputs from this function are returned to the command line.
function varargout = wave_clus_aux_OutputFcn(hObject, eventdata, handles)
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

for i = 4:8
	si = num2str(i);
	set(eval(['handles.isi' si '_accept_button']),'value',1);
	set(eval(['handles.isi' si '_reject_button']),'value',0);
	set(eval(['handles.fix' si '_button']),'value',0);
	
	eval(['set(handles.isi' si '_nbins,''string'',par.nbins' si ');']);
	eval(['set(handles.isi' si '_bin_step,''string'',par.bin_step' si ');']);
	
end

plot_spikes_aux(handles,USER_DATA,0)


% FIX buttons
% --------------------------------------------------------
function fix_button_Callback(hObject, eventdata, handles, cn)
main_fig = findobj( 0, 'type', 'figure', 'tag', 'wave_clus_figure');

USER_DATA = get(main_fig,'userdata');
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
set(main_fig,'userdata',USER_DATA);



