function varargout = wave_clus_aux5(varargin)
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

% Last Modified by GUIDE v2.5 16-Dec-2004 18:37:28

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @wave_clus_aux5_OpeningFcn, ...
                   'gui_OutputFcn',  @wave_clus_aux5_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
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


% --- Executes just before wave_clus_aux is made visible.
function wave_clus_aux5_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to wave_clus_aux (see VARARGIN)

% Choose default command line output for wave_clus_aux
handles.output = hObject;
set(handles.isi29_accept_button,'value',1);
set(handles.isi30_accept_button,'value',1);
set(handles.isi31_accept_button,'value',1);
set(handles.isi32_accept_button,'value',1);
set(handles.isi33_accept_button,'value',1);
set(handles.fix29_button,'value',0);
set(handles.fix30_button,'value',0);
set(handles.fix31_button,'value',0);
set(handles.fix32_button,'value',0);
set(handles.fix33_button,'value',0);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes wave_clus_aux wait for user response (see UIRESUME)
% uiwait(handles.wave_clus_aux);


% --- Outputs from this function are returned to the command line.
function varargout = wave_clus_aux5_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

set(handles.isi29_accept_button,'value',1);
set(handles.isi30_accept_button,'value',1);
set(handles.isi31_accept_button,'value',1);
set(handles.isi32_accept_button,'value',1);
set(handles.isi33_accept_button,'value',1);
set(handles.isi29_reject_button,'value',0);
set(handles.isi30_reject_button,'value',0);
set(handles.isi31_reject_button,'value',0);
set(handles.isi32_reject_button,'value',0);
set(handles.isi33_reject_button,'value',0);
set(handles.fix29_button,'value',0);
set(handles.fix30_button,'value',0);
set(handles.fix31_button,'value',0);
set(handles.fix32_button,'value',0);
set(handles.fix33_button,'value',0);

h_figs=get(0,'children');
h_fig = findobj(h_figs,'tag','wave_clus_figure');
h_fig1 = findobj(h_figs,'tag','wave_clus_aux');
h_fig2 = findobj(h_figs,'tag','wave_clus_aux1');
h_fig3 = findobj(h_figs,'tag','wave_clus_aux2');
h_fig4 = findobj(h_figs,'tag','wave_clus_aux3');
h_fig5 = findobj(h_figs,'tag','wave_clus_aux4');
USER_DATA = get(h_fig,'UserData');
par = USER_DATA{1};

set(handles.isi29_nbins,'string',par.nbins29);
set(handles.isi30_nbins,'string',par.nbins30);
set(handles.isi31_nbins,'string',par.nbins31);
set(handles.isi32_nbins,'string',par.nbins32);
set(handles.isi33_nbins,'string',par.nbins33);
set(handles.isi29_bin_step,'string',par.bin_step29);
set(handles.isi30_bin_step,'string',par.bin_step30);
set(handles.isi31_bin_step,'string',par.bin_step31);
set(handles.isi32_bin_step,'string',par.bin_step32);
set(handles.isi33_bin_step,'string',par.bin_step33);

% That's for passing the fix button settings to plot_spikes.
if get(handles.fix29_button,'value') ==1     
    par.fix29 = 1;
else
    par.fix29 = 0;
end
if get(handles.fix30_button,'value') ==1     
    par.fix30 = 1;
else
    par.fix30 = 0;
end
if get(handles.fix31_button,'value') ==1     
    par.fix31 = 1;
else
    par.fix31 = 0;
end
if get(handles.fix32_button,'value') ==1     
    par.fix32 = 1;
else
    par.fix32 = 0;
end
if get(handles.fix33_button,'value') ==1     
    par.fix33 = 1;
else
    par.fix33 = 0;
end
USER_DATA{1} = par;
set(handles.wave_clus_aux5,'userdata',USER_DATA)
set(h_fig,'userdata',USER_DATA)
set(h_fig1,'userdata',USER_DATA)
set(h_fig2,'userdata',USER_DATA)
set(h_fig3,'userdata',USER_DATA)
set(h_fig4,'userdata',USER_DATA)
set(h_fig5,'userdata',USER_DATA)

plot_spikes_aux5(handles)




% Change nbins
% -------------------------------------------------------------
function isi29_nbins_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_aux5,'userdata');
par = USER_DATA{1};
par.nbins29 = str2num(get(hObject, 'String'));
par.axes_nr = 30;
classes = USER_DATA{6};
par.class_to_plot = find(classes==25);
USER_DATA{1} = par;
USER_DATA{6} = classes;
set(handles.wave_clus_aux5,'userdata',USER_DATA);
plot_spikes_aux5(handles)
% --------------------------------------------------------------------
function isi30_nbins_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_aux5,'userdata');
par = USER_DATA{1};
par.nbins30 = str2num(get(hObject, 'String'));
par.axes_nr = 31;
classes = USER_DATA{6};
par.class_to_plot = find(classes==30);
USER_DATA{1} = par;
USER_DATA{6} = classes;
set(handles.wave_clus_aux5,'userdata',USER_DATA);
plot_spikes_aux5(handles)
% --------------------------------------------------------------------
function isi31_nbins_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_aux5,'userdata');
par = USER_DATA{1};
par.nbins31 = str2num(get(hObject, 'String'));
par.axes_nr = 32;
classes = USER_DATA{6};
par.class_to_plot = find(classes==31);
USER_DATA{1} = par;
USER_DATA{6} = classes;
set(handles.wave_clus_aux5,'userdata',USER_DATA);
plot_spikes_aux5(handles)
% --------------------------------------------------------------------
function isi32_nbins_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_aux5,'userdata');
par = USER_DATA{1};
par.nbins32 = str2num(get(hObject, 'String'));
par.axes_nr = 33;
classes = USER_DATA{6};
par.class_to_plot = find(classes==32);
USER_DATA{1} = par;
USER_DATA{6} = classes;
set(handles.wave_clus_aux5,'userdata',USER_DATA);
plot_spikes_aux5(handles)
% --------------------------------------------------------------------
function isi33_nbins_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_aux5,'userdata');
par = USER_DATA{1};
par.nbins33 = str2num(get(hObject, 'String'));
par.axes_nr = 34;
classes = USER_DATA{6};
par.class_to_plot = find(classes==33);
USER_DATA{1} = par;
USER_DATA{6} = classes;
set(handles.wave_clus_aux5,'userdata',USER_DATA);
plot_spikes_aux5(handles)
% --------------------------------------------------------------------


% Change bin steps
% --------------------------------------------------------
function isi29_bin_step_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_aux5,'userdata');
par = USER_DATA{1};
par.bin_step29 = str2num(get(hObject, 'String'));
par.axes_nr = 30;
classes = USER_DATA{6};
par.class_to_plot = find(classes==29);
USER_DATA{1} = par;
USER_DATA{6} = classes;
set(handles.wave_clus_aux5,'userdata',USER_DATA);
plot_spikes_aux5(handles)
% --------------------------------------------------------------------
function isi30_bin_step_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_aux5,'userdata');
par = USER_DATA{1};
par.bin_step30 = str2num(get(hObject, 'String'));
par.axes_nr = 31;
classes = USER_DATA{6};
par.class_to_plot = find(classes==30);
USER_DATA{1} = par;
USER_DATA{6} = classes;
set(handles.wave_clus_aux5,'userdata',USER_DATA);
plot_spikes_aux5(handles)
% --------------------------------------------------------------------
function isi31_bin_step_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_aux5,'userdata');
par = USER_DATA{1};
par.bin_step31 = str2num(get(hObject, 'String'));
par.axes_nr = 32;
classes = USER_DATA{6};
par.class_to_plot = find(classes==31);
USER_DATA{1} = par;
USER_DATA{6} = classes;
set(handles.wave_clus_aux5,'userdata',USER_DATA);
plot_spikes_aux5(handles)
% --------------------------------------------------------------------
function isi32_bin_step_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_aux5,'userdata');
par = USER_DATA{1};
par.bin_step32 = str2num(get(hObject, 'String'));
par.axes_nr = 33;
classes = USER_DATA{6};
par.class_to_plot = find(classes==32);
USER_DATA{1} = par;
USER_DATA{6} = classes;
set(handles.wave_clus_aux5,'userdata',USER_DATA);
plot_spikes_aux5(handles)
% --------------------------------------------------------------------
function isi33_bin_step_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_aux5,'userdata');
par = USER_DATA{1};
par.bin_step33 = str2num(get(hObject, 'String'));
par.axes_nr = 34;
classes = USER_DATA{6};
par.class_to_plot = find(classes==33);
USER_DATA{1} = par;
USER_DATA{6} = classes;
set(handles.wave_clus_aux5,'userdata',USER_DATA);
plot_spikes_aux5(handles)
% --------------------------------------------------------------------


% Accept and Reject buttons
% --------------------------------------------------------
function isi29_accept_button_Callback(hObject, eventdata, handles)
set(gcbo,'value',1);
set(handles.isi29_reject_button,'value',0);
% --------------------------------------------------------------------
function isi29_reject_button_Callback(hObject, eventdata, handles)
set(gcbo,'value',1);
set(handles.isi29_accept_button,'value',0);
USER_DATA = get(handles.wave_clus_aux5,'userdata');
classes = USER_DATA{6};
classes(find(classes==29))=0;
USER_DATA{6} = classes;
USER_DATA{9} = classes;
h_figs=get(0,'children');
h_fig = findobj(h_figs,'tag','wave_clus_figure');
h_fig1 = findobj(h_figs,'tag','wave_clus_aux');
h_fig2 = findobj(h_figs,'tag','wave_clus_aux1');
h_fig3 = findobj(h_figs,'tag','wave_clus_aux2');
h_fig4 = findobj(h_figs,'tag','wave_clus_aux3');
h_fig5 = findobj(h_figs,'tag','wave_clus_aux4');
set(handles.wave_clus_aux5,'userdata',USER_DATA);
set(h_fig,'userdata',USER_DATA)
set(h_fig1,'userdata',USER_DATA)
set(h_fig2,'userdata',USER_DATA)
set(h_fig3,'userdata',USER_DATA)
set(h_fig4,'userdata',USER_DATA)
set(h_fig5,'userdata',USER_DATA)
axes(handles.spikes29); 
cla reset
axes(handles.isi29); 
cla reset
set(gcbo,'value',0);
set(handles.isi29_accept_button,'value',1);

% --------------------------------------------------------
function isi30_accept_button_Callback(hObject, eventdata, handles)
set(gcbo,'value',1);
set(handles.isi30_reject_button,'value',0);
% --------------------------------------------------------------------
function isi30_reject_button_Callback(hObject, eventdata, handles)
set(gcbo,'value',1);
set(handles.isi30_accept_button,'value',0);
USER_DATA = get(handles.wave_clus_aux5,'userdata');
classes = USER_DATA{6};
classes(find(classes==30))=0;
USER_DATA{6} = classes;
USER_DATA{9} = classes;
h_figs=get(0,'children');
h_fig = findobj(h_figs,'tag','wave_clus_figure');
h_fig1 = findobj(h_figs,'tag','wave_clus_aux');
h_fig2 = findobj(h_figs,'tag','wave_clus_aux1');
h_fig3 = findobj(h_figs,'tag','wave_clus_aux2');
h_fig4 = findobj(h_figs,'tag','wave_clus_aux3');
h_fig5 = findobj(h_figs,'tag','wave_clus_aux4');
set(handles.wave_clus_aux5,'userdata',USER_DATA);
set(h_fig,'userdata',USER_DATA)
set(h_fig1,'userdata',USER_DATA)
set(h_fig2,'userdata',USER_DATA)
set(h_fig3,'userdata',USER_DATA)
set(h_fig4,'userdata',USER_DATA)
set(h_fig5,'userdata',USER_DATA)
axes(handles.spikes30); 
cla reset
axes(handles.isi30); 
cla reset
set(gcbo,'value',0);
set(handles.isi30_accept_button,'value',1);

% --------------------------------------------------------
function isi31_accept_button_Callback(hObject, eventdata, handles)
set(gcbo,'value',1);
set(handles.isi31_reject_button,'value',0);
% --------------------------------------------------------------------
function isi31_reject_button_Callback(hObject, eventdata, handles)
set(gcbo,'value',1);
set(handles.isi31_accept_button,'value',0);
USER_DATA = get(handles.wave_clus_aux5,'userdata');
classes = USER_DATA{6};
classes(find(classes==31))=0;
USER_DATA{6} = classes;
USER_DATA{9} = classes;
h_figs=get(0,'children');
h_fig = findobj(h_figs,'tag','wave_clus_figure');
h_fig1 = findobj(h_figs,'tag','wave_clus_aux');
h_fig2 = findobj(h_figs,'tag','wave_clus_aux1');
h_fig3 = findobj(h_figs,'tag','wave_clus_aux2');
h_fig4 = findobj(h_figs,'tag','wave_clus_aux3');
h_fig5 = findobj(h_figs,'tag','wave_clus_aux4');
set(handles.wave_clus_aux5,'userdata',USER_DATA);
set(h_fig,'userdata',USER_DATA)
set(h_fig1,'userdata',USER_DATA)
set(h_fig2,'userdata',USER_DATA)
set(h_fig3,'userdata',USER_DATA)
set(h_fig4,'userdata',USER_DATA)
set(h_fig5,'userdata',USER_DATA)
axes(handles.spikes31); 
cla reset
axes(handles.isi31); 
cla reset
set(gcbo,'value',0);
set(handles.isi31_accept_button,'value',1);

% --------------------------------------------------------
function isi32_accept_button_Callback(hObject, eventdata, handles)
set(gcbo,'value',1);
set(handles.isi32_reject_button,'value',0);
% --------------------------------------------------------------------
function isi32_reject_button_Callback(hObject, eventdata, handles)
set(gcbo,'value',1);
set(handles.isi32_accept_button,'value',0);
USER_DATA = get(handles.wave_clus_aux5,'userdata');
classes = USER_DATA{6};
classes(find(classes==32))=0;
USER_DATA{6} = classes;
USER_DATA{9} = classes;
h_figs=get(0,'children');
h_fig = findobj(h_figs,'tag','wave_clus_figure');
h_fig1 = findobj(h_figs,'tag','wave_clus_aux');
h_fig2 = findobj(h_figs,'tag','wave_clus_aux1');
h_fig3 = findobj(h_figs,'tag','wave_clus_aux2');
h_fig4 = findobj(h_figs,'tag','wave_clus_aux3');
h_fig5 = findobj(h_figs,'tag','wave_clus_aux4');
set(handles.wave_clus_aux5,'userdata',USER_DATA);
set(h_fig,'userdata',USER_DATA)
set(h_fig1,'userdata',USER_DATA)
set(h_fig2,'userdata',USER_DATA)
set(h_fig3,'userdata',USER_DATA)
set(h_fig4,'userdata',USER_DATA)
set(h_fig5,'userdata',USER_DATA)
axes(handles.spikes32); 
cla reset
axes(handles.isi32); 
cla reset
set(gcbo,'value',0);
set(handles.isi32_accept_button,'value',1);

% --------------------------------------------------------
function isi33_accept_button_Callback(hObject, eventdata, handles)
set(gcbo,'value',1);
set(handles.isi33_reject_button,'value',0);
% --------------------------------------------------------------------
function isi33_reject_button_Callback(hObject, eventdata, handles)
set(gcbo,'value',1);
set(handles.isi33_accept_button,'value',0);
USER_DATA = get(handles.wave_clus_aux5,'userdata');
classes = USER_DATA{6};
classes(find(classes==33))=0;
USER_DATA{6} = classes;
USER_DATA{9} = classes;
h_figs=get(0,'children');
h_fig = findobj(h_figs,'tag','wave_clus_figure');
h_fig1 = findobj(h_figs,'tag','wave_clus_aux');
h_fig2 = findobj(h_figs,'tag','wave_clus_aux1');
h_fig3 = findobj(h_figs,'tag','wave_clus_aux2');
h_fig4 = findobj(h_figs,'tag','wave_clus_aux3');
h_fig5 = findobj(h_figs,'tag','wave_clus_aux4');
set(handles.wave_clus_aux5,'userdata',USER_DATA);
set(h_fig,'userdata',USER_DATA)
set(h_fig1,'userdata',USER_DATA)
set(h_fig2,'userdata',USER_DATA)
set(h_fig3,'userdata',USER_DATA)
set(h_fig4,'userdata',USER_DATA)
set(h_fig5,'userdata',USER_DATA)
axes(handles.spikes33); 
cla reset
axes(handles.isi33); 
cla reset
set(gcbo,'value',0);
set(handles.isi33_accept_button,'value',1);



% FIX buttons
% --------------------------------------------------------
function fix29_button_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_aux5,'userdata');
par = USER_DATA{1};
classes = USER_DATA{6};
fix_class = find(classes==29);
if get(handles.fix29_button,'value') ==1
    USER_DATA{48} = fix_class;
    par.fix29 = 1;
else
    USER_DATA{48} = [];
    par.fix29 = 0
end
USER_DATA{1} = par;
h_figs=get(0,'children');
h_fig = findobj(h_figs,'tag','wave_clus_figure');
h_fig1 = findobj(h_figs,'tag','wave_clus_aux');
h_fig2 = findobj(h_figs,'tag','wave_clus_aux1');
h_fig3 = findobj(h_figs,'tag','wave_clus_aux2');
h_fig4 = findobj(h_figs,'tag','wave_clus_aux3');
h_fig5 = findobj(h_figs,'tag','wave_clus_aux4');
set(handles.wave_clus_aux5,'userdata',USER_DATA);
set(h_fig,'userdata',USER_DATA)
set(h_fig1,'userdata',USER_DATA)
set(h_fig2,'userdata',USER_DATA)
set(h_fig3,'userdata',USER_DATA)
set(h_fig4,'userdata',USER_DATA)
set(h_fig5,'userdata',USER_DATA)
% --------------------------------------------------------
function fix30_button_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_aux4,'userdata');
par = USER_DATA{1};
classes = USER_DATA{6};
fix_class = find(classes==30);
if get(handles.fix30_button,'value') ==1
    USER_DATA{49} = fix_class;
    par.fix30 = 1;
else
    USER_DATA{49} = [];
    par.fix30 = 0;
end
USER_DATA{1} = par;
h_figs=get(0,'children');
h_fig = findobj(h_figs,'tag','wave_clus_figure');
h_fig1 = findobj(h_figs,'tag','wave_clus_aux');
h_fig2 = findobj(h_figs,'tag','wave_clus_aux1');
h_fig3 = findobj(h_figs,'tag','wave_clus_aux2');
h_fig4 = findobj(h_figs,'tag','wave_clus_aux3');
h_fig5 = findobj(h_figs,'tag','wave_clus_aux4');
set(handles.wave_clus_aux5,'userdata',USER_DATA);
set(h_fig,'userdata',USER_DATA)
set(h_fig1,'userdata',USER_DATA)
set(h_fig2,'userdata',USER_DATA)
set(h_fig3,'userdata',USER_DATA)
set(h_fig4,'userdata',USER_DATA)
set(h_fig5,'userdata',USER_DATA)
% --------------------------------------------------------
function fix31_button_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_aux4,'userdata');
par = USER_DATA{1};
classes = USER_DATA{6};
fix_class = find(classes==31);
if get(handles.fix31_button,'value') ==1
    USER_DATA{50} = fix_class;
    par.fix31 = 1;
else
    USER_DATA{50} = [];
    par.fix31 = 0;
end
USER_DATA{1} = par;
h_figs=get(0,'children');
h_fig = findobj(h_figs,'tag','wave_clus_figure');
h_fig1 = findobj(h_figs,'tag','wave_clus_aux');
h_fig2 = findobj(h_figs,'tag','wave_clus_aux1');
h_fig3 = findobj(h_figs,'tag','wave_clus_aux2');
h_fig4 = findobj(h_figs,'tag','wave_clus_aux3');
h_fig5 = findobj(h_figs,'tag','wave_clus_aux4');
set(handles.wave_clus_aux5,'userdata',USER_DATA);
set(h_fig,'userdata',USER_DATA)
set(h_fig1,'userdata',USER_DATA)
set(h_fig2,'userdata',USER_DATA)
set(h_fig3,'userdata',USER_DATA)
set(h_fig4,'userdata',USER_DATA)
set(h_fig5,'userdata',USER_DATA)
% --------------------------------------------------------
function fix32_button_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_aux4,'userdata');
par = USER_DATA{1};
classes = USER_DATA{6};
fix_class = find(classes==32);
if get(handles.fix32_button,'value') ==1
    USER_DATA{51} = fix_class;
    par.fix32 = 1;
else
    USER_DATA{51} = [];
    par.fix32 = 0;
end
USER_DATA{1} = par;
h_figs=get(0,'children');
h_fig = findobj(h_figs,'tag','wave_clus_figure');
h_fig1 = findobj(h_figs,'tag','wave_clus_aux');
h_fig2 = findobj(h_figs,'tag','wave_clus_aux1');
h_fig3 = findobj(h_figs,'tag','wave_clus_aux2');
h_fig4 = findobj(h_figs,'tag','wave_clus_aux3');
h_fig5 = findobj(h_figs,'tag','wave_clus_aux4');
set(handles.wave_clus_aux5,'userdata',USER_DATA);
set(h_fig,'userdata',USER_DATA)
set(h_fig1,'userdata',USER_DATA)
set(h_fig2,'userdata',USER_DATA)
set(h_fig3,'userdata',USER_DATA)
set(h_fig4,'userdata',USER_DATA)
set(h_fig5,'userdata',USER_DATA)
% --------------------------------------------------------
function fix33_button_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_aux4,'userdata');
par = USER_DATA{1};
classes = USER_DATA{6};
fix_class = find(classes==33);
if get(handles.fix33_button,'value') ==1
    USER_DATA{52} = fix_class;
    par.fix33 = 1;
else
    USER_DATA{52} = [];
    par.fix33 = 0;
end
USER_DATA{1} = par;
h_figs=get(0,'children');
h_fig = findobj(h_figs,'tag','wave_clus_figure');
h_fig1 = findobj(h_figs,'tag','wave_clus_aux');
h_fig2 = findobj(h_figs,'tag','wave_clus_aux1');
h_fig3 = findobj(h_figs,'tag','wave_clus_aux2');
h_fig4 = findobj(h_figs,'tag','wave_clus_aux3');
h_fig5 = findobj(h_figs,'tag','wave_clus_aux4');
set(handles.wave_clus_aux5,'userdata',USER_DATA);
set(h_fig,'userdata',USER_DATA)
set(h_fig1,'userdata',USER_DATA)
set(h_fig2,'userdata',USER_DATA)
set(h_fig3,'userdata',USER_DATA)
set(h_fig4,'userdata',USER_DATA)
set(h_fig5,'userdata',USER_DATA)



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% --- Executes during object creation, after setting all properties.
function isi29_nbins_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes during object creation, after setting all properties.
function isi29_bin_step_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes during object creation, after setting all properties.
function isi30_nbins_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes during object creation, after setting all properties.
function isi30_bin_step_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes during object creation, after setting all properties.
function isi31_nbins_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes during object creation, after setting all properties.
function isi31_bin_step_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes during object creation, after setting all properties.
function isi32_nbins_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes during object creation, after setting all properties.
function isi32_bin_step_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes during object creation, after setting all properties.
function isi33_nbins_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes during object creation, after setting all properties.
function isi33_bin_step_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
