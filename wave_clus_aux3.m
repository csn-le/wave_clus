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

% Last Modified by GUIDE v2.5 16-Dec-2004 18:37:28

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @wave_clus_aux3_OpeningFcn, ...
                   'gui_OutputFcn',  @wave_clus_aux3_OutputFcn, ...
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

set(handles.isi19_accept_button,'value',1);
set(handles.isi20_accept_button,'value',1);
set(handles.isi21_accept_button,'value',1);
set(handles.isi22_accept_button,'value',1);
set(handles.isi23_accept_button,'value',1);
set(handles.isi19_reject_button,'value',0);
set(handles.isi20_reject_button,'value',0);
set(handles.isi21_reject_button,'value',0);
set(handles.isi22_reject_button,'value',0);
set(handles.isi23_reject_button,'value',0);
set(handles.fix19_button,'value',0);
set(handles.fix20_button,'value',0);
set(handles.fix21_button,'value',0);
set(handles.fix22_button,'value',0);
set(handles.fix23_button,'value',0);

h_figs=get(0,'children');
h_fig = findobj(h_figs,'tag','wave_clus_figure');
h_fig1 = findobj(h_figs,'tag','wave_clus_aux');
h_fig2 = findobj(h_figs,'tag','wave_clus_aux1');
h_fig3 = findobj(h_figs,'tag','wave_clus_aux2');
h_fig4 = findobj(h_figs,'tag','wave_clus_aux4');
h_fig5 = findobj(h_figs,'tag','wave_clus_aux5');
USER_DATA = get(h_fig,'UserData');
par = USER_DATA{1};

set(handles.isi19_nbins,'string',par.nbins9);
set(handles.isi20_nbins,'string',par.nbins10);
set(handles.isi21_nbins,'string',par.nbins11);
set(handles.isi22_nbins,'string',par.nbins12);
set(handles.isi23_nbins,'string',par.nbins13);
set(handles.isi19_bin_step,'string',par.bin_step9);
set(handles.isi20_bin_step,'string',par.bin_step10);
set(handles.isi21_bin_step,'string',par.bin_step11);
set(handles.isi22_bin_step,'string',par.bin_step12);
set(handles.isi23_bin_step,'string',par.bin_step13);

% That's for passing the fix button settings to plot_spikes.
if get(handles.fix19_button,'value') ==1     
    par.fix19 = 1;
else
    par.fix19 = 0;
end
if get(handles.fix20_button,'value') ==1     
    par.fix20 = 1;
else
    par.fix20 = 0;
end
if get(handles.fix21_button,'value') ==1     
    par.fix21 = 1;
else
    par.fix21 = 0;
end
if get(handles.fix22_button,'value') ==1     
    par.fix22 = 1;
else
    par.fix22 = 0;
end
if get(handles.fix23_button,'value') ==1     
    par.fix23 = 1;
else
    par.fix23 = 0;
end
USER_DATA{1} = par;
set(handles.wave_clus_aux3,'userdata',USER_DATA)
set(h_fig,'userdata',USER_DATA)
set(h_fig1,'userdata',USER_DATA)
set(h_fig2,'userdata',USER_DATA)
set(h_fig3,'userdata',USER_DATA)

plot_spikes_aux3(handles)




% Change nbins
% -------------------------------------------------------------
function isi19_nbins_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_aux3,'userdata');
par = USER_DATA{1};
par.nbins19 = str2num(get(hObject, 'String'));
par.axes_nr = 20;
classes = USER_DATA{6};
par.class_to_plot = find(classes==19);
USER_DATA{1} = par;
USER_DATA{6} = classes;
set(handles.wave_clus_aux3,'userdata',USER_DATA);
plot_spikes_aux3(handles)
% --------------------------------------------------------------------
function isi20_nbins_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_aux3,'userdata');
par = USER_DATA{1};
par.nbins20 = str2num(get(hObject, 'String'));
par.axes_nr = 21;
classes = USER_DATA{6};
par.class_to_plot = find(classes==20);
USER_DATA{1} = par;
USER_DATA{6} = classes;
set(handles.wave_clus_aux3,'userdata',USER_DATA);
plot_spikes_aux3(handles)
% --------------------------------------------------------------------
function isi21_nbins_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_aux3,'userdata');
par = USER_DATA{1};
par.nbins21 = str2num(get(hObject, 'String'));
par.axes_nr = 22;
classes = USER_DATA{6};
par.class_to_plot = find(classes==21);
USER_DATA{1} = par;
USER_DATA{6} = classes;
set(handles.wave_clus_aux3,'userdata',USER_DATA);
plot_spikes_aux3(handles)
% --------------------------------------------------------------------
function isi22_nbins_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_aux3,'userdata');
par = USER_DATA{1};
par.nbins22 = str2num(get(hObject, 'String'));
par.axes_nr = 23;
classes = USER_DATA{6};
par.class_to_plot = find(classes==22);
USER_DATA{1} = par;
USER_DATA{6} = classes;
set(handles.wave_clus_aux3,'userdata',USER_DATA);
plot_spikes_aux3(handles)
% --------------------------------------------------------------------
function isi23_nbins_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_aux3,'userdata');
par = USER_DATA{1};
par.nbins23 = str2num(get(hObject, 'String'));
par.axes_nr = 24;
classes = USER_DATA{6};
par.class_to_plot = find(classes==23);
USER_DATA{1} = par;
USER_DATA{6} = classes;
set(handles.wave_clus_aux3,'userdata',USER_DATA);
plot_spikes_aux3(handles)
% --------------------------------------------------------------------


% Change bin steps
% --------------------------------------------------------
function isi19_bin_step_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_aux3,'userdata');
par = USER_DATA{1};
par.bin_step19 = str2num(get(hObject, 'String'));
par.axes_nr = 20;
classes = USER_DATA{6};
par.class_to_plot = find(classes==19);
USER_DATA{1} = par;
USER_DATA{6} = classes;
set(handles.wave_clus_aux3,'userdata',USER_DATA);
plot_spikes_aux3(handles)
% --------------------------------------------------------------------
function isi20_bin_step_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_aux3,'userdata');
par = USER_DATA{1};
par.bin_step20 = str2num(get(hObject, 'String'));
par.axes_nr = 21;
classes = USER_DATA{6};
par.class_to_plot = find(classes==20);
USER_DATA{1} = par;
USER_DATA{6} = classes;
set(handles.wave_clus_aux3,'userdata',USER_DATA);
plot_spikes_aux3(handles)
% --------------------------------------------------------------------
function isi21_bin_step_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_aux3,'userdata');
par = USER_DATA{1};
par.bin_step21 = str2num(get(hObject, 'String'));
par.axes_nr = 22;
classes = USER_DATA{6};
par.class_to_plot = find(classes==21);
USER_DATA{1} = par;
USER_DATA{6} = classes;
set(handles.wave_clus_aux3,'userdata',USER_DATA);
plot_spikes_aux3(handles)
% --------------------------------------------------------------------
function isi22_bin_step_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_aux3,'userdata');
par = USER_DATA{1};
par.bin_step22 = str2num(get(hObject, 'String'));
par.axes_nr = 23;
classes = USER_DATA{6};
par.class_to_plot = find(classes==22);
USER_DATA{1} = par;
USER_DATA{6} = classes;
set(handles.wave_clus_aux3,'userdata',USER_DATA);
plot_spikes_aux3(handles)
% --------------------------------------------------------------------
function isi23_bin_step_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_aux3,'userdata');
par = USER_DATA{1};
par.bin_step23 = str2num(get(hObject, 'String'));
par.axes_nr = 24;
classes = USER_DATA{6};
par.class_to_plot = find(classes==23);
USER_DATA{1} = par;
USER_DATA{6} = classes;
set(handles.wave_clus_aux3,'userdata',USER_DATA);
plot_spikes_aux3(handles)
% --------------------------------------------------------------------


% Accept and Reject buttons
% --------------------------------------------------------
function isi19_accept_button_Callback(hObject, eventdata, handles)
set(gcbo,'value',1);
set(handles.isi19_reject_button,'value',0);
% --------------------------------------------------------------------
function isi19_reject_button_Callback(hObject, eventdata, handles)
set(gcbo,'value',1);
set(handles.isi19_accept_button,'value',0);
USER_DATA = get(handles.wave_clus_aux3,'userdata');
classes = USER_DATA{6};
classes(find(classes==19))=0;
USER_DATA{6} = classes;
USER_DATA{9} = classes;
h_figs=get(0,'children');
h_fig = findobj(h_figs,'tag','wave_clus_figure');
h_fig1 = findobj(h_figs,'tag','wave_clus_aux');
h_fig2 = findobj(h_figs,'tag','wave_clus_aux1');
h_fig3 = findobj(h_figs,'tag','wave_clus_aux2');
h_fig4 = findobj(h_figs,'tag','wave_clus_aux4');
h_fig5 = findobj(h_figs,'tag','wave_clus_aux5');
set(handles.wave_clus_aux3,'userdata',USER_DATA);
set(h_fig,'userdata',USER_DATA)
set(h_fig1,'userdata',USER_DATA)
set(h_fig2,'userdata',USER_DATA)
set(h_fig3,'userdata',USER_DATA)
set(h_fig4,'userdata',USER_DATA)
set(h_fig5,'userdata',USER_DATA)
axes(handles.spikes19); 
cla reset
axes(handles.isi19); 
cla reset
set(gcbo,'value',0);
set(handles.isi19_accept_button,'value',1);

% --------------------------------------------------------
function isi20_accept_button_Callback(hObject, eventdata, handles)
set(gcbo,'value',1);
set(handles.isi20_reject_button,'value',0);
% --------------------------------------------------------------------
function isi20_reject_button_Callback(hObject, eventdata, handles)
set(gcbo,'value',1);
set(handles.isi20_accept_button,'value',0);
USER_DATA = get(handles.wave_clus_aux3,'userdata');
classes = USER_DATA{6};
classes(find(classes==20))=0;
USER_DATA{6} = classes;
USER_DATA{9} = classes;
h_figs=get(0,'children');
h_fig = findobj(h_figs,'tag','wave_clus_figure');
h_fig1 = findobj(h_figs,'tag','wave_clus_aux');
h_fig2 = findobj(h_figs,'tag','wave_clus_aux1');
h_fig3 = findobj(h_figs,'tag','wave_clus_aux2');
h_fig4 = findobj(h_figs,'tag','wave_clus_aux4');
h_fig5 = findobj(h_figs,'tag','wave_clus_aux5');
set(handles.wave_clus_aux3,'userdata',USER_DATA);
set(h_fig,'userdata',USER_DATA)
set(h_fig1,'userdata',USER_DATA)
set(h_fig2,'userdata',USER_DATA)
set(h_fig3,'userdata',USER_DATA)
set(h_fig4,'userdata',USER_DATA)
set(h_fig5,'userdata',USER_DATA)
axes(handles.spikes20); 
cla reset
axes(handles.isi20); 
cla reset
set(gcbo,'value',0);
set(handles.isi20_accept_button,'value',1);

% --------------------------------------------------------
function isi21_accept_button_Callback(hObject, eventdata, handles)
set(gcbo,'value',1);
set(handles.isi21_reject_button,'value',0);
% --------------------------------------------------------------------
function isi21_reject_button_Callback(hObject, eventdata, handles)
set(gcbo,'value',1);
set(handles.isi21_accept_button,'value',0);
USER_DATA = get(handles.wave_clus_aux3,'userdata');
classes = USER_DATA{6};
classes(find(classes==16))=0;
USER_DATA{6} = classes;
USER_DATA{9} = classes;
h_figs=get(0,'children');
h_fig = findobj(h_figs,'tag','wave_clus_figure');
h_fig1 = findobj(h_figs,'tag','wave_clus_aux');
h_fig2 = findobj(h_figs,'tag','wave_clus_aux1');
h_fig3 = findobj(h_figs,'tag','wave_clus_aux2');
h_fig4 = findobj(h_figs,'tag','wave_clus_aux4');
h_fig5 = findobj(h_figs,'tag','wave_clus_aux5');
set(handles.wave_clus_aux3,'userdata',USER_DATA);
set(h_fig,'userdata',USER_DATA)
set(h_fig1,'userdata',USER_DATA)
set(h_fig2,'userdata',USER_DATA)
set(h_fig3,'userdata',USER_DATA)
set(h_fig4,'userdata',USER_DATA)
set(h_fig5,'userdata',USER_DATA)
axes(handles.spikes21); 
cla reset
axes(handles.isi21); 
cla reset
set(gcbo,'value',0);
set(handles.isi21_accept_button,'value',1);

% --------------------------------------------------------
function isi22_accept_button_Callback(hObject, eventdata, handles)
set(gcbo,'value',1);
set(handles.isi22_reject_button,'value',0);
% --------------------------------------------------------------------
function isi22_reject_button_Callback(hObject, eventdata, handles)
set(gcbo,'value',1);
set(handles.isi22_accept_button,'value',0);
USER_DATA = get(handles.wave_clus_aux3,'userdata');
classes = USER_DATA{6};
classes(find(classes==22))=0;
USER_DATA{6} = classes;
USER_DATA{9} = classes;
h_figs=get(0,'children');
h_fig = findobj(h_figs,'tag','wave_clus_figure');
h_fig1 = findobj(h_figs,'tag','wave_clus_aux');
h_fig2 = findobj(h_figs,'tag','wave_clus_aux1');
h_fig3 = findobj(h_figs,'tag','wave_clus_aux2');
set(handles.wave_clus_aux3,'userdata',USER_DATA);
set(h_fig,'userdata',USER_DATA)
set(h_fig1,'userdata',USER_DATA)
set(h_fig2,'userdata',USER_DATA)
set(h_fig3,'userdata',USER_DATA)
axes(handles.spikes22); 
cla reset
axes(handles.isi22); 
cla reset
set(gcbo,'value',0);
set(handles.isi22_accept_button,'value',1);

% --------------------------------------------------------
function isi23_accept_button_Callback(hObject, eventdata, handles)
set(gcbo,'value',1);
set(handles.isi23_reject_button,'value',0);
% --------------------------------------------------------------------
function isi23_reject_button_Callback(hObject, eventdata, handles)
set(gcbo,'value',1);
set(handles.isi23_accept_button,'value',0);
USER_DATA = get(handles.wave_clus_aux3,'userdata');
classes = USER_DATA{6};
classes(find(classes==23))=0;
USER_DATA{6} = classes;
USER_DATA{9} = classes;
h_figs=get(0,'children');
h_fig = findobj(h_figs,'tag','wave_clus_figure');
h_fig1 = findobj(h_figs,'tag','wave_clus_aux');
h_fig2 = findobj(h_figs,'tag','wave_clus_aux1');
h_fig3 = findobj(h_figs,'tag','wave_clus_aux2');
h_fig4 = findobj(h_figs,'tag','wave_clus_aux4');
h_fig5 = findobj(h_figs,'tag','wave_clus_aux5');
set(handles.wave_clus_aux3,'userdata',USER_DATA);
set(h_fig,'userdata',USER_DATA)
set(h_fig1,'userdata',USER_DATA)
set(h_fig2,'userdata',USER_DATA)
set(h_fig3,'userdata',USER_DATA)
set(h_fig4,'userdata',USER_DATA)
set(h_fig5,'userdata',USER_DATA)
axes(handles.spikes23); 
cla reset
axes(handles.isi23); 
cla reset
set(gcbo,'value',0);
set(handles.isi23_accept_button,'value',1);



% FIX buttons
% --------------------------------------------------------
function fix19_button_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_aux3,'userdata');
par = USER_DATA{1};
classes = USER_DATA{6};
fix_class = find(classes==19);
if get(handles.fix19_button,'value') ==1
    USER_DATA{38} = fix_class;
    par.fix19 = 1;
else
    USER_DATA{38} = [];
    par.fix19 = 0
end
USER_DATA{1} = par;
h_figs=get(0,'children');
h_fig = findobj(h_figs,'tag','wave_clus_figure');
h_fig1 = findobj(h_figs,'tag','wave_clus_aux');
h_fig2 = findobj(h_figs,'tag','wave_clus_aux1');
h_fig3 = findobj(h_figs,'tag','wave_clus_aux2');
set(handles.wave_clus_aux3,'userdata',USER_DATA);
set(h_fig,'userdata',USER_DATA)
set(h_fig1,'userdata',USER_DATA)
set(h_fig2,'userdata',USER_DATA)
set(h_fig3,'userdata',USER_DATA)
% --------------------------------------------------------
function fix20_button_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_aux3,'userdata');
par = USER_DATA{1};
classes = USER_DATA{6};
fix_class = find(classes==20);
if get(handles.fix20_button,'value') ==1
    USER_DATA{39} = fix_class;
    par.fix20 = 1;
else
    USER_DATA{39} = [];
    par.fix20 = 0;
end
USER_DATA{1} = par;
h_figs=get(0,'children');
h_fig = findobj(h_figs,'tag','wave_clus_figure');
h_fig1 = findobj(h_figs,'tag','wave_clus_aux');
h_fig2 = findobj(h_figs,'tag','wave_clus_aux1');
h_fig3 = findobj(h_figs,'tag','wave_clus_aux2');
h_fig4 = findobj(h_figs,'tag','wave_clus_aux4');
h_fig5 = findobj(h_figs,'tag','wave_clus_aux5');
set(handles.wave_clus_aux3,'userdata',USER_DATA);
set(h_fig,'userdata',USER_DATA)
set(h_fig1,'userdata',USER_DATA)
set(h_fig2,'userdata',USER_DATA)
set(h_fig3,'userdata',USER_DATA)
set(h_fig4,'userdata',USER_DATA)
set(h_fig5,'userdata',USER_DATA)
% --------------------------------------------------------
function fix21_button_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_aux3,'userdata');
par = USER_DATA{1};
classes = USER_DATA{6};
fix_class = find(classes==21);
if get(handles.fix21_button,'value') ==1
    USER_DATA{40} = fix_class;
    par.fix21 = 1;
else
    USER_DATA{40} = [];
    par.fix21 = 0;
end
USER_DATA{1} = par;
h_figs=get(0,'children');
h_fig = findobj(h_figs,'tag','wave_clus_figure');
h_fig1 = findobj(h_figs,'tag','wave_clus_aux');
h_fig2 = findobj(h_figs,'tag','wave_clus_aux1');
h_fig3 = findobj(h_figs,'tag','wave_clus_aux2');
h_fig4 = findobj(h_figs,'tag','wave_clus_aux4');
h_fig5 = findobj(h_figs,'tag','wave_clus_aux5');
set(handles.wave_clus_aux3,'userdata',USER_DATA);
set(h_fig,'userdata',USER_DATA)
set(h_fig1,'userdata',USER_DATA)
set(h_fig2,'userdata',USER_DATA)
set(h_fig3,'userdata',USER_DATA)
set(h_fig4,'userdata',USER_DATA)
set(h_fig5,'userdata',USER_DATA)
% --------------------------------------------------------
function fix22_button_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_aux3,'userdata');
par = USER_DATA{1};
classes = USER_DATA{6};
fix_class = find(classes==22);
if get(handles.fix22_button,'value') ==1
    USER_DATA{41} = fix_class;
    par.fix22 = 1;
else
    USER_DATA{41} = [];
    par.fix22 = 0;
end
USER_DATA{1} = par;
h_figs=get(0,'children');
h_fig = findobj(h_figs,'tag','wave_clus_figure');
h_fig1 = findobj(h_figs,'tag','wave_clus_aux');
h_fig2 = findobj(h_figs,'tag','wave_clus_aux1');
h_fig3 = findobj(h_figs,'tag','wave_clus_aux2');
h_fig4 = findobj(h_figs,'tag','wave_clus_aux4');
h_fig5 = findobj(h_figs,'tag','wave_clus_aux5');
set(handles.wave_clus_aux3,'userdata',USER_DATA);
set(h_fig,'userdata',USER_DATA)
set(h_fig1,'userdata',USER_DATA)
set(h_fig2,'userdata',USER_DATA)
set(h_fig3,'userdata',USER_DATA)
set(h_fig4,'userdata',USER_DATA)
set(h_fig5,'userdata',USER_DATA)
% --------------------------------------------------------
function fix23_button_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_aux3,'userdata');
par = USER_DATA{1};
classes = USER_DATA{6};
fix_class = find(classes==23);
if get(handles.fix23_button,'value') ==1
    USER_DATA{42} = fix_class;
    par.fix23 = 1;
else
    USER_DATA{42} = [];
    par.fix23 = 0;
end
USER_DATA{1} = par;
h_figs=get(0,'children');
h_fig = findobj(h_figs,'tag','wave_clus_figure');
h_fig1 = findobj(h_figs,'tag','wave_clus_aux');
h_fig2 = findobj(h_figs,'tag','wave_clus_aux1');
h_fig3 = findobj(h_figs,'tag','wave_clus_aux2');
h_fig4 = findobj(h_figs,'tag','wave_clus_aux4');
h_fig5 = findobj(h_figs,'tag','wave_clus_aux5');
set(handles.wave_clus_aux3,'userdata',USER_DATA);
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
function isi19_nbins_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes during object creation, after setting all properties.
function isi19_bin_step_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes during object creation, after setting all properties.
function isi20_nbins_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes during object creation, after setting all properties.
function isi20_bin_step_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes during object creation, after setting all properties.
function isi21_nbins_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes during object creation, after setting all properties.
function isi21_bin_step_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes during object creation, after setting all properties.
function isi22_nbins_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes during object creation, after setting all properties.
function isi22_bin_step_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes during object creation, after setting all properties.
function isi23_nbins_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes during object creation, after setting all properties.
function isi23_bin_step_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
