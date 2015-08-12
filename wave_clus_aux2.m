function varargout = wave_clus_aux2(varargin)
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
                   'gui_OpeningFcn', @wave_clus_aux2_OpeningFcn, ...
                   'gui_OutputFcn',  @wave_clus_aux2_OutputFcn, ...
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
function wave_clus_aux2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to wave_clus_aux (see VARARGIN)

% Choose default command line output for wave_clus_aux
handles.output = hObject;
set(handles.isi14_accept_button,'value',1);
set(handles.isi15_accept_button,'value',1);
set(handles.isi16_accept_button,'value',1);
set(handles.isi17_accept_button,'value',1);
set(handles.isi18_accept_button,'value',1);
set(handles.fix14_button,'value',0);
set(handles.fix15_button,'value',0);
set(handles.fix16_button,'value',0);
set(handles.fix17_button,'value',0);
set(handles.fix18_button,'value',0);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes wave_clus_aux wait for user response (see UIRESUME)
% uiwait(handles.wave_clus_aux);


% --- Outputs from this function are returned to the command line.
function varargout = wave_clus_aux2_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

set(handles.isi14_accept_button,'value',1);
set(handles.isi15_accept_button,'value',1);
set(handles.isi16_accept_button,'value',1);
set(handles.isi17_accept_button,'value',1);
set(handles.isi18_accept_button,'value',1);
set(handles.isi14_reject_button,'value',0);
set(handles.isi15_reject_button,'value',0);
set(handles.isi16_reject_button,'value',0);
set(handles.isi17_reject_button,'value',0);
set(handles.isi18_reject_button,'value',0);
set(handles.fix14_button,'value',0);
set(handles.fix15_button,'value',0);
set(handles.fix16_button,'value',0);
set(handles.fix17_button,'value',0);
set(handles.fix18_button,'value',0);

h_figs=get(0,'children');
h_fig = findobj(h_figs,'tag','wave_clus_figure');
h_fig1 = findobj(h_figs,'tag','wave_clus_aux');
h_fig2 = findobj(h_figs,'tag','wave_clus_aux1');
h_fig3 = findobj(h_figs,'tag','wave_clus_aux3');
h_fig4 = findobj(h_figs,'tag','wave_clus_aux4');
h_fig5 = findobj(h_figs,'tag','wave_clus_aux5');
USER_DATA = get(h_fig,'UserData');
par = USER_DATA{1};

set(handles.isi14_nbins,'string',par.nbins14);
set(handles.isi15_nbins,'string',par.nbins15);
set(handles.isi16_nbins,'string',par.nbins16);
set(handles.isi17_nbins,'string',par.nbins17);
set(handles.isi18_nbins,'string',par.nbins18);
set(handles.isi14_bin_step,'string',par.bin_step14);
set(handles.isi15_bin_step,'string',par.bin_step15);
set(handles.isi16_bin_step,'string',par.bin_step16);
set(handles.isi17_bin_step,'string',par.bin_step17);
set(handles.isi18_bin_step,'string',par.bin_step18);

% That's for passing the fix button settings to plot_spikes.
if get(handles.fix14_button,'value') ==1     
    par.fix14 = 1;
else
    par.fix14 = 0;
end
if get(handles.fix15_button,'value') ==1     
    par.fix15 = 1;
else
    par.fix15 = 0;
end
if get(handles.fix16_button,'value') ==1     
    par.fix16 = 1;
else
    par.fix16 = 0;
end
if get(handles.fix17_button,'value') ==1     
    par.fix17 = 1;
else
    par.fix17 = 0;
end
if get(handles.fix18_button,'value') ==1     
    par.fix18 = 1;
else
    par.fix18 = 0;
end
USER_DATA{1} = par;
set(handles.wave_clus_aux2,'userdata',USER_DATA)
set(h_fig,'userdata',USER_DATA)
set(h_fig1,'userdata',USER_DATA)
set(h_fig2,'userdata',USER_DATA)
set(h_fig3,'userdata',USER_DATA)

plot_spikes_aux2(handles)




% Change nbins
% -------------------------------------------------------------
function isi14_nbins_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_aux2,'userdata');
par = USER_DATA{1};
par.nbins14 = str2num(get(hObject, 'String'));
par.axes_nr = 15;
classes = USER_DATA{6};
par.class_to_plot = find(classes==14);
USER_DATA{1} = par;
USER_DATA{6} = classes;
set(handles.wave_clus_aux2,'userdata',USER_DATA);
plot_spikes_aux2(handles)
% --------------------------------------------------------------------
function isi15_nbins_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_aux2,'userdata');
par = USER_DATA{1};
par.nbins15 = str2num(get(hObject, 'String'));
par.axes_nr = 16;
classes = USER_DATA{6};
par.class_to_plot = find(classes==15);
USER_DATA{1} = par;
USER_DATA{6} = classes;
set(handles.wave_clus_aux2,'userdata',USER_DATA);
plot_spikes_aux2(handles)
% --------------------------------------------------------------------
function isi16_nbins_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_aux2,'userdata');
par = USER_DATA{1};
par.nbins16 = str2num(get(hObject, 'String'));
par.axes_nr = 17;
classes = USER_DATA{6};
par.class_to_plot = find(classes==16);
USER_DATA{1} = par;
USER_DATA{6} = classes;
set(handles.wave_clus_aux2,'userdata',USER_DATA);
plot_spikes_aux2(handles)
% --------------------------------------------------------------------
function isi17_nbins_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_aux2,'userdata');
par = USER_DATA{1};
par.nbins17 = str2num(get(hObject, 'String'));
par.axes_nr = 18;
classes = USER_DATA{6};
par.class_to_plot = find(classes==17);
USER_DATA{1} = par;
USER_DATA{6} = classes;
set(handles.wave_clus_aux2,'userdata',USER_DATA);
plot_spikes_aux2(handles)
% --------------------------------------------------------------------
function isi18_nbins_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_aux2,'userdata');
par = USER_DATA{1};
par.nbins18 = str2num(get(hObject, 'String'));
par.axes_nr = 19;
classes = USER_DATA{6};
par.class_to_plot = find(classes==18);
USER_DATA{1} = par;
USER_DATA{6} = classes;
set(handles.wave_clus_aux2,'userdata',USER_DATA);
plot_spikes_aux2(handles)
% --------------------------------------------------------------------


% Change bin steps
% --------------------------------------------------------
function isi14_bin_step_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_aux2,'userdata');
par = USER_DATA{1};
par.bin_step14 = str2num(get(hObject, 'String'));
par.axes_nr = 15;
classes = USER_DATA{6};
par.class_to_plot = find(classes==14);
USER_DATA{1} = par;
USER_DATA{6} = classes;
set(handles.wave_clus_aux2,'userdata',USER_DATA);
plot_spikes_aux2(handles)
% --------------------------------------------------------------------
function isi15_bin_step_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_aux2,'userdata');
par = USER_DATA{1};
par.bin_step15 = str2num(get(hObject, 'String'));
par.axes_nr = 16;
classes = USER_DATA{6};
par.class_to_plot = find(classes==15);
USER_DATA{1} = par;
USER_DATA{6} = classes;
set(handles.wave_clus_aux2,'userdata',USER_DATA);
plot_spikes_aux2(handles)
% --------------------------------------------------------------------
function isi16_bin_step_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_aux2,'userdata');
par = USER_DATA{1};
par.bin_step16 = str2num(get(hObject, 'String'));
par.axes_nr = 17;
classes = USER_DATA{6};
par.class_to_plot = find(classes==16);
USER_DATA{1} = par;
USER_DATA{6} = classes;
set(handles.wave_clus_aux2,'userdata',USER_DATA);
plot_spikes_aux2(handles)
% --------------------------------------------------------------------
function isi17_bin_step_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_aux2,'userdata');
par = USER_DATA{1};
par.bin_step17 = str2num(get(hObject, 'String'));
par.axes_nr = 18;
classes = USER_DATA{6};
par.class_to_plot = find(classes==17);
USER_DATA{1} = par;
USER_DATA{6} = classes;
set(handles.wave_clus_aux2,'userdata',USER_DATA);
plot_spikes_aux2(handles)
% --------------------------------------------------------------------
function isi18_bin_step_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_aux2,'userdata');
par = USER_DATA{1};
par.bin_step18 = str2num(get(hObject, 'String'));
par.axes_nr = 19;
classes = USER_DATA{6};
par.class_to_plot = find(classes==18);
USER_DATA{1} = par;
USER_DATA{6} = classes;
set(handles.wave_clus_aux2,'userdata',USER_DATA);
plot_spikes_aux2(handles)
% --------------------------------------------------------------------


% Accept and Reject buttons
% --------------------------------------------------------
function isi14_accept_button_Callback(hObject, eventdata, handles)
set(gcbo,'value',1);
set(handles.isi14_reject_button,'value',0);
% --------------------------------------------------------------------
function isi14_reject_button_Callback(hObject, eventdata, handles)
set(gcbo,'value',1);
set(handles.isi14_accept_button,'value',0);
USER_DATA = get(handles.wave_clus_aux2,'userdata');
classes = USER_DATA{6};
classes(find(classes==14))=0;
USER_DATA{6} = classes;
USER_DATA{9} = classes;
h_figs=get(0,'children');
h_fig = findobj(h_figs,'tag','wave_clus_figure');
h_fig1 = findobj(h_figs,'tag','wave_clus_aux');
h_fig2 = findobj(h_figs,'tag','wave_clus_aux1');
h_fig3 = findobj(h_figs,'tag','wave_clus_aux3');
set(handles.wave_clus_aux2,'userdata',USER_DATA);
set(h_fig,'userdata',USER_DATA)
set(h_fig1,'userdata',USER_DATA)
set(h_fig2,'userdata',USER_DATA)
set(h_fig3,'userdata',USER_DATA)
axes(handles.spikes14); 
cla reset
axes(handles.isi14); 
cla reset
set(gcbo,'value',0);
set(handles.isi14_accept_button,'value',1);

% --------------------------------------------------------
function isi15_accept_button_Callback(hObject, eventdata, handles)
set(gcbo,'value',1);
set(handles.isi15_reject_button,'value',0);
% --------------------------------------------------------------------
function isi15_reject_button_Callback(hObject, eventdata, handles)
set(gcbo,'value',1);
set(handles.isi15_accept_button,'value',0);
USER_DATA = get(handles.wave_clus_aux2,'userdata');
classes = USER_DATA{6};
classes(find(classes==15))=0;
USER_DATA{6} = classes;
USER_DATA{9} = classes;
h_figs=get(0,'children');
h_fig = findobj(h_figs,'tag','wave_clus_figure');
h_fig1 = findobj(h_figs,'tag','wave_clus_aux');
h_fig2 = findobj(h_figs,'tag','wave_clus_aux1');
h_fig3 = findobj(h_figs,'tag','wave_clus_aux3');
set(handles.wave_clus_aux2,'userdata',USER_DATA);
set(h_fig,'userdata',USER_DATA)
set(h_fig1,'userdata',USER_DATA)
set(h_fig2,'userdata',USER_DATA)
set(h_fig3,'userdata',USER_DATA)
axes(handles.spikes15); 
cla reset
axes(handles.isi15); 
cla reset
set(gcbo,'value',0);
set(handles.isi15_accept_button,'value',1);

% --------------------------------------------------------
function isi16_accept_button_Callback(hObject, eventdata, handles)
set(gcbo,'value',1);
set(handles.isi16_reject_button,'value',0);
% --------------------------------------------------------------------
function isi16_reject_button_Callback(hObject, eventdata, handles)
set(gcbo,'value',1);
set(handles.isi16_accept_button,'value',0);
USER_DATA = get(handles.wave_clus_aux2,'userdata');
classes = USER_DATA{6};
classes(find(classes==16))=0;
USER_DATA{6} = classes;
USER_DATA{9} = classes;
h_figs=get(0,'children');
h_fig = findobj(h_figs,'tag','wave_clus_figure');
h_fig1 = findobj(h_figs,'tag','wave_clus_aux');
h_fig2 = findobj(h_figs,'tag','wave_clus_aux1');
h_fig3 = findobj(h_figs,'tag','wave_clus_aux3');
set(handles.wave_clus_aux2,'userdata',USER_DATA);
set(h_fig,'userdata',USER_DATA)
set(h_fig1,'userdata',USER_DATA)
set(h_fig2,'userdata',USER_DATA)
set(h_fig3,'userdata',USER_DATA)
axes(handles.spikes16); 
cla reset
axes(handles.isi16); 
cla reset
set(gcbo,'value',0);
set(handles.isi16_accept_button,'value',1);

% --------------------------------------------------------
function isi17_accept_button_Callback(hObject, eventdata, handles)
set(gcbo,'value',1);
set(handles.isi17_reject_button,'value',0);
% --------------------------------------------------------------------
function isi17_reject_button_Callback(hObject, eventdata, handles)
set(gcbo,'value',1);
set(handles.isi17_accept_button,'value',0);
USER_DATA = get(handles.wave_clus_aux2,'userdata');
classes = USER_DATA{6};
classes(find(classes==17))=0;
USER_DATA{6} = classes;
USER_DATA{9} = classes;
h_figs=get(0,'children');
h_fig = findobj(h_figs,'tag','wave_clus_figure');
h_fig1 = findobj(h_figs,'tag','wave_clus_aux');
h_fig2 = findobj(h_figs,'tag','wave_clus_aux1');
h_fig3 = findobj(h_figs,'tag','wave_clus_aux3');
set(handles.wave_clus_aux2,'userdata',USER_DATA);
set(h_fig,'userdata',USER_DATA)
set(h_fig1,'userdata',USER_DATA)
set(h_fig2,'userdata',USER_DATA)
set(h_fig3,'userdata',USER_DATA)
axes(handles.spikes17); 
cla reset
axes(handles.isi17); 
cla reset
set(gcbo,'value',0);
set(handles.isi17_accept_button,'value',1);

% --------------------------------------------------------
function isi18_accept_button_Callback(hObject, eventdata, handles)
set(gcbo,'value',1);
set(handles.isi18_reject_button,'value',0);
% --------------------------------------------------------------------
function isi18_reject_button_Callback(hObject, eventdata, handles)
set(gcbo,'value',1);
set(handles.isi18_accept_button,'value',0);
USER_DATA = get(handles.wave_clus_aux2,'userdata');
classes = USER_DATA{6};
classes(find(classes==18))=0;
USER_DATA{6} = classes;
USER_DATA{9} = classes;
h_figs=get(0,'children');
h_fig = findobj(h_figs,'tag','wave_clus_figure');
h_fig1 = findobj(h_figs,'tag','wave_clus_aux');
h_fig2 = findobj(h_figs,'tag','wave_clus_aux1');
h_fig3 = findobj(h_figs,'tag','wave_clus_aux3');
set(handles.wave_clus_aux2,'userdata',USER_DATA);
set(h_fig,'userdata',USER_DATA)
set(h_fig1,'userdata',USER_DATA)
set(h_fig2,'userdata',USER_DATA)
set(h_fig3,'userdata',USER_DATA)
axes(handles.spikes18); 
cla reset
axes(handles.isi18); 
cla reset
set(gcbo,'value',0);
set(handles.isi18_accept_button,'value',1);



% FIX buttons
% --------------------------------------------------------
function fix14_button_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_aux2,'userdata');
par = USER_DATA{1};
classes = USER_DATA{6};
fix_class = find(classes==14);
if get(handles.fix14_button,'value') ==1
    USER_DATA{33} = fix_class;
    par.fix14 = 1;
else
    USER_DATA{33} = [];
    par.fix14 = 0
end
USER_DATA{1} = par;
h_figs=get(0,'children');
h_fig = findobj(h_figs,'tag','wave_clus_figure');
h_fig1 = findobj(h_figs,'tag','wave_clus_aux');
h_fig2 = findobj(h_figs,'tag','wave_clus_aux1');
h_fig3 = findobj(h_figs,'tag','wave_clus_aux3');
set(handles.wave_clus_aux2,'userdata',USER_DATA);
set(h_fig,'userdata',USER_DATA)
set(h_fig1,'userdata',USER_DATA)
set(h_fig2,'userdata',USER_DATA)
set(h_fig3,'userdata',USER_DATA)
% --------------------------------------------------------
function fix15_button_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_aux2,'userdata');
par = USER_DATA{1};
classes = USER_DATA{6};
fix_class = find(classes==15);
if get(handles.fix15_button,'value') ==1
    USER_DATA{34} = fix_class;
    par.fix15 = 1;
else
    USER_DATA{34} = [];
    par.fix15 = 0;
end
USER_DATA{1} = par;
h_figs=get(0,'children');
h_fig = findobj(h_figs,'tag','wave_clus_figure');
h_fig1 = findobj(h_figs,'tag','wave_clus_aux');
h_fig2 = findobj(h_figs,'tag','wave_clus_aux1');
h_fig3 = findobj(h_figs,'tag','wave_clus_aux3');
set(handles.wave_clus_aux2,'userdata',USER_DATA);
set(h_fig,'userdata',USER_DATA)
set(h_fig1,'userdata',USER_DATA)
set(h_fig2,'userdata',USER_DATA)
set(h_fig3,'userdata',USER_DATA)
% --------------------------------------------------------
function fix16_button_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_aux2,'userdata');
par = USER_DATA{1};
classes = USER_DATA{6};
fix_class = find(classes==16);
if get(handles.fix16_button,'value') ==1
    USER_DATA{35} = fix_class;
    par.fix16 = 1;
else
    USER_DATA{35} = [];
    par.fix16 = 0;
end
USER_DATA{1} = par;
h_figs=get(0,'children');
h_fig = findobj(h_figs,'tag','wave_clus_figure');
h_fig1 = findobj(h_figs,'tag','wave_clus_aux');
h_fig2 = findobj(h_figs,'tag','wave_clus_aux1');
h_fig3 = findobj(h_figs,'tag','wave_clus_aux3');
set(handles.wave_clus_aux2,'userdata',USER_DATA);
set(h_fig,'userdata',USER_DATA)
set(h_fig1,'userdata',USER_DATA)
set(h_fig2,'userdata',USER_DATA)
set(h_fig3,'userdata',USER_DATA)
% --------------------------------------------------------
function fix17_button_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_aux2,'userdata');
par = USER_DATA{1};
classes = USER_DATA{6};
fix_class = find(classes==17);
if get(handles.fix17_button,'value') ==1
    USER_DATA{36} = fix_class;
    par.fix17 = 1;
else
    USER_DATA{36} = [];
    par.fix17 = 0;
end
USER_DATA{1} = par;
h_figs=get(0,'children');
h_fig = findobj(h_figs,'tag','wave_clus_figure');
h_fig1 = findobj(h_figs,'tag','wave_clus_aux');
h_fig2 = findobj(h_figs,'tag','wave_clus_aux1');
h_fig3 = findobj(h_figs,'tag','wave_clus_aux3');
set(handles.wave_clus_aux2,'userdata',USER_DATA);
set(h_fig,'userdata',USER_DATA)
set(h_fig1,'userdata',USER_DATA)
set(h_fig2,'userdata',USER_DATA)
set(h_fig3,'userdata',USER_DATA)
% --------------------------------------------------------
function fix18_button_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_aux2,'userdata');
par = USER_DATA{1};
classes = USER_DATA{6};
fix_class = find(classes==18);
if get(handles.fix18_button,'value') ==1
    USER_DATA{37} = fix_class;
    par.fix18 = 1;
else
    USER_DATA{37} = [];
    par.fix18 = 0;
end
USER_DATA{1} = par;
h_figs=get(0,'children');
h_fig = findobj(h_figs,'tag','wave_clus_figure');
h_fig1 = findobj(h_figs,'tag','wave_clus_aux');
h_fig2 = findobj(h_figs,'tag','wave_clus_aux1');
h_fig3 = findobj(h_figs,'tag','wave_clus_aux3');
set(handles.wave_clus_aux2,'userdata',USER_DATA);
set(h_fig,'userdata',USER_DATA)
set(h_fig1,'userdata',USER_DATA)
set(h_fig2,'userdata',USER_DATA)
set(h_fig3,'userdata',USER_DATA)



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% --- Executes during object creation, after setting all properties.
function isi14_nbins_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes during object creation, after setting all properties.
function isi14_bin_step_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes during object creation, after setting all properties.
function isi15_nbins_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes during object creation, after setting all properties.
function isi15_bin_step_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes during object creation, after setting all properties.
function isi16_nbins_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes during object creation, after setting all properties.
function isi16_bin_step_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes during object creation, after setting all properties.
function isi17_nbins_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes during object creation, after setting all properties.
function isi17_bin_step_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes during object creation, after setting all properties.
function isi18_nbins_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes during object creation, after setting all properties.
function isi18_bin_step_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
