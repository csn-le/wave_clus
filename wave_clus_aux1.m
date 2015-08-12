function varargout = wave_clus_aux1(varargin)
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
                   'gui_OpeningFcn', @wave_clus_aux1_OpeningFcn, ...
                   'gui_OutputFcn',  @wave_clus_aux1_OutputFcn, ...
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
function wave_clus_aux1_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to wave_clus_aux (see VARARGIN)

% Choose default command line output for wave_clus_aux
handles.output = hObject;
set(handles.isi9_accept_button,'value',1);
set(handles.isi10_accept_button,'value',1);
set(handles.isi11_accept_button,'value',1);
set(handles.isi12_accept_button,'value',1);
set(handles.isi13_accept_button,'value',1);
set(handles.fix9_button,'value',0);
set(handles.fix10_button,'value',0);
set(handles.fix11_button,'value',0);
set(handles.fix12_button,'value',0);
set(handles.fix13_button,'value',0);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes wave_clus_aux wait for user response (see UIRESUME)
% uiwait(handles.wave_clus_aux);


% --- Outputs from this function are returned to the command line.
function varargout = wave_clus_aux1_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

set(handles.isi9_accept_button,'value',1);
set(handles.isi10_accept_button,'value',1);
set(handles.isi11_accept_button,'value',1);
set(handles.isi12_accept_button,'value',1);
set(handles.isi13_accept_button,'value',1);
set(handles.isi9_reject_button,'value',0);
set(handles.isi10_reject_button,'value',0);
set(handles.isi11_reject_button,'value',0);
set(handles.isi12_reject_button,'value',0);
set(handles.isi13_reject_button,'value',0);
set(handles.fix9_button,'value',0);
set(handles.fix10_button,'value',0);
set(handles.fix11_button,'value',0);
set(handles.fix12_button,'value',0);
set(handles.fix13_button,'value',0);

h_figs=get(0,'children');
h_fig = findobj(h_figs,'tag','wave_clus_figure');
h_fig1 = findobj(h_figs,'tag','wave_clus_aux');
h_fig2 = findobj(h_figs,'tag','wave_clus_aux2');
h_fig3 = findobj(h_figs,'tag','wave_clus_aux3');
h_fig4 = findobj(h_figs,'tag','wave_clus_aux4');
h_fig5 = findobj(h_figs,'tag','wave_clus_aux5');
USER_DATA = get(h_fig,'UserData');
par = USER_DATA{1};

set(handles.isi9_nbins,'string',par.nbins9);
set(handles.isi10_nbins,'string',par.nbins10);
set(handles.isi11_nbins,'string',par.nbins11);
set(handles.isi12_nbins,'string',par.nbins12);
set(handles.isi13_nbins,'string',par.nbins13);
set(handles.isi9_bin_step,'string',par.bin_step9);
set(handles.isi10_bin_step,'string',par.bin_step10);
set(handles.isi11_bin_step,'string',par.bin_step11);
set(handles.isi12_bin_step,'string',par.bin_step12);
set(handles.isi13_bin_step,'string',par.bin_step13);

% That's for passing the fix button settings to plot_spikes.
if get(handles.fix9_button,'value') ==1     
    par.fix9 = 1;
else
    par.fix9 = 0;
end
if get(handles.fix10_button,'value') ==1     
    par.fix10 = 1;
else
    par.fix10 = 0;
end
if get(handles.fix11_button,'value') ==1     
    par.fix11 = 1;
else
    par.fix11 = 0;
end
if get(handles.fix12_button,'value') ==1     
    par.fix12 = 1;
else
    par.fix12 = 0;
end
if get(handles.fix13_button,'value') ==1     
    par.fix13 = 1;
else
    par.fix13 = 0;
end
USER_DATA{1} = par;
set(handles.wave_clus_aux1,'userdata',USER_DATA)
set(h_fig,'userdata',USER_DATA)
set(h_fig1,'userdata',USER_DATA)
set(h_fig2,'userdata',USER_DATA)
set(h_fig3,'userdata',USER_DATA)

plot_spikes_aux1(handles)




% Change nbins
% -------------------------------------------------------------
function isi9_nbins_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_aux1,'userdata');
par = USER_DATA{1};
par.nbins9 = str2num(get(hObject, 'String'));
par.axes_nr = 10;
classes = USER_DATA{6};
par.class_to_plot = find(classes==9);
USER_DATA{1} = par;
USER_DATA{6} = classes;
set(handles.wave_clus_aux1,'userdata',USER_DATA);
plot_spikes_aux1(handles)
% --------------------------------------------------------------------
function isi10_nbins_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_aux1,'userdata');
par = USER_DATA{1};
par.nbins10 = str2num(get(hObject, 'String'));
par.axes_nr = 11;
classes = USER_DATA{6};
par.class_to_plot = find(classes==10);
USER_DATA{1} = par;
USER_DATA{6} = classes;
set(handles.wave_clus_aux1,'userdata',USER_DATA);
plot_spikes_aux1(handles)
% --------------------------------------------------------------------
function isi11_nbins_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_aux1,'userdata');
par = USER_DATA{1};
par.nbins11 = str2num(get(hObject, 'String'));
par.axes_nr = 12;
classes = USER_DATA{6};
par.class_to_plot = find(classes==11);
USER_DATA{1} = par;
USER_DATA{6} = classes;
set(handles.wave_clus_aux1,'userdata',USER_DATA);
plot_spikes_aux1(handles)
% --------------------------------------------------------------------
function isi12_nbins_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_aux1,'userdata');
par = USER_DATA{1};
par.nbins12 = str2num(get(hObject, 'String'));
par.axes_nr = 13;
classes = USER_DATA{6};
par.class_to_plot = find(classes==12);
USER_DATA{1} = par;
USER_DATA{6} = classes;
set(handles.wave_clus_aux1,'userdata',USER_DATA);
plot_spikes_aux1(handles)
% --------------------------------------------------------------------
function isi13_nbins_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_aux1,'userdata');
par = USER_DATA{1};
par.nbins13 = str2num(get(hObject, 'String'));
par.axes_nr = 14;
classes = USER_DATA{6};
par.class_to_plot = find(classes==13);
USER_DATA{1} = par;
USER_DATA{6} = classes;
set(handles.wave_clus_aux1,'userdata',USER_DATA);
plot_spikes_aux1(handles)
% --------------------------------------------------------------------


% Change bin steps
% --------------------------------------------------------
function isi9_bin_step_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_aux1,'userdata');
par = USER_DATA{1};
par.bin_step9 = str2num(get(hObject, 'String'));
par.axes_nr = 10;
classes = USER_DATA{6};
par.class_to_plot = find(classes==9);
USER_DATA{1} = par;
USER_DATA{6} = classes;
set(handles.wave_clus_aux1,'userdata',USER_DATA);
plot_spikes_aux1(handles)
% --------------------------------------------------------------------
function isi10_bin_step_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_aux1,'userdata');
par = USER_DATA{1};
par.bin_step10 = str2num(get(hObject, 'String'));
par.axes_nr = 11;
classes = USER_DATA{6};
par.class_to_plot = find(classes==10);
USER_DATA{1} = par;
USER_DATA{6} = classes;
set(handles.wave_clus_aux1,'userdata',USER_DATA);
plot_spikes_aux1(handles)
% --------------------------------------------------------------------
function isi11_bin_step_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_aux1,'userdata');
par = USER_DATA{1};
par.bin_step11 = str2num(get(hObject, 'String'));
par.axes_nr = 12;
classes = USER_DATA{6};
par.class_to_plot = find(classes==11);
USER_DATA{1} = par;
USER_DATA{6} = classes;
set(handles.wave_clus_aux1,'userdata',USER_DATA);
plot_spikes_aux1(handles)
% --------------------------------------------------------------------
function isi12_bin_step_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_aux1,'userdata');
par = USER_DATA{1};
par.bin_step12 = str2num(get(hObject, 'String'));
par.axes_nr = 13;
classes = USER_DATA{6};
par.class_to_plot = find(classes==12);
USER_DATA{1} = par;
USER_DATA{6} = classes;
set(handles.wave_clus_aux1,'userdata',USER_DATA);
plot_spikes_aux1(handles)
% --------------------------------------------------------------------
function isi13_bin_step_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_aux1,'userdata');
par = USER_DATA{1};
par.bin_step13 = str2num(get(hObject, 'String'));
par.axes_nr = 14;
classes = USER_DATA{6};
par.class_to_plot = find(classes==13);
USER_DATA{1} = par;
USER_DATA{6} = classes;
set(handles.wave_clus_aux1,'userdata',USER_DATA);
plot_spikes_aux1(handles)
% --------------------------------------------------------------------


% Accept and Reject buttons
% --------------------------------------------------------
function isi9_accept_button_Callback(hObject, eventdata, handles)
set(gcbo,'value',1);
set(handles.isi9_reject_button,'value',0);
% --------------------------------------------------------------------
function isi9_reject_button_Callback(hObject, eventdata, handles)
set(gcbo,'value',1);
set(handles.isi9_accept_button,'value',0);
USER_DATA = get(handles.wave_clus_aux1,'userdata');
classes = USER_DATA{6};
classes(find(classes==9))=0;
USER_DATA{6} = classes;
USER_DATA{9} = classes;
h_figs=get(0,'children');
h_fig = findobj(h_figs,'tag','wave_clus_figure');
h_fig1 = findobj(h_figs,'tag','wave_clus_aux');
h_fig2 = findobj(h_figs,'tag','wave_clus_aux2');
h_fig3 = findobj(h_figs,'tag','wave_clus_aux3');
set(handles.wave_clus_aux1,'userdata',USER_DATA);
set(h_fig,'userdata',USER_DATA)
set(h_fig1,'userdata',USER_DATA)
set(h_fig2,'userdata',USER_DATA)
set(h_fig3,'userdata',USER_DATA)
axes(handles.spikes9); 
cla reset
axes(handles.isi9); 
cla reset
set(gcbo,'value',0);
set(handles.isi9_accept_button,'value',1);

% --------------------------------------------------------
function isi10_accept_button_Callback(hObject, eventdata, handles)
set(gcbo,'value',1);
set(handles.isi10_reject_button,'value',0);
% --------------------------------------------------------------------
function isi10_reject_button_Callback(hObject, eventdata, handles)
set(gcbo,'value',1);
set(handles.isi10_accept_button,'value',0);
USER_DATA = get(handles.wave_clus_aux1,'userdata');
classes = USER_DATA{6};
classes(find(classes==10))=0;
USER_DATA{6} = classes;
USER_DATA{9} = classes;
h_figs=get(0,'children');
h_fig = findobj(h_figs,'tag','wave_clus_figure');
h_fig1 = findobj(h_figs,'tag','wave_clus_aux');
h_fig2 = findobj(h_figs,'tag','wave_clus_aux2');
h_fig3 = findobj(h_figs,'tag','wave_clus_aux3');
set(handles.wave_clus_aux1,'userdata',USER_DATA);
set(h_fig,'userdata',USER_DATA)
set(h_fig1,'userdata',USER_DATA)
set(h_fig2,'userdata',USER_DATA)
set(h_fig3,'userdata',USER_DATA)
axes(handles.spikes10); 
cla reset
axes(handles.isi10); 
cla reset
set(gcbo,'value',0);
set(handles.isi10_accept_button,'value',1);

% --------------------------------------------------------
function isi11_accept_button_Callback(hObject, eventdata, handles)
set(gcbo,'value',1);
set(handles.isi11_reject_button,'value',0);
% --------------------------------------------------------------------
function isi11_reject_button_Callback(hObject, eventdata, handles)
set(gcbo,'value',1);
set(handles.isi11_accept_button,'value',0);
USER_DATA = get(handles.wave_clus_aux1,'userdata');
classes = USER_DATA{6};
classes(find(classes==11))=0;
USER_DATA{6} = classes;
USER_DATA{9} = classes;
h_figs=get(0,'children');
h_fig = findobj(h_figs,'tag','wave_clus_figure');
h_fig1 = findobj(h_figs,'tag','wave_clus_aux');
h_fig2 = findobj(h_figs,'tag','wave_clus_aux2');
h_fig3 = findobj(h_figs,'tag','wave_clus_aux3');
set(handles.wave_clus_aux1,'userdata',USER_DATA);
set(h_fig,'userdata',USER_DATA)
set(h_fig1,'userdata',USER_DATA)
set(h_fig2,'userdata',USER_DATA)
set(h_fig3,'userdata',USER_DATA)
axes(handles.spikes11); 
cla reset
axes(handles.isi11); 
cla reset
set(gcbo,'value',0);
set(handles.isi11_accept_button,'value',1);

% --------------------------------------------------------
function isi12_accept_button_Callback(hObject, eventdata, handles)
set(gcbo,'value',1);
set(handles.isi12_reject_button,'value',0);
% --------------------------------------------------------------------
function isi12_reject_button_Callback(hObject, eventdata, handles)
set(gcbo,'value',1);
set(handles.isi12_accept_button,'value',0);
USER_DATA = get(handles.wave_clus_aux1,'userdata');
classes = USER_DATA{6};
classes(find(classes==12))=0;
USER_DATA{6} = classes;
USER_DATA{9} = classes;
h_figs=get(0,'children');
h_fig = findobj(h_figs,'tag','wave_clus_figure');
h_fig1 = findobj(h_figs,'tag','wave_clus_aux');
h_fig2 = findobj(h_figs,'tag','wave_clus_aux2');
h_fig3 = findobj(h_figs,'tag','wave_clus_aux3');
set(handles.wave_clus_aux1,'userdata',USER_DATA);
set(h_fig,'userdata',USER_DATA)
set(h_fig1,'userdata',USER_DATA)
set(h_fig2,'userdata',USER_DATA)
set(h_fig3,'userdata',USER_DATA)
axes(handles.spikes12); 
cla reset
axes(handles.isi12); 
cla reset
set(gcbo,'value',0);
set(handles.isi12_accept_button,'value',1);

% --------------------------------------------------------
function isi13_accept_button_Callback(hObject, eventdata, handles)
set(gcbo,'value',1);
set(handles.isi13_reject_button,'value',0);
% --------------------------------------------------------------------
function isi13_reject_button_Callback(hObject, eventdata, handles)
set(gcbo,'value',1);
set(handles.isi13_accept_button,'value',0);
USER_DATA = get(handles.wave_clus_aux1,'userdata');
classes = USER_DATA{6};
classes(find(classes==13))=0;
USER_DATA{6} = classes;
USER_DATA{9} = classes;
h_figs=get(0,'children');
h_fig = findobj(h_figs,'tag','wave_clus_figure');
h_fig1 = findobj(h_figs,'tag','wave_clus_aux');
h_fig2 = findobj(h_figs,'tag','wave_clus_aux2');
h_fig3 = findobj(h_figs,'tag','wave_clus_aux3');
set(handles.wave_clus_aux1,'userdata',USER_DATA);
set(h_fig,'userdata',USER_DATA)
set(h_fig1,'userdata',USER_DATA)
set(h_fig2,'userdata',USER_DATA)
set(h_fig3,'userdata',USER_DATA)
axes(handles.spikes13); 
cla reset
axes(handles.isi13); 
cla reset
set(gcbo,'value',0);
set(handles.isi13_accept_button,'value',1);



% FIX buttons
% --------------------------------------------------------
function fix9_button_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_aux1,'userdata');
par = USER_DATA{1};
classes = USER_DATA{6};
fix_class = find(classes==9);
if get(handles.fix9_button,'value') ==1
    USER_DATA{28} = fix_class;
    par.fix9 = 1;
else
    USER_DATA{28} = [];
    par.fix9 = 0
end
USER_DATA{1} = par;
h_figs=get(0,'children');
h_fig = findobj(h_figs,'tag','wave_clus_figure');
h_fig1 = findobj(h_figs,'tag','wave_clus_aux');
h_fig2 = findobj(h_figs,'tag','wave_clus_aux2');
h_fig3 = findobj(h_figs,'tag','wave_clus_aux3');
set(handles.wave_clus_aux1,'userdata',USER_DATA);
set(h_fig,'userdata',USER_DATA)
set(h_fig1,'userdata',USER_DATA)
set(h_fig2,'userdata',USER_DATA)
set(h_fig3,'userdata',USER_DATA)
% --------------------------------------------------------
function fix10_button_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_aux1,'userdata');
par = USER_DATA{1};
classes = USER_DATA{6};
fix_class = find(classes==10);
if get(handles.fix10_button,'value') ==1
    USER_DATA{29} = fix_class;
    par.fix10 = 1;
else
    USER_DATA{29} = [];
    par.fix10 = 0;
end
USER_DATA{1} = par;
h_figs=get(0,'children');
h_fig = findobj(h_figs,'tag','wave_clus_figure');
h_fig1 = findobj(h_figs,'tag','wave_clus_aux');
h_fig2 = findobj(h_figs,'tag','wave_clus_aux2');
h_fig3 = findobj(h_figs,'tag','wave_clus_aux3');
set(handles.wave_clus_aux1,'userdata',USER_DATA);
set(h_fig,'userdata',USER_DATA)
set(h_fig1,'userdata',USER_DATA)
set(h_fig2,'userdata',USER_DATA)
set(h_fig3,'userdata',USER_DATA)
% --------------------------------------------------------
function fix11_button_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_aux1,'userdata');
par = USER_DATA{1};
classes = USER_DATA{6};
fix_class = find(classes==11);
if get(handles.fix11_button,'value') ==1
    USER_DATA{30} = fix_class;
    par.fix11 = 1;
else
    USER_DATA{30} = [];
    par.fix11 = 0;
end
USER_DATA{1} = par;
h_figs=get(0,'children');
h_fig = findobj(h_figs,'tag','wave_clus_figure');
h_fig1 = findobj(h_figs,'tag','wave_clus_aux');
h_fig2 = findobj(h_figs,'tag','wave_clus_aux2');
h_fig3 = findobj(h_figs,'tag','wave_clus_aux3');
set(handles.wave_clus_aux1,'userdata',USER_DATA);
set(h_fig,'userdata',USER_DATA)
set(h_fig1,'userdata',USER_DATA)
set(h_fig2,'userdata',USER_DATA)
set(h_fig3,'userdata',USER_DATA)
% --------------------------------------------------------
function fix12_button_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_aux1,'userdata');
par = USER_DATA{1};
classes = USER_DATA{6};
fix_class = find(classes==12);
if get(handles.fix12_button,'value') ==1
    USER_DATA{31} = fix_class;
    par.fix12 = 1;
else
    USER_DATA{31} = [];
    par.fix12 = 0;
end
USER_DATA{1} = par;
h_figs=get(0,'children');
h_fig = findobj(h_figs,'tag','wave_clus_figure');
h_fig1 = findobj(h_figs,'tag','wave_clus_aux');
h_fig2 = findobj(h_figs,'tag','wave_clus_aux2');
h_fig3 = findobj(h_figs,'tag','wave_clus_aux3');
set(handles.wave_clus_aux1,'userdata',USER_DATA);
set(h_fig,'userdata',USER_DATA)
set(h_fig1,'userdata',USER_DATA)
set(h_fig2,'userdata',USER_DATA)
set(h_fig3,'userdata',USER_DATA)
% --------------------------------------------------------
function fix13_button_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_aux1,'userdata');
par = USER_DATA{1};
classes = USER_DATA{6};
fix_class = find(classes==13);
if get(handles.fix13_button,'value') ==1
    USER_DATA{32} = fix_class;
    par.fix13 = 1;
else
    USER_DATA{32} = [];
    par.fix13 = 0;
end
USER_DATA{1} = par;
h_figs=get(0,'children');
h_fig = findobj(h_figs,'tag','wave_clus_figure');
h_fig1 = findobj(h_figs,'tag','wave_clus_aux');
h_fig2 = findobj(h_figs,'tag','wave_clus_aux2');
h_fig3 = findobj(h_figs,'tag','wave_clus_aux3');
set(handles.wave_clus_aux1,'userdata',USER_DATA);
set(h_fig,'userdata',USER_DATA)
set(h_fig1,'userdata',USER_DATA)
set(h_fig2,'userdata',USER_DATA)
set(h_fig3,'userdata',USER_DATA)



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% --- Executes during object creation, after setting all properties.
function isi9_nbins_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes during object creation, after setting all properties.
function isi9_bin_step_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes during object creation, after setting all properties.
function isi10_nbins_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes during object creation, after setting all properties.
function isi10_bin_step_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes during object creation, after setting all properties.
function isi11_nbins_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes during object creation, after setting all properties.
function isi11_bin_step_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes during object creation, after setting all properties.
function isi12_nbins_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes during object creation, after setting all properties.
function isi12_bin_step_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes during object creation, after setting all properties.
function isi13_nbins_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes during object creation, after setting all properties.
function isi13_bin_step_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
