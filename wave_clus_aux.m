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

% Edit the above text to modify the response to help wave_clus_aux

% Last Modified by GUIDE v2.5 16-Dec-2004 18:37:28

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @wave_clus_aux_OpeningFcn, ...
                   'gui_OutputFcn',  @wave_clus_aux_OutputFcn, ...
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

set(handles.isi4_accept_button,'value',1);
set(handles.isi5_accept_button,'value',1);
set(handles.isi6_accept_button,'value',1);
set(handles.isi7_accept_button,'value',1);
set(handles.isi8_accept_button,'value',1);
set(handles.isi4_reject_button,'value',0);
set(handles.isi5_reject_button,'value',0);
set(handles.isi6_reject_button,'value',0);
set(handles.isi7_reject_button,'value',0);
set(handles.isi8_reject_button,'value',0);
set(handles.fix4_button,'value',0);
set(handles.fix5_button,'value',0);
set(handles.fix6_button,'value',0);
set(handles.fix7_button,'value',0);
set(handles.fix8_button,'value',0);

h_figs=get(0,'children');
h_fig = findobj(h_figs,'tag','wave_clus_figure');
h_fig1 = findobj(h_figs,'tag','wave_clus_aux1');
h_fig2 = findobj(h_figs,'tag','wave_clus_aux2');
h_fig3 = findobj(h_figs,'tag','wave_clus_aux3');
h_fig4 = findobj(h_figs,'tag','wave_clus_aux4');
h_fig5 = findobj(h_figs,'tag','wave_clus_aux5');
USER_DATA = get(h_fig,'UserData');
par = USER_DATA{1};

set(handles.isi4_nbins,'string',par.nbins4);
set(handles.isi5_nbins,'string',par.nbins5);
set(handles.isi6_nbins,'string',par.nbins6);
set(handles.isi7_nbins,'string',par.nbins7);
set(handles.isi8_nbins,'string',par.nbins8);
set(handles.isi4_bin_step,'string',par.bin_step4);
set(handles.isi5_bin_step,'string',par.bin_step5);
set(handles.isi6_bin_step,'string',par.bin_step6);
set(handles.isi7_bin_step,'string',par.bin_step7);
set(handles.isi8_bin_step,'string',par.bin_step8);

% That's for passing the fix button settings to plot_spikes.
if get(handles.fix4_button,'value') ==1     
    par.fix4 = 1;
else
    par.fix4 = 0;
end
if get(handles.fix5_button,'value') ==1     
    par.fix5 = 1;
else
    par.fix5 = 0;
end
if get(handles.fix6_button,'value') ==1     
    par.fix6 = 1;
else
    par.fix6 = 0;
end
if get(handles.fix7_button,'value') ==1     
    par.fix7 = 1;
else
    par.fix7 = 0;
end
if get(handles.fix8_button,'value') ==1     
    par.fix8 = 1;
else
    par.fix8 = 0;
end
USER_DATA{1} = par;
set(handles.wave_clus_aux,'userdata',USER_DATA)
set(h_fig,'userdata',USER_DATA)
set(h_fig1,'userdata',USER_DATA)
set(h_fig,'userdata',USER_DATA)
set(h_fig1,'userdata',USER_DATA)

plot_spikes_aux(handles)




% Change nbins
% -------------------------------------------------------------
function isi4_nbins_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_aux,'userdata');
par = USER_DATA{1};
par.nbins4 = str2num(get(hObject, 'String'));
par.axes_nr = 5;
classes = USER_DATA{6};
par.class_to_plot = find(classes==4);
USER_DATA{1} = par;
USER_DATA{6} = classes;
set(handles.wave_clus_aux,'userdata',USER_DATA);
plot_spikes_aux(handles)
% --------------------------------------------------------------------
function isi5_nbins_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_aux,'userdata');
par = USER_DATA{1};
par.nbins5 = str2num(get(hObject, 'String'));
par.axes_nr = 6;
classes = USER_DATA{6};
par.class_to_plot = find(classes==5);
USER_DATA{1} = par;
USER_DATA{6} = classes;
set(handles.wave_clus_aux,'userdata',USER_DATA);
plot_spikes_aux(handles)
% --------------------------------------------------------------------
function isi6_nbins_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_aux,'userdata');
par = USER_DATA{1};
par.nbins6 = str2num(get(hObject, 'String'));
par.axes_nr = 7;
classes = USER_DATA{6};
par.class_to_plot = find(classes==6);
USER_DATA{1} = par;
USER_DATA{6} = classes;
set(handles.wave_clus_aux,'userdata',USER_DATA);
plot_spikes_aux(handles)
% --------------------------------------------------------------------
function isi7_nbins_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_aux,'userdata');
par = USER_DATA{1};
par.nbins7 = str2num(get(hObject, 'String'));
par.axes_nr = 8;
classes = USER_DATA{6};
par.class_to_plot = find(classes==7);
USER_DATA{1} = par;
USER_DATA{6} = classes;
set(handles.wave_clus_aux,'userdata',USER_DATA);
plot_spikes_aux(handles)
% --------------------------------------------------------------------
function isi8_nbins_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_aux,'userdata');
par = USER_DATA{1};
par.nbins8 = str2num(get(hObject, 'String'));
par.axes_nr = 9;
classes = USER_DATA{6};
par.class_to_plot = find(classes==8);
USER_DATA{1} = par;
USER_DATA{6} = classes;
set(handles.wave_clus_aux,'userdata',USER_DATA);
plot_spikes_aux(handles)
% --------------------------------------------------------------------


% Change bin steps
% --------------------------------------------------------
function isi4_bin_step_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_aux,'userdata');
par = USER_DATA{1};
par.bin_step4 = str2num(get(hObject, 'String'));
par.axes_nr = 5;
classes = USER_DATA{6};
par.class_to_plot = find(classes==4);
USER_DATA{1} = par;
USER_DATA{6} = classes;
set(handles.wave_clus_aux,'userdata',USER_DATA);
plot_spikes_aux(handles)
% --------------------------------------------------------------------
function isi5_bin_step_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_aux,'userdata');
par = USER_DATA{1};
par.bin_step5 = str2num(get(hObject, 'String'));
par.axes_nr = 6;
classes = USER_DATA{6};
par.class_to_plot = find(classes==5);
USER_DATA{1} = par;
USER_DATA{6} = classes;
set(handles.wave_clus_aux,'userdata',USER_DATA);
plot_spikes_aux(handles)
% --------------------------------------------------------------------
function isi6_bin_step_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_aux,'userdata');
par = USER_DATA{1};
par.bin_step6 = str2num(get(hObject, 'String'));
par.axes_nr = 7;
classes = USER_DATA{6};
par.class_to_plot = find(classes==6);
USER_DATA{1} = par;
USER_DATA{6} = classes;
set(handles.wave_clus_aux,'userdata',USER_DATA);
plot_spikes_aux(handles)
% --------------------------------------------------------------------
function isi7_bin_step_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_aux,'userdata');
par = USER_DATA{1};
par.bin_step7 = str2num(get(hObject, 'String'));
par.axes_nr = 8;
classes = USER_DATA{6};
par.class_to_plot = find(classes==7);
USER_DATA{1} = par;
USER_DATA{6} = classes;
set(handles.wave_clus_aux,'userdata',USER_DATA);
plot_spikes_aux(handles)
% --------------------------------------------------------------------
function isi8_bin_step_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_aux,'userdata');
par = USER_DATA{1};
par.bin_step8 = str2num(get(hObject, 'String'));
par.axes_nr = 9;
classes = USER_DATA{6};
par.class_to_plot = find(classes==8);
USER_DATA{1} = par;
USER_DATA{6} = classes;
set(handles.wave_clus_aux,'userdata',USER_DATA);
plot_spikes_aux(handles)
% --------------------------------------------------------------------


% Accept and Reject buttons
% --------------------------------------------------------
function isi4_accept_button_Callback(hObject, eventdata, handles)
set(gcbo,'value',1);
set(handles.isi4_reject_button,'value',0);
% --------------------------------------------------------------------
function isi4_reject_button_Callback(hObject, eventdata, handles)
set(gcbo,'value',1);
set(handles.isi4_accept_button,'value',0);
USER_DATA = get(handles.wave_clus_aux,'userdata');
classes = USER_DATA{6};
classes(find(classes==4))=0;
USER_DATA{6} = classes;
USER_DATA{9} = classes;
h_figs=get(0,'children');
h_fig = findobj(h_figs,'tag','wave_clus_figure');
h_fig1 = findobj(h_figs,'tag','wave_clus_aux1');
h_fig2 = findobj(h_figs,'tag','wave_clus_aux2');
h_fig3 = findobj(h_figs,'tag','wave_clus_aux3');
set(handles.wave_clus_aux,'userdata',USER_DATA);
set(h_fig,'userdata',USER_DATA)
set(h_fig1,'userdata',USER_DATA)
set(h_fig2,'userdata',USER_DATA)
set(h_fig3,'userdata',USER_DATA)
axes(handles.spikes4); 
cla reset
axes(handles.isi4); 
cla reset
set(gcbo,'value',0);
set(handles.isi4_accept_button,'value',1);

% --------------------------------------------------------
function isi5_accept_button_Callback(hObject, eventdata, handles)
set(gcbo,'value',1);
set(handles.isi5_reject_button,'value',0);
% --------------------------------------------------------------------
function isi5_reject_button_Callback(hObject, eventdata, handles)
set(gcbo,'value',1);
set(handles.isi5_accept_button,'value',0);
USER_DATA = get(handles.wave_clus_aux,'userdata');
classes = USER_DATA{6};
classes(find(classes==5))=0;
USER_DATA{6} = classes;
USER_DATA{9} = classes;
h_figs=get(0,'children');
h_fig = findobj(h_figs,'tag','wave_clus_figure');
h_fig1 = findobj(h_figs,'tag','wave_clus_aux1');
h_fig2 = findobj(h_figs,'tag','wave_clus_aux2');
h_fig3 = findobj(h_figs,'tag','wave_clus_aux3');
set(handles.wave_clus_aux,'userdata',USER_DATA);
set(h_fig,'userdata',USER_DATA)
set(h_fig1,'userdata',USER_DATA)
set(h_fig2,'userdata',USER_DATA)
set(h_fig3,'userdata',USER_DATA)
axes(handles.spikes5); 
cla reset
axes(handles.isi5); 
cla reset
set(gcbo,'value',0);
set(handles.isi5_accept_button,'value',1);

% --------------------------------------------------------
function isi6_accept_button_Callback(hObject, eventdata, handles)
set(gcbo,'value',1);
set(handles.isi6_reject_button,'value',0);
% --------------------------------------------------------------------
function isi6_reject_button_Callback(hObject, eventdata, handles)
set(gcbo,'value',1);
set(handles.isi6_accept_button,'value',0);
USER_DATA = get(handles.wave_clus_aux,'userdata');
classes = USER_DATA{6};
classes(find(classes==6))=0;
USER_DATA{6} = classes;
USER_DATA{9} = classes;
h_figs=get(0,'children');
h_fig = findobj(h_figs,'tag','wave_clus_figure');
h_fig1 = findobj(h_figs,'tag','wave_clus_aux1');
h_fig2 = findobj(h_figs,'tag','wave_clus_aux2');
h_fig3 = findobj(h_figs,'tag','wave_clus_aux3');
set(handles.wave_clus_aux,'userdata',USER_DATA);
set(h_fig,'userdata',USER_DATA)
set(h_fig1,'userdata',USER_DATA)
set(h_fig2,'userdata',USER_DATA)
set(h_fig3,'userdata',USER_DATA)
axes(handles.spikes6); 
cla reset
axes(handles.isi6); 
cla reset
set(gcbo,'value',0);
set(handles.isi6_accept_button,'value',1);

% --------------------------------------------------------
function isi7_accept_button_Callback(hObject, eventdata, handles)
set(gcbo,'value',1);
set(handles.isi7_reject_button,'value',0);
% --------------------------------------------------------------------
function isi7_reject_button_Callback(hObject, eventdata, handles)
set(gcbo,'value',1);
set(handles.isi7_accept_button,'value',0);
USER_DATA = get(handles.wave_clus_aux,'userdata');
classes = USER_DATA{6};
classes(find(classes==7))=0;
USER_DATA{6} = classes;
USER_DATA{9} = classes;
h_figs=get(0,'children');
h_fig = findobj(h_figs,'tag','wave_clus_figure');
h_fig1 = findobj(h_figs,'tag','wave_clus_aux1');
h_fig2 = findobj(h_figs,'tag','wave_clus_aux2');
h_fig3 = findobj(h_figs,'tag','wave_clus_aux3');
set(handles.wave_clus_aux,'userdata',USER_DATA);
set(h_fig,'userdata',USER_DATA)
set(h_fig1,'userdata',USER_DATA)
set(h_fig2,'userdata',USER_DATA)
set(h_fig3,'userdata',USER_DATA)
axes(handles.spikes7); 
cla reset
axes(handles.isi7); 
cla reset
set(gcbo,'value',0);
set(handles.isi7_accept_button,'value',1);

% --------------------------------------------------------
function isi8_accept_button_Callback(hObject, eventdata, handles)
set(gcbo,'value',1);
set(handles.isi8_reject_button,'value',0);
% --------------------------------------------------------------------
function isi8_reject_button_Callback(hObject, eventdata, handles)
set(gcbo,'value',1);
set(handles.isi8_accept_button,'value',0);
USER_DATA = get(handles.wave_clus_aux,'userdata');
classes = USER_DATA{6};
classes(find(classes==8))=0;
USER_DATA{6} = classes;
USER_DATA{9} = classes;
h_figs=get(0,'children');
h_fig = findobj(h_figs,'tag','wave_clus_figure');
h_fig1 = findobj(h_figs,'tag','wave_clus_aux1');
h_fig2 = findobj(h_figs,'tag','wave_clus_aux2');
h_fig3 = findobj(h_figs,'tag','wave_clus_aux3');
set(handles.wave_clus_aux,'userdata',USER_DATA);
set(h_fig,'userdata',USER_DATA)
set(h_fig1,'userdata',USER_DATA)
set(h_fig2,'userdata',USER_DATA)
set(h_fig3,'userdata',USER_DATA)
axes(handles.spikes8); 
cla reset
axes(handles.isi8); 
cla reset
set(gcbo,'value',0);
set(handles.isi8_accept_button,'value',1);



% FIX buttons
% --------------------------------------------------------
function fix4_button_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_aux,'userdata');
par = USER_DATA{1};
classes = USER_DATA{6};
fix_class = find(classes==4);
if get(handles.fix4_button,'value') ==1
    USER_DATA{23} = fix_class;
    par.fix4 = 1;
else
    USER_DATA{23} = [];
    par.fix4 = 0
end
USER_DATA{1} = par;
h_figs=get(0,'children');
h_fig = findobj(h_figs,'tag','wave_clus_figure');
h_fig1 = findobj(h_figs,'tag','wave_clus_aux1');
h_fig2 = findobj(h_figs,'tag','wave_clus_aux2');
h_fig3 = findobj(h_figs,'tag','wave_clus_aux3');
set(handles.wave_clus_aux,'userdata',USER_DATA);
set(h_fig,'userdata',USER_DATA)
set(h_fig1,'userdata',USER_DATA)
set(h_fig2,'userdata',USER_DATA)
set(h_fig3,'userdata',USER_DATA)
% --------------------------------------------------------
function fix5_button_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_aux,'userdata');
par = USER_DATA{1};
classes = USER_DATA{6};
fix_class = find(classes==5);
if get(handles.fix5_button,'value') ==1
    USER_DATA{24} = fix_class;
    par.fix5 = 1;
else
    USER_DATA{24} = [];
    par.fix5 = 0;
end
USER_DATA{1} = par;
h_figs=get(0,'children');
h_fig = findobj(h_figs,'tag','wave_clus_figure');
h_fig1 = findobj(h_figs,'tag','wave_clus_aux1');
h_fig2 = findobj(h_figs,'tag','wave_clus_aux2');
h_fig3 = findobj(h_figs,'tag','wave_clus_aux3');
set(handles.wave_clus_aux,'userdata',USER_DATA);
set(h_fig,'userdata',USER_DATA)
set(h_fig1,'userdata',USER_DATA)
set(h_fig2,'userdata',USER_DATA)
set(h_fig3,'userdata',USER_DATA)
% --------------------------------------------------------
function fix6_button_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_aux,'userdata');
par = USER_DATA{1};
classes = USER_DATA{6};
fix_class = find(classes==6);
if get(handles.fix6_button,'value') ==1
    USER_DATA{25} = fix_class;
    par.fix6 = 1;
else
    USER_DATA{25} = [];
    par.fix6 = 0;
end
USER_DATA{1} = par;
h_figs=get(0,'children');
h_fig = findobj(h_figs,'tag','wave_clus_figure');
h_fig1 = findobj(h_figs,'tag','wave_clus_aux1');
h_fig2 = findobj(h_figs,'tag','wave_clus_aux2');
h_fig3 = findobj(h_figs,'tag','wave_clus_aux3');
set(handles.wave_clus_aux,'userdata',USER_DATA);
set(h_fig,'userdata',USER_DATA)
set(h_fig1,'userdata',USER_DATA)
set(h_fig2,'userdata',USER_DATA)
set(h_fig3,'userdata',USER_DATA)
% --------------------------------------------------------
function fix7_button_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_aux,'userdata');
par = USER_DATA{1};
classes = USER_DATA{6};
fix_class = find(classes==7);
if get(handles.fix7_button,'value') ==1
    USER_DATA{26} = fix_class;
    par.fix7 = 1;
else
    USER_DATA{26} = [];
    par.fix7 = 0;
end
USER_DATA{1} = par;
h_figs=get(0,'children');
h_fig = findobj(h_figs,'tag','wave_clus_figure');
h_fig1 = findobj(h_figs,'tag','wave_clus_aux1');
h_fig2 = findobj(h_figs,'tag','wave_clus_aux2');
h_fig3 = findobj(h_figs,'tag','wave_clus_aux3');
set(handles.wave_clus_aux,'userdata',USER_DATA);
set(h_fig,'userdata',USER_DATA)
set(h_fig1,'userdata',USER_DATA)
set(h_fig2,'userdata',USER_DATA)
set(h_fig3,'userdata',USER_DATA)
% --------------------------------------------------------
function fix8_button_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_aux,'userdata');
par = USER_DATA{1};
classes = USER_DATA{6};
fix_class = find(classes==8);
if get(handles.fix8_button,'value') ==1
    USER_DATA{27} = fix_class;
    par.fix8 = 1;
else
    USER_DATA{27} = [];
    par.fix8 = 0;
end
USER_DATA{1} = par;
h_figs=get(0,'children');
h_fig = findobj(h_figs,'tag','wave_clus_figure');
h_fig1 = findobj(h_figs,'tag','wave_clus_aux1');
h_fig2 = findobj(h_figs,'tag','wave_clus_aux2');
h_fig3 = findobj(h_figs,'tag','wave_clus_aux3');
set(handles.wave_clus_aux,'userdata',USER_DATA);
set(h_fig,'userdata',USER_DATA)
set(h_fig1,'userdata',USER_DATA)
set(h_fig2,'userdata',USER_DATA)
set(h_fig3,'userdata',USER_DATA)



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% --- Executes during object creation, after setting all properties.
function isi4_nbins_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes during object creation, after setting all properties.
function isi4_bin_step_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes during object creation, after setting all properties.
function isi5_nbins_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes during object creation, after setting all properties.
function isi5_bin_step_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes during object creation, after setting all properties.
function isi6_nbins_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes during object creation, after setting all properties.
function isi6_bin_step_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes during object creation, after setting all properties.
function isi7_nbins_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes during object creation, after setting all properties.
function isi7_bin_step_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes during object creation, after setting all properties.
function isi8_nbins_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes during object creation, after setting all properties.
function isi8_bin_step_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
