function varargout = set_parameters_ui(varargin)
% SET_PARAMETERS_UI MATLAB code for set_parameters_ui.fig
%      SET_PARAMETERS_UI, by itself, creates a new SET_PARAMETERS_UI or raises the existing
%      singleton*.
%
%      H = SET_PARAMETERS_UI returns the handle to a new SET_PARAMETERS_UI or the handle to
%      the existing singleton*.
%
%      SET_PARAMETERS_UI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SET_PARAMETERS_UI.M with the given input arguments.
%
%      SET_PARAMETERS_UI('Property','Value',...) creates a new SET_PARAMETERS_UI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before set_parameters_ui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to set_parameters_ui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help set_parameters_ui

% Last Modified by GUIDE v2.5 05-Oct-2015 18:05:52

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @set_parameters_ui_OpeningFcn, ...
                   'gui_OutputFcn',  @set_parameters_ui_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before set_parameters_ui is made visible.
function set_parameters_ui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to set_parameters_ui (see VARARGIN)

% Choose default command line output for set_parameters_ui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes set_parameters_ui wait for user response (see UIRESUME)
% uiwait(handles.set_parameters_ui);
% this_folder = dir();
% if ~ismember('set_parameters.m',{this_folder.name})
%      copyfile(, [pwd filesep 'set_parameters.m']);
% end
%edit([pwd filesep 'set_parameters.m'])
text = fileread([fileparts(mfilename('fullpath')) filesep 'set_parameters.m']);

new_lines = strfind(text, sprintf('\n'));
text = text(new_lines(1)+1:end);
set(handles.text_editor,'string',text);


% --- Outputs from this function are returned to the command line.
function varargout = set_parameters_ui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



% --- Executes on button press in default.
function default_Callback(hObject, eventdata, handles)
% hObject    handle to default (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
text = fileread([fileparts(mfilename('fullpath')) filesep 'set_parameters_DEFAULT.m']);
new_lines = strfind(text, sprintf('\n'));
text = text(new_lines(1)+1:end);
set(handles.text_editor,'string',text);


function text_editor_Callback(hObject, eventdata, handles)
% hObject    handle to text_editor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_editor as text
%        str2double(get(hObject,'String')) returns contents of text_editor as a double


% --- Executes during object creation, after setting all properties.
function text_editor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_editor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject, 'Max', 2); %// Enable multi-line string input to the editbox



% --- Executes on button press in save.
function save_Callback(hObject, eventdata, handles)
% hObject    handle to save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
text = get(handles.text_editor, 'String');
fout = fopen([fileparts(mfilename('fullpath')) filesep 'set_parameters.m'],'w');

fprintf(fout, 'function par = set_parameters() \n');
for row = 1:size(text,1)
    fprintf(fout, '%s \n', strtrim(text(row,1:end)));
end

fclose(fout);
