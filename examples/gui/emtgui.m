function varargout = emtgui(varargin)
% EMTGUI MATLAB code for emtgui.fig
%      EMTGUI, by itself, creates a new EMTGUI or raises the existing
%      singleton*.
%
%      H = EMTGUI returns the handle to a new EMTGUI or the handle to
%      the existing singleton*.
%
%      EMTGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EMTGUI.M with the given input arguments.
%
%      EMTGUI('Property','Value',...) creates a new EMTGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before emtgui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to emtgui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help emtgui

% Last Modified by GUIDE v2.5 15-Jul-2022 11:43:58

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @emtgui_OpeningFcn, ...
                   'gui_OutputFcn',  @emtgui_OutputFcn, ...
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


% --- Executes just before emtgui is made visible.
function emtgui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to emtgui (see VARARGIN)

% Choose default command line output for emtgui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes emtgui wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% Clear axes
set(handles.axes,'YTickLabel',[]);
set(handles.axes,'XTickLabel',[]);
set(handles.axes,'XTick',[])
set(handles.axes,'YTick',[])

% --- Outputs from this function are returned to the command line.
function varargout = emtgui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in btn_loadnetwork.
function btn_loadnetwork_Callback(hObject, eventdata, handles)
% hObject    handle to btn_loadnetwork (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, filepath] = uigetfile({'*.inp'},'Load EPANET Inp File.');
if ~isnumeric(filepath)
    cla( handles.axes )

    d=epanet([filepath, filename]);
    h = d.plot('axes', handles.axes, 'legend', 'hide');

    set(handles.networkname,'String', ['File name: "', filename, '"']);
    set(handles.linkcount,'String', ['Link count: ', num2str(d.LinkCount)]);
    set(handles.nodecount,'String', ['Node count: ', num2str(d.NodeCount)]);

    d.unload % Unload network and libraries.
end


% --- Executes when figure1 is resized.
function figure1_SizeChangedFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
