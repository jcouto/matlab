function varargout = browseKlustersDAT(varargin)
% BROWSEKLUSTERSDAT MATLAB code for browseKlustersDAT.fig
%      BROWSEKLUSTERSDAT, by itself, creates a new BROWSEKLUSTERSDAT or raises the existing
%      singleton*.
%
%      H = BROWSEKLUSTERSDAT returns the handle to a new BROWSEKLUSTERSDAT or the handle to
%      the existing singleton*.
%
%      BROWSEKLUSTERSDAT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BROWSEKLUSTERSDAT.M with the given input arguments.
%
%      BROWSEKLUSTERSDAT('Property','Value',...) creates a new BROWSEKLUSTERSDAT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before browseKlustersDAT_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to browseKlustersDAT_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help browseKlustersDAT

% Last Modified by GUIDE v2.5 19-Feb-2015 13:10:58

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @browseKlustersDAT_OpeningFcn, ...
    'gui_OutputFcn',  @browseKlustersDAT_OutputFcn, ...
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

% --- Executes just before browseKlustersDAT is made visible.
function browseKlustersDAT_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to browseKlustersDAT (see VARARGIN)

% Choose default command line output for browseKlustersDAT
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
% Load file and information

handles.info = load('01_run_gratings.info.mat');
info = handles.info;
handles.mfile = memmapfile('01_run_gratings.dat','Format', ...
    {'int16' [info.nchannels info.nsamples ] 'data'}, 'Repeat', 1, 'Writable', false);
handles.timebin = 0.25; % seconds
handles.voltageOffset = 200; % seconds

handles.timeoffset = 0; %seconds
handles.timemax = info.nsamples./info.srate;
handles.timebin_samples = handles.timebin.*handles.info.srate;
guidata(hObject,handles);
set(handles.timeSlider,'Min',0,'Max',handles.timemax,'Value',handles.timeoffset);
% UIWAIT makes browseKlustersDAT wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = browseKlustersDAT_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function timeSlider_Callback(hObject, eventdata, handles)
% hObject    handle to timeSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles = guidata(hObject);
get(hObject,'Value')
handles.timeoffset = ceil(get(hObject,'Value'));
handles = updatePlot(handles);


% --- Executes during object creation, after setting all properties.
function timeSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to timeSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function voltageSlider_Callback(hObject, eventdata, handles)
% hObject    handle to voltageSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function voltageSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to voltageSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function voltageOffset_Callback(hObject, eventdata, handles)
% hObject    handle to voltageOffset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of voltageOffset as text
%        str2double(get(hObject,'String')) returns contents of voltageOffset as a double


% --- Executes during object creation, after setting all properties.
function voltageOffset_CreateFcn(hObject, eventdata, handles)
% hObject    handle to voltageOffset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in resetPlot.
function resetPlot_Callback(hObject, eventdata, handles)
% hObject    handle to resetPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(gcf,'toolbar','figure');
handles = initialPlot(handles);
guidata(hObject,handles);

function handles = initialPlot(handles)
set(handles.timeSlider,'Min',0,'Max',handles.timemax,'Value',handles.timeoffset);

axes(handles.mainPlot);
cla;hold on
info = handles.info;
idx = handles.timeoffset*info.srate +[1:handles.timebin_samples];
handles.trace = {};
set(handles.timetext,'String',sprintf('%.2fs',handles.timeoffset))
for ch = 1:info.nchannels
    handles.trace{ch} = plot(handles.info.convertToFloat(handles.mfile.Data.data(ch,idx))+((ch-1)*handles.voltageOffset),'k');
end
axis tight

function handles = updatePlot(handles)
info = handles.info;
idx = handles.timeoffset*info.srate +[1:handles.timebin_samples];
if sum(idx>1 & idx<handles.timemax.*handles.info.srate) == length(idx)
    set(handles.timetext,'String',sprintf('%.2fs',handles.timeoffset))
    for ch = 1:info.nchannels
        
        set(handles.trace{ch},'ydata',handles.info.convertToFloat(handles.mfile.Data.data(ch,idx))+((ch-1)*handles.voltageOffset));
    end
else
    disp('Window moved outside of plot.')
end

% --- Executes during object creation, after setting all properties.
function mainPlot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mainPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate mainPlot



function timetext_Callback(hObject, eventdata, handles)
% hObject    handle to timetext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of timetext as text
%        str2double(get(hObject,'String')) returns contents of timetext as a double


% --- Executes during object creation, after setting all properties.
function timetext_CreateFcn(hObject, eventdata, handles)
% hObject    handle to timetext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when figure1 is resized.
function figure1_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
