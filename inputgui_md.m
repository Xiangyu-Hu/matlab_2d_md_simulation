function varargout = inputgui_md(varargin)
% INPUTGUI_MD M-file for inputgui_md.fig
%      INPUTGUI_MD, by itself, creates a new INPUTGUI_MD or raises the existing
%      singleton*.
%
%      H = INPUTGUI_MD returns the handle to a new INPUTGUI_MD or the handle to
%      the existing singleton*.
%
%      INPUTGUI_MD('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in INPUTGUI_MD.M with the given input arguments.
%
%      INPUTGUI_MD('Property','Value',...) creates a new INPUTGUI_MD or raises
%      the existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before inputgui_md_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to inputgui_md_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help inputgui_md

% Last Modified by GUIDE v2.5 30-Oct-2008 14:39:21

% Begin initialization code - DO NOT EDIT

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @inputgui_md_OpeningFcn, ...
                   'gui_OutputFcn',  @inputgui_md_OutputFcn, ...
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

% --- Executes just before inputgui_md is made visible.
function inputgui_md_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to inputgui_md (see VARARGIN)

% Choose default command line output for inputgui_md
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

initialize_gui(hObject, handles, false);
scrsz = get(0,'ScreenSize');
cursize=get(handles.figure1,'Position');
cursize(1)=0.05*scrsz(3);cursize(2)=.5*scrsz(4);
set(handles.figure1,'Position',cursize)

% UIWAIT makes inputgui_md wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = inputgui_md_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function np_CreateFcn(hObject, eventdata, handles)
% hObject    handle to np (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function np_Callback(hObject, eventdata, handles)
% hObject    handle to np (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of np as text
%        str2double(get(hObject,'String')) returns contents of np as a double
density = str2double(get(hObject, 'String'));
if isnan(density)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end

% Save the new np value
handles.metricdata.density = density;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function volume_CreateFcn(hObject, eventdata, handles)
% hObject    handle to volume (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
density = str2double(get(hObject, 'String'));
if isnan(density)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function volume_Callback(hObject, eventdata, handles)
% hObject    handle to volume (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of volume as text
%        str2double(get(hObject,'String')) returns contents of volume as a double
volume = str2double(get(hObject, 'String'));
if isnan(volume)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end

% Save the new volume value
guidata(hObject,handles);

% --- Executes on button press in calculate.
function calculate_Callback(hObject, eventdata, handles)
% hObject    handle to calculate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Disable inputs
set(handles.calculate,'Enable','off');
set(handles.nsteps,'Enable','off');
set(handles.temp,'Enable','off');
set(handles.density,'Enable','off');
set(handles.dt,'Enable','off');
set(handles.np,'Enable','off');
set(handles.saveloc,'Enable','off');

%Perform Monte Carlo Simulation!
%Perform Molecular Dynamics Simulation!
numpart=str2num(get(handles.np,'String'));
len=linspace(.1,0.9,ceil(numpart^(.5)));
[x,y]=meshgrid(len,len);
pos=[];
pos(:,1)=x(1:numpart);
pos(:,2)=y(1:numpart);
%pos(:,3)=z(1:numpart);
set(handles.stop,'Enable','on');
set(handles.pause,'Enable','on');

if get(handles.saveloc,'Value')==0
    hand=figure;
    %set(handles.saveloc,'Value',hand);
else
    hand=get(handles.saveloc,'Value');
    close(hand)
    hand=figure;
    set(handles.saveloc,'Value',hand);
end
set(handles.status,'String','MD Working','BackgroundColor',get(handles.uipanel1,'BackgroundColor'));
set(handles.statuspanel,'BackgroundColor',get(handles.uipanel1,'BackgroundColor'));
ff=fopen('simdata','w+');
fwrite(ff,0);
fclose(ff);
[rdfx rdfy flag]=md1(pos,[1 1],str2num(get(handles.nsteps,'String')),str2num(get(handles.dt,'String')),200,str2num(get(handles.density,'String')),str2num(get(handles.temp,'String')),hand,get(handles.saveloc,'String'));
%[rdfx rdfy flag]=md1(pos,[1 1 1],str2num(get(handles.nsteps,'String')),str2num(get(handles.dt,'String')),200,str2num(get(handles.density,'String')),str2num(get(handles.temp,'String')),hand,get(handles.saveloc,'String'));
set(handles.stop,'Enable','off');
set(handles.pause,'Enable','off');
set(handles.calculate,'Enable','on');
set(handles.nsteps,'Enable','on');
set(handles.temp,'Enable','on');
set(handles.density,'Enable','on');
set(handles.dt,'Enable','on');
set(handles.np,'Enable','on');
set(handles.saveloc,'Enable','on');
if flag
set(handles.status,'String','Stopped','BackgroundColor','r');
set(handles.statuspanel,'BackgroundColor','r');
else
    set(handles.status,'String','Finished','BackgroundColor','g');
    set(handles.statuspanel,'BackgroundColor','g');
end
initialize_gui(gcbf, handles, true);

% --- Executes on button press in stop.
function stop_Callback(hObject, eventdata, handles);
% hObject    handle to stop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.stop,'Enable','off');
set(handles.pause,'Enable','off');
set(handles.calculate,'Enable','on');
set(handles.nsteps,'Enable','on');
set(handles.temp,'Enable','on');
set(handles.density,'Enable','on');
set(handles.dt,'Enable','on');
set(handles.np,'Enable','on');
set(handles.saveloc,'Enable','on');
set(handles.pause,'String','Pause');
fop=fopen('simdata','w+');
fwrite(fop,2);
fclose(fop);
initialize_gui(gcbf, handles, true);

% --- Executes when selected object changed in unitgroup.
function unitgroup_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in unitgroup 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --------------------------------------------------------------------
function initialize_gui(fig_handle, handles, isreset)
% If the metricdata field is present and the stop flag is false, it means
% we are we are just re-initializing a GUI by calling it from the cmd line
% while it is up. So, bail out as we dont want to stop the data.
if isfield(handles, 'metricdata') && ~isreset
    return;
end

% Update handles structure
guidata(handles.figure1, handles);



function saveloc_Callback(hObject, eventdata, handles)
% hObject    handle to saveloc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of saveloc as text
%        str2double(get(hObject,'String')) returns contents of saveloc as a double


% --- Executes during object creation, after setting all properties.
function saveloc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to saveloc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function temp_Callback(hObject, eventdata, handles)
% hObject    handle to temp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of temp as text
%        str2double(get(hObject,'String')) returns contents of temp as a double
density = str2double(get(hObject, 'String'));
if isnan(density)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end

% --- Executes during object creation, after setting all properties.
function temp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to temp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function density_Callback(hObject, eventdata, handles)
% hObject    handle to density (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of density as text
%        str2double(get(hObject,'String')) returns contents of density as a double
density = str2double(get(hObject, 'String'));
if isnan(density)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end

% --- Executes during object creation, after setting all properties.
function density_CreateFcn(hObject, eventdata, handles)
% hObject    handle to density (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dt_Callback(hObject, eventdata, handles)
% hObject    handle to dt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dt as text
%        str2double(get(hObject,'String')) returns contents of dt as a double
density = str2double(get(hObject, 'String'));
if isnan(density)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end

% --- Executes during object creation, after setting all properties.
function dt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nsteps_Callback(hObject, eventdata, handles)
% hObject    handle to nsteps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nsteps as text
%        str2double(get(hObject,'String')) returns contents of nsteps as a double
density = str2double(get(hObject, 'String'));
if isnan(density)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end

% --- Executes during object creation, after setting all properties.
function nsteps_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nsteps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pause.
function pause_Callback(hObject, eventdata, handles)
% hObject    handle to pause (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


fop=fopen('simdata','r');
stat=fread(fop);
fclose(fop);

if stat==0
    set(handles.pause,'String','Continue');
else
        set(handles.pause,'String','Pause');
end
fop=fopen('simdata','w+');
fwrite(fop,1-stat);
fclose(fop);



function edit22_Callback(hObject, eventdata, handles)
% hObject    handle to ncalc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ncalc as text
%        str2double(get(hObject,'String')) returns contents of ncalc as a double


% --- Executes during object creation, after setting all properties.
function edit22_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ncalc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ncalc_Callback(hObject, eventdata, handles)
% hObject    handle to ncalc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ncalc as text
%        str2double(get(hObject,'String')) returns contents of ncalc as a double


% --- Executes during object creation, after setting all properties.
function ncalc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ncalc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function savelocmc_Callback(hObject, eventdata, handles)
% hObject    handle to savelocmc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of savelocmc as text
%        str2double(get(hObject,'String')) returns contents of savelocmc as a double


% --- Executes during object creation, after setting all properties.
function savelocmc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to savelocmc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dr_Callback(hObject, eventdata, handles)
% hObject    handle to dr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dr as text
%        str2double(get(hObject,'String')) returns contents of dr as a double


% --- Executes during object creation, after setting all properties.
function dr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nequil_Callback(hObject, eventdata, handles)
% hObject    handle to nequil (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nequil as text
%        str2double(get(hObject,'String')) returns contents of nequil as a double


% --- Executes during object creation, after setting all properties.
function nequil_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nequil (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function npartmc_Callback(hObject, eventdata, handles)
% hObject    handle to npartmc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of npartmc as text
%        str2double(get(hObject,'String')) returns contents of npartmc as a double


% --- Executes during object creation, after setting all properties.
function npartmc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to npartmc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function movepercycle_Callback(hObject, eventdata, handles)
% hObject    handle to movepercycle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of movepercycle as text
%        str2double(get(hObject,'String')) returns contents of movepercycle as a double


% --- Executes during object creation, after setting all properties.
function movepercycle_CreateFcn(hObject, eventdata, handles)
% hObject    handle to movepercycle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


