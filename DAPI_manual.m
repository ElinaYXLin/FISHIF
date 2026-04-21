function varargout = DAPI_manual(varargin)
% DAPI_MANUAL M-file for DAPI_manual.fig
%      DAPI_MANUAL, by itself, creates a new DAPI_MANUAL or raises the existing
%      singleton*.
%
%      H = DAPI_MANUAL returns the handle to a new DAPI_MANUAL or the handle to
%      the existing singleton*.
%
%      DAPI_MANUAL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DAPI_MANUAL.M with the given input arguments.
%
%      DAPI_MANUAL('Property','Value',...) creates a new DAPI_MANUAL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before DAPI_manual_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to DAPI_manual_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help DAPI_manual

% Last Modified by GUIDE v2.5 20-Apr-2011 16:10:17

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @DAPI_manual_OpeningFcn, ...
                   'gui_OutputFcn',  @DAPI_manual_OutputFcn, ...
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


% --- Executes just before DAPI_manual is made visible.
function DAPI_manual_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to DAPI_manual (see VARARGIN)

global imstack imstack0 immask dlayer dcontrast dcontrmin dtpsize dfsize dimthresh dlogsize dlogth dmask resolution0 resolution0Z dclose dopen limtpsize limgsize limthresh limcontrast limclose limopen limlog limlogth

%% Default value setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dlayer = 1;
dcontrast = 1;
dcontrmin = 0;
dtpsize = 10;
dfsize = 1;
dimthresh = 1;
dlogsize = 1.5;
dlogth = 1;
dmask = true;
resolution0 = 0.13;
resolution0Z = 0.38;
dclose = 5;
dopen = 5;
limtpsize = [1,10];
limgsize = [1,10];
limthresh = [0,1];
limcontrast = [0,1];
limclose = [1,20];
limopen = [1,20];
limlog = [0.1,10];
limlogth = [0,1];
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Choose default command line output for stack_th
handles.output = hObject;
set(hObject,'Resize','On');
load_file_Callback(hObject, eventdata, handles)

% Update handles structure
guidata(hObject, handles);



% --- Outputs from this function are returned to the command line.
function varargout = DAPI_manual_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





% --- Executes on button press in load_file.
function load_file_Callback(hObject, eventdata, handles)
% hObject    handle to load_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global imstack imstack0 immask dlayer dcontrast dcontrmin dtpsize dfsize dimthresh dlogsize dlogth dmask resolution0 resolution0Z dclose dopen mstatus L_ratio L_ratioZ vlayer vcontrast vcontrmin vtpsize vfsize vimthresh vlogsize vlogth vmask vopen vclose limtpsize limgsize limthresh limcontrast imname imfolder limlog limlogth
set(handles.message_show,'String','Wait');
% Update handles structure
guidata(hObject, handles);

imstack = zeros(0);
imstack0 = zeros(0);
immask = false(0);
Idim = 1;

%% Image stack loading: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%if isempty(varargin)
    [imname,imfolder,~] = uigetfile('*.*','Select a file to load');   %%% Pick up a lsm file
%end

set(handles.load_file,'String','Busy');
% Update handles structure
guidata(hObject, handles);

lsm_stack = tiffread([imfolder,imname]); %%% Lsm file loading
immax = length(lsm_stack);
tempmax = zeros(1,immax);
if isa(lsm_stack(1).data,'cell')
    for I_layer = 1:immax
        temp = lsm_stack(I_layer).data;
        imstack0(:,:,I_layer) = temp{Idim};
        tempmax(I_layer) = max(max(imstack0(:,:,I_layer)));
    end
else
    for I_layer = 1:immax
        temp = lsm_stack(I_layer).data;
        imstack0(:,:,I_layer) = temp;
        tempmax(I_layer) = max(max(imstack(:,:,I_layer)));
    end
end

imstack = imstack0;
imstack = double(imstack)/65535;   %%% Normalization
imstack0 = double(imstack);
immask = false(size(imstack));
resolution = lsm_stack(1).lsm.VoxelSizeX/(1e-6);
resolutionZ = lsm_stack(1).lsm.VoxelSizeZ/(1e-6);
L_ratio = (resolution/resolution0);
L_ratioZ = (resolutionZ/resolution0Z);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% initialization: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%vlayer = dlayer;
[~,vlayer] = max(tempmax);
vcontrast = dcontrast;
vcontrmin = dcontrmin;
vtpsize = dtpsize;
vfsize = dfsize;
vimthresh = dimthresh;
vlogsize = dlogsize;
vlogth = dlogth;
vmask = dmask;
vopen = dopen;
vclose = dclose;

set(handles.file_name,'String',[imfolder,imname]);
set(handles.layer_slide,'Max',immax);
set(handles.layer_slide,'Min',1);
set(handles.layer_slide,'SliderStep',[1/immax,1/immax]);
set(handles.layer_slide,'Value',vlayer);
set(handles.layer_value,'String',num2str(vlayer));
set(handles.max_layer,'String',num2str(immax));
set(handles.contrast_min,'String',num2str(vcontrmin));
set(handles.contrast_max,'String',num2str(vcontrast));
set(handles.tp_value,'String',num2str(vtpsize));
set(handles.gau_value,'String',num2str(vfsize));
set(handles.imthresh,'String',num2str(vimthresh));
set(handles.log_th,'String',num2str(vlogth));
set(handles.mask_on,'Value',vmask);
set(handles.close_value,'String',num2str(vclose));
set(handles.open_value,'String',num2str(vopen));
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Image output: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
overlay = DAPI_show(imstack,immask,vlayer,vcontrast,vcontrmin,vmask);
axes(handles.stack_show);
imshow(overlay)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(handles.message_show,'String','Done');
% Update handles structure
guidata(hObject, handles);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% --- Executes on slider movement.
function layer_slide_Callback(hObject, eventdata, handles)
% hObject    handle to layer_slide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global imstack imstack0 immask dlayer dcontrast dcontrmin dtpsize dfsize dimthresh dmask resolution0 resolution0Z dclose dopen mstatus L_ratio L_ratioZ vlayer vcontrast vcontrmin vtpsize vfsize vimthresh vmask vopen vclose limtpsize limgsize limthresh limcontrast imname imfolder
vlayer = floor(get(handles.layer_slide,'Value'));
set(handles.layer_value,'String', num2str(vlayer));

overlay = DAPI_show(imstack,immask,vlayer,vcontrast,vcontrmin,vmask);
axes(handles.stack_show);
v = axis; 
imshow(overlay)
axis(v);
% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function layer_slide_CreateFcn(hObject, eventdata, handles)
% hObject    handle to layer_slide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





function layer_value_Callback(hObject, eventdata, handles)
% hObject    handle to layer_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global imstack imstack0 immask dlayer dcontrast dcontrmin dtpsize dfsize dimthresh dmask resolution0 resolution0Z dclose dopen mstatus L_ratio L_ratioZ vlayer vcontrast vcontrmin vtpsize vfsize vimthresh vmask vopen vclose limtpsize limgsize limthresh limcontrast imname imfolder
minT = get(handles.layer_slide,'Min');
maxT = get(handles.layer_slide,'Max');
vlayer = floor(str2num(get(handles.layer_value,'String')));
if (isempty(vlayer) || vlayer < minT)
    vlayer = minT;
elseif vlayer > maxT
    vlayer = maxT;
end
set(handles.layer_slide,'Value',vlayer);
set(handles.layer_value,'String',num2str(vlayer));

overlay = DAPI_show(imstack,immask,vlayer,vcontrast,vcontrmin,vmask);
axes(handles.stack_show);
v = axis; 
imshow(overlay)
axis(v);
% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function layer_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to layer_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






function contrast_min_Callback(hObject, eventdata, handles)
% hObject    handle to contrast_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global imstack imstack0 immask dlayer dcontrast dcontrmin dtpsize dfsize dimthresh dmask resolution0 resolution0Z dclose dopen mstatus L_ratio L_ratioZ vlayer vcontrast vcontrmin vtpsize vfsize vimthresh vmask vopen vclose limtpsize limgsize limthresh limcontrast imname imfolder
minT = limcontrast(1);
maxT = vcontrast;
vcontrmin = str2num(get(handles.contrast_min,'String'));
if (isempty(vcontrmin) || vcontrmin < minT)
    vcontrmin = minT;
elseif vcontrmin > maxT
    vcontrmin = maxT;
end
set(handles.contrast_min,'String',num2str(vcontrmin));

overlay = DAPI_show(imstack,immask,vlayer,vcontrast,vcontrmin,vmask);
axes(handles.stack_show);
v = axis; 
imshow(overlay)
axis(v);
% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function contrast_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to contrast_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







function contrast_max_Callback(hObject, eventdata, handles)
% hObject    handle to contrast_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global imstack imstack0 immask dlayer dcontrast dcontrmin dtpsize dfsize dimthresh dmask resolution0 resolution0Z dclose dopen mstatus L_ratio L_ratioZ vlayer vcontrast vcontrmin vtpsize vfsize vimthresh vmask vopen vclose limtpsize limgsize limthresh limcontrast imname imfolder
minT = vcontrmin;
maxT = limcontrast(2);
vcontrast = str2num(get(handles.contrast_max,'String'));
if (isempty(vcontrast) || vcontrast < minT)
    vcontrast = minT;
elseif vcontrast > maxT
    vcontrast = maxT;
end
set(handles.contrast_max,'String',num2str(vcontrast));

overlay = DAPI_show(imstack,immask,vlayer,vcontrast,vcontrmin,vmask);
axes(handles.stack_show);
v = axis; 
imshow(overlay)
axis(v);
% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function contrast_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to contrast_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








function tp_value_Callback(hObject, eventdata, handles)
% hObject    handle to tp_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global imstack imstack0 immask dlayer dcontrast dcontrmin dtpsize dfsize dimthresh dmask resolution0 resolution0Z dclose dopen mstatus L_ratio L_ratioZ vlayer vcontrast vcontrmin vtpsize vfsize vimthresh vmask vopen vclose limtpsize limgsize limthresh limcontrast imname imfolder
minT = limtpsize(1);
maxT = limtpsize(2);
vtpsize = round(str2num(get(handles.tp_value,'String')));
if (isempty(vtpsize) || vtpsize < minT)
    vtpsize = minT;
elseif vtpsize > maxT
    vtpsize = maxT;
end
set(handles.tp_value,'String',num2str(vtpsize));
% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function tp_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tp_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








% --- Executes on button press in tp_on.
function tp_on_Callback(hObject, eventdata, handles)
% hObject    handle to tp_on (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global imstack imstack0 immask dlayer dcontrast dcontrmin dtpsize dfsize dimthresh dmask resolution0 resolution0Z dclose dopen mstatus L_ratio L_ratioZ vlayer vcontrast vcontrmin vtpsize vfsize vimthresh vmask vopen vclose limtpsize limgsize limthresh limcontrast imname imfolder
set(handles.message_show,'String','Wait');
% Update handles structure
guidata(hObject, handles);

% se = ball_gen(vtpsize);
% imstack = imtophat(imstack,se);
se = strel('disk',vtpsize);
for kk = 1:size(imstack,3)
    imstack(:,:,kk) = imtophat(imstack(:,:,kk),se);
end
imstack = imstack./max(max(max(imstack)));
overlay = DAPI_show(imstack,immask,vlayer,vcontrast,vcontrmin,vmask);
axes(handles.stack_show);
v = axis; 
imshow(overlay)
axis(v);
set(handles.message_show,'String','Done');
% Update handles structure
guidata(hObject, handles);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%









function gau_value_Callback(hObject, eventdata, handles)
% hObject    handle to gau_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global imstack imstack0 immask dlayer dcontrast dcontrmin dtpsize dfsize dimthresh dmask resolution0 resolution0Z dclose dopen mstatus L_ratio L_ratioZ vlayer vcontrast vcontrmin vtpsize vfsize vimthresh vmask vopen vclose limtpsize limgsize limthresh limcontrast imname imfolder
minT = limgsize(1);
maxT = limgsize(2);
vfsize = round(str2num(get(handles.gau_value,'String')));
if (isempty(vfsize) || vfsize < minT)
    vfsize = minT;
elseif vfsize > maxT
    vfsize = maxT;
end
set(handles.gau_value,'String',num2str(vfsize));
% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function gau_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gau_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







% --- Executes on button press in gau_on.
function gau_on_Callback(hObject, eventdata, handles)
% hObject    handle to gau_on (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global imstack imstack0 immask dlayer dcontrast dcontrmin dtpsize dfsize dimthresh dmask resolution0 resolution0Z dclose dopen mstatus L_ratio L_ratioZ vlayer vcontrast vcontrmin vtpsize vfsize vimthresh vmask vopen vclose limtpsize limgsize limthresh limcontrast imname imfolder
% gau2 = fspecial('gaussian',5,vfsize);
% gau1 = gau2(3,:);
% se = repmat(gau2,[1,1,5]).*repmat(permute(gau1,[1,3,2]),[5,5,1]);
% imstack = imfilter(imstack,se,'symmetric','conv');
set(handles.message_show,'String','Wait');
% Update handles structure
guidata(hObject, handles);

se = fspecial('gaussian',10,vfsize);
for kk = 1:size(imstack,3)
    imstack(:,:,kk) = imfilter(imstack(:,:,kk),se,'symmetric','conv');
end
imstack = imstack./max(max(max(imstack)));
overlay = DAPI_show(imstack,immask,vlayer,vcontrast,vcontrmin,vmask);
axes(handles.stack_show);
v = axis; 
imshow(overlay)
axis(v);

set(handles.message_show,'String','Done');
% Update handles structure
guidata(hObject, handles);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








function imthresh_Callback(hObject, eventdata, handles)
% hObject    handle to imthresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global imstack imstack0 immask dlayer dcontrast dcontrmin dtpsize dfsize dimthresh dmask resolution0 resolution0Z dclose dopen mstatus L_ratio L_ratioZ vlayer vcontrast vcontrmin vtpsize vfsize vimthresh vmask vopen vclose limtpsize limgsize limthresh limcontrast imname imfolder
minT = limthresh(1);
maxT = limthresh(2);
vimthresh = str2num(get(handles.imthresh,'String'));
if (isempty(vimthresh) || vimthresh < minT)
    vimthresh = minT;
elseif vimthresh > maxT
    vimthresh = maxT;
end
set(handles.imthresh,'String',num2str(vimthresh));
% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function imthresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to imthresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








% --- Executes on button press in Thresh_default_on.
function Thresh_default_on_Callback(hObject, eventdata, handles)
% hObject    handle to Thresh_default_on (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global imstack imstack0 immask dlayer dcontrast dcontrmin dtpsize dfsize dimthresh dmask resolution0 resolution0Z dclose dopen mstatus L_ratio L_ratioZ vlayer vcontrast vcontrmin vtpsize vfsize vimthresh vmask vopen vclose limtpsize limgsize limthresh limcontrast imname imfolder
set(handles.message_show,'String','Done');
% Update handles structure
guidata(hObject, handles);

vimthresh = graythresh(imstack);
set(handles.imthresh,'String',num2str(vimthresh));
set(handles.message_show,'String','Done');
% Update handles structure
guidata(hObject, handles);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








% --- Executes on button press in thresh_on.
function thresh_on_Callback(hObject, eventdata, handles)
% hObject    handle to thresh_on (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global imstack imstack0 immask dlayer dcontrast dcontrmin dtpsize dfsize dimthresh dmask resolution0 resolution0Z dclose dopen mstatus L_ratio L_ratioZ vlayer vcontrast vcontrmin vtpsize vfsize vimthresh vmask vopen vclose limtpsize limgsize limthresh limcontrast imname imfolder
set(handles.message_show,'String','Wait');
% Update handles structure
guidata(hObject, handles);

for kk = 1:size(imstack,3)
    immask(:,:,kk) = im2bw(imstack(:,:,kk),vimthresh);
end
overlay = DAPI_show(imstack,immask,vlayer,vcontrast,vcontrmin,vmask);
axes(handles.stack_show);
v = axis; 
imshow(overlay)
axis(v);
set(handles.message_show,'String','Done');
% Update handles structure
guidata(hObject, handles);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








% --- Executes on button press in im_lock.
function im_lock_Callback(hObject, eventdata, handles)
% hObject    handle to im_lock (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








% --- Executes on button press in mask_on.
function mask_on_Callback(hObject, eventdata, handles)
% hObject    handle to mask_on (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global imstack imstack0 immask dlayer dcontrast dcontrmin dtpsize dfsize dimthresh dmask resolution0 resolution0Z dclose dopen mstatus L_ratio L_ratioZ vlayer vcontrast vcontrmin vtpsize vfsize vimthresh vmask vopen vclose limtpsize limgsize limthresh limcontrast imname imfolder
vmask = get(handles.mask_on,'Value');
overlay = DAPI_show(imstack,immask,vlayer,vcontrast,vcontrmin,vmask);
axes(handles.stack_show);
v = axis; 
imshow(overlay)
axis(v);
% Update handles structure
guidata(hObject, handles);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








% --- Executes on button press in im_cancel.
function im_cancel_Callback(hObject, eventdata, handles)
% hObject    handle to im_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global imstack imstack0 immask dlayer dcontrast dcontrmin dtpsize dfsize dimthresh dmask resolution0 resolution0Z dclose dopen mstatus L_ratio L_ratioZ vlayer vcontrast vcontrmin vtpsize vfsize vimthresh vmask vopen vclose limtpsize limgsize limthresh limcontrast imname imfolder
imstack = imstack0;
immask(:,:,:) = 0;
overlay = DAPI_show(imstack,immask,vlayer,vcontrast,vcontrmin,vmask);
axes(handles.stack_show);
v = axis; 
imshow(overlay)
axis(v);
% Update handles structure
guidata(hObject, handles);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








function close_value_Callback(hObject, eventdata, handles)
% hObject    handle to close_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global imstack imstack0 immask dlayer dcontrast dcontrmin dtpsize dfsize dimthresh dmask resolution0 resolution0Z dclose dopen mstatus L_ratio L_ratioZ vlayer vcontrast vcontrmin vtpsize vfsize vimthresh vmask vopen vclose limtpsize limgsize limthresh limcontrast limclose limopen imname imfolder
minT = limclose(1);
maxT = limclose(2);
vclose = round(str2num(get(handles.close_value,'String')));
if (isempty(vclose) || vclose < minT)
    vclose = minT;
elseif vclose > maxT
    vclose = maxT;
end
set(handles.close_value,'String',num2str(vclose));
% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function close_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to close_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%









% --- Executes on button press in close_on.
function close_on_Callback(hObject, eventdata, handles)
% hObject    handle to close_on (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global imstack imstack0 immask dlayer dcontrast dcontrmin dtpsize dfsize dimthresh dmask resolution0 resolution0Z dclose dopen mstatus L_ratio L_ratioZ vlayer vcontrast vcontrmin vtpsize vfsize vimthresh vmask vopen vclose limtpsize limgsize limthresh limcontrast imname imfolder
set(handles.message_show,'String','Wait');
% Update handles structure
guidata(hObject, handles);

% se = ball_gen(vclose);
% immask = imclose(immask,se);
se = strel('disk',vclose);
for kk = 1:size(imstack,3)
    immask(:,:,kk) = imclose(immask(:,:,kk),se);
end
overlay = DAPI_show(imstack,immask,vlayer,vcontrast,vcontrmin,vmask);
axes(handles.stack_show);
v = axis; 
imshow(overlay)
axis(v);
set(handles.message_show,'String','Done');
% Update handles structure
guidata(hObject, handles);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%









function open_value_Callback(hObject, eventdata, handles)
% hObject    handle to open_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global imstack imstack0 immask dlayer dcontrast dcontrmin dtpsize dfsize dimthresh dmask resolution0 resolution0Z dclose dopen mstatus L_ratio L_ratioZ vlayer vcontrast vcontrmin vtpsize vfsize vimthresh vmask vopen vclose limtpsize limgsize limthresh limcontrast limclose limopen imname imfolder
minT = limopen(1);
maxT = limopen(2);
vopen = round(str2num(get(handles.open_value,'String')));
if (isempty(vopen) || vopen < minT)
    vopen = minT;
elseif vopen > maxT
    vopen = maxT;
end
set(handles.open_value,'String',num2str(vopen));
% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function open_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to open_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%









% --- Executes on button press in open_on.
function open_on_Callback(hObject, eventdata, handles)
% hObject    handle to open_on (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global imstack imstack0 immask dlayer dcontrast dcontrmin dtpsize dfsize dimthresh dmask resolution0 resolution0Z dclose dopen mstatus L_ratio L_ratioZ vlayer vcontrast vcontrmin vtpsize vfsize vimthresh vmask vopen vclose limtpsize limgsize limthresh limcontrast imname imfolder
set(handles.message_show,'String','Wait');
% Update handles structure
guidata(hObject, handles);

% se = ball_gen(vopen);
% immask = imopen(immask,se);
se = strel('disk',vopen);
for kk = 1:size(imstack,3)
    immask(:,:,kk) = imopen(immask(:,:,kk),se);
end
overlay = DAPI_show(imstack,immask,vlayer,vcontrast,vcontrmin,vmask);
axes(handles.stack_show);
v = axis; 
imshow(overlay)
axis(v);
set(handles.message_show,'String','Done');
% Update handles structure
guidata(hObject, handles);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%









% --- Executes on button press in fill_on.
function fill_on_Callback(hObject, eventdata, handles)
% hObject    handle to fill_on (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global imstack imstack0 immask dlayer dcontrast dcontrmin dtpsize dfsize dimthresh dmask resolution0 resolution0Z dclose dopen mstatus L_ratio L_ratioZ vlayer vcontrast vcontrmin vtpsize vfsize vimthresh vmask vopen vclose limtpsize limgsize limthresh limcontrast imname imfolder
set(handles.message_show,'String','Wait');
% Update handles structure
guidata(hObject, handles);

for kk = 1:size(imstack,3)
    immask(:,:,kk) = imfill(immask(:,:,kk),'holes');
end
overlay = DAPI_show(imstack,immask,vlayer,vcontrast,vcontrmin,vmask);
axes(handles.stack_show);
v = axis; 
imshow(overlay)
axis(v);
set(handles.message_show,'String','Done');
% Update handles structure
guidata(hObject, handles);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%









% --- Executes on button press in convex_on.
function convex_on_Callback(hObject, eventdata, handles)
% hObject    handle to convex_on (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%










function log_th_Callback(hObject, eventdata, handles)
% hObject    handle to log_th (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global imstack imstack0 immask dlayer dcontrast dcontrmin dtpsize dfsize dimthresh dlogsize dlogth dmask resolution0 resolution0Z dclose dopen mstatus L_ratio L_ratioZ vlayer vcontrast vcontrmin vtpsize vfsize vimthresh vlogsize vlogth vmask vopen vclose limtpsize limgsize limthresh limcontrast imname imfolder limlog limlogth
minT = limlogth(1);
maxT = limlogth(2);
vlogth = str2num(get(handles.log_th,'String'));
if (isempty(vlogth) || vlogth < minT)
    vlogth = minT;
elseif vlogth > maxT
    vlogth = maxT;
end
set(handles.log_th,'String',num2str(vlogth));
% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function log_th_CreateFcn(hObject, eventdata, handles)
% hObject    handle to log_th (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%









% --- Executes on button press in default_log.
function default_log_Callback(hObject, eventdata, handles)
% hObject    handle to default_log (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global imstack imstack0 immask dlayer dcontrast dcontrmin dtpsize dfsize dimthresh dlogsize dlogth dmask resolution0 resolution0Z dclose dopen mstatus L_ratio L_ratioZ vlayer vcontrast vcontrmin vtpsize vfsize vimthresh vlogsize vlogth vmask vopen vclose limtpsize limgsize limthresh limcontrast imname imfolder limlog limlogth g
set(handles.message_show,'String','Wait');
% Update handles structure
guidata(hObject, handles);

thtemp = zeros(1,size(imstack,3));
for kk = 1:size(imstack,3)
    [~,thtemp(kk)] = edge(imstack(:,:,kk),'sobel');
end
vlogth = max(thtemp);
set(handles.log_th,'String',num2str(vlogth));
set(handles.message_show,'String','Done');
% Update handles structure
guidata(hObject, handles);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








% --- Executes on button press in log_on.
function log_on_Callback(hObject, eventdata, handles)
% hObject    handle to log_on (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global imstack imstack0 immask dlayer dcontrast dcontrmin dtpsize dfsize dimthresh dlogsize dlogth dmask resolution0 resolution0Z dclose dopen mstatus L_ratio L_ratioZ vlayer vcontrast vcontrmin vtpsize vfsize vimthresh vlogsize vlogth vmask vopen vclose limtpsize limgsize limthresh limcontrast imname imfolder limlog limlogth g
% gau2 = fspecial('gaussian',5,vfsize);
% gau1 = gau2(3,:);
% se = repmat(gau2,[1,1,5]).*repmat(permute(gau1,[1,3,2]),[5,5,1]);
% imstack = imfilter(imstack,se,'symmetric','conv');
set(handles.message_show,'String','Wait');
% Update handles structure
guidata(hObject, handles);

for kk = 1:size(imstack,3)
    immask(:,:,kk) = edge(imstack(:,:,kk),'sobel',vlogth);
end
overlay = DAPI_show(imstack,immask,vlayer,vcontrast,vcontrmin,vmask);
axes(handles.stack_show);
v = axis; 
imshow(overlay)
axis(v);

set(handles.message_show,'String','Done');
% Update handles structure
guidata(hObject, handles);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
