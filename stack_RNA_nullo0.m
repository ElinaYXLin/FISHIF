function varargout = stack_RNA_nullo0(varargin)
% STACK_RNA_NULLO0 M-file for stack_RNA_nullo0.fig
%      STACK_RNA_NULLO0, by itself, creates a new STACK_RNA_NULLO0 or raises the existing
%      singleton*.
%
%      H = STACK_RNA_NULLO0 returns the handle to a new STACK_RNA_NULLO0 or the handle to
%      the existing singleton*.
%
%      STACK_RNA_NULLO0('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in STACK_RNA_NULLO0.M with the given input arguments.
%
%      STACK_RNA_NULLO0('Property','Value',...) creates a new STACK_RNA_NULLO0 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before stack_RNA_nullo0_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to stack_RNA_nullo0_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help stack_RNA_nullo0

% Last Modified by GUIDE v2.5 07-Apr-2012 21:37:53

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @stack_RNA_nullo0_OpeningFcn, ...
                   'gui_OutputFcn',  @stack_RNA_nullo0_OutputFcn, ...
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


% --- Executes just before stack_RNA_nullo0 is made visible.
function stack_RNA_nullo0_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to stack_RNA_nullo0 (see VARARGIN)
global imstack immask dlayer dcontrast dcontrmin dfsize dimthresh dmask dpeak mstatus resolution0 L_ratio g vlayer vcontrast vcontrmin vfsize vimthresh vmask vpeak limsize limthresh limcontrast bk_label sub_folder bk_name sg_name stack_th stack_th0

%% Default value setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dlayer = 1;
dcontrast = 1;
dcontrmin = 0;
dfsize = 1.5;
dimthresh = 0.01;
dmask = true;
dpeak = true;
stack_th0 = 0;
stack_th = stack_th0;
bk_label = 'OreR';
mstatus = {'Mask on','Mask off'};
resolution0 = 0.13;
limsize = [0,10];
limthresh = [0,1];
limcontrast = [0,1];
sub_folder = 'Histogram_nullo/';
bk_name = 'bk_record.xls';
sg_name = 'signal_record.xls';
set(handles.show_bk,'Value',false);
set(handles.presult_on,'Value',true);
set(handles.Fitted_on,'Value',false);
% set(handles.sigma_show,'Color','none');
% set(handles.sigma_zoom,'Color','none');
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Choose default command line output for stack_RNA_nullo0
handles.output = hObject;
set(hObject,'Resize','On');
%load_file_Callback(hObject, eventdata, handles)

set(gcf,'WindowButtonDownFcn',{@ButttonDownFcn2});
% Update handles structure
guidata(hObject, handles);

% if ~isempty(varargin)
%     hist_auto_Callback(hObject, eventdata, handles, varargin{1});
% end
% UIWAIT makes stack_RNA_nullo0 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = stack_RNA_nullo0_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% --- Executes on slider movement.
function layer_slide_Callback(hObject, eventdata, handles)
% hObject    handle to layer_slide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global imstack immask dlayer dcontrast dcontrmin dfsize dimthresh dmask dpeak mstatus resolution0 L_ratio g vlayer vcontrast vcontrmin vfsize vimthresh vmask vpeak limsize limthresh limcontrast
vlayer = floor(get(handles.layer_slide,'Value'));
set(handles.layer_value,'String', num2str(vlayer));

overlay = cf_show1(imstack,g,vlayer,vimthresh,vpeak,vcontrast,vcontrmin,vmask,get(handles.Fitted_on,'Value'));
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function layer_value_Callback(hObject, eventdata, handles)
% hObject    handle to layer_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global imstack immask dlayer dcontrast dcontrmin dfsize dimthresh dmask dpeak mstatus resolution0 L_ratio g vlayer vcontrast vcontrmin vfsize vimthresh vmask vpeak limsize limthresh limcontrast
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

overlay = cf_show1(imstack,g,vlayer,vimthresh,vpeak,vcontrast,vcontrmin,vmask,get(handles.Fitted_on,'Value'));
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function fsize_Callback(hObject, eventdata, handles)
% hObject    handle to fsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global imstack immask dlayer dcontrast dcontrmin dfsize dimthresh dmask dpeak mstatus resolution0 L_ratio g vlayer vcontrast vcontrmin vfsize vimthresh vmask vpeak limsize limthresh limcontrast
minT = limsize(1);
maxT = limsize(2);
vfsize = str2num(get(handles.fsize,'String'));
if (isempty(vfsize) || vfsize < minT)
    vfsize = minT;
elseif vfsize > maxT
    vfsize = maxT;
end
set(handles.fsize,'String',num2str(vfsize));

[overlay,g] = cf_show(imstack,vlayer,vfsize,L_ratio,vimthresh,vpeak,vcontrast,vcontrmin,vmask,get(handles.Fitted_on,'Value'));
axes(handles.stack_show);
v = axis; 
imshow(overlay)
axis(v);
% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function fsize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function imthresh_Callback(hObject, eventdata, handles)
% hObject    handle to imthresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global imstack immask dlayer dcontrast dcontrmin dfsize dimthresh dmask dpeak mstatus resolution0 L_ratio g vlayer vcontrast vcontrmin vfsize vimthresh vmask vpeak limsize limthresh limcontrast
minT = limthresh(1);
maxT = limthresh(2);
vimthresh = str2num(get(handles.imthresh,'String'));
if (isempty(vimthresh) || vimthresh < minT)
    vimthresh = minT;
elseif vimthresh > maxT
    vimthresh = maxT;
end
set(handles.imthresh,'String',num2str(vimthresh));

overlay = cf_show1(imstack,g,vlayer,vimthresh,vpeak,vcontrast,vcontrmin,vmask,get(handles.Fitted_on,'Value'));
axes(handles.stack_show);
v = axis; 
imshow(overlay)
axis(v);
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% --- Executes on button press in load_file.
function load_file_Callback(hObject, eventdata, handles, varargin)
% hObject    handle to load_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global imstack imstack0 immask dlayer dcontrast dcontrmin dfsize dimthresh dmask dpeak mstatus resolution0 L_ratio g vlayer vcontrast vcontrmin vfsize vimthresh vmask vpeak limsize limthresh limcontrast max_spot n0 xout0 n xout imname imfolder bk_label bk_use mask_out stack_th
imstack = zeros(0);
imstack0 = zeros(0);
immask = false(0);
max_spot = zeros(0);
n0 = zeros(0);
xout0 = zeros(0);
n = zeros(0);
xout = zeros(0);
bk_use = cell(0);
%Idim = 2;
standard_record = 'Calibration/Results/standard.mat';
global standard_data
standard_data = load(standard_record);
list_name = 'matchlist.xls';
image_tail = '*.tif';

%% Image stack loading: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(varargin)
    imfolder = uigetdir('','Select a stack to load');   %%% Pick up a image stack folder
    imname = imfolder(find((imfolder == '/') | (imfolder == '\'),1,'last')+1:end);
    imfolder = imfolder(1:find((imfolder == '/') | (imfolder == '\'),1,'last'));
end
% load([imfolder(1:find((imfolder(1:end-1) == '/') | (imfolder(1:end-1) == '\'),1,'last')),'masks/',imname,'/mask.mat']);
[num_list, file_list] = xlsread([imfolder,list_name]);
I_file = strcmp([imname,'/'],file_list(:,3));
I_file = I_file | strcmp([imname,'\'],file_list(:,3));

set(handles.load_file,'String','Busy');
% Update handles structure
guidata(hObject, handles);

lsm_stack = dir([imfolder,imname,'/',image_tail]); %%% image file list loading
resolution = num_list(I_file,9);
%RNA_channel = num_list(I_file,10);
RNA_channel = num_list(I_file,12);
L_ratio = (resolution/resolution0);

immax = length(lsm_stack);
tempmax = zeros(1,immax);

for I_layer = 1:immax
    temp0 = imread([imfolder,imname,'/',lsm_stack(I_layer).name]);
    temp = temp0(:,:,RNA_channel);
    imstack0(:,:,I_layer) = temp;
    imstack(:,:,I_layer) = imfilter(temp,fspecial('gaussian',5,1),'symmetric','conv');
    tempmax(I_layer) = max(max(imstack(:,:,I_layer)));
end

% %imstack0 = imstack0(ceil(size(imstack0,1)*1/4):ceil(size(imstack0,1)*3/4),ceil(size(imstack0,2)*1/4):ceil(size(imstack0,2)*3/4),:);
% %imstack  = imstack(ceil(size(imstack,1)*1/4):ceil(size(imstack,1)*3/4),ceil(size(imstack,2)*1/4):ceil(size(imstack,2)*3/4),:);
% 
% mask1D = max(max(logical(mask_stack),[],3),[],1);
% ELmin = find(mask1D,1);
% ELmax = find(mask1D,1,'last');
% 
% rxmin1 = 4/16;
% rxmax1 = 6/16;
% rymin1 = 3/8;
% rymax1 = 5/8;
% 
% xmin1 = round(ELmin+(ELmax-ELmin)*rxmin1);
% xmax1 = round(ELmin+(ELmax-ELmin)*rxmax1);
% ymin1 = round(rymin1*size(imstack,1));
% ymax1 = round(rymax1*size(imstack,1));
% 
% 
% rxmin2 = 1-6/16;
% rxmax2 = 1-4/16;
% rymin2 = 3/8;
% rymax2 = 5/8;
% 
% xmin2 = round(ELmin+(ELmax-ELmin)*rxmin2);
% xmax2 = round(ELmin+(ELmax-ELmin)*rxmax2);
% ymin2 = round(rymin2*size(imstack,1));
% ymax2 = round(rymax2*size(imstack,1));
% 
% if mean(mean(mean(imstack0(ymin1:ymax1,xmin1:xmax1,:)))) >= mean(mean(mean(imstack0(ymin2:ymax2,xmin2:xmax2,:))))
%     imstack0 = imstack0(ymin1:ymax1,xmin1:xmax1,:);
%     imstack  = imstack(ymin1:ymax1,xmin1:xmax1,:);
% else
%     imstack0 = imstack0(ymin2:ymax2,xmin2:xmax2,:);
%     imstack  = imstack(ymin2:ymax2,xmin2:xmax2,:);
% end



temp0 = fspecial('gaussian',5,1);
kernel0 = zeros(1,1,5);
kernel0(1,1,:) = temp0(3,:)./sum(temp0(3,:));
imstack = imfilter(imstack,kernel0,'symmetric','conv');

imstack = double(imstack)/65535;   %%% Normalization
imstack0 = double(imstack0);
immask = false(size(immask));
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% initialization: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%vlayer = dlayer;
[~,vlayer] = max(tempmax);
vcontrast = dcontrast;
vcontrmin = dcontrmin;
vfsize = dfsize;
vimthresh = dimthresh;
vmask = dmask;
vpeak = dpeak;
set(handles.file_name,'String',[imfolder,imname]);
set(handles.layer_slide,'Max',immax);
set(handles.layer_slide,'Min',1);
set(handles.layer_slide,'SliderStep',[1/(immax-1),1/(immax-1)]);
set(handles.layer_slide,'Value',vlayer);
set(handles.layer_value,'String',num2str(vlayer));
set(handles.max_layer,'String',num2str(immax));
set(handles.contrast_min,'String',num2str(vcontrmin));
set(handles.contrast_value,'String',num2str(vcontrast));
set(handles.fsize,'String',num2str(vfsize));
set(handles.imthresh,'String',num2str(vimthresh));
set(handles.Mask_status,'Value',vmask);
set(handles.Mask_status,'String',mstatus{vmask+1});
set(handles.Peak_on,'Value',vpeak);
% if ~isempty(strfind(imname,bk_label))
%     set(handles.background_on,'Value',true);
% else
%     set(handles.background_on,'Value',false);
% end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Image output: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mask_out = spmask(imstack,1,stack_th);   %%% 3D local maximal mask generation (Primary filter)
%mask_out = mask_out & (imstack*65535 >= 7000);   %%% Prefilter using peak intensity threshold
[overlay,g] = cf_show(imstack,vlayer,vfsize,L_ratio,vimthresh,vpeak,vcontrast,vcontrmin,vmask,get(handles.Fitted_on,'Value'));
axes(handles.stack_show);
imshow(overlay)

axes(handles.hist_show);
cla
axes(handles.hist_zoom);
cla
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(handles.load_file,'String','Load');
% Update handles structure
guidata(hObject, handles);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% --- Executes on button press in Mask_status.
function Mask_status_Callback(hObject, eventdata, handles)
% hObject    handle to Mask_status (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global imstack immask dlayer dcontrast dcontrmin dfsize dimthresh dmask dpeak mstatus resolution0 L_ratio g vlayer vcontrast vcontrmin vfsize vimthresh vmask vpeak limsize limthresh limcontrast
vmask = get(handles.Mask_status,'Value');
set(handles.Mask_status,'String',mstatus{vmask+1});

overlay = cf_show1(imstack,g,vlayer,vimthresh,vpeak,vcontrast,vcontrmin,vmask,get(handles.Fitted_on,'Value'));
axes(handles.stack_show);
v = axis; 
imshow(overlay)
axis(v);
% Update handles structure
guidata(hObject, handles);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




function contrast_value_Callback(hObject, eventdata, handles)
% hObject    handle to contrast_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global imstack immask dlayer dcontrast dcontrmin dfsize dimthresh dmask dpeak mstatus resolution0 L_ratio g vlayer vcontrast vcontrmin vfsize vimthresh vmask vpeak limsize limthresh limcontrast
minT = limcontrast(1);
maxT = limcontrast(2);
vcontrast = str2num(get(handles.contrast_value,'String'));
if (isempty(vcontrast) || vcontrast < minT)
    vcontrast = minT;
elseif vcontrast > maxT
    vcontrast = maxT;
end
set(handles.contrast_value,'String',num2str(vcontrast));

overlay = cf_show1(imstack,g,vlayer,vimthresh,vpeak,vcontrast,vcontrmin,vmask,get(handles.Fitted_on,'Value'));
axes(handles.stack_show);
v = axis; 
imshow(overlay)
axis(v);
% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function contrast_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to contrast_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function contrast_min_Callback(hObject, eventdata, handles)
% hObject    handle to contrast_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global imstack immask dlayer dcontrast dcontrmin dfsize dimthresh dmask dpeak mstatus resolution0 L_ratio g vlayer vcontrast vcontrmin vfsize vimthresh vmask vpeak limsize limthresh limcontrast
minT = limcontrast(1);
maxT = vcontrast;
vcontrmin = str2num(get(handles.contrast_min,'String'));
if (isempty(vcontrmin) || vcontrmin < minT)
    vcontrmin = minT;
elseif vcontrmin > maxT
    vcontrmin = maxT;
end
set(handles.contrast_min,'String',num2str(vcontrmin));

overlay = cf_show1(imstack,g,vlayer,vimthresh,vpeak,vcontrast,vcontrmin,vmask,get(handles.Fitted_on,'Value'));
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






% --- Executes on button press in hist_on.
function hist_on_Callback(hObject, eventdata, handles)
% hObject    handle to hist_on (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global imstack imstack0 immask dlayer dcontrast dcontrmin dfsize dimthresh dmask dpeak mstatus resolution0 L_ratio g vlayer vcontrast vcontrmin vfsize vimthresh vmask vpeak limsize limthresh limcontrast max_spot n0 xout0 n xout sg_name bk_use imname imfolder sub_folder mask_out m_data stack_th stack_th0 RNA_mask lim t0 dt0
n = zeros(0);
xout = zeros(0);
lim = [100:20:800];
low_th = 250;
%lim = [0:100:8000];
% decrease_th = 0.95^((lim(2)-lim(1))/100);
decrease_th = 0.85;
set(handles.hist_on,'String','Wait')
guidata(hObject, handles);

xls_tail = '.xls';
xls_name2 = [imname,'_raw',xls_tail];
k0 = strfind(imfolder,'stacks');
imfolder0 = imfolder(1:(k0(end)-1));
[stack_th,t0,dt0] = th_find(spmask(imstack,1),lim,decrease_th,low_th);
% [~,temp_I] = findpeaks(dt0);
% stack_th = lim(temp_I(1));

if (stack_th > lim(end)) ||(stack_th < lim(1)) || isempty(stack_th) || isnan(stack_th)
    stack_th = stack_th0;
end
%stack_th = 300;
%stack_th = 0;

if get(handles.presult_on,'Value') && exist([imfolder0,sub_folder,xls_name2],'file')
    [max_spot,~,~] = xlsread([imfolder0,sub_folder,xls_name2]);
    max_spot = spfilter(max_spot,10,10,1);   %%% Fine filter to get rid of noise
else
    % for I_layer = 1:get(handles.layer_slide,'Max')
    %     immask(:,:,I_layer) = im2bw(g(:,:,I_layer),vimthresh);
    % end
    % STATS = regionprops(immask,imstack*65535,'MaxIntensity');
    % max_spot = imstack(imregionalmax(g)&immask)*65535;
    mask_out = spmask(imstack,1,stack_th);   %%% 3D local maximal mask generation (Primary filter)
    %mask_out = mask_out & (imstack*65535 >= 4000);   %%% Prefilter using peak intensity threshold
    % max_spot = imstack(mask_out)*65535;
    peak_parameter = ptrack(imstack0,imstack,mask_out);   %%% Track/fit mRNA spots
    max_spot = spfilter(peak_parameter,10,10,1);   %%% Fine filter to get rid of noise
end
M_spot = max_spot(:,1).*max_spot(:,2).*max_spot(:,3)*2*pi;

RNA_mask = false(size(imstack,1),size(imstack,2),size(imstack,3));
for ii = 1:size(max_spot,1)
    x_spot = min(max(round(max_spot(ii,6)),1),size(imstack,1));
    y_spot = min(max(round(max_spot(ii,7)),1),size(imstack,2));
    RNA_mask(x_spot,y_spot,round(max_spot(ii,8))) = true;
end

%% Background calculation: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if get(handles.show_bk,'Value') && ~get(handles.background_on,'Value')
    if isempty(bk_use)
        if exist([imfolder,sg_name],'file')
            [~,sg_record,~] = xlsread([imfolder0,sg_name]);
            I_record = find(strcmp(imname,sg_record(:,1)),1);
            if (~isempty(I_record)) && (size(sg_record,2) > 1)
                bk_use = sg_record(I_record,2:end);
            end
        end
    end
    if isempty(bk_use)
        sel_bk_Callback(hObject, eventdata, handles)
    end
end

n_bk = zeros(0);   %%% Initialization of background distribution
xout_bk = zeros(0);

if ~isempty(bk_use)
    for I_bk = 1:length(bk_use)
        if ~isempty(bk_use{I_bk})
            [bk_temp,~,~] = xlsread([imfolder0,sub_folder,bk_use{I_bk}]);
            n_bk = cat(2,n_bk,bk_temp(:,1));
        end
    end
    xout_bk = bk_temp(:,2);
    n_bk = mean(n_bk,2);
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Histogram plot: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% ---------------------------------------------------------------------
axes(handles.hist_show);
%[n0,xout0] = hist([STATS.MaxIntensity],60);
[n0,xout0] = hist(M_spot,60);
%n0 = n0/sum(n0)*100;
bar(xout0,n0)
if get(handles.show_bk,'Value') && ~isempty(n_bk)
    [~,Imax] = max(n_bk);
    %[~,Imin] = min(abs(xout0-xout_bk(Imax)));
    [~,Imin] = max(n0);
    nlevel = n0(Imin);
    hold on
    plot(xout_bk,n_bk/max(n_bk)*nlevel,'r')
    hold off
end
set(gca,'YScale','log')
ylabel('#')
ylim0 = ylim;
ylim0(1) = 1;
ylim(ylim0);
xlim0 = xlim;
%%%%% ---------------------------------------------------------------------
%axes(handles.sigma_show);
% plot(max_spot(:,1),max_spot(:,2),'ro',max_spot(:,1),max_spot(:,3),'gx',max_spot(:,1),max_spot(:,4),'k+');
% ylabel('Sigma (pixel)')
%set(gca,'YAxisLocation','right','Color','none')
% xlim(xlim0);
% legend('Sigma x','Sigma y','Sigma z')
%%% =======================================================================
axes(handles.hist_zoom);
% low_th = 60000;
% low_bin = 500;
% low_th2 = 50000;
% low_bin2 = 5000;
low_th = 1.01e5;
low_bin = 1e3;
low_th2 = 1e5;
low_bin2 = 5e3;

[n,xout] = hist(M_spot,[0:low_bin:low_th]);
%[n,xout] = hist(M_spot(:,1),[0:1000:60000]);
%n = n/sum(n)*100;
bar(xout,n)
if get(handles.show_bk,'Value') && ~isempty(n_bk)
    [~,Imax] = max(n_bk);
    %[~,Imin] = min(abs(xout-xout_bk(Imax)));
    [~,Imin] = max(n);
    nlevel = n(Imin);
    hold on
    plot(xout_bk,n_bk/max(n_bk)*nlevel,'r')
    hold off
end
xlim([0,low_th2]);
%xlim([0,10000]);
set(gca,'XTick',[0:low_bin2:low_th2])%,'YScale','log')
%set(gca,'XTick',[0:5000:50000],'YScale','log')

ylabel('#')
ylim0 = ylim;
if ylim0(1) > 1 
    ylim0(1) = 1;
end
ylim(ylim0);
xlim0 = xlim;

%%%%% Statistics:
m_temp = M_spot(M_spot(:,1)<low_th,1);
m_mean = mean(m_temp);
m_std = std(m_temp);
m_median = median(m_temp);
m_peak = xout(find(n == max(n(1:(end-1))),1));
m_p =  1./(48/(m_mean/m_std)^2+1);
m_data = [m_mean,m_std,m_median,m_peak,m_p];
set(handles.text_input,'String',['Low intensity zoom in (mean = ',num2str(m_mean),', std = ',num2str(m_std),', median = ',num2str(m_median),', peak = ',num2str(m_peak),', p0 = ',num2str(m_p)])

%%%%% ---------------------------------------------------------------------
%axes(handles.sigma_zoom);
% plot(max_spot(:,1),max_spot(:,2),'ro',max_spot(:,1),max_spot(:,3),'gx',max_spot(:,1),max_spot(:,4),'k+');
% ylabel('Sigma (pixel)')
%set(gca,'YAxisLocation','right','Color','none')
% xlim(xlim0);
% legend('Sigma x','Sigma y','Sigma z')
%%%%% ---------------------------------------------------------------------
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(handles.hist_on,'String','Histogram plot')
% Update handles structure
guidata(hObject, handles);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% --- Executes on button press in Peak_on.
function Peak_on_Callback(hObject, eventdata, handles)
% hObject    handle to Peak_on (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global imstack immask dlayer dcontrast dcontrmin dfsize dimthresh dmask dpeak mstatus resolution0 L_ratio g vlayer vcontrast vcontrmin vfsize vimthresh vmask vpeak limsize limthresh limcontrast max_spot
vpeak = get(handles.Peak_on,'Value');

overlay = cf_show1(imstack,g,vlayer,vimthresh,vpeak,vcontrast,vcontrmin,vmask,get(handles.Fitted_on,'Value'));
axes(handles.stack_show);
v = axis; 
imshow(overlay)
axis(v);
% Update handles structure
guidata(hObject, handles);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% --- Executes on button press in hist_save.
function hist_save_Callback(hObject, eventdata, handles)
% hObject    handle to hist_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global imstack0 imstack immask dlayer dcontrast dcontrmin dfsize dimthresh dmask dpeak mstatus resolution0 L_ratio g vlayer vcontrast vcontrmin vfsize vimthresh vmask vpeak limsize limthresh limcontrast imname imfolder n0 xout0 n xout bk_name sg_name bk_use sub_folder max_spot m_data stack_th lim t0 dt0
hist_tail = '.fig';
tif_tail = '.tif';
xls_tail = '.xls';
hist_name = [imname,hist_tail];
tif_name = [imname,tif_tail];
xls_name0 = [imname,'_all',xls_tail];
xls_name1 = [imname,'_low',xls_tail];
xls_name2 = [imname,'_raw',xls_tail];
k0 = strfind(imfolder,'stacks');
imfolder0 = imfolder(1:(k0(end)-1));

%% Save images and data: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist([imfolder0,sub_folder],'dir')
    mkdir([imfolder0,sub_folder]);
end
hgsave([imfolder0,sub_folder,hist_name]);
saveas(gcf,[imfolder0,sub_folder,tif_name]);
if ~isempty(n0)
    xlswrite([imfolder0,sub_folder,xls_name0],[n0',xout0']);
    xlswrite([imfolder0,sub_folder,xls_name1],[n',xout']);
    xlswrite([imfolder0,sub_folder,xls_name2],max_spot);
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Save the record: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if get(handles.background_on,'Value')
        %%% Save background record: %%% =======================================
        if exist([imfolder0,bk_name],'file')
            [bk_data,bk_record,~] = xlsread([imfolder0,bk_name]);
        else
            bk_record = {};
            bk_data = [];
        end
        if ~any(strcmp(xls_name1,bk_record))
            I_bk = size(bk_record,1)+1;
        else
            I_bk = find(strcmp(xls_name1,bk_record),1);
        end
        bk_record{I_bk,1} = xls_name1;
        bk_data(I_bk,:) = m_data;
        
        if size(bk_data,1) < size(bk_record,1)
            bk_data([(size(bk_data,1)+1):size(bk_record,1)],:) = zeros(size(bk_record,1)-size(bk_data,1),size(bk_data,2));
        end
        xlswrite([imfolder0,bk_name],cat(2,bk_record,num2cell(bk_data)));
        %%% ===================================================================
    else
        %%% Save data record: %%% =============================================
        if exist([imfolder0,sg_name],'file')
            [sg_data,sg_record,~] = xlsread([imfolder0,sg_name]);
        else
            sg_record = {};
            sg_data = [];
        end
        if ~isempty(sg_record)
            I_record = find(strcmp(imname,sg_record(:,1)),1);
        else
            I_record = [];
        end
        if isempty(I_record)
            I_record = size(sg_record,1)+1;
        end
        sg_record{I_record,1} = imname;
        sg_data(I_record,:) = m_data;
        if ~isempty(bk_use)
            sg_record(I_record,2:(length(bk_use)+1)) = reshape(bk_use,1,length(bk_use));
        end
        if (length(bk_use)+1) < size(sg_record,2)
            sg_record(I_record,(length(bk_use)+2):size(sg_record,2)) = cell(1,size(sg_record,2)-(length(bk_use)+1));
        end
        
        if size(sg_data,1) < size(sg_record,1)
            sg_data([(size(sg_data,1)+1):size(sg_record,1)],:) = zeros(size(sg_record,1)-size(sg_data,1),size(sg_data,2));
        end
        xlswrite([imfolder0,sg_name],cat(2,sg_record,num2cell(sg_data)));
        %%% ===================================================================
    end
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


%%% Foci recognition output: %%% ==========================================
foci_add = '_foci';
seg_add = '_seg';
figure_tail = '.fig';
out_folder = 'Results/';
result_folder = [imfolder(1:find((imfolder(1:(end-1)) == '/')|(imfolder(1:(end-1)) == '\'),1,'last')),out_folder];
foci_bw = false(size(imstack0,1),size(imstack0,2));

for ii = 1:size(max_spot,1)
    x_spot = min(max(round(max_spot(ii,6)),1),size(imstack,1));
    y_spot = min(max(round(max_spot(ii,7)),1),size(imstack,2));
    foci_bw(x_spot,y_spot) = true;
end

foci_bw0 = imdilate(foci_bw,strel('disk',4));
bw_perim_g = bwperim(foci_bw0);
%bw_perim_WGA = bwperim(seg_bw);
new_image(:,:,1) = max(imstack0,[],3)./max(max(max(imstack0,[],3)));
new_image(:,:,2) = 0;
new_image(:,:,3) = 0;
overlay = imoverlay(new_image, bw_perim_g, [1,1,1]);
%overlay = imoverlay(overlay,bw_perim_WGA,[0,0,1]);

figure(4)
imshow(overlay)
title([result_folder,imname,', white: transcription foci recognition'],'Interpreter','none')
saveas(4,[imfolder0,sub_folder,imname,foci_add,seg_add,figure_tail]);
close(4)
%%% =======================================================================


%%% Foci recognition mask selection: %%% ==================================
th_add = '_th';
figure(3)
[AX,H1,H2] = plotyy(lim,t0,lim,dt0);
hold on
plot(stack_th*[1,1],ylim,'--')
title(['Threshold - # profile and the selected threshold value (',imfolder0,') threshold = ',num2str(stack_th)],'Interpreter','none')
xlabel('Threshold value (A.U.)')
set(get(AX(1),'Ylabel'),'String','#') 
set(get(AX(2),'Ylabel'),'String','#/#') 
%legend('# of recognized regions','#(n)/#(n-1)','Threshold value')
saveas(3,[imfolder0,sub_folder,imname,foci_add,th_add,figure_tail]);
close(3)
%%% =======================================================================


% Update handles structure
guidata(hObject, handles);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% --- Executes on button press in hist_auto.
function hist_auto_Callback(hObject, eventdata, handles)%, varargin)
% hObject    handle to hist_auto (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global imstack immask dlayer dcontrast dcontrmin dfsize dimthresh dmask dpeak mstatus resolution0 L_ratio g vlayer vcontrast vcontrmin vfsize vimthresh vmask vpeak limsize limthresh limcontrast imname imfolder
list_name = 'matchlist.xls';

set(handles.hist_auto,'String','Busy');
% Update handles structure
guidata(hObject, handles);

%imfolder = uigetdir;
% if isempty(varargin)
     imfolder_list = uigetfile_n_dir;
% else
%     imfolder_list = varargin{1};
% end

for I_folder = 1:length(imfolder_list)
    imfolder = imfolder_list{I_folder};
    if ~(imfolder(end) == '/' || imfolder(end) == '\')
        imfolder = [imfolder,'/'];
    end
    
    I_cut = strfind(imfolder,'stacks');
    if I_cut
        imname = imfolder(I_cut+7:end-1);
        imfolder = imfolder(1:I_cut+6);
        set(handles.Fitted_on,'Value',false);
        load_file_Callback(hObject, eventdata, handles, true);
        hist_on_Callback(hObject, eventdata, handles);
        hist_save_Callback(hObject, eventdata, handles);
    else           
        imfolder = [imfolder,'stacks/'];
        [~,file_list] = xlsread([imfolder,list_name]);
        if size(file_list,1) > 0
            for I_file = 1:size(file_list,1)
                imname = file_list{I_file,3};
                imname = imname(1:(end-1));
                set(handles.Fitted_on,'Value',false);
                load_file_Callback(hObject, eventdata, handles, true);
                hist_on_Callback(hObject, eventdata, handles);
                hist_save_Callback(hObject, eventdata, handles);
            end
        end
    end
    
end

set(handles.hist_auto,'String','Auto');
% Update handles structure
guidata(hObject, handles);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





% --- Executes on button press in background_on.
function background_on_Callback(hObject, eventdata, handles)
% hObject    handle to background_on (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





% --- Executes on button press in sel_bk.
function sel_bk_Callback(hObject, eventdata, handles)
% hObject    handle to sel_bk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global imstack immask dlayer dcontrast dcontrmin dfsize dimthresh dmask dpeak mstatus resolution0 L_ratio g vlayer vcontrast vcontrmin vfsize vimthresh vmask vpeak limsize limthresh limcontrast imname imfolder bk_name bk_use
k0 = strfind(imfolder,'stacks');
imfolder0 = imfolder(1:(k0(end)-1));
if exist([imfolder0,bk_name],'file')
    [~,bk_record,~] = xlsread([imfolder0,bk_name]);
    [Selection,ok] = listdlg('ListString',bk_record,'Name','Select your background profile');
    if ok
        bk_use = bk_record(Selection);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





% --- Executes on button press in show_bk.
function show_bk_Callback(hObject, eventdata, handles)
% hObject    handle to show_bk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%hist_on_Callback(hObject, eventdata, handles)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





% --- Executes on button press in presult_on.
function presult_on_Callback(hObject, eventdata, handles)
% hObject    handle to presult_on (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






% --- Executes on button press in Fitted_on.
function Fitted_on_Callback(hObject, eventdata, handles)
% hObject    handle to Fitted_on (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







function ButttonDownFcn2(src,event)
global imstack0 imstack immask dlayer dcontrast dcontrmin dfsize dimthresh dmask dpeak mstatus resolution0 L_ratio g vlayer vcontrast vcontrmin vfsize vimthresh vmask vpeak limsize limthresh limcontrast mask_out max_spot
pt = get(gca,'CurrentPoint');
x = floor(pt(1,1)+0.5);
y = floor(pt(1,2)+0.5);
xy_radius = 15;
foci_sub = imstack((y-xy_radius):(y+xy_radius),(x-xy_radius):(x+xy_radius),vlayer)*65535;

    temp_mask = false(size(mask_out));
    temp_mask((y-xy_radius):(y+xy_radius),(x-xy_radius):(x+xy_radius),vlayer) = mask_out((y-xy_radius):(y+xy_radius),(x-xy_radius):(x+xy_radius),vlayer);
    peak_parameter = ptrack(imstack0,imstack,temp_mask)   %%% Track/fit mRNA spots
    temp_spot = spfilter(peak_parameter,10,10)   %%% Fine filter to get rid of noise
    
%[X,Y] = meshgrid([-xy_radius:xy_radius],[-xy_radius:xy_radius]);
figure
subplot(1,2,1)
    surf(foci_sub);
subplot(1,2,2)
    plot(reshape(imstack(y,x,:),1,size(imstack,3))*65535,'ob')
    hold on
    plot(vlayer*[1,1],imstack(y,x,vlayer)*65535*[0,1],'--r')
    hold off
    
    
    
    
function [stack_th,t0,dt0] = th_find(mask_out,lim,decrease_th,varargin)
global imstack

if ~isempty(varargin)
    low_th = varargin{1};
else
    low_th = 0;
end
t0 = zeros(size(lim));
all_spot = imstack(mask_out)*65535;
for I_th = 1:length(lim)
    t0(I_th) = nnz(all_spot >=lim(I_th));
end
dt0 = exp(diff(log(t0)));
dt0(end+1) = nan;
stack_th1 = lim(find((dt0 == max(dt0(lim >= low_th)) & (lim >= low_th)),1));

[~,I0] = min(dt0(t0 > 10));
I1 = find((dt0(I0:end) >= decrease_th) & (lim(I0:end) >= low_th),1)+I0-1;
stack_th2 = lim(I1);

if ~isempty(stack_th2)
    stack_th = min(stack_th1,stack_th2);
else
    stack_th = stack_th1;
end
