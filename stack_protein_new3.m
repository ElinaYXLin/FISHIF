function varargout = stack_protein_new3(varargin)
% STACK_PROTEIN_NEW3 M-file for stack_protein_new3.fig
%      STACK_PROTEIN_NEW3, by itself, creates a new STACK_PROTEIN_NEW3 or raises the existing
%      singleton*.
%
%      H = STACK_PROTEIN_NEW3 returns the handle to a new STACK_PROTEIN_NEW3 or the handle to
%      the existing singleton*.
%
%      STACK_PROTEIN_NEW3('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in STACK_PROTEIN_NEW3.M with the given input arguments.
%
%      STACK_PROTEIN_NEW3('Property','Value',...) creates a new STACK_PROTEIN_NEW3 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before stack_protein_new3_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to stack_protein_new3_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help stack_protein_new3

% Last Modified by GUIDE v2.5 24-Nov-2018 11:15:42

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @stack_protein_new3_OpeningFcn, ...
                   'gui_OutputFcn',  @stack_protein_new3_OutputFcn, ...
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


% --- Executes just before stack_protein_new3 is made visible.
function stack_protein_new3_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to stack_protein_new3 (see VARARGIN)
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
if get(handles.type_Anterior,'Value')
    sub_folder = 'Histogram_protein2_A/';
else
    sub_folder = 'Histogram_protein2_P/';
end
bk_name = 'bk_record.xls';
sg_name = 'signal_record.xls';
set(handles.show_bk,'Value',false);
set(handles.presult_on,'Value',true);
set(handles.Fitted_on,'Value',false);
set(handles.ana_type,'SelectionChangeFcn',@ana_type_SelectionChangeFcn);
%set(handles.type_Anterior,'Value',false);
% set(handles.sigma_show,'Color','none');
% set(handles.sigma_zoom,'Color','none');
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Choose default command line output for stack_protein_new3
handles.output = hObject;
set(hObject,'Resize','On');
%load_file_Callback(hObject, eventdata, handles)

set(gcf,'WindowButtonDownFcn',{@ButttonDownFcn2});
% Update handles structure
guidata(hObject, handles);

% if ~isempty(varargin)
%     if length(varargin) > 1
%         set(handles.type_Posterior,'Value',varargin{2});
%     end
%     if get(handles.type_Posterior,'Value')
%         sub_folder = 'Histogram_A/';
%     else
%         sub_folder = 'Histogram/';
%     end
%     guidata(hObject, handles);
%     hist_auto_Callback(hObject, eventdata, handles, varargin{1});
% end
% UIWAIT makes stack_protein_new3 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = stack_protein_new3_OutputFcn(hObject, eventdata, handles) 
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
global imstack imstack0 immask dlayer dcontrast dcontrmin dfsize dimthresh dmask dpeak mstatus resolution0 L_ratio g vlayer vcontrast vcontrmin vfsize vimthresh vmask vpeak limsize limthresh limcontrast max_spot n0 xout0 n xout imname imfolder bk_label bk_use mask_out stack_th mask_stack
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
[num_list, file_list] = xlsread([imfolder,list_name]);
I_file = strcmp([imname,'/'],file_list(:,3));
I_file = I_file | strcmp([imname,'\'],file_list(:,3));

set(handles.load_file,'String','Busy');
% Update handles structure
guidata(hObject, handles);

lsm_stack = dir([imfolder,imname,'/',image_tail]); %%% image file list loading
resolution = num_list(I_file,9);
% RNA_channel = num_list(I_file,10);
% RNA_channel = num_list(I_file,8);
%RNA_channel = num_list(I_file,12);
RNA_channel = num_list(I_file,14);
L_ratio = (resolution/resolution0);

immax = length(lsm_stack);
tempmax = zeros(1,immax);

for I_layer = 1:immax
    temp0 = imread([imfolder,imname,'/',lsm_stack(I_layer).name]);
    temp = temp0(:,:,RNA_channel);
    imstack0(:,:,I_layer) = temp;
    imstack(:,:,I_layer) = imfilter(temp,fspecial('gaussian',3,1),'symmetric','conv');
    tempmax(I_layer) = max(max(imstack(:,:,I_layer)));
end

%imstack0 = imstack0(ceil(size(imstack0,1)*1/4):ceil(size(imstack0,1)*3/4),ceil(size(imstack0,2)*1/4):ceil(size(imstack0,2)*3/4),:);
%imstack  = imstack(ceil(size(imstack,1)*1/4):ceil(size(imstack,1)*3/4),ceil(size(imstack,2)*1/4):ceil(size(imstack,2)*3/4),:);

try
    load([imfolder(1:find((imfolder(1:end-1) == '/') | (imfolder(1:end-1) == '\'),1,'last')),'masks/',imname,'/mask.mat']);
catch
    mask_stack = true(size(imstack));
end
%     mask_stack = true(size(imstack));
mask1D = max(max(logical(mask_stack),[],3),[],1);
ELmin = find(mask1D,1);
ELmax = find(mask1D,1,'last');

rxmin1 = 3/16;
rxmax1 = 5/16;
rymin1 = 7/16;
rymax1 = 9/16;

xmin1 = round(ELmin+(ELmax-ELmin)*rxmin1);
xmax1 = round(ELmin+(ELmax-ELmin)*rxmax1);
ymin1 = round(rymin1*size(imstack,1));
ymax1 = round(rymax1*size(imstack,1));


rxmin2 = 1-rxmax1;%6/16;
rxmax2 = 1-rxmin1;%2/16;
rymin2 = rymin1;%2/8;
rymax2 = rymax1;%6/8;

xmin2 = round(ELmin+(ELmax-ELmin)*rxmin2);
xmax2 = round(ELmin+(ELmax-ELmin)*rxmax2);
ymin2 = round(rymin2*size(imstack,1));
ymax2 = round(rymax2*size(imstack,1));

Itrue = mean(mean(mean(imstack0(ymin1:ymax1,xmin1:xmax1,:)))) >= mean(mean(mean(imstack0(ymin2:ymax2,xmin2:xmax2,:))));
if ~xor(Itrue , get(handles.type_Anterior,'Value'))
    imstack0 = imstack0(ymin1:ymax1,xmin1:xmax1,:);
    imstack  = imstack(ymin1:ymax1,xmin1:xmax1,:);
    mask_stack  = mask_stack(ymin1:ymax1,xmin1:xmax1,:);
else
    imstack0 = imstack0(ymin2:ymax2,xmin2:xmax2,:);
    imstack  = imstack(ymin2:ymax2,xmin2:xmax2,:);
    mask_stack  = mask_stack(ymin2:ymax2,xmin2:xmax2,:);
end
% mask_stack = true(size(imstack));


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
global imstack imstack0 immask dlayer dcontrast dcontrmin dfsize dimthresh dmask dpeak mstatus resolution0 L_ratio g vlayer vcontrast vcontrmin vfsize vimthresh vmask vpeak limsize limthresh limcontrast max_spot n0 xout0 n xout sg_name bk_use imname imfolder sub_folder mask_out m_data stack_th stack_th0 RNA_mask lim t0 dt0 mask_stack
n = zeros(0);
xout = zeros(0);
low_th = 500;
% lim = [0:10:500];
rxy = 10;
set(handles.hist_on,'String','Wait')
guidata(hObject, handles);
xls_tail = '.xls';
xls_name2 = [imname,'_raw',xls_tail];
k0 = strfind(imfolder,'stacks');
imfolder0 = imfolder(1:(k0(end)-1));

if get(handles.type_Posterior,'Value')
    lim = [0:20:2000];
    decrease_th = 0.9^((lim(2)-lim(1))/20);
    [stack_th,t0,dt0] = th_find_low(spmask1(imstack,1),lim,decrease_th,low_th);
else
	lim = [100:100:8000];
    decrease_th = 0.85^((lim(2)-lim(1))/200);
    [stack_th,t0,dt0] = th_find(spmask1(imstack,1),lim,decrease_th,low_th);
end

% [~,temp_I] = findpeaks(dt0);
% stack_th = lim(temp_I(1));
% % [n_hist,~] = hist(imstack(mask_out)*65535,lim);
% % [~,temp_I] = findpeaks(-n_hist);
% % if isempty(temp_I)
% %     temp_I = 1;
% % end
% % stack_th = lim(temp_I(1));

if (stack_th > lim(end)) ||(stack_th < lim(1)) || isempty(stack_th) || isnan(stack_th)
    stack_th = stack_th0;
end
% stack_th = 2000;
% stack_th = 1100;
% stack_th = 340;
stack_th = 0;

if get(handles.presult_on,'Value') && exist([imfolder0,sub_folder,xls_name2],'file')
    [max_spot,~,~] = xlsread([imfolder0,sub_folder,xls_name2]);
%     max_spot = spfilter(max_spot,rxy,rxy,1);   %%% Fine filter to get rid of noise
else
    % for I_layer = 1:get(handles.layer_slide,'Max')
    %     immask(:,:,I_layer) = im2bw(g(:,:,I_layer),vimthresh);
    % end
    % STATS = regionprops(immask,imstack*65535,'MaxIntensity');
    % max_spot = imstack(imregionalmax(g)&immask)*65535;
    [mask_out,mask_out2D] = spmask1(imstack,1,stack_th);   %%% 3D local maximal mask generation (Primary filter)
    %mask_out = mask_out & (imstack*65535 >= 4000);   %%% Prefilter using peak intensity threshold
    % max_spot = imstack(mask_out)*65535;

    embryo_mask = repmat(imdilate(bwconvhull(max(logical(mask_stack),[],3)),strel('disk',20)),[1,1,size(mask_out,3)]);
    embryo_mask = embryo_mask & (~imdilate(mask_stack,strel('disk',2)));
    mask_out = mask_out & embryo_mask;
    mask_out2D = false(size(mask_out2D)) & embryo_mask;

    peak_parameter = ptrack(imstack0,imstack,mask_out,mask_out2D);   %%% Track/fit mRNA spots
    max_spot = spfilter(peak_parameter,rxy,rxy,1);   %%% Fine filter to get rid of noise
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
xy_radius = 10;
foci_sub = imstack((y-xy_radius):(y+xy_radius),(x-xy_radius):(x+xy_radius),vlayer)*65535;

    temp_mask = false(size(mask_out));
    temp_mask((y-xy_radius):(y+xy_radius),(x-xy_radius):(x+xy_radius),vlayer) = mask_out((y-xy_radius):(y+xy_radius),(x-xy_radius):(x+xy_radius),vlayer);
    [peak_parameter,xf_all] = ptrack3D(imstack0,imstack,temp_mask);   %%% Track/fit mRNA spots
    peak_parameter
    peak_parameter(:,2:4)
    temp_spot = spfilter(peak_parameter,10,10)   %%% Fine filter to get rid of noise
    
%[X,Y] = meshgrid([-xy_radius:xy_radius],[-xy_radius:xy_radius]);
figure
subplot(1,2,1)
    surf(foci_sub);
subplot(1,2,2)
    plot(reshape(imstack0(y,x,:),1,size(imstack0,3)),'ob')
%     plot(reshape(imstack(y,x,:),1,size(imstack,3))*65535,'ob')
    hold on
    plot(vlayer*[1,1],imstack0(y,x,vlayer)*[0,1],'--r')
%     plot(vlayer*[1,1],imstack(y,x,vlayer)*65535*[0,1],'--r')
    hold off
    
% % % fit_plot_cmp3D(xf_all,[y,x,vlayer],xy_radius,imstack0(:,:,vlayer));
% fit_plot_cmp(xf_all,y,x,xy_radius,imstack0(:,:,vlayer));
    
    
    

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





function [stack_th,t0,dt0] = th_find_low(mask_out,lim,decrease_th,varargin)
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
dt0 = diff(t0);
dt0(end+1) = nan;

H = fspecial('gaussian',10,4);
[~,I0] = findpeaks(conv(dt0,H(5,:)/sum(H(5,:)),'same'));
[~,I1] = findpeaks(conv(diff(conv(dt0,H(5,:)/sum(H(5,:)),'same')),H(5,:)/sum(H(5,:)),'same'));
[~,I2] = findpeaks(-conv(diff(conv(dt0,H(5,:)/sum(H(5,:)),'same')),H(5,:)/sum(H(5,:)),'same'));
lim1 = lim(I1(find(lim(I1) >= low_th,1))+2);
if isempty(lim1)
    lim1 = low_th;
end
lim2 = lim(I2(find(lim(I2) >= lim1,1))+2);
if ~isempty(lim2)
    lim0 = lim(I0(find(lim(I0) >= lim1 & lim(I0) <= lim2,1))+1);
else
    lim0 = zeros(0);
end

if ~isempty(lim0)
    stack_th = lim0;
else
    stack_th = lim1;
end





function ana_type_SelectionChangeFcn(hObject, eventdata)
 
%retrieve GUI data, i.e. the handles structure
handles = guidata(hObject); 
global imstack immask dlayer dcontrast dcontrmin dfsize dimthresh dmask dpeak mstatus resolution0 L_ratio g vlayer vcontrast vcontrmin vfsize vimthresh vmask vpeak limsize limthresh limcontrast bk_label sub_folder bk_name sg_name stack_th stack_th0

switch get(eventdata.NewValue,'Tag')   % Get Tag of selected object
    case 'type_Anterior'
        sub_folder = 'Histogram_protein_A/';
    case 'type_Posterior'
        sub_folder = 'Histogram_protein_P/';
end

type_name = get(eventdata.NewValue,'Tag');
switch_reply = questdlg(['Switch analysis type to ',type_name(find(type_name == '_',1)+1:end),'?']);

if strcmp(switch_reply,'Yes')
    stack_protein_new3_OpeningFcn(hObject, eventdata, handles)
else
    children_handle = get(handles.ana_type,'Children');
    current_handle = get(handles.ana_type,'SelectedObject');
    next_handle = children_handle(children_handle ~= current_handle);
    set(handles.ana_type,'SelectedObject',next_handle)
end
%updates the handles structure
guidata(hObject, handles);

        



function Inten_th_Callback(hObject, eventdata, handles)
% hObject    handle to Inten_th (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Inten_th as text
%        str2double(get(hObject,'String')) returns contents of Inten_th as a double


% --- Executes during object creation, after setting all properties.
function Inten_th_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Inten_th (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
