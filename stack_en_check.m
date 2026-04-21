function varargout = stack_en_check(varargin)
% STACK_EN_CHECK M-file for stack_en_check.fig
%      STACK_EN_CHECK, by itself, creates a new STACK_EN_CHECK or raises the existing
%      singleton*.
%
%      H = STACK_EN_CHECK returns the handle to a new STACK_EN_CHECK or the handle to
%      the existing singleton*.
%
%      STACK_EN_CHECK('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in STACK_EN_CHECK.M with the given input arguments.
%
%      STACK_EN_CHECK('Property','Value',...) creates a new STACK_EN_CHECK or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before stack_en_check_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to stack_en_check_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help stack_en_check

% Last Modified by GUIDE v2.5 18-Nov-2018 00:03:06

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @stack_en_check_OpeningFcn, ...
                   'gui_OutputFcn',  @stack_en_check_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1}) && ~isempty(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before stack_en_check is made visible.
function stack_en_check_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to stack_en_check (see VARARGIN)

% Parameter setting:
set(gcf,'WindowButtonDownFcn',{}); set(gcf,'KeyPressFcn',{@KeyPressFcn0,hObject,eventdata,handles});
global varargin0 scale0 in_folder input_name out_folder hist_folder hist_folder_single mask_folder mask_name channel0_name channel_name channel2_name mat_tail hist_tail fit_tail hist_add d0 ccode0 current_handle h00 ev00 resolution0
h00 = hObject;
ev00 = eventdata;
in_folder = 'stacks/';
input_name = 'matchlist.xls';
out_folder = 'Results/';
channel2_add = '_protein2';

if length(varargin) >= 2 && ~isempty(varargin{2})
    scale0 = varargin{2};
else
    scale0 = 1;
end

varargin0 = varargin;

if get(handles.ch1,'Value')
    hist_folder = 'Histogram_en/';
    hist_folder_single = 'Histogram_protein_A/';
    channel0_name = 'protein_channel';
    channel_name = 'RNA_channel';
    channel2_name = 'signal2_channel';
else
    hist_folder = ['Histogram_en4',channel2_add,'/'];
    hist_folder_single = ['Histogram_protein_A',channel2_add,'/'];
    channel0_name = 'protein2_channel';
    channel_name = 'RNA_channel';
    channel2_name = 'signal2_channel';
end
current_handle = get(handles.channel_sel,'SelectedObject');

mask_folder = 'masks/';
mask_name = 'mask';
mat_tail = '.mat';
hist_tail = '_raw.xls';
fit_tail = '_spot_fit.mat';
hist_add = '_add';
d0 = 3;
ccode0 = [0,1,0];
resolution0 = 0.083;
% Choose default command line output for stack_en_check
handles.output = hObject;
%set(hObject,'toolbar','figure');
set(hObject,'Resize','On');

%%% Add button groups:
set(handles.seg_group,'SelectionChangeFcn',@seg_group_SelectionChangeFcn);
set(gcf,'CloseRequestFcn',@th_closereq);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes stack_en_check wait for user response (see UIRESUME)
% uiwait(handles.figure1);



% --- Outputs from this function are returned to the command line.
function varargout = stack_en_check_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% --- Executes on button press in load_folder.
function load_folder_Callback(hObject, eventdata, handles)
% hObject    handle to load_folder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(gcf,'WindowButtonDownFcn',{}); set(gcf,'KeyPressFcn',{@KeyPressFcn0,hObject,eventdata,handles});
global scale0 sub_list sub_num in_folder input_name open_folder

J1 = 1;

set(handles.load_folder,'String','Wait');
% File loading:
open_folder = uigetdir;
open_folder = [open_folder,'\'];
[sub_num, sub_list] = xlsread([open_folder,in_folder,input_name]);
[M1,~] = size(sub_list);

J0 = str2num(get(handles.list_J,'String'));
M10 = str2num(get(handles.all_J,'String'));
if J0 > 0 && J0 <= M1 && M10 == 0
    J1 = J0;
end

set(handles.list_J,'String',num2str(J1));
set(handles.all_J,'String',num2str(M1));
set(handles.folder_name,'String',open_folder);
set(handles.file_name,'String',sub_list{J1,3}(1:end-1));

list_J_core(hObject, eventdata, handles);

set(handles.load_folder,'String','Load');
% Update handles structure
guidata(hObject, handles);




% --- Executes on button press in next_file.
function next_file_Callback(hObject, eventdata, handles)
% hObject    handle to next_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(gcf,'WindowButtonDownFcn',{}); set(gcf,'KeyPressFcn',{@KeyPressFcn0,hObject,eventdata,handles});
fig = gcf;
JJ = str2num(get(handles.list_J,'String'))+1;
M1 = str2num(get(handles.all_J,'String'));

if JJ <= M1
    quit_reply = questdlg('Move to the next file (Make sure you have saved the result)?');
    if strcmp(quit_reply,'Yes')
        set(handles.next_file,'String','Wait');
        % File loading:
        set(handles.list_J,'String',num2str(JJ));
        % Update handles structure
        guidata(hObject, handles);

        list_J_core(hObject, eventdata, handles);
        set(handles.next_file,'String','Next');
    end
end
% Update handles structure
guidata(hObject, handles);




function list_J_Callback(hObject, eventdata, handles)
% hObject    handle to list_J (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of list_J as text
%        str2double(get(hObject,'String')) returns contents of list_J as a double
set(gcf,'WindowButtonDownFcn',{}); set(gcf,'KeyPressFcn',{@KeyPressFcn0,hObject,eventdata,handles});
global sub_list 

fig = gcf;
JJ = str2num(get(handles.list_J,'String'));
M1 = str2num(get(handles.all_J,'String'));

if JJ <= M1
    quit_reply = questdlg(['Move to file #',sub_list{JJ,3}(1:(length(sub_list{JJ,3})-1)),' (Make sure you have saved the result)?']);
    if strcmp(quit_reply,'Yes')
        set(handles.next_file,'String','Wait');
        % Update handles structure
        guidata(hObject, handles);
        list_J_core(hObject, eventdata, handles);
        set(handles.next_file,'String','Next');
    end
end
    
% Update handles structure
guidata(hObject, handles);




% --- Executes during object creation, after setting all properties.
function list_J_CreateFcn(hObject, eventdata, handles)
% hObject    handle to list_J (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function list_J_core(hObject, eventdata, handles)
% hObject    handle to list_J (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of list_J as text
%        str2double(get(hObject,'String')) returns contents of list_J as a double
global scale0 size0 bin0 sub_list sub_num in_folder open_folder out_folder hist_folder hist_folder_single mask_folder mask_name mat_tail hist_tail fit_tail channel0_name channel_name channel2_name channel0_value channel_value channel2_value foci_list foci_list0 foci_im foci_id del_im add_im add_id b d0 ccode0 im_protein_3D im_RNA_3D im_RNA2_3D im0 im1 im2 layer_N layer_I mask_stack resolution0 L_ratio xrange
fig = gcf;
JJ = str2num(get(handles.list_J,'String'));
M1 = str2num(get(handles.all_J,'String'));

% Parameter initialization:
set(handles.contra_min,'String','0');
set(handles.contra_max,'String','1');
set(handles.nu_on,'Value',false);
set(handles.ref_on,'Value',false);

resolution = sub_num(JJ,9);
L_ratio = (resolution/resolution0);
xrange = max(5,round(3/L_ratio));

%%% Load image stack: %%% =================================================
image_type = '*.tif';
load([open_folder,out_folder,sub_list{JJ,3}(1:end-1),mat_tail],'max_image',channel0_name,channel_name,channel2_name);
size0 = size(max_image);
bin0 = round(1/scale0);
if bin0 > 1
    temp0 = zeros(0,'uint16');
    for ii = 1:size0(3)
        temp0 = cat(3,temp0,blkproc(max_image(:,:,ii), [bin0,bin0], 'mean2'));
    end
    max_image = temp0;
end

if exist(channel0_name)
    channel0_value = eval(channel0_name);
else
    channel0_value = 0;
end
if exist(channel_name)
    channel_value = eval(channel_name);
else
    channel_value = 0;
end
if exist(channel2_name)
    channel2_value = eval(channel2_name);
else
    channel2_value = 0;
end
imlist = dir([open_folder,in_folder,sub_list{JJ,3}(1:end-1),'\',image_type]); %%% get the image list from image folder

for image_I = 1:length(imlist)
    raw_im = imread([open_folder,in_folder,sub_list{JJ,3}(1:end-1),'\',imlist(image_I).name]);
    if bin0 > 1
        raw_im = imresize(raw_im,1/bin0,'nearest');
    end
    if image_I == 1
        im_protein_3D = zeros(size(raw_im,1),size(raw_im,2),length(imlist));
        im_RNA_3D = zeros(size(raw_im,1),size(raw_im,2),length(imlist));
        im_RNA2_3D = zeros(size(raw_im,1),size(raw_im,2),length(imlist));
    end
    im_protein_3D(:,:,image_I) = imfilter(double(raw_im(:,:,channel0_value)),fspecial('gaussian',3,1),'symmetric','conv');
    im_RNA_3D(:,:,image_I) = imfilter(double(raw_im(:,:,channel_value)),fspecial('gaussian',3,1),'symmetric','conv');
    im_RNA2_3D(:,:,image_I) = imfilter(double(raw_im(:,:,channel2_value)),fspecial('gaussian',3,1),'symmetric','conv');
% % %      new_DAPI3D(:,:,image_I) = imfilter(double(raw_im(:,:,DAPI_channel)),fspecial('gaussian',10,0.75),'same','conv');
end
im_protein_3D = im_protein_3D/max(im_protein_3D(:));
im_RNA_3D = im_RNA_3D/max(im_RNA_3D(:));
im_RNA2_3D = im_RNA2_3D/max(im_RNA2_3D(:));
im0 = im_protein_3D;
im1 = im_RNA_3D;
im2 = im_RNA2_3D;
%%% =======================================================================


%%% Load nuclear mask: %%% ================================================
load([open_folder,mask_folder,sub_list{JJ,3},mask_name,mat_tail],'mask_stack');
if bin0 > 1
    mask_stack = imresize(mask_stack,1/bin0,'nearest');
end
%%% =======================================================================


%%% Load foci_list: %%% ===================================================
if exist([open_folder,hist_folder,sub_list{JJ,3}(1:(end-1)),hist_tail])
    [foci_list,~,~] = xlsread([open_folder,hist_folder,sub_list{JJ,3}(1:(end-1)),hist_tail]);
else
    [foci_list,~,~] = xlsread([open_folder,hist_folder,sub_list{JJ,3}(1:(end-1)),hist_tail,'x']);
end
if exist([open_folder,hist_folder_single,sub_list{JJ,3}(1:(end-1)),fit_tail])
    load([open_folder,hist_folder_single,sub_list{JJ,3}(1:(end-1)),fit_tail],'b');
else
    b = 1;
end
foci_list0 = foci_list;
foci_list = foci_combine(foci_list);
if bin0 > 1
    foci_list(:,6:7) = foci_list(:,6:7)/bin0;
end
foci_xyz = max(round(foci_list(:,6:8)),1);
id_list = sub2ind(size(im_protein_3D),foci_xyz(:,1),foci_xyz(:,2),foci_xyz(:,3));
%%% foci image generation:
foci_im = false(size(im_protein_3D));
foci_im(id_list) = true;
for ii = 1:size(foci_im,3)
    foci_im(:,:,ii) = bwthicken(foci_im(:,:,ii),2);
end
%%% foci indices matrix generation:
foci_id = double(foci_im);
foci_id(id_list) = 1:numel(id_list);
for ii = size(foci_id,3)
    foci_st = regionprops(foci_im(:,:,ii),foci_id(:,:,ii),'MaxIntensity');
    id_list = [0,[foci_st.MaxIntensity]];
    foci_id(:,:,ii) = id_list(bwlabel(foci_im(:,:,ii))+1);
end
% Recording matrices initialization:
del_im = false(size(im_protein_3D));
add_im = false(size(im_protein_3D));
% del_id = false(size(im0));
add_id = false(size(im_protein_3D));
%%% =======================================================================


%%% Image output: %%% =====================================================
layer_N = length(imlist);
layer_I = 1;
set(handles.N_layer,'String',num2str(layer_N))
set(handles.I_layer,'String',num2str(layer_I))
set(handles.slider_layer,'Max',layer_N);
set(handles.slider_layer,'Min',1);
set(handles.slider_layer,'SliderStep',[1./(layer_N-1+(layer_N == 1)),1./(layer_N-1+(layer_N == 1))]);
set(handles.slider_layer,'Value',layer_I);

% image output
axes(handles.im_show);
overlay = im_overlay(im0,im1,im2,foci_im,del_im,add_im,mask_stack,layer_I,handles);
imshow(overlay)

if get(handles.nu_on,'Value')
hold on
text(foci_list(:,7)+d0,foci_list(:,6)+d0,num2str(prod(foci_list(:,1:3),2)*2*pi/b),'Color',ccode0)
end
%%% =======================================================================


% Update handles structure
guidata(hObject, handles);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




function contra_min_Callback(hObject, eventdata, handles)
% hObject    handle to contra_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of contra_min as text
%        str2double(get(hObject,'String')) returns contents of contra_min as a double
global im_protein_3D im_RNA_3D im_RNA2_3D im0 im1 im2 foci_im del_im add_im mask_stack foci_list b d0 ccode0 layer_I
cmin = str2num(get(handles.contra_min,'String'));
cmax = str2num(get(handles.contra_max,'String'));
if isempty(cmin) || cmin < 0 || cmin > cmax
    cmin = 0;
end
set(handles.contra_min,'String', num2str(cmin));
im0 = (im_protein_3D-cmin)/(cmax-cmin);
im1 = (im_RNA_3D-cmin)/(cmax-cmin);
im2 = (im_RNA2_3D-cmin)/(cmax-cmin);

% image output
axes(handles.im_show);
overlay = im_overlay(im0,im1,im2,foci_im,del_im,add_im,mask_stack,layer_I,handles);
v = axis; 
imshow(overlay)

if get(handles.nu_on,'Value')
hold on
text(foci_list(:,7)+d0,foci_list(:,6)+d0,num2str(prod(foci_list(:,1:3),2)*2*pi/b),'Color',ccode0)
end

axis(v);



% --- Executes during object creation, after setting all properties.
function contra_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to contra_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function contra_max_Callback(hObject, eventdata, handles)
% hObject    handle to contra_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of contra_max as text
%        str2double(get(hObject,'String')) returns contents of contra_max as a double
global im_protein_3D im_RNA_3D im_RNA2_3D im0 im1 im2 foci_im del_im add_im mask_stack foci_list b d0 ccode0 layer_I
cmin = str2num(get(handles.contra_min,'String'));
cmax = str2num(get(handles.contra_max,'String'));
if isempty(cmax) || cmax < cmin || cmax > 1
    cmax = 1;
end
set(handles.contra_max,'String', num2str(cmax));
im0 = (im_protein_3D-cmin)/(cmax-cmin);
im1 = (im_RNA_3D-cmin)/(cmax-cmin);
im2 = (im_RNA2_3D-cmin)/(cmax-cmin);

% image output
axes(handles.im_show);
overlay = im_overlay(im0,im1,im2,foci_im,del_im,add_im,mask_stack,layer_I,handles);
v = axis; 
imshow(overlay)

if get(handles.nu_on,'Value')
hold on
text(foci_list(:,7)+d0,foci_list(:,6)+d0,num2str(prod(foci_list(:,1:3),2)*2*pi/b),'Color',ccode0)
end

axis(v);



% --- Executes during object creation, after setting all properties.
function contra_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to contra_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





function seg_group_SelectionChangeFcn(hObject, eventdata)
 
%retrieve GUI data, i.e. the handles structure
handles = guidata(hObject); 
switch get(eventdata.NewValue,'Tag')   % Get Tag of selected object
    case 'sel_cls'
      %execute this code when "select/clear" button is selected
        axes(handles.im_show);
        zoom off
        pan off
        set(gcf,'WindowButtonDownFcn',{@ButttonDownFcn1,handles});

    case 'del_cls'
      %execute this code when "delete/clear" button is selected
        axes(handles.im_show);
        zoom off
        pan off
        set(gcf,'WindowButtonDownFcn',{@ButttonDownFcn2,handles});
end
%updates the handles structure
guidata(hObject, handles);




function ButttonDownFcn1(src,event,handles)
global sub_list open_folder out_folder hist_folder mat_tail hist_tail channel_name channel2_name im0 im1 im2 foci_list foci_im foci_id del_im add_im add_id b d0 ccode0 layer_I mask_stack
b_type = get(gcf,'SELECTIONTYPE');
pt = get(gca,'CurrentPoint');
x = floor(pt(1,1)+0.5);
y = floor(pt(1,2)+0.5);
axes(handles.im_show);
v = axis; 

if strcmpi(b_type,'normal')
    if ~foci_im(y,x,layer_I) && ~add_im(y,x,layer_I) && x >= v(1) && x <= v(1)+v(2) && y >= v(3) && y <= v(3)+v(4)
        add_id(y,x,layer_I) = true;
    %     add_im = bwthicken(add_id,4);
        rr = 50;
        imin = max(1,y-rr);
        imax = min(size(im0,1),y+rr);
        jmin = max(1,x-rr);
        jmax = min(size(im0,2),x+rr);
        add_im(imin:imax,jmin:jmax,layer_I) = bwthicken(add_id(imin:imax,jmin:jmax,layer_I),1);
    elseif ~foci_im(y,x,layer_I) && add_im(y,x,layer_I) && x >= v(1) && x <= v(1)+v(2) && y >= v(3) && y <= v(3)+v(4)
        add_label = bwlabel(add_im(:,:,layer_I));
        add_id(:,:,layer_I) = add_id(:,:,layer_I) & (add_label ~= add_label(y,x));
        add_im(:,:,layer_I) = add_im(:,:,layer_I) & (add_label ~= add_label(y,x));
    end
elseif strcmpi(b_type,'alt') && x >= v(1) && x <= v(1)+v(2) && y >= v(3) && y <= v(3)+v(4)
    set(handles.ref_on,'Value',~get(handles.ref_on,'Value'))
    ref_on_Callback(handles.ref_on,event,handles);
    set(handles.ref2_on,'Value',~get(handles.ref2_on,'Value'))
    ref2_on_Callback(handles.ref2_on,event,handles);
end    

%%% Image display
axes(handles.im_show);
overlay = im_overlay(im0,im1,im2,foci_im,del_im,add_im,mask_stack,layer_I,handles);
imshow(overlay)

if get(handles.nu_on,'Value')
hold on
text(foci_list(:,7)+d0,foci_list(:,6)+d0,num2str(prod(foci_list(:,1:3),2)*2*pi/b),'Color',ccode0)
end

axis(v);




function ButttonDownFcn2(src,event,handles)
global sub_list open_folder out_folder hist_folder mat_tail hist_tail channel_name channel2_name im0 im1 im2 foci_list foci_im foci_id del_im add_im add_id b d0 ccode0 layer_I mask_stack
b_type = get(gcf,'SELECTIONTYPE');
pt = get(gca,'CurrentPoint');
x = floor(pt(1,1)+0.5);
y = floor(pt(1,2)+0.5);
axes(handles.im_show);
v = axis; 

if strcmpi(b_type,'normal')
    if foci_im(y,x,layer_I) && ~del_im(y,x,layer_I) && x >= v(1) && x <= v(1)+v(2) && y >= v(3) && y <= v(3)+v(4)
        del_im(:,:,layer_I) = del_im(:,:,layer_I) | (foci_id(:,:,layer_I) == foci_id(y,x,layer_I));
    elseif del_im(y,x,layer_I) && x >= v(1) && x >= v(1) && x <= v(1)+v(2) && y >= v(3) && y <= v(3)+v(4)
        del_im(:,:,layer_I) = del_im(:,:,layer_I) & (foci_id(:,:,layer_I) ~= foci_id(y,x,layer_I));
    end
elseif strcmpi(b_type,'alt') && x >= v(1) && x <= v(1)+v(2) && y >= v(3) && y <= v(3)+v(4)
    set(handles.ref_on,'Value',~get(handles.ref_on,'Value'))
    ref_on_Callback(handles.ref_on,event,handles);
    set(handles.ref2_on,'Value',~get(handles.ref2_on,'Value'))
    ref2_on_Callback(handles.ref2_on,event,handles);
end    

%%% Image display
axes(handles.im_show);
overlay = im_overlay(im0,im1,im2,foci_im,del_im,add_im,mask_stack,layer_I,handles);
imshow(overlay)

if get(handles.nu_on,'Value')
hold on
text(foci_list(:,7)+d0,foci_list(:,6)+d0,num2str(prod(foci_list(:,1:3),2)*2*pi/b),'Color',ccode0)
end

axis(v);




% --- Executes on button press in auto_sel.
function auto_sel_Callback(hObject, eventdata, handles)
% hObject    handle to auto_sel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global sub_list sub_num open_folder out_folder hist_folder mat_tail hist_tail channel_name channel2_name im0 im1 im2 foci_list foci_im foci_id del_im add_im add_id b d0 ccode0 layer_I mask_stack
cmin = str2num(get(handles.contra_min,'String'));
cmax = str2num(get(handles.contra_max,'String'));
cth = cmax;

JJ = str2num(get(handles.list_J,'String'));
d0 = 10;
Nbin1 = sub_num(JJ,2);
Mdim = sub_num(JJ,3);

if (Nbin1-floor(Nbin1))*10 >= 1
    Nbin1(2) = (Nbin1-floor(Nbin1))*10;
    Nbin1(1) = floor(Nbin1(1));
    Nbin1 = round(Nbin1);
    Mdim(2) = (Mdim-floor(Mdim))*10;
    Mdim(1) = floor(Mdim(1));
    Mdim = round(Mdim);
end

size0 = size(im0);
im_manual = imdilate(bwconvhull(logical(mask_stack)),strel('disk',15));
for tt = 1:length(Nbin1)
    for ss = 1:(Nbin1(1)-1)
        im_manual(:,(ss*size0(Mdim(tt))/Nbin1(tt)-d0):(ss*size0(Mdim(tt))/Nbin1(tt)+d0),:) = false;
    end
end

% % im1 = im0 >= cth;
% % prop0 = regionprops(im1,foci_im,'MaxIntensity');
% % im2 = ismember(bwlabel(im1),find([prop0.MaxIntensity] == 0));
im0f = imfilter(im0,fspecial('gaussian',3,1),'symmetric','conv');
add_id = imregionalmax(im0f) & im0f >= cth & ~foci_im;% & im_manual;
% add_id = imregionalmax(imfilter(im0,fspecial('gaussian',3,1),'symmetric','conv')) & im0 >= cth & ~foci_im;% & im_manual;
add_im = bwthicken(add_id,1);

%%% Image display
overlay = im_overlay(im0,im1,im2,foci_im,del_im,add_im,mask_stack,layer_I,handles);
axes(handles.im_show);
v = axis; 
imshow(overlay)

if get(handles.nu_on,'Value')
hold on
text(foci_list(:,7)+d0,foci_list(:,6)+d0,num2str(prod(foci_list(:,1:3),2)*2*pi/b),'Color',ccode0)
end

axis(v);





% --- Executes on button press in auto_desel.
function auto_desel_Callback(hObject, eventdata, handles)
% hObject    handle to auto_desel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global sub_list open_folder out_folder hist_folder mat_tail hist_tail channel_name channel2_name im0 im1 im2 foci_list foci_im foci_id del_im add_im add_id b d0 ccode0 layer_I mask_stack
cmin = str2num(get(handles.contra_min,'String'));
cmax = str2num(get(handles.contra_max,'String'));
cth = cmin;

% del_id = imregionalmax(imfilter(im0,fspecial('gaussian',3,1),'symmetric','conv')) & im0 <= cth;
% prop0 = regionprops(foci_im,del_id,'MaxIntensity');
% del_im = ismember(bwlabel(foci_im),find([prop0.MaxIntensity]));
for ii = 1:size(foci_im,3)
    prop0 = regionprops(foci_im(:,:,layer_I),imfilter(im0(:,:,layer_I),fspecial('gaussian',3,1),'symmetric','conv'),'MaxIntensity');
    del_im(:,:,layer_I) = ismember(bwlabel(foci_im(:,:,layer_I)),find([prop0.MaxIntensity] <= cth));
end

%%% Image display
overlay = im_overlay(im0,im1,im2,foci_im,del_im,add_im,mask_stack,layer_I,handles);
axes(handles.im_show);
v = axis; 
imshow(overlay)

if get(handles.nu_on,'Value')
hold on
text(foci_list(:,7)+d0,foci_list(:,6)+d0,num2str(prod(foci_list(:,1:3),2)*2*pi/b),'Color',ccode0)
end

axis(v);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





% --- Executes on button press in save_seg.
function save_seg_Callback(hObject, eventdata, handles)
% hObject    handle to save_seg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global sub_list open_folder out_folder hist_folder mat_tail hist_tail hist_add channel_name channel2_name channel_value channel2_value im0 im1 im2 foci_list foci_im foci_id del_im add_im add_id b d0 ccode0 layer_I mask_stack
outim_tail = '_foci_seg.fig';

quit_reply = questdlg('Save the spot recognition result?');
if strcmp(quit_reply,'Yes')
    set(handles.save_seg,'String','Wait');
    % Update handles structure
    guidata(hObject, handles);
    
    %%% Data storage
    JJ = str2num(get(handles.list_J,'String'));
    del_list = unique(foci_id(del_im));
    id_old = 1:size(foci_list,1);
    id_new = setdiff(id_old, del_list);
    foci_list = foci_list(id_new,:);
    delete([open_folder,hist_folder,sub_list{JJ,3}(1:(end-1)),hist_tail]);
    xlswrite([open_folder,hist_folder,sub_list{JJ,3}(1:(end-1)),hist_tail],foci_list);
    
    if nnz(add_id)
        save([open_folder,hist_folder,sub_list{JJ,3}(1:(end-1)),hist_add,mat_tail],'add_id','channel_value')
    end
    
    % Foci new id assignment:
    id_list = [0,zeros(size(id_old))];
    id_list(id_new+1) = 1:length(id_new);
    foci_id = id_list(foci_id+1);
    foci_im = logical(foci_id);

    % Recording matrices reinitialization:
    del_im = false(size(im0));

    % image output
    axes(handles.im_show);
    overlay = im_overlay(im0,im1,im2,foci_im,del_im,add_im,mask_stack,layer_I,handles);
    v = axis; 
    imshow(overlay)

    if get(handles.nu_on,'Value')
    hold on
    text(foci_list(:,7)+d0,foci_list(:,6)+d0,num2str(prod(foci_list(:,1:3),2)*2*pi/b),'Color',ccode0)
    end

    axis(v);
    
    h = figure;
    overlay = im_overlay(im0,im1,im2,foci_im,del_im,add_im,mask_stack,0,handles);
    imshow(overlay)

    if get(handles.nu_on,'Value')
    hold on
    text(foci_list(:,7)+d0,foci_list(:,6)+d0,num2str(prod(foci_list(:,1:3),2)*2*pi/b),'Color',ccode0)
    end

    title([open_folder,out_folder,sub_list{JJ,3}(1:(end-1)),', white: protein spots recognition'])
    saveas(h,[open_folder,hist_folder,sub_list{JJ,3}(1:(end-1)),outim_tail])
    close(h)
    
    set(handles.save_seg,'String','Save');
    % Update handles structure
    guidata(hObject, handles);
end
% Update handles structure
guidata(hObject, handles);





% --- Executes on button press in cancel_seg.
function cancel_seg_Callback(hObject, eventdata, handles)
% hObject    handle to cancel_seg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(gcf,'WindowButtonDownFcn',{}); set(gcf,'KeyPressFcn',{@KeyPressFcn0,hObject,eventdata,handles});
global sub_list open_folder out_folder hist_folder mat_tail hist_tail hist_add channel_name channel2_name im0 im1 im2 foci_list foci_im foci_id del_im add_im add_id b d0 ccode0 layer_I mask_stack
quit_reply = questdlg('Cancel selection/deletion?');
if strcmp(quit_reply,'Yes')
    set(handles.cancel_seg,'String','Wait');
    % Update handles structure
    guidata(hObject, handles);
    
    % Recording matrices initialization:
    del_im = false(size(im0));
    add_im = false(size(im0));
    % del_id = false(size(im0));
    add_id = false(size(im0));

    % image output
    axes(handles.im_show);
    overlay = im_overlay(im0,im1,im2,foci_im,del_im,add_im,mask_stack,layer_I,handles);
    v = axis; 
    imshow(overlay)

    if get(handles.nu_on,'Value')
    hold on
    text(foci_list(:,7)+d0,foci_list(:,6)+d0,num2str(prod(foci_list(:,1:3),2)*2*pi/b),'Color',ccode0)
    end

    axis(v);
    
    set(handles.cancel_seg,'String','Cancel');
end
%%% Update handles structure
guidata(hObject, handles);





% --- Executes on button press in apply_recog.
function apply_recog_Callback(hObject, eventdata, handles)
% hObject    handle to apply_recog (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global size0 bin0 sub_list open_folder in_folder out_folder hist_folder mat_tail hist_tail hist_add channel_name channel2_name channel_value channel2_value im0 im1 im2 foci_list foci_list0 foci_im foci_id del_im add_im add_id b d0 ccode0 layer_I mask_stack L_ratio xrange
quit_reply = questdlg('Apply recognition/deletion?');
if strcmp(quit_reply,'Yes')
    set(handles.apply_recog,'String','Wait');
    % Update handles structure
    guidata(hObject, handles);
    
    % Foci deletion:
    del_list = unique(foci_id(del_im));
    id_old = 1:size(foci_list,1);
    id_new = setdiff(id_old, del_list);
    foci_list = foci_list(id_new,:);
    foci_list0 = foci_list;
    del_im = false(size(im0));
    % Foci new id assignment:
    id_list = [0,zeros(size(id_old))];
    id_list(id_new+1) = 1:length(id_new);
    foci_id = id_list(foci_id+1);
    foci_im = logical(foci_id);
    
    
    % Foci recognition:
    if nnz(add_id)
        JJ = str2num(get(handles.list_J,'String'));
        
        foci_list1 = foci_list;
        if bin0 > 1
            foci_list1(:,6:7) = foci_list(:,6:7)*bin0;
            add_id1 = false(size0(1:2));
            [x0,y0] = find(add_id);
            x0 = round(x0*bin0-bin0/2);
            y0 = round(y0*bin0-bin0/2);
            i00 = sub2ind(size0(1:2),x0,y0);
            add_id1(i00) = true;
        else
            add_id1 = add_id;
        end
        
        new_list = manual_foci3D([open_folder,in_folder,sub_list{JJ,3}],channel_value,add_id1,foci_list1,[],xrange);
        
        if bin0 > 1
            new_list(:,6:7) = new_list(:,6:7)/bin0;
        end
        
        foci_list = [foci_list;new_list];
        foci_list0 = foci_list;
        foci_list = foci_combine(foci_list);
        foci_xyz = max(round(foci_list(:,6:8)),1);
        id_list = sub2ind(size(im0),foci_xyz(:,1),foci_xyz(:,2),foci_xyz(:,3));
        %%% foci image generation:
        foci_im = false(size(im0));
        foci_im(id_list) = true;
        for ii = 1:size(foci_im,3)
            foci_im(:,:,ii) = bwthicken(foci_im(:,:,ii),4);
        end
        %%% foci indices matrix generation:
        foci_id = double(foci_im);
        foci_id(id_list) = 1:numel(id_list);
        for ii = size(foci_id,3)
            foci_st = regionprops(foci_im(:,:,ii),foci_id(:,:,ii),'MaxIntensity');
            id_list = [0,[foci_st.MaxIntensity]];
            foci_id(:,:,ii) = id_list(bwlabel(foci_im(:,:,ii))+1);
        end
        %%% Recording matrices reinitialization:
        add_im = false(size(im0));
        add_id = false(size(im0));
    end

    
    % image output
    axes(handles.im_show);
    overlay = im_overlay(im0,im1,im2,foci_im,del_im,add_im,mask_stack,layer_I,handles);
    v = axis; 
    imshow(overlay)

    if get(handles.nu_on,'Value')
    hold on
    text(foci_list(:,7)+d0,foci_list(:,6)+d0,num2str(prod(foci_list(:,1:3),2)*2*pi/b),'Color',ccode0)
    end

    axis(v);
    
    set(handles.apply_recog,'String','Apply');
end
%%% Update handles structure
guidata(hObject, handles);






function overlay = im_overlay(im0,im1,im2,foci_im,del_im,add_im,mask_stack,layer_I,handles)
if layer_I > 0
    foci_per = bwperim(foci_im(:,:,layer_I));
    del_per = bwperim(del_im(:,:,layer_I));
    add_per = bwperim(add_im(:,:,layer_I));
    mask_per = bwperim(mask_stack(:,:,layer_I));

    overlay(:,:,2) = im0(:,:,layer_I)+foci_per-del_per+add_per+0.5*mask_per;
    overlay(:,:,1) = get(handles.ref_on,'Value')*im1(:,:,layer_I)+foci_per-del_per;
    overlay(:,:,3) = get(handles.ref2_on,'Value')*im2(:,:,layer_I)+foci_per+del_per+0.5*mask_per;
else
    foci_per = bwperim(max(foci_im,[],3));
    del_per = bwperim(max(del_im,[],3));
    add_per = bwperim(max(add_im,[],3));
    mask_per = bwperim(max(mask_stack,[],3));

    overlay(:,:,2) = max(im0,[],3)+foci_per-del_per+add_per+0.5*mask_per;
    overlay(:,:,1) = get(handles.ref_on,'Value')*max(im1,[],3)+foci_per-del_per;
    overlay(:,:,3) = get(handles.ref2_on,'Value')*max(im2,[],3)+foci_per+del_per+0.5*mask_per;
end




function th_closereq(hObject, eventdata, handles)
quit_reply = questdlg('Quit the program (Make sure you have saved the result)?');
if strcmp(quit_reply,'Yes')
    clear global
    delete(gcf)
else
    return
end




function out_list = foci_combine(in_list)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A function to refine foci list and combine adjacent foci spots: %%%%%%%%
%% out_list: refined foci list.
%% in_list: raw foci list.
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Lth = 1.5;
isrun = true;
z_th = 0;
out_list = in_list;

while isrun
    foci_xyz = max(round(in_list(:,6:8)),1);
    foci_xyz = foci_xyz.*repmat([1,1,0.35/0.083],size(foci_xyz,1),1);
    dxy = squareform(pdist(foci_xyz));
    id_pair = sparse(tril(dxy <= Lth,-1));

    if nnz(id_pair)
        good_id = true(1,size(out_list,1));
        [N_comp, id_comp] = graphconncomp(id_pair,'Weak',true);
        Nid = hist(id_comp,1:N_comp);
        for I_comp = find(Nid > 1)
            id_list0 = find(id_comp == I_comp);
            I00 = prod(in_list(id_list0,1:3),2);
            z00 = in_list(id_list0,8);
            [~,IX0] = max(I00);
            z0max = z00(IX0);
            id_list = id_list0((z00 >= z0max-z_th) & (z00 <= z0max+z_th));
            good_id(setdiff(id_list0,id_list)) = false;
            
            if length(id_list) > 1
                sigmax = sqrt(sum(in_list(id_list,2).^2));
                sigmay = sqrt(sum(in_list(id_list,3).^2));
                I0 = prod(in_list(id_list,1:3),2);
                Iall = sum(I0)./sigmax./sigmay;
                xmean = sum(in_list(id_list,6).*I0)./sum(I0);
                ymean = sum(in_list(id_list,7).*I0)./sum(I0);
                z0 = in_list(id_list,8);
                [~,IX] = max(I0);
                zmax = z0(IX);
                out_list(id_list(1),[1:3,6:8]) = [Iall,sigmax,sigmay,xmean,ymean,zmax];
                good_id(id_list(2:end)) = false;
            end
        end
        out_list = out_list(good_id,:);
        in_list = out_list;
    else
        isrun = false;
    end
end


% --- Executes when selected object is changed in channel_sel.
function channel_sel_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in channel_sel 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
 
%retrieve GUI data, i.e. the handles structure
% handles = guidata(hObject); 
global current_handle h00 ev00 varargin0

type_name = get(eventdata.NewValue,'Tag');
switch_reply = questdlg(['Switch to ',type_name,'?']);

if strcmp(switch_reply,'Yes')
    stack_en_check_OpeningFcn(h00, ev00, handles, varargin0{:})
else
% %     children_handle = get(handles.ana_type,'Children');
% %     current_handle = get(handles.ana_type,'SelectedObject');
% %     next_handle = children_handle(children_handle ~= current_handle);
    set(handles.channel_sel,'SelectedObject',current_handle)
end
%updates the handles structure
guidata(hObject, handles);


% --- Executes on button press in nu_on.
function nu_on_Callback(hObject, eventdata, handles)
% hObject    handle to nu_on (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of nu_on
global im0 im1 im2 foci_im del_im add_im foci_list b d0 ccode0 text0 mask_stack layer_I

% image output
axes(handles.im_show);
overlay = im_overlay(im0,im1,im2,foci_im,del_im,add_im,mask_stack,layer_I,handles);
v = axis; 
imshow(overlay)

if get(handles.nu_on,'Value')
    hold on
    text0 = text(foci_list(:,7)+d0,foci_list(:,6)+d0,num2str(prod(foci_list(:,1:3),2)*2*pi/b),'Color',ccode0);
end

if ~get(handles.nu_on,'Value')
    delete(text0)
end

axis(v);


% --- Executes on button press in ref_on.
function ref_on_Callback(hObject, eventdata, handles)
% hObject    handle to ref_on (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ref_on

global im0 im1 im2 foci_im del_im add_im foci_list b d0 ccode0 mask_stack layer_I

% image output
axes(handles.im_show);
overlay = im_overlay(im0,im1,im2,foci_im,del_im,add_im,mask_stack,layer_I,handles);
v = axis; 
imshow(overlay)

if get(handles.nu_on,'Value')
hold on
text(foci_list(:,7)+d0,foci_list(:,6)+d0,num2str(prod(foci_list(:,1:3),2)*2*pi/b),'Color',ccode0)
end

axis(v);



% --- Executes on button press in ref2_on.
function ref2_on_Callback(hObject, eventdata, handles)
% hObject    handle to ref2_on (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ref2_on

global im0 im1 im2 foci_im del_im add_im foci_list b d0 ccode0 mask_stack layer_I

% image output
axes(handles.im_show);
overlay = im_overlay(im0,im1,im2,foci_im,del_im,add_im,mask_stack,layer_I,handles);
v = axis; 
imshow(overlay)

if get(handles.nu_on,'Value')
hold on
text(foci_list(:,7)+d0,foci_list(:,6)+d0,num2str(prod(foci_list(:,1:3),2)*2*pi/b),'Color',ccode0)
end

axis(v);



% --- Executes on slider movement.
function slider_layer_Callback(hObject, eventdata, handles)
% hObject    handle to slider_layer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

temp_I = round(get(handles.slider_layer,'Value'));
set(handles.I_layer,'String', num2str(temp_I));

%%% Update handles structure
guidata(hObject, handles);

I_layer_Callback(hObject, eventdata, handles);



% --- Executes during object creation, after setting all properties.
function slider_layer_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_layer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function I_layer_Callback(hObject, eventdata, handles)
% hObject    handle to I_layer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of I_layer as text
%        str2double(get(hObject,'String')) returns contents of I_layer as a double

global im0 im1 im2 foci_im del_im add_im foci_list b d0 ccode0 mask_stack layer_I
layer_I = round(str2num(get(handles.I_layer,'String')));
minT = get(handles.slider_layer,'Min');
maxT = get(handles.slider_layer,'Max');
if (isempty(layer_I) || layer_I < minT)
    layer_I = minT;
elseif layer_I > maxT
    layer_I = maxT;
end
set(handles.slider_layer,'Value', layer_I);
set(handles.I_layer,'String', num2str(layer_I));

% image output
axes(handles.im_show);
overlay = im_overlay(im0,im1,im2,foci_im,del_im,add_im,mask_stack,layer_I,handles);
v = axis; 
imshow(overlay)

if get(handles.nu_on,'Value')
hold on
text(foci_list(:,7)+d0,foci_list(:,6)+d0,num2str(prod(foci_list(:,1:3),2)*2*pi/b),'Color',ccode0)
end

axis(v);

%%% Update handles structure
guidata(hObject, handles);




% --- Executes during object creation, after setting all properties.
function I_layer_CreateFcn(hObject, eventdata, handles)
% hObject    handle to I_layer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
