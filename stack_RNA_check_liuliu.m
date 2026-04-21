function varargout = stack_RNA_check_liuliu(varargin)
% STACK_RNA_CHECK_LIULIU M-file for stack_RNA_check_liuliu.fig
%      STACK_RNA_CHECK_LIULIU, by itself, creates a new STACK_RNA_CHECK_LIULIU or raises the existing
%      singleton*.
%
%      H = STACK_RNA_CHECK_LIULIU returns the handle to a new STACK_RNA_CHECK_LIULIU or the handle to
%      the existing singleton*.
%
%      STACK_RNA_CHECK_LIULIU('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in STACK_RNA_CHECK_LIULIU.M with the given input arguments.
%
%      STACK_RNA_CHECK_LIULIU('Property','Value',...) creates a new STACK_RNA_CHECK_LIULIU or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before stack_RNA_check_liuliu_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to stack_RNA_check_liuliu_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help stack_RNA_check_liuliu

% Last Modified by GUIDE v2.5 12-Feb-2015 16:42:29

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @stack_RNA_check_liuliu_OpeningFcn, ...
                   'gui_OutputFcn',  @stack_RNA_check_liuliu_OutputFcn, ...
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


% --- Executes just before stack_RNA_check_liuliu is made visible.
function stack_RNA_check_liuliu_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to stack_RNA_check_liuliu (see VARARGIN)

% Parameter setting:
set(gcf,'WindowButtonDownFcn',{}); set(gcf,'KeyPressFcn',{@KeyPressFcn0,hObject,eventdata,handles});
global in_folder input_name out_folder hist_folder mask_folder mask_name channel_name mat_tail hist_tail hist_add
in_folder = 'stacks/';
input_name = 'matchlist.xls';
out_folder = 'Results/';
hist_folder = 'Histogram/';
mask_folder = 'masks/';
mask_name = 'mask';
channel_name = 'RNA_channel';
mat_tail = '.mat';
hist_tail = '_raw.xls';
hist_add = '_add';

% Choose default command line output for stack_RNA_check_liuliu
handles.output = hObject;
%set(hObject,'toolbar','figure');
set(hObject,'Resize','On');

%%% Add button groups:
set(handles.seg_group,'SelectionChangeFcn',@seg_group_SelectionChangeFcn);
set(gcf,'CloseRequestFcn',@th_closereq);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes stack_RNA_check_liuliu wait for user response (see UIRESUME)
% uiwait(handles.figure1);



% --- Outputs from this function are returned to the command line.
function varargout = stack_RNA_check_liuliu_OutputFcn(hObject, eventdata, handles) 
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
global sub_list sub_num in_folder input_name open_folder
J1 = 1;

set(handles.load_folder,'String','Wait');
% File loading:
open_folder = uigetdir;
open_folder = [open_folder,'\'];
[sub_num, sub_list] = xlsread([open_folder,in_folder,input_name]);
[M1,~] = size(sub_list);
set(handles.list_J,'String',num2str(J1));
set(handles.all_J,'String',num2str(M1));
set(handles.folder_name,'String',open_folder);

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
global zmin zmax sub_list sub_num open_folder in_folder out_folder hist_folder mask_folder mask_name mat_tail hist_tail channel_name channel_value im0 im_now foci_list foci_list0 foci_im foci_id del_im add_im add_id mask2D
fig = gcf;
JJ = str2num(get(handles.list_J,'String'));
M1 = str2num(get(handles.all_J,'String'));

% Parameter initialization:
set(handles.contra_min,'String','0');
set(handles.contra_max,'String','1');
set(handles.file_name,'String',sub_list{JJ,3}(1:end-1));

% Load max_image:
% load([open_folder,out_folder,sub_list{JJ,3}(1:end-1),mat_tail],'max_image',channel_name);
% channel_value = eval(channel_name);

channel_value = sub_num(JJ,10);
fname = dir([open_folder,in_folder,sub_list{JJ,3},'*.tif']);

if size(sub_num,2) >= 13 && ~isempty(sub_num(JJ,13)) && ~isnan(sub_num(JJ,13))
    zmin = sub_num(JJ,13);
else
    zmin = 1;
end
if size(sub_num,2) >= 14 && ~isempty(sub_num(JJ,14)) && ~isnan(sub_num(JJ,14))
    zmax = sub_num(JJ,14);
else
    zmax = length(fname);
end

for ii = zmin:zmax
    temp = imread([open_folder,in_folder,sub_list{JJ,3},fname(ii).name]);
    if ii == zmin
        max_image = temp(:,:,channel_value);
    else
        max_image = max(max_image,temp(:,:,channel_value));
    end
end

%%%%% Reload image and exclude the top layers:
im_raw = zeros(0);
fname0 = dir([open_folder,in_folder,sub_list{JJ,3},'*.tif']);

if size(sub_num,2) >= 13 && ~isempty(sub_num(JJ,13)) && ~isnan(sub_num(JJ,13))
    zmin = sub_num(JJ,13);
else
    zmin = 1;
end
if size(sub_num,2) >= 14 && ~isempty(sub_num(JJ,14)) && ~isnan(sub_num(JJ,14))
    zmax = sub_num(JJ,14);
else
    zmax = length(fname0);
end

for ii = zmin:zmax
    raw0 = imread([open_folder,in_folder,sub_list{JJ,3},fname0(ii).name]);
    im_raw = cat(3,im_raw,raw0(:,:,channel_value));
end
%%%%% ****************************************
im00 = max(im_raw,[],3);
% im00 = max_image(:,:,channel_value);
im0 = double(im00)/double(max(im00(:)));
im_now = im0;   %%% Current FISH image with certain contrast

% Load nuclear mask:
try
    load([open_folder,mask_folder,sub_list{JJ,3},mask_name,mat_tail],'mask_stack');
    mask2D = max(mask_stack,[],3);
catch
    mask_stack = zeros(size(max_image,1),size(max_image,2));
    mask2D = mask_stack;
end    
clear mask_stack

% Load foci_list:
if exist([open_folder,hist_folder,sub_list{JJ,3}(1:(end-1)),hist_tail])
    [foci_list,~,~] = xlsread([open_folder,hist_folder,sub_list{JJ,3}(1:(end-1)),hist_tail]);
else
    foci_list = zeros(0,10);
end
if isempty(foci_list)
    foci_list = 10*ones(1,10);
end

foci_list0 = foci_list;
foci_list = foci_combine(foci_list);
foci_xy = max(round(foci_list(:,6:7)),1);
id_list = sub2ind(size(im0),foci_xy(:,1),foci_xy(:,2));
%%% foci image generation:
foci_im = false(size(im0));
foci_im(id_list) = true;
foci_im = bwthicken(foci_im,4);
%%% foci indices matrix generation:
foci_id = double(foci_im);
foci_id(id_list) = 1:numel(id_list);
foci_st = regionprops(foci_im,foci_id,'MaxIntensity');
id_list = [0,[foci_st.MaxIntensity]];
foci_id = id_list(bwlabel(foci_im)+1);

% Recording matrices initialization:
del_im = false(size(im0));
add_im = false(size(im0));
% del_id = false(size(im0));
add_id = false(size(im0));

% image output
axes(handles.im_show);
overlay = im_overlay(im_now,foci_im,del_im,add_im,mask2D);
imshow(overlay)

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
global im0 im_now foci_im del_im add_im mask2D
cmin = str2num(get(handles.contra_min,'String'));
cmax = str2num(get(handles.contra_max,'String'));
if isempty(cmin) || cmin < 0 || cmin > cmax
    cmin = 0;
end
set(handles.contra_min,'String', num2str(cmin));
im_now = (im0-cmin)/(cmax-cmin);

% image output
axes(handles.im_show);
overlay = im_overlay(im_now,foci_im,del_im,add_im,mask2D);
v = axis; 
imshow(overlay)
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
global im0 im_now foci_im del_im add_im mask2D
cmin = str2num(get(handles.contra_min,'String'));
cmax = str2num(get(handles.contra_max,'String'));
if isempty(cmax) || cmax < cmin || cmax > 1
    cmax = 1;
end
set(handles.contra_max,'String', num2str(cmax));
im_now = (im0-cmin)/(cmax-cmin);

% image output
axes(handles.im_show);
overlay = im_overlay(im_now,foci_im,del_im,add_im,mask2D);
v = axis; 
imshow(overlay)
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
global sub_list open_folder out_folder hist_folder mat_tail hist_tail channel_name im0 im_now foci_list foci_im foci_id del_im add_im add_id mask2D
pt = get(gca,'CurrentPoint');
x = floor(pt(1,1)+0.5);
y = floor(pt(1,2)+0.5);

if ~foci_im(y,x) && ~add_im(y,x)
    add_id(y,x) = true;
%     add_im = bwthicken(add_id,4);
    rr = 50;
    imin = max(1,y-rr);
    imax = min(size(im0,1),y+rr);
    jmin = max(1,x-rr);
    jmax = min(size(im0,2),x+rr);
    add_im(imin:imax,jmin:jmax) = bwthicken(add_id(imin:imax,jmin:jmax),1);
elseif ~foci_im(y,x) && add_im(y,x)
    add_label = bwlabel(add_im);
    add_id(add_label == add_label(y,x)) = false;
    add_im(add_label == add_label(y,x)) = false;
end

%%% Image display
overlay = im_overlay(im_now,foci_im,del_im,add_im,mask2D);
axes(handles.im_show);
v = axis; 
imshow(overlay)
axis(v);




function ButttonDownFcn2(src,event,handles)
global sub_list open_folder out_folder hist_folder mat_tail hist_tail channel_name im0 im_now foci_list foci_im foci_id del_im add_im add_id mask2D
pt = get(gca,'CurrentPoint');
x = floor(pt(1,1)+0.5);
y = floor(pt(1,2)+0.5);

if foci_im(y,x) && ~del_im(y,x)
    del_im(foci_id == foci_id(y,x)) = true;
elseif del_im(y,x)
    del_im(foci_id == foci_id(y,x)) = false;
end

%%% Image display
overlay = im_overlay(im_now,foci_im,del_im,add_im,mask2D);
axes(handles.im_show);
v = axis; 
imshow(overlay)
axis(v);





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





% --- Executes on button press in save_seg.
function save_seg_Callback(hObject, eventdata, handles)
% hObject    handle to save_seg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global sub_list open_folder out_folder hist_folder mat_tail hist_tail hist_add channel_name channel_value im0 im_now foci_list foci_im foci_id del_im add_im add_id mask2D
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
    overlay = im_overlay(im_now,foci_im,del_im,add_im,mask2D);
    v = axis; 
    imshow(overlay)
    axis(v);
    
    h = figure;
    overlay = im_overlay(im0,foci_im,del_im,add_im,mask2D);
    imshow(overlay)
    title([open_folder,out_folder,sub_list{JJ,3}(1:(end-1)),', white: transcription foci recognition'])
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
global sub_list open_folder out_folder hist_folder mat_tail hist_tail hist_add channel_name im0 im_now foci_list foci_im foci_id del_im add_im add_id mask2D
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
    overlay = im_overlay(im_now,foci_im,del_im,add_im,mask2D);
    v = axis; 
    imshow(overlay)
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
global zmin zmax sub_list open_folder in_folder out_folder hist_folder mat_tail hist_tail hist_add channel_name channel_value im0 im_now foci_list foci_list0 foci_im foci_id del_im add_im add_id mask2D
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
        new_list = manual_foci([open_folder,in_folder,sub_list{JJ,3}],channel_value,add_id,foci_list,[zmin,zmax]);
        foci_list = [foci_list;new_list];
        foci_list0 = foci_list;
        foci_list = foci_combine(foci_list);
        foci_xy = max(round(foci_list(:,6:7)),1);
        id_list = sub2ind(size(im0),foci_xy(:,1),foci_xy(:,2));
        %%% foci image generation:
        foci_im = false(size(im0));
        foci_im(id_list) = true;
        foci_im = bwthicken(foci_im,4);
        %%% foci indices matrix generation:
        foci_id = double(foci_im);
        foci_id(id_list) = 1:numel(id_list);
        foci_st = regionprops(foci_im,foci_id,'MaxIntensity');
        id_list = [0,[foci_st.MaxIntensity]];
        foci_id = id_list(bwlabel(foci_im)+1);
        %%% Recording matrices reinitialization:
        add_im = false(size(im0));
        add_id = false(size(im0));
    end
    
    
    % image output
    axes(handles.im_show);
    overlay = im_overlay(im_now,foci_im,del_im,add_im,mask2D);
    v = axis; 
    imshow(overlay)
    axis(v);
    
    set(handles.apply_recog,'String','Apply');
end
%%% Update handles structure
guidata(hObject, handles);






function overlay = im_overlay(im_now,foci_im,del_im,add_im,mask2D)
foci_per = bwperim(foci_im);
del_per = bwperim(del_im);
add_per = bwperim(add_im);
mask_per = bwperim(mask2D);

overlay(:,:,1) = im_now+foci_per-del_per;
overlay(:,:,2) = foci_per-del_per+add_per+0.5*mask_per;
overlay(:,:,3) = foci_per+del_per+0.5*mask_per;





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
z_th = 2;
out_list = in_list;

while isrun
    foci_xy = max(round(in_list(:,6:7)),1);
    dxy = squareform(pdist(foci_xy));
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




% --- Executes on button press in auto_desel.
function auto_desel_Callback(hObject, eventdata, handles)
% hObject    handle to auto_desel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global sub_list open_folder out_folder hist_folder mat_tail hist_tail channel_name im0 im_now foci_list foci_im foci_id del_im add_im add_id mask2D
cmin = str2num(get(handles.contra_min,'String'));
cmax = str2num(get(handles.contra_max,'String'));
cth = cmin;

% del_id = imregionalmax(imfilter(im0,fspecial('gaussian',3,1),'symmetric','conv')) & im0 <= cth;
% prop0 = regionprops(foci_im,del_id,'MaxIntensity');
% del_im = ismember(bwlabel(foci_im),find([prop0.MaxIntensity]));
del_id = im0 <= cth;
prop0 = regionprops(foci_im,del_id,'MinIntensity');
del_im = ismember(bwlabel(foci_im),find([prop0.MinIntensity]));

%%% Image display
overlay = im_overlay(im_now,foci_im,del_im,add_im,mask2D);
axes(handles.im_show);
v = axis; 
imshow(overlay)
axis(v);






% --- Executes on button press in auto_sel.
function auto_sel_Callback(hObject, eventdata, handles)
% hObject    handle to auto_sel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global sub_list open_folder out_folder hist_folder mat_tail hist_tail channel_name im0 im_now foci_list foci_im foci_id del_im add_im add_id mask2D
cmin = str2num(get(handles.contra_min,'String'));
cmax = str2num(get(handles.contra_max,'String'));
cth = cmax;

d0 = 10;
size0 = size(im0);
% im_manual = imdilate(bwconvhull(logical(mask2D)),strel('disk',15));
% im_manual(:,[(size0(2)/3-d0):(size0(2)/3+d0),(size0(2)*2/3-d0):(size0(2)*2/3+d0)]) = false;

% % im1 = im0 >= cth;
% % prop0 = regionprops(im1,foci_im,'MaxIntensity');
% % im2 = ismember(bwlabel(im1),find([prop0.MaxIntensity] == 0));
add_id = imregionalmax(imfilter(im0,fspecial('gaussian',3,1),'symmetric','conv')) & im0 >= cth & ~foci_im;% & im_manual;
add_im = bwthicken(add_id,1);

%%% Image display
overlay = im_overlay(im_now,foci_im,del_im,add_im,mask2D);
axes(handles.im_show);
v = axis; 
imshow(overlay)
axis(v);
