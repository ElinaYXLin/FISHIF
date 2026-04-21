function varargout = stack_RNA_check_live(varargin)
% STACK_RNA_CHECK_LIVE M-file for stack_RNA_check_live.fig
%      STACK_RNA_CHECK_LIVE, by itself, creates a new STACK_RNA_CHECK_LIVE or raises the existing
%      singleton*.
%
%      H = STACK_RNA_CHECK_LIVE returns the handle to a new STACK_RNA_CHECK_LIVE or the handle to
%      the existing singleton*.
%
%      STACK_RNA_CHECK_LIVE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in STACK_RNA_CHECK_LIVE.M with the given input arguments.
%
%      STACK_RNA_CHECK_LIVE('Property','Value',...) creates a new STACK_RNA_CHECK_LIVE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before stack_RNA_check_live_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to stack_RNA_check_live_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help stack_RNA_check_live

% Last Modified by GUIDE v2.5 15-Jun-2022 14:55:22

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @stack_RNA_check_live_OpeningFcn, ...
                   'gui_OutputFcn',  @stack_RNA_check_live_OutputFcn, ...
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


% --- Executes just before stack_RNA_check_live is made visible.
function stack_RNA_check_live_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to stack_RNA_check_live (see VARARGIN)

% Parameter setting:
set(gcf,'WindowButtonDownFcn',{}); set(gcf,'KeyPressFcn',{@KeyPressFcn0,hObject,eventdata,handles});
global varargin0 scale0 in_folder input_name out_folder hist_folder hist_folder_single mask_folder mask_name channel_name channel2_name mat_tail hist_tail fit_tail hist_add d0 ccode0 current_handle h00 ev00 resolution0 prompt dlgtitle dims definput opts
h00 = hObject;
ev00 = eventdata;
in_folder = 'stacks/';
input_name = 'matchlist.xls';
out_folder = 'Results/';
channel2_add = '_RNA2';

if length(varargin) >= 2 && ~isempty(varargin{2})
    scale0 = varargin{2};
else
    scale0 = 1;
end

varargin0 = varargin;

if get(handles.ch1,'Value')
%     hist_folder = 'Histogram_alignment/';
    hist_folder = 'Histogram/';
    hist_folder_single = 'Histogram_A/';
    channel_name = 'RNA_channel';
    %channel2_name = 'signal2_channel';
     channel2_name = 'protein_channel';
else
%     hist_folder = ['Histogram_alignment',channel2_add,'/'];
    hist_folder = ['Histogram',channel2_add,'/'];
    hist_folder_single = ['Histogram_A',channel2_add,'/'];
    channel_name = 'signal2_channel';
    channel2_name = 'RNA_channel';
end
current_handle = get(handles.channel_sel,'SelectedObject');

% hist_folder = 'Histogram/';
% hist_folder = 'Histogram_RNA2/';
% hist_folder = 'Histogram_det/';
% hist_folder = 'Histogram_gal4/';
% hist_folder_single = 'Histogram_A_RNA2/';
mask_folder = 'masks/';
mask_name = 'mask';
% channel_name = 'signal2_channel';
% channel_name = 'signal_channel';
% channel_name = 'RNA_channel';
mat_tail = '.mat';
hist_tail = '_raw.xls';
fit_tail = '_spot_fit.mat';
hist_add = '_add';
d0 = 3;
ccode0 = [0,1,0];
resolution0 = 0.083;
%time
prompt = {'Enter a time to check'};
dlgtitle = 'Live Time';
definput = {'0001'};
dims = [1 40];
opts.Interpreter = 'tex';
% answer = inputdlg(prompt,dlgtitle,dims,definput,opts);
% Choose default command line output for stack_RNA_check_live
handles.output = hObject;
%set(hObject,'toolbar','figure');
set(hObject,'Resize','On');

%%% Add button groups:
set(handles.seg_group,'SelectionChangeFcn',@seg_group_SelectionChangeFcn);
set(gcf,'CloseRequestFcn',@th_closereq);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes stack_RNA_check_live wait for user response (see UIRESUME)
% uiwait(handles.figure1);



% --- Outputs from this function are returned to the command line.
function varargout = stack_RNA_check_live_OutputFcn(hObject, eventdata, handles) 
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
global scale0 sub_list sub_num in_folder input_name open_folder prompt dlgtitle dims definput opts Times

J1 = 1; 
TimeAnswer = inputdlg(prompt,dlgtitle,dims,definput,opts);
Times=TimeAnswer{1};
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
global scale0 size0 bin0 sub_list sub_num open_folder out_folder in_folder hist_folder hist_folder_single mask_folder mask_name mat_tail hist_tail fit_tail channel_name channel2_name channel_value channel2_value im0 im1 im_now foci_list foci_list0 foci_im foci_id del_im add_im add_id mask2D b d0 ccode0 resolution0 L_ratio xrange Times TimeAdd TimeSheet TimeAddMask
image_type = '*.tif';
Ex = fspecial('gaussian',10,2);
% Ix = fspecial('gaussian',30,10);
Ix = getnhood(strel('disk',50)); Ix = Ix/sum(Ix(:));
fig = gcf;
JJ = str2num(get(handles.list_J,'String'));
M1 = str2num(get(handles.all_J,'String'));

% Parameter initialization:
set(handles.contra_min,'String','0');
set(handles.contra_max,'String','1');
set(handles.Ref_max,'String','1');
set(handles.nu_on,'Value',false);
set(handles.ref_on,'Value',false);

resolution = sub_num(JJ,9);
L_ratio = (resolution/resolution0);
xrange = max(5,round(3/L_ratio));
TimeAdd=['_time',Times];
% Load max_image:
load([open_folder,out_folder,sub_list{JJ,3}(1:end-1),TimeAdd,mat_tail],'max_image',channel_name,channel2_name);
size0 = size(max_image);
bin0 = round(1/scale0);

max_image = imfilter(max_image,Ex,'replicate'); 

% % % % if exist(channel2_name)
% % % %     imlist = dir([open_folder,in_folder,sub_list{JJ,3},image_type]); %%% get the image list from image folder1
% % % %     
% % % %     temp0 = zeros(0,'uint16');
% % % %     for image_I = 1:length(imlist)
% % % %         raw_im = imread([open_folder,in_folder,sub_list{JJ,3},imlist(image_I).name]);
% % % %         temp0 = cat(4,temp0,raw_im);
% % % %     end
% % % %     max_image = max(temp0,4);
% % % %     clear temp0
% % % % end

% % if exist(channel2_name)
% %     imlist = dir([open_folder,in_folder,sub_list{JJ,3},image_type]); %%% get the image list from image folder1
% %     
% %     temp0 = zeros(0,'uint16');
% %     for image_I = 1:length(imlist)
% %         raw_im = imread([open_folder,in_folder,sub_list{JJ,3},imlist(image_I).name]);
% %         outE = imfilter(raw_im,Ex,'replicate'); 
% %         outI = imfilter(raw_im,Ix,'replicate'); 
% % %        outI = uint16(repmat(reshape(mean(mean(raw_im,1),2),1,1,size(raw_im,3)),size(raw_im,1),size(raw_im,2),1)); 
% %         temp1 = outE-outI;
% %         temp1(temp1 < 0) = 0;
% %         temp0 = cat(4,temp0,temp1);
% %     end
% %     max_image = max(temp0,4);
% %     clear temp0
% % end

if bin0 > 1
    temp0 = zeros(0,'uint16');
    for ii = 1:size0(3)
        temp0 = cat(3,temp0,blkproc(max_image(:,:,ii), [bin0,bin0], 'mean2'));
    end
    max_image = temp0;
end

if exist(channel_name)
    channel_value = eval(channel_name);
    im00 = max_image(:,:,channel_value);
    im0 = double(im00)/double(max(im00(:)));
else
    im0 = zeros(size(max_image,1),size(max_image,2));
end

if exist(channel2_name)
    channel2_value = eval(channel2_name);
    im11 = max_image(:,:,channel2_value);
    im1 = double(im11)/double(max(im11(:)));
else
    im1 = zeros(size(max_image,1),size(max_image,2));
end

im_now = cat(3,im0,im1);   %%% Current FISH image with certain contrast

% Load nuclear mask:
TimeAddMask=['time',Times];
if exist([open_folder,mask_folder,sub_list{JJ,3},TimeAddMask,mask_name,mat_tail])
    load([open_folder,mask_folder,sub_list{JJ,3},TimeAddMask,mask_name,mat_tail],'mask_stack');
    mask2D = max(mask_stack,[],3);
    if bin0 > 1
        mask2D = blkproc(mask2D, [bin0,bin0], 'mean2');
    end
    clear mask_stack
else
    mask2D = zeros(size(im0));
end

% Load foci_list:
TimeSheet=str2double(Times);
fit_tail_time='_spot_fit';
[foci_list,~,~] = xlsread([open_folder,hist_folder,sub_list{JJ,3}(1:(end-1)),hist_tail],TimeSheet);
if exist([open_folder,hist_folder_single,sub_list{JJ,3}(1:(end-1)),fit_tail_time,TimeAdd,'.mat'])
    load([open_folder,hist_folder_single,sub_list{JJ,3}(1:(end-1)),fit_tail_time,TimeAdd,'.mat'],'b');
else
    b = 1;
end
if ~isempty(foci_list)
    foci_list0 = foci_list;
    foci_list = foci_combine(foci_list);
    if bin0 > 1
        foci_list(:,6:7) = foci_list(:,6:7)/bin0;
    end
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
    overlay = im_overlay(im_now,foci_im,del_im,add_im,mask2D,handles);
    imshow(overlay)
else
    % image output
    NanFoci=logical(zeros(size(mask2D)));
    axes(handles.im_show);
    overlay = im_overlay(im_now,NanFoci,NanFoci,NanFoci,mask2D,handles);
    imshow(overlay)
end


if get(handles.nu_on,'Value')
hold on
text(foci_list(:,7)+d0,foci_list(:,6)+d0,num2str(prod(foci_list(:,1:3),2)*2*pi/b),'Color',ccode0)
end

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
global im0 im1 im_now foci_im del_im add_im mask2D foci_list b d0 ccode0
cmin = str2num(get(handles.contra_min,'String'));
cmax = str2num(get(handles.contra_max,'String'));
if isempty(cmin) || cmin < 0 || cmin > cmax
    cmin = 0;
end
set(handles.contra_min,'String', num2str(cmin));
im_now = (cat(3,im0,im1)-cmin)/(cmax-cmin);

% image output
axes(handles.im_show);
overlay = im_overlay(im_now,foci_im,del_im,add_im,mask2D,handles);
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
global im0 im1 im_now foci_im del_im add_im mask2D foci_list b d0 ccode0
cmin = str2num(get(handles.contra_min,'String'));
cmax = str2num(get(handles.contra_max,'String'));
if isempty(cmax) || cmax < cmin || cmax > 1
    cmax = 1;
end
set(handles.contra_max,'String', num2str(cmax));
im_now = (cat(3,im0,im1)-cmin)/(cmax-cmin);

% image output
axes(handles.im_show);
overlay = im_overlay(im_now,foci_im,del_im,add_im,mask2D,handles);
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
global im_now foci_list foci_im foci_id del_im add_id add_im mask2D d0

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
        
    case 'sel_cls_group'
      %execute this code when "Group clear" button is selected
        set(handles.sel_cls_group,'ForegroundColor',[1,0,0]);
        axes(handles.im_show);
        zoom off
        pan off
        v = axis; 
        set(gcf,'WindowButtonDownFcn',{}); set(gcf,'KeyPressFcn',{});

        %%% group de-select
        bw_temp = roipoly;
        add_label = bwlabel(add_im);
        id0 = setdiff(unique(add_label(bw_temp)),0);
        add_id(ismember(add_label,id0)) = false;
        add_im(ismember(add_label,id0)) = false;
        
        %%% Image display
        overlay = im_overlay(im_now,foci_im,del_im,add_im,mask2D,handles);
        imshow(overlay)
        if get(handles.nu_on,'Value')
            hold on
            text(foci_list(:,7)+d0,foci_list(:,6)+d0,num2str(prod(foci_list(:,1:3),2)*2*pi/b),'Color',ccode0)
        end
        axis(v);
        
        set(handles.sel_cls_group,'Value',0);
        set(handles.sel_cls_group,'ForegroundColor',[0,0.5,0]);
        set(handles.sel_cls,'Value',1);
        set(gcf,'WindowButtonDownFcn',{@ButttonDownFcn1,handles}); set(gcf,'KeyPressFcn',{@KeyPressFcn0,hObject,eventdata,handles});

    case 'del_cls_group'
      %execute this code when "Group delete" button is selected
        set(handles.del_cls_group,'ForegroundColor',[1,0,0]);
        axes(handles.im_show);
        zoom off
        pan off
        v = axis; 
        set(gcf,'WindowButtonDownFcn',{}); set(gcf,'KeyPressFcn',{});

        %%% Group de-select
        bw_temp = roipoly;
        id0 = setdiff(unique(foci_id(bw_temp)),0);
        del_im(ismember(foci_id,id0)) = ~del_im(ismember(foci_id,id0));
        
        %%% Image display
        overlay = im_overlay(im_now,foci_im,del_im,add_im,mask2D,handles);
        imshow(overlay)
        if get(handles.nu_on,'Value')
            hold on
            text(foci_list(:,7)+d0,foci_list(:,6)+d0,num2str(prod(foci_list(:,1:3),2)*2*pi/b),'Color',ccode0)
        end
        axis(v);

        set(handles.del_cls_group,'Value',0);
        set(handles.del_cls_group,'ForegroundColor',[0,0,1]);
        set(handles.del_cls,'Value',1);
        set(gcf,'WindowButtonDownFcn',{@ButttonDownFcn2,handles}); set(gcf,'KeyPressFcn',{@KeyPressFcn0,hObject,eventdata,handles});
        
    case 'del_cls_out'
      %execute this code when "Del out" button is selected
        set(handles.del_cls_group,'ForegroundColor',[1,0,0]);
        axes(handles.im_show);
        zoom off
        pan off
        v = axis; 
        set(gcf,'WindowButtonDownFcn',{}); set(gcf,'KeyPressFcn',{});

        %%% Group de-select
        bw_temp = ~mask2D;
        id0 = setdiff(unique(foci_id(bw_temp)),0);
        del_im(ismember(foci_id,id0)) = ~del_im(ismember(foci_id,id0));
        
        %%% Image display
        overlay = im_overlay(im_now,foci_im,del_im,add_im,mask2D,handles);
        imshow(overlay)
        if get(handles.nu_on,'Value')
            hold on
            text(foci_list(:,7)+d0,foci_list(:,6)+d0,num2str(prod(foci_list(:,1:3),2)*2*pi/b),'Color',ccode0)
        end
        axis(v);

        set(handles.del_cls_group,'Value',0);
        set(handles.del_cls_group,'ForegroundColor',[0,0,1]);
        set(handles.del_cls,'Value',1);
        set(gcf,'WindowButtonDownFcn',{@ButttonDownFcn2,handles}); set(gcf,'KeyPressFcn',{@KeyPressFcn0,hObject,eventdata,handles});
        
end
%updates the handles structure
guidata(hObject, handles);




function ButttonDownFcn1(src,event,handles)
global sub_list open_folder out_folder hist_folder mat_tail hist_tail channel_name channel2_name im0 im1 im_now foci_list foci_im foci_id del_im add_im add_id mask2D b d0 ccode0
b_type = get(gcf,'SELECTIONTYPE');
pt = get(gca,'CurrentPoint');
x = floor(pt(1,1)+0.5);
y = floor(pt(1,2)+0.5);
axes(handles.im_show);
v = axis; 

if strcmpi(b_type,'normal')
    if ~foci_im(y,x) && ~add_im(y,x) && x >= v(1) && x <= v(2) && y >= v(3) && y <= v(4)
        rr = 3;
        imin = max(1,y-rr);
        imax = min(size(im0,1),y+rr);
        jmin = max(1,x-rr);
        jmax = min(size(im0,2),x+rr);
        [max0,Imax] = max(~foci_im(imin:imax,jmin:jmax).*~add_im(imin:imax,jmin:jmax).*im0(imin:imax,jmin:jmax));
        [~,j00] = max(max0);
        i00 = Imax(j00);
%         [i00,j00] = ind2sub([imax-imin+1,jmax-jmin+1],Imax);
        add_id(imin+i00-1,jmin+j00-1) = true;
% %         add_id(y,x) = true;
% %     %     add_im = bwthicken(add_id,4);
        rr = 50;
        imin = max(1,y-rr);
        imax = min(size(im0,1),y+rr);
        jmin = max(1,x-rr);
        jmax = min(size(im0,2),x+rr);
        add_im(imin:imax,jmin:jmax) = bwthicken(add_id(imin:imax,jmin:jmax),1);
    elseif ~foci_im(y,x) && add_im(y,x) && x >= v(1) && x <= v(2) && y >= v(3) && y <= v(4)
        add_label = bwlabel(add_im);
        add_id(add_label == add_label(y,x)) = false;
        add_im(add_label == add_label(y,x)) = false;
    end
elseif strcmpi(b_type,'alt') && x >= v(1) && x <= v(2) && y >= v(3) && y <= v(4)
    set(handles.ref_on,'Value',~get(handles.ref_on,'Value'))
    ref_on_Callback(handles.ref_on,event,handles);
end    

%%% Image display
overlay = im_overlay(im_now,foci_im,del_im,add_im,mask2D,handles);
imshow(overlay)

if get(handles.nu_on,'Value')
hold on
text(foci_list(:,7)+d0,foci_list(:,6)+d0,num2str(prod(foci_list(:,1:3),2)*2*pi/b),'Color',ccode0)
end

axis(v);




function ButttonDownFcn2(src,event,handles)
global sub_list open_folder out_folder hist_folder mat_tail hist_tail channel_name channel2_name im0 im1 im_now foci_list foci_im foci_id del_im add_im add_id mask2D b d0 ccode0
b_type = get(gcf,'SELECTIONTYPE');
pt = get(gca,'CurrentPoint');
x = floor(pt(1,1)+0.5);
y = floor(pt(1,2)+0.5);
axes(handles.im_show);
v = axis; 

if strcmpi(b_type,'normal')
    if foci_im(y,x) && ~del_im(y,x) && x >= v(1) && x <= v(2) && y >= v(3) && y <= v(4)
        del_im(foci_id == foci_id(y,x)) = true;
    elseif del_im(y,x) && x >= v(1) && x >= v(1) && x <= v(2) && y >= v(3) && y <= v(4)
        del_im(foci_id == foci_id(y,x)) = false;
    end
elseif strcmpi(b_type,'alt') && x >= v(1) && x <= v(2) && y >= v(3) && y <= v(4)
    set(handles.ref_on,'Value',~get(handles.ref_on,'Value'))
    ref_on_Callback(handles.ref_on,event,handles);
end    

%%% Image display
overlay = im_overlay(im_now,foci_im,del_im,add_im,mask2D,handles);
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
global sub_list sub_num open_folder out_folder hist_folder mat_tail hist_tail channel_name channel2_name im0 im1 im_now foci_list foci_im foci_id del_im add_im add_id mask2D b d0 ccode0
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
% % im_manual = imdilate(bwconvhull(logical(mask2D)),strel('disk',15));
im_manual = imdilate(logical(mask2D),strel('disk',5));
% % % % % for tt = 1:length(Nbin1)
% % % % %     for ss = 1:(Nbin1(1)-1)
% % % % %         im_manual(:,(ss*size0(Mdim(tt))/Nbin1(tt)-d0):(ss*size0(Mdim(tt))/Nbin1(tt)+d0)) = false;
% % % % %     end
% % % % % end

% % im1 = im0 >= cth;
% % prop0 = regionprops(im1,foci_im,'MaxIntensity');
% % im2 = ismember(bwlabel(im1),find([prop0.MaxIntensity] == 0));
im0f = imfilter(im0,fspecial('gaussian',15,1.5)-fspecial('gaussian',15,4),'symmetric','conv');
% % add_id = imregionalmax(im0f) & im0f >= cth & ~foci_im;% & im_manual;
add_id = imregionalmax(im0f) & im0f >= cth/8 & ~imdilate(foci_im,strel('disk',5)) & im_manual;
% add_id = imregionalmax(imfilter(im0,fspecial('gaussian',3,1),'symmetric','conv')) & im0 >= cth & ~foci_im;% & im_manual;
add_im = bwthicken(add_id,1);

%%% Image display
overlay = im_overlay(im_now,foci_im,del_im,add_im,mask2D,handles);
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
global sub_list open_folder out_folder hist_folder mat_tail hist_tail channel_name channel2_name im0 im1 im_now foci_list foci_im foci_id del_im add_im add_id mask2D b d0 ccode0
cmin = str2num(get(handles.contra_min,'String'));
cmax = str2num(get(handles.contra_max,'String'));
cth = cmin;

% del_id = imregionalmax(imfilter(im0,fspecial('gaussian',3,1),'symmetric','conv')) & im0 <= cth;
% prop0 = regionprops(foci_im,del_id,'MaxIntensity');
% del_im = ismember(bwlabel(foci_im),find([prop0.MaxIntensity]));
prop0 = regionprops(foci_im,imfilter(im0,fspecial('gaussian',3,1),'symmetric','conv'),'MaxIntensity');
del_im = ismember(bwlabel(foci_im),find([prop0.MaxIntensity] <= cth));

%%% Image display
overlay = im_overlay(im_now,foci_im,del_im,add_im,mask2D,handles);
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
global bin0 sub_list open_folder out_folder hist_folder mat_tail hist_tail hist_add channel_name channel2_name channel_value channel2_value im0 im1 im_now foci_list foci_im foci_id del_im add_im add_id mask2D b d0 ccode0 Times TimeAdd TimeSheet TimeAddMask
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
    
    foci_list1 = foci_list;
    if bin0 > 1
        foci_list1(:,6:7) = foci_list(:,6:7)*bin0;
    end
    
%     delete([open_folder,hist_folder,sub_list{JJ,3}(1:(end-1)),hist_tail]);
    if isempty(foci_list1)
        foci_list_save=zeros(1,10);
    else
        foci_list_save=foci_list1;
    end
    xlswrite([open_folder,hist_folder,sub_list{JJ,3}(1:(end-1)),hist_tail],' ',TimeSheet,'A1:J500');
    xlswrite([open_folder,hist_folder,sub_list{JJ,3}(1:(end-1)),hist_tail],foci_list_save,TimeSheet);
    if nnz(add_id)
        save([open_folder,hist_folder,sub_list{JJ,3}(1:(end-1)),hist_add,TimeAdd,mat_tail],'add_id','channel_value')
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
    overlay = im_overlay(im_now,foci_im,del_im,add_im,mask2D,handles);
    v = axis; 
    imshow(overlay)

    if get(handles.nu_on,'Value')
    hold on
    text(foci_list(:,7)+d0,foci_list(:,6)+d0,num2str(prod(foci_list(:,1:3),2)*2*pi/b),'Color',ccode0)
    end

    axis(v);
    
    h = figure;
    overlay = im_overlay(im0,foci_im,del_im,add_im,mask2D,handles);
    imshow(overlay)

    if get(handles.nu_on,'Value')
    hold on
    text(foci_list(:,7)+d0,foci_list(:,6)+d0,num2str(prod(foci_list(:,1:3),2)*2*pi/b),'Color',ccode0)
    end

    title([open_folder,out_folder,sub_list{JJ,3}(1:(end-1)),', white: transcription foci recognition'])
    outim_tail_time=['_foci_seg',TimeAdd,'.fig'];
    saveas(h,[open_folder,hist_folder,sub_list{JJ,3}(1:(end-1)),outim_tail_time])
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
global sub_list open_folder out_folder hist_folder mat_tail hist_tail hist_add channel_name channel2_name im0 im1 im_now foci_list foci_im foci_id del_im add_im add_id mask2D b d0 ccode0
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
    overlay = im_overlay(im_now,foci_im,del_im,add_im,mask2D,handles);
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
global size0 bin0 sub_list open_folder in_folder out_folder hist_folder mat_tail hist_tail hist_add channel_name channel2_name channel_value channel2_value im0 im1 im_now foci_list foci_list0 foci_im foci_id del_im add_im add_id mask2D b d0 ccode0 L_ratio xrange
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
        
        new_list = manual_foci([open_folder,in_folder,sub_list{JJ,3}],channel_value,add_id1,foci_list1,[],xrange);
        
        if bin0 > 1
            new_list(:,6:7) = new_list(:,6:7)/bin0;
        end
        
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
    overlay = im_overlay(im_now,foci_im,del_im,add_im,mask2D,handles);
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






function overlay = im_overlay(im_now,foci_im,del_im,add_im,mask2D,handles)
foci_per = bwperim(foci_im);
del_per = bwperim(del_im);
add_per = bwperim(add_im);
mask_per = bwperim(mask2D);

overlay(:,:,1) = im_now(:,:,1)+foci_per-del_per;
% if size(im_now,3) > 1
%     overlay(:,:,2) = 0.7*get(handles.ref_on,'Value')*im_now(:,:,2)+foci_per-del_per+add_per+0.5*mask_per;
% else
    overlay(:,:,2) = foci_per-del_per+add_per+0.5*mask_per;
% end
if size(im_now,3) > 1
    overlay(:,:,3) = str2num(get(handles.Ref_max,'String'))*get(handles.ref_on,'Value')*im_now(:,:,2)+foci_per+del_per+0.5*mask_per;
else
    overlay(:,:,3) = foci_per+del_per+0.5*mask_per;
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
Lth = 1;%old:1.5
isrun = true;
z_th = 0;
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
    stack_RNA_check_live_OpeningFcn(h00, ev00, handles, varargin0{:})
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
global im0 im1 im_now foci_im del_im add_im mask2D foci_list b d0 ccode0 text0

% image output
axes(handles.im_show);
overlay = im_overlay(im_now,foci_im,del_im,add_im,mask2D,handles);
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

global im0 im1 im_now foci_im del_im add_im mask2D foci_list b d0 ccode0

% image output
axes(handles.im_show);
overlay = im_overlay(im_now,foci_im,del_im,add_im,mask2D,handles);
v = axis; 
imshow(overlay)

if get(handles.nu_on,'Value')
hold on
text(foci_list(:,7)+d0,foci_list(:,6)+d0,num2str(prod(foci_list(:,1:3),2)*2*pi/b),'Color',ccode0)
end

axis(v);



function Ref_max_Callback(hObject, eventdata, handles)
% hObject    handle to Ref_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Ref_max as text
%        str2double(get(hObject,'String')) returns contents of Ref_max as a double

global im0 im1 im_now foci_im del_im add_im mask2D foci_list b d0 ccode0
rmax = str2num(get(handles.Ref_max,'String'));
if isempty(rmax) || rmax <= 0
    rmax = 1;
end
set(handles.Ref_max,'String', num2str(rmax));

% image output
axes(handles.im_show);
overlay = im_overlay(im_now,foci_im,del_im,add_im,mask2D,handles);
v = axis; 
imshow(overlay)

if get(handles.nu_on,'Value')
hold on
text(foci_list(:,7)+d0,foci_list(:,6)+d0,num2str(prod(foci_list(:,1:3),2)*2*pi/b),'Color',ccode0)
end

axis(v);



% --- Executes during object creation, after setting all properties.
function Ref_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Ref_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on key press with focus on list_J and none of its controls.
function list_J_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to list_J (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)



function TimeCheck_Callback(hObject, eventdata, handles)
% hObject    handle to TimeCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(gcf,'WindowButtonDownFcn',{}); set(gcf,'KeyPressFcn',{@KeyPressFcn0,hObject,eventdata,handles});
global sub_list Times

fig = gcf;
Times = get(handles.TimeCheck,'String');
JJ = str2num(get(handles.list_J,'String'));
M1 = str2num(get(handles.all_J,'String'));
TimesQ=[' Time ',Times];
if JJ <= M1
    quit_reply = questdlg(['Move to file #',sub_list{JJ,3}(1:(length(sub_list{JJ,3})-1)),TimesQ,' (Make sure you have saved the result)?']);
    if strcmp(quit_reply,'Yes')
        set(handles.next_file,'String','Wait');
        % Update handles structure
        guidata(hObject, handles);
        list_J_core(hObject, eventdata, handles);
        set(handles.next_file,'String','Next');
    end
end
% Hints: get(hObject,'String') returns contents of TimeCheck as text
%        str2double(get(hObject,'String')) returns contents of TimeCheck as a double


% --- Executes during object creation, after setting all properties.
function TimeCheck_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TimeCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
