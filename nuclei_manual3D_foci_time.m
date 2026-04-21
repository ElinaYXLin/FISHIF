function varargout = nuclei_manual3D_foci_time(T,varargin)
% NUCLEI_MANUAL3D_FOCI M-file for nuclei_manual3D_foci.fig
%      NUCLEI_MANUAL3D_FOCI, by itself, creates a new NUCLEI_MANUAL3D_FOCI or raises the existing
%      singleton*.
%
%      H = NUCLEI_MANUAL3D_FOCI returns the handle to a new NUCLEI_MANUAL3D_FOCI or the handle to
%      the existing singleton*.
%
%      NUCLEI_MANUAL3D_FOCI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in NUCLEI_MANUAL3D_FOCI.M with the given input arguments.
%
%      NUCLEI_MANUAL3D_FOCI('Property','Value',...) creates a new NUCLEI_MANUAL3D_FOCI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before nuclei_manual3D_foci_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to nuclei_manual3D_foci_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help nuclei_manual3D_foci

% Last Modified by GUIDE v2.5 18-Aug-2021 17:00:07

% Begin initialization code - DO NOT EDIT
global time_index
time_index=T;
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @nuclei_manual3D_foci_OpeningFcn, ...
                   'gui_OutputFcn',  @nuclei_manual3D_foci_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin>1 && ischar(varargin{1}) && ~isempty(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before nuclei_manual3D_foci is made visible.
function nuclei_manual3D_foci_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to nuclei_manual3D_foci (see VARARGIN)

% Parameter setting:
set(gcf,'WindowButtonDownFcn',{}); set(gcf,'KeyPressFcn',{@KeyPressFcn0,hObject,eventdata,handles});
global time_index varargin0 scale0 bin0 size0 sub_list sub_num in_folder input_name out_folder mask_name quick_mask_name hist_folder hist_tail hist_tail2 mat_tail resolution0 z_range image_folder DAPI_channel DAPI_th0 new_DAPI new_DAPI3D bw_DAPI bw_applied bw_applied3D sel_DAPI sel_applied L_ratio DAPI_on0 applied_on0 foci_on0 new_DAPI3D_save bw_applied3D_save v_save mask_stack layer_I layer_num  bw_cut area_thresh

if length(varargin) >= 2 && ~isempty(varargin{2})
    scale0 = varargin{2};
else
    scale0 = 1;
end

bin0 = round(1/scale0);
varargin0 = varargin;
in_folder = 'stacks/';
input_name = 'matchlist.xls';
out_folder = 'masks/';
mask_name = ['time',num2str(time_index,'%04u'),'mask'];
quick_mask_name =['time',num2str(time_index,'%04u'),'quick_mask'];
mat_tail = '.mat';
hist_folder = 'Histogram_alignment/';
hist_tail = '_raw.xls';
hist_tail2 = '_raw.xlsx';
resolution0 = 0.091;
DAPI_th0 = zeros(0);
area_thresh = 400;
z_range = 3;
% File loading:
open_folder = uigetdir;
open_folder = [open_folder,'\'];
[sub_num, sub_list] = xlsread([open_folder,in_folder,input_name]);
[M1,~] = size(sub_list);

if length(varargin) >= 3 && ~isempty(varargin{3}) && M1 >= varargin{3} && varargin{3} > 0
    JJ = varargin{3};
else
    JJ = 1;
end

set(handles.list_J,'String',num2str(JJ));
set(handles.all_J,'String',num2str(M1));
set(handles.folder_name,'String',open_folder);
set(handles.file_name,'String',sub_list{JJ,3}(1:(length(sub_list{JJ,3})-1)));
% Choose default command line output for nuclei_manual3D
handles.output = hObject;
set(hObject,'Resize','On');

list_J_Callback(hObject, eventdata, handles);

%%% Add button groups:
set(handles.seg_group,'SelectionChangeFcn',@seg_group_SelectionChangeFcn);
set(gcf,'CloseRequestFcn',@th_closereq);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes nuclei_manual3D_foci wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = nuclei_manual3D_foci_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function list_J_Callback(hObject, eventdata, handles)
% hObject    handle to list_J (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of list_J as text
%        str2double(get(hObject,'String')) returns contents of list_J as a double
set(gcf,'WindowButtonDownFcn',{}); set(gcf,'KeyPressFcn',{@KeyPressFcn0,hObject,eventdata,handles});
global time_index varargin0  bin0 size0 sub_list sub_num in_folder input_name out_folder mask_name quick_mask_name hist_folder hist_tail hist_tail2 mat_tail resolution0 z_range image_folder DAPI_channel DAPI_th0 new_DAPI new_DAPI3D bw_DAPI bw_applied bw_applied3D sel_DAPI sel_applied L_ratio DAPI_on0 applied_on0 foci_on0 new_DAPI3D_save bw_applied3D_save v_save mask_stack im_foci im_foci3D im_foci3D_save layer_I layer_num  bw_cut area_thresh

fig = gcf;
JJ = str2num(get(handles.list_J,'String'));
M1 = str2num(get(handles.all_J,'String'));
layer_I = 1;

if JJ == 1
    quit_reply = 'Yes';
elseif JJ <= M1
    quit_reply = questdlg(['Move to file #',sub_list{JJ,3}(1:(length(sub_list{JJ,3})-1)),' (Make sure you have saved the result)?']);
end
    
if strcmp(quit_reply,'Yes')
    set(handles.next_file,'String','Wait');
    % Update handles structure
    guidata(hObject, handles);
    % File loading:

    if JJ <= M1
        open_folder = get(handles.folder_name,'String');
        set(handles.list_J,'String',num2str(JJ));
        set(handles.file_name,'String',sub_list{JJ,3}(1:(length(sub_list{JJ,3})-1)));
        DAPI_th0 = zeros(0);
        
        %%% Load 3D mask: %%% =====================================================
        image_folder = [open_folder,in_folder,sub_list{JJ,3}];
        if ~exist([open_folder,out_folder,sub_list{JJ,3},mask_name,mat_tail],'file')
            if ~exist([open_folder,out_folder,sub_list{JJ,3}],'dir')
                mkdir([open_folder,out_folder,sub_list{JJ,3}]);
            end
            DAPI_seg3D2(image_folder);
        end
        
        if ~exist([open_folder,out_folder,sub_list{JJ,3},quick_mask_name,mat_tail])
            load([open_folder,out_folder,sub_list{JJ,3},mask_name,mat_tail]);
            mask_stack = logical(mask_stack);
        else
            load([open_folder,out_folder,sub_list{JJ,3},quick_mask_name,mat_tail]);
            mask_stack = logical(bw_applied3D_save);
        end
        
        size0 = size(mask_stack);
        if bin0 > 1
%             temp0 = false(0);
%             for ii = 1:size0(3)
%                 temp0 = cat(3,temp0,blkproc(mask_stack(:,:,ii), [bin0,bin0], 'mean2') >= 0.5);
%             end
%             mask_stack = temp0;
            mask_stack = imresize(mask_stack,1/bin0,'nearest');
        end
        %%% =======================================================================

        
        %%% Load foci: %%% ========================================================
        if exist([open_folder,hist_folder,sub_list{JJ,3}(1:(end-1)),hist_tail],'file')
            [foci_list,~,~] = xlsread([open_folder,hist_folder,sub_list{JJ,3}(1:(end-1)),hist_tail]);
        elseif exist([open_folder,hist_folder,sub_list{JJ,3}(1:(end-1)),hist_tail2],'file')
            [foci_list,~,~] = xlsread([open_folder,hist_folder,sub_list{JJ,3}(1:(end-1)),hist_tail2]);
        else
            foci_list = [];
        end
        im_foci3D = zeros(size(mask_stack));
        if ~isempty(foci_list)
            ind_foci = sub2ind(size(mask_stack),min(max(round(foci_list(:,6)/bin0),1),size(mask_stack,1)),min(max(round(foci_list(:,7)/bin0),1),size(mask_stack,2)),max(min(round(foci_list(:,8)),size(mask_stack,3)),1));
        else
            ind_foci = [];
        end
        im_foci3D(ind_foci) = 1;
        temp0 = double(getnhood(strel('disk',4)));
        kernel0 = repmat(temp0,[1,1,2*z_range+1]);
        im_foci3D = min(imfilter(im_foci3D,kernel0,'symmetric','conv'),1);
        im_foci = im_foci3D(:,:,layer_I);
        %%% =======================================================================

        
        %%% Load image stack: %%% =================================================
        image_type = ['time',num2str(time_index,'%04u'),'*.tif'];
        DAPI_channel = sub_num(JJ,7);
        imlist = dir([image_folder,image_type]); %%% get the image list from image folder1

%         sigmaE = 40;
%         sigmaI = 60;
%         FiltSize = round(1.3*sigmaI);
%         H = fspecial('disk',30); 
%         Ex = fspecial('gaussian',FiltSize,sigmaE);
%         Ix = fspecial('gaussian',FiltSize,sigmaI);
        
        for image_I = 1:length(imlist)
            raw_im = imread([image_folder,imlist(image_I).name]);
            if bin0 > 1
%                 temp1 = blkproc(raw_im(:,:,DAPI_channel), [bin0,bin0], 'mean2');
                temp1 = imresize(raw_im(:,:,DAPI_channel),1/bin0,'nearest');
            else
                temp1 = raw_im(:,:,DAPI_channel);
            end
            if image_I == 1
                new_DAPI3D = zeros(size(temp1,1),size(temp1,2),length(imlist));
            end
            new_DAPI3D(:,:,image_I) = imfilter(double(temp1),fspecial('gaussian',10,0.75),'same','conv');
            
% %              new_DAPI3D(:,:,image_I) = imfilter((double(raw_im(:,:,2))+double(raw_im(:,:,5)))/2,fspecial('gaussian',10,1.5),'same','conv');
% % %              new_DAPI3D(:,:,image_I) = imfilter(double(raw_im(:,:,DAPI_channel)),fspecial('gaussian',10,0.75),'same','conv');
%              new_DAPI3D(:,:,image_I) = imfilter(double(raw_im(:,:,DAPI_channel)),fspecial('gaussian',10,1.5),'same','conv');
%              new_DAPI3D(:,:,image_I) = imfilter(double(raw_im(:,:,DAPI_channel)),fspecial('gaussian',10,3),'same','conv');
%             new_DAPI3D(:,:,image_I) = imfilter(double(raw_im(:,:,DAPI_channel)),fspecial('disk',10,5),'same','conv');
            
%             I = imfilter(double(raw_im(:,:,DAPI_channel)),H,'same','conv');
% %             I =  adapthisteq(I);
%             outE = imfilter(single(I),Ex,'replicate');
%             outI = imfilter(single(I),Ix,'replicate');
%             new_DAPI3D(:,:,image_I) = outE - outI;

            
        end
%         new_DAPI3D(new_DAPI3D < 0) = 0;
        new_DAPI3D = new_DAPI3D/max(new_DAPI3D(:));

        layer_num = length(imlist);
%         [~,layer_I] = max(mean(mean(new_DAPI)));
        new_DAPI = new_DAPI3D(:,:,layer_I);
        %%% =======================================================================

        resolution = sub_num(JJ,9);
        L_ratio = (resolution/resolution0);
        set(handles.list_J,'UserData', L_ratio);

        % Choose default command line output for nuclei_manual3D_foci
        handles.output = hObject;
        %set(hObject,'toolbar','figure');
%        set(hObject,'Resize','On');

        % Initialization of GUI:
        if isempty(DAPI_th0)
            DAPI_th0 = 0.1;
        end
        set(handles.DAPI_th,'String', num2str(DAPI_th0));
        set(handles.applied_on,'Value',true);
        set(handles.DAPI_on,'Value',true);
        set(handles.current_seg,'Value',true);
        set(handles.im_lock,'Value',false);
        set(handles.im_lock,'String','Lock');
        set(handles.slider_layer,'Max',layer_num);
        set(handles.slider_layer,'Min',1);
        set(handles.slider_layer,'SliderStep',[1./(layer_num-1+(layer_num == 1)),1./(layer_num-1+(layer_num == 1))]);
        set(handles.slider_layer,'Value',layer_I);
        set(handles.I_layer,'String', num2str(layer_I));
        set(handles.N_layer,'String', num2str(layer_num));
        set(handles.sel_empty,'value',true);

        DAPI_on0 = 1;
        applied_on0 = 1;
        foci_on0 = 0;

        % Show image:
        bw_DAPI = im2bw(new_DAPI,DAPI_th0);
        bw_DAPI = imfill(bw_DAPI,'holes');
        bw_DAPI = bwareaopen(bw_DAPI, area_thresh);

        bw_applied = logical(mask_stack(:,:,layer_I));
        bw_applied3D = logical(mask_stack);
        sel_DAPI = false(size(bw_applied));
        sel_applied = sel_DAPI;

        im_foci3D_save = im_foci3D;
        new_DAPI3D_save = new_DAPI3D;
        bw_applied3D_save = bw_applied3D;
        bw_cut = false(size(bw_DAPI));

        %%% Image display
        overlay = im_overlay(new_DAPI,sel_DAPI,sel_applied,bw_cut,bw_DAPI,bw_applied,im_foci,handles);
        axes(handles.im_show);
        imshow(overlay)
        v_save = [1,size(new_DAPI,2),1,size(new_DAPI,1)];

    else
        close(gcf)
    end
    set(handles.next_file,'String','Next');
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



% --- Executes on button press in next_file.
function next_file_Callback(hObject, eventdata, handles)
% hObject    handle to next_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(gcf,'WindowButtonDownFcn',{}); set(gcf,'KeyPressFcn',{@KeyPressFcn0,hObject,eventdata,handles});
fig = gcf;
quit_reply = questdlg('Move to the next file (Make sure you have saved the result)?');
if strcmp(quit_reply,'Yes')
    set(handles.next_file,'String','Wait');
    % File loading:
    JJ = str2num(get(handles.list_J,'String'))+1;
    set(handles.list_J,'String',num2str(JJ));
    % Update handles structure
    guidata(hObject, handles);
    
    list_J_Callback(hObject, eventdata, handles);
    set(handles.next_file,'String','Next');
end
% Update handles structure
guidata(hObject, handles);



% --- Executes on button press in save_seg.
function save_seg_Callback(hObject, eventdata, handles)
% hObject    handle to save_seg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global mask_stack0 bin0 size0 sub_list sub_num in_folder input_name out_folder mask_name mat_tail resolution0 image_folder DAPI_channel DAPI_th0 new_DAPI new_DAPI3D bw_DAPI bw_applied bw_applied3D sel_DAPI sel_applied L_ratio DAPI_on0 applied_on0 foci_on0 new_DAPI3D_save bw_applied3D_save v_save mask_stack im_foci im_foci3D im_foci3D_save layer_I layer_num  bw_cut
quit_reply = questdlg('Save the segmentation result?');
if strcmp(quit_reply,'Yes')
    set(handles.save_seg,'String','Wait');
    % Update handles structure
    guidata(hObject, handles);
    
    %%% Data loading
    label_layer_Callback(hObject, eventdata, handles)    
    mask_stack0 = mask_stack;
    if bin0 > 1
        mask_stack = imresize(mask_stack,bin0,'nearest');
        mask_stack = mask_stack(1:size0(1),1:size0(2),:);
    end
    
    %%% Data storage
    open_folder = get(handles.folder_name,'String');
    JJ = str2num(get(handles.list_J,'String'));
    if exist([open_folder,out_folder,sub_list{JJ,3},mask_name,mat_tail])
        movefile([open_folder,out_folder,sub_list{JJ,3},mask_name,mat_tail],[open_folder,out_folder,sub_list{JJ,3},mask_name,'_old',mat_tail])
    end
    save([open_folder,out_folder,sub_list{JJ,3},mask_name,mat_tail],'mask_stack','-v7.3');
%     save([open_folder,out_folder,sub_list{JJ,3},mask_name,mat_tail],'mask_stack','-append','-v7.3');
    mask_stack = mask_stack0;
    
    set(handles.save_seg,'String','Save');
    guidata(hObject, handles);
end
% Update handles structure
guidata(hObject, handles);





% --- Executes on button press in quick_save_seg.
function quick_save_seg_Callback(hObject, eventdata, handles)
% hObject    handle to quick_save_seg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global temp00 bin0 size0 sub_list sub_num in_folder input_name out_folder mask_name quick_mask_name mat_tail resolution0 image_folder  DAPI_th0 new_DAPI new_DAPI3D bw_DAPI bw_applied bw_applied3D sel_DAPI sel_applied L_ratio DAPI_on0 applied_on0 foci_on0 new_DAPI3D_save bw_applied3D_save v_save mask_stack im_foci im_foci3D im_foci3D_save layer_I layer_num  bw_cut
quit_reply = questdlg('Quick save the segmentation result?');
if strcmp(quit_reply,'Yes')
    set(handles.quick_save_seg,'String','Wait');
    % Update handles structure
    guidata(hObject, handles);
    
% % %     %%% Data loading
% % %     label_layer_Callback(hObject, eventdata, handles)    
    
    %%% Data storage
    open_folder = get(handles.folder_name,'String');
    JJ = str2num(get(handles.list_J,'String'));
    
    bw_applied3D(:,:,layer_I) = bw_applied;
    bw_applied3D_save(v_save(3):v_save(4),v_save(1):v_save(2),:) = bw_applied3D;

    temp00 = bw_applied3D_save;
    if bin0 > 1
        bw_applied3D_save = imresize(bw_applied3D_save,bin0,'nearest');
        bw_applied3D_save = bw_applied3D_save(1:size0(1),1:size0(2),:);
    end
    
    if exist([open_folder,out_folder,sub_list{JJ,3},quick_mask_name,'_old',mat_tail])
        movefile([open_folder,out_folder,sub_list{JJ,3},quick_mask_name,'_old',mat_tail],[open_folder,out_folder,sub_list{JJ,3},quick_mask_name,'_old2',mat_tail])
    end
    if exist([open_folder,out_folder,sub_list{JJ,3},quick_mask_name,mat_tail])
        movefile([open_folder,out_folder,sub_list{JJ,3},quick_mask_name,mat_tail],[open_folder,out_folder,sub_list{JJ,3},quick_mask_name,'_old',mat_tail])
    end
    save([open_folder,out_folder,sub_list{JJ,3},quick_mask_name,mat_tail],'bw_applied3D_save','-v7.3');
    bw_applied3D_save = temp00;
    
    set(handles.quick_save_seg,'String','Quick save');
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
quit_reply = questdlg('All parameters return to default values?');
if strcmp(quit_reply,'Yes')
    bw_cut = false(size(bw_DAPI));
    set(handles.cancel_seg,'String','Wait');
    % Update handles structure
    guidata(hObject, handles);
    list_J_Callback(hObject, eventdata, handles);
    set(handles.cancel_seg,'String','Cancel');
end
%%% Update handles structure
guidata(hObject, handles);




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
global sub_list sub_num in_folder input_name out_folder mask_name mat_tail resolution0 image_folder DAPI_channel DAPI_th0 new_DAPI new_DAPI3D bw_DAPI bw_applied bw_applied3D sel_DAPI sel_applied L_ratio DAPI_on0 applied_on0 foci_on0 new_DAPI3D_save bw_applied3D_save v_save mask_stack im_foci im_foci3D im_foci3D_save layer_I layer_num  bw_cut area_thresh
bw_applied3D(:,:,layer_I) = bw_applied;

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

bw_cut = false(size(bw_DAPI));
im_foci = im_foci3D(:,:,layer_I);
new_DAPI = new_DAPI3D(:,:,layer_I);
bw_DAPI = im2bw(new_DAPI,DAPI_th0);
bw_DAPI = imfill(bw_DAPI,'holes');
bw_DAPI = bwareaopen(bw_DAPI, area_thresh);

bw_applied = bw_applied3D(:,:,layer_I);
sel_DAPI = false(size(bw_applied));
sel_applied = sel_DAPI;

%%% Image display
overlay = im_overlay(new_DAPI,sel_DAPI,sel_applied,bw_cut,bw_DAPI,bw_applied,im_foci,handles);
axes(handles.im_show);
v = axis; 
imshow(overlay)
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




function DAPI_th_Callback(hObject, eventdata, handles)
% hObject    handle to DAPI_th (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DAPI_th as text
%        str2double(get(hObject,'String')) returns contents of DAPI_th as a double
global sub_list sub_num in_folder input_name out_folder mask_name mat_tail resolution0 image_folder DAPI_channel DAPI_th0 new_DAPI new_DAPI3D bw_DAPI bw_applied bw_applied3D sel_DAPI sel_applied L_ratio DAPI_on0 applied_on0 foci_on0 new_DAPI3D_save bw_applied3D_save v_save mask_stack im_foci im_foci3D im_foci3D_save layer_I layer_num  bw_cut area_thresh

minT = get(handles.DAPI_th,'Min');
maxT = get(handles.DAPI_th,'Max');

%get the string for the editText component
DAPI_th0 = str2num(get(handles.DAPI_th,'String'));
 
%if user inputs something is not a number, or if the input is less than 0
%or greater than 100, then the slider value defaults to 0
if (isempty(DAPI_th0) || DAPI_th0 < minT)
    DAPI_th0 = minT;
    set(handles.DAPI_th,'String',num2str(DAPI_th0));
elseif DAPI_th0 > maxT
    DAPI_th0 = maxT;
    set(handles.DAPI_th,'String',num2str(DAPI_th0));
end

bw_DAPI = im2bw(new_DAPI,DAPI_th0);
bw_DAPI = imfill(bw_DAPI,'holes');
bw_DAPI = bwareaopen(bw_DAPI, area_thresh);

%%% Image display
overlay = im_overlay(new_DAPI,sel_DAPI,sel_applied,bw_cut,bw_DAPI,bw_applied,im_foci,handles);
axes(handles.im_show);
v = axis; 
imshow(overlay)
axis(v);




% --- Executes during object creation, after setting all properties.
function DAPI_th_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DAPI_th (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in DAPI_on.
function DAPI_on_Callback(hObject, eventdata, handles)
% hObject    handle to DAPI_on (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of DAPI_on
global sub_list sub_num in_folder input_name out_folder mask_name mat_tail resolution0 image_folder DAPI_channel DAPI_th0 new_DAPI new_DAPI3D bw_DAPI bw_applied bw_applied3D sel_DAPI sel_applied L_ratio DAPI_on0 applied_on0 foci_on0 new_DAPI3D_save bw_applied3D_save v_save mask_stack im_foci im_foci3D im_foci3D_save layer_I layer_num  bw_cut
DAPI_on0 = get(handles.DAPI_on,'Value');

%%% Image display
overlay = im_overlay(new_DAPI,sel_DAPI,sel_applied,bw_cut,bw_DAPI,bw_applied,im_foci,handles);
axes(handles.im_show);
v = axis; 
imshow(overlay)
axis(v);




% --- Executes on button press in seg_convex.
function seg_convex_Callback(hObject, eventdata, handles)
% hObject    handle to seg_convex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of seg_convex
global sub_list sub_num in_folder input_name out_folder mask_name mat_tail resolution0 image_folder DAPI_channel DAPI_th0 new_DAPI new_DAPI3D bw_DAPI bw_applied bw_applied3D sel_DAPI sel_applied L_ratio DAPI_on0 applied_on0 foci_on0 new_DAPI3D_save bw_applied3D_save v_save mask_stack im_foci im_foci3D im_foci3D_save layer_I layer_num  bw_cut
if get(handles.current_seg,'Value')
    sel_DAPI = bwconvhull(sel_DAPI,'objects');
    bw_DAPI = bw_DAPI | sel_DAPI;
else
    bw_cut0 = sel_applied | bw_applied;
    
    sel_applied = bwconvhull(sel_applied,'objects');
    bw_applied = bw_applied | sel_applied;
    
%     temp_sel = ismember(temp_label,temp_label(bw_cut0));
    D= bwdist(~bw_applied);
    g2 = imimposemin(-D,bw_cut0);
    bw_cut0 = imdilate(watershed(g2) == 0,strel('disk',1));% & (temp_sel);
    bw_applied = bw_applied & (~ bw_cut0);
    sel_applied = sel_applied & (~ bw_cut0);
    
end
%%% Image display
overlay = im_overlay(new_DAPI,sel_DAPI,sel_applied,bw_cut,bw_DAPI,bw_applied,im_foci,handles);
axes(handles.im_show);
v = axis; 
imshow(overlay)
axis(v);




% --- Executes on button press in dilation_r.
function dilation_r_Callback(hObject, eventdata, handles)
% hObject    handle to dilation_r (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of dilation_r
global sub_list sub_num in_folder input_name out_folder mask_name mat_tail resolution0 image_folder DAPI_channel DAPI_th0 new_DAPI new_DAPI3D bw_DAPI bw_applied bw_applied3D sel_DAPI sel_applied L_ratio DAPI_on0 applied_on0 foci_on0 new_DAPI3D_save bw_applied3D_save v_save mask_stack im_foci im_foci3D im_foci3D_save layer_I layer_num  bw_cut
if get(handles.current_seg,'Value')
    label_DAPI = bwlabel(bwmorph(bw_DAPI,'thicken',1));
    sel_DAPI = ismember(label_DAPI,unique(label_DAPI(sel_DAPI)));
    bw_DAPI = bw_DAPI | sel_DAPI;
else
    label_applied = bwlabel(bwmorph(bw_applied,'thicken',1));
    sel_applied = ismember(label_applied,unique(label_applied(sel_applied)));
    bw_applied = bw_applied | sel_applied;
end
%%% Image display
overlay = im_overlay(new_DAPI,sel_DAPI,sel_applied,bw_cut,bw_DAPI,bw_applied,im_foci,handles);
axes(handles.im_show);
v = axis; 
imshow(overlay)
axis(v);




% --- Executes on button press in erosion_r.
function erosion_r_Callback(hObject, eventdata, handles)
% hObject    handle to erosion_r (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of erosion_r
global sub_list sub_num in_folder input_name out_folder mask_name mat_tail resolution0 image_folder DAPI_channel DAPI_th0 new_DAPI new_DAPI3D bw_DAPI bw_applied bw_applied3D sel_DAPI sel_applied L_ratio DAPI_on0 applied_on0 foci_on0 new_DAPI3D_save bw_applied3D_save v_save mask_stack im_foci im_foci3D im_foci3D_save layer_I layer_num  bw_cut
if get(handles.current_seg,'Value')
    bw_DAPI = bw_DAPI & (~sel_DAPI);
    sel_DAPI = imerode(sel_DAPI, strel('disk',2));
    bw_DAPI = bw_DAPI | sel_DAPI;
else
    bw_applied = bw_applied & (~sel_applied);
    sel_applied = imerode(sel_applied, strel('disk',2));
    bw_applied = bw_applied | sel_applied;
end
%%% Image display
overlay = im_overlay(new_DAPI,sel_DAPI,sel_applied,bw_cut,bw_DAPI,bw_applied,im_foci,handles);
axes(handles.im_show);
v = axis; 
imshow(overlay)
axis(v);




% --- Executes on button press in applied_on.
function applied_on_Callback(hObject, eventdata, handles)
% hObject    handle to applied_on (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of applied_on
global sub_list sub_num in_folder input_name out_folder mask_name mat_tail resolution0 image_folder DAPI_channel DAPI_th0 new_DAPI new_DAPI3D bw_DAPI bw_applied bw_applied3D sel_DAPI sel_applied L_ratio DAPI_on0 applied_on0 foci_on0 new_DAPI3D_save bw_applied3D_save v_save mask_stack im_foci im_foci3D im_foci3D_save layer_I layer_num  bw_cut
applied_on0 = get(handles.applied_on,'Value');
%%% Image display
overlay = im_overlay(new_DAPI,sel_DAPI,sel_applied,bw_cut,bw_DAPI,bw_applied,im_foci,handles);
axes(handles.im_show);
v = axis; 
imshow(overlay)
axis(v);




% --- Executes on button press in delete_seg.
function delete_seg_Callback(hObject, eventdata, handles)
% hObject    handle to delete_seg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of delete_seg
global sub_list sub_num in_folder input_name out_folder mask_name mat_tail resolution0 image_folder DAPI_channel DAPI_th0 new_DAPI new_DAPI3D bw_DAPI bw_applied bw_applied3D sel_DAPI sel_applied L_ratio DAPI_on0 applied_on0 foci_on0 new_DAPI3D_save bw_applied3D_save v_save mask_stack im_foci im_foci3D im_foci3D_save layer_I layer_num  bw_cut
if get(handles.current_seg,'Value')
    bw_DAPI(sel_DAPI) = false;
    sel_DAPI(sel_DAPI) = false;
else
    bw_applied(sel_applied) = false;
    sel_applied(sel_applied) = false;
end
%%% Image display
overlay = im_overlay(new_DAPI,sel_DAPI,sel_applied,bw_cut,bw_DAPI,bw_applied,im_foci,handles);
axes(handles.im_show);
v = axis; 
imshow(overlay)
axis(v);




function seg_group_SelectionChangeFcn(hObject, eventdata)
 
%retrieve GUI data, i.e. the handles structure
handles = guidata(hObject); 
global sub_list sub_num in_folder input_name out_folder mask_name mat_tail resolution0 image_folder DAPI_channel DAPI_th0 new_DAPI new_DAPI3D bw_DAPI bw_applied bw_applied3D sel_DAPI sel_applied L_ratio DAPI_on0 applied_on0 foci_on0 new_DAPI3D_save bw_applied3D_save v_save mask_stack im_foci im_foci3D im_foci3D_save layer_I layer_num bw_cut
 
switch get(eventdata.NewValue,'Tag')   % Get Tag of selected object
    case 'sel_cls'
      %execute this code when "select/clear" button is selected
        if ~get(handles.manual_draw,'Value') && ~get(handles.manual_erase,'Value')
            if get(handles.sel_cls,'Value')
                axes(handles.im_show);
                zoom off
                pan off
                bw_cut = false(size(bw_DAPI));
                set(gcf,'WindowButtonDownFcn',{@ButttonDownFcn1,handles}); set(gcf,'KeyPressFcn',{@KeyPressFcn0,hObject,eventdata,handles});
            else
                bw_cut = false(size(bw_DAPI));
                axes(handles.im_show);
                set(gcf,'WindowButtonDownFcn',{}); set(gcf,'KeyPressFcn',{@KeyPressFcn0,hObject,eventdata,handles});
            end
        else
            bw_cut = false(size(bw_DAPI));
            set(handles.sel_cls,'Value',0);
        end
        
    case 'manual_draw'
      %execute this code when "manual draw" button is selected        
        if get(handles.manual_draw,'Value')
            set(handles.manual_draw,'ForegroundColor',[1,0,0]);
            axes(handles.im_show);
            zoom off
            pan off
            set(gcf,'WindowButtonDownFcn',{}); set(gcf,'KeyPressFcn',{});
            bw_cut = false(size(bw_DAPI));
            add_switch = true;
            while add_switch
                %%% Image display
                overlay = im_overlay(new_DAPI,sel_DAPI,sel_applied,bw_cut,bw_DAPI,bw_applied,im_foci,handles);
                axes(handles.im_show);
                v = axis; 
                imshow(overlay)
                axis(v);
                
                bw_add = roipoly;
                add_switch = any(any(bw_add));
                if add_switch
                    if get(handles.current_seg,'Value')
                        bw_DAPI = bw_DAPI | bw_add;
                        label_DAPI = bwlabel(bw_DAPI);
                        sel_DAPI = sel_DAPI | ismember(label_DAPI,unique(label_DAPI(logical(bw_add))));
                    else
                        bw_applied = bw_applied | bw_add;
                        label_applied = bwlabel(bw_applied);
                        sel_applied = sel_applied | ismember(label_applied,unique(label_applied(logical(bw_add))));
                    end
                end
                add_switch = false;
            end
            %%% Image display
            overlay = im_overlay(new_DAPI,sel_DAPI,sel_applied,bw_cut,bw_DAPI,bw_applied,im_foci,handles);
            axes(handles.im_show);
            v = axis; 
            imshow(overlay)
            axis(v);
            
            set(handles.manual_draw,'Value',0);
            set(handles.manual_draw,'ForegroundColor',[0,0,0]);
            set(gcf,'WindowButtonDownFcn',{}); set(gcf,'KeyPressFcn',{@KeyPressFcn0,hObject,eventdata,handles});
       end
 
    case 'manual_erase'
      %execute this code when "manual erase" button is selected        
        if get(handles.manual_erase,'Value')
            set(handles.manual_erase,'ForegroundColor',[1,0,0]);
            axes(handles.im_show);
            zoom off
            pan off
            set(gcf,'WindowButtonDownFcn',{}); set(gcf,'KeyPressFcn',{});
            add_switch = true;
            bw_cut = false(size(bw_DAPI));
            while add_switch
                %%% Image display
                overlay = im_overlay(new_DAPI,sel_DAPI,sel_applied,bw_cut,bw_DAPI,bw_applied,im_foci,handles);
                axes(handles.im_show);
                v = axis; 
                imshow(overlay)
                axis(v);
                
                bw_add = roipoly;
                add_switch = any(any(bw_add));
                if add_switch
                    if get(handles.current_seg,'Value')
                        bw_DAPI = bw_DAPI & (~bw_add);
                        sel_DAPI = sel_DAPI & (~bw_add);
                    else
                        bw_applied = bw_applied & (~bw_add);
                        sel_applied = sel_applied & (~bw_add);
                    end
                end
                add_switch = false;
            end
            %%% Image display
            overlay = im_overlay(new_DAPI,sel_DAPI,sel_applied,bw_cut,bw_DAPI,bw_applied,im_foci,handles);
            axes(handles.im_show);
            v = axis; 
            imshow(overlay)
            axis(v);
            
            set(handles.manual_erase,'Value',0);
            set(handles.manual_erase,'ForegroundColor',[0,0,0]);
            set(gcf,'WindowButtonDownFcn',{}); set(gcf,'KeyPressFcn',{@KeyPressFcn0,hObject,eventdata,handles});
        end

    case 'auto_cut'
      %execute this code when "select/clear" button is selected
        if ~get(handles.manual_draw,'Value') && ~get(handles.manual_erase,'Value')
            if get(handles.auto_cut,'Value')
                axes(handles.im_show);
                zoom off
                pan off
                set(gcf,'WindowButtonDownFcn',{@ButttonDownFcn2,handles});
                set(gcf,'KeyPressFcn',{@KeyPressFcn2,handles});
                
                
            else
                axes(handles.im_show);
                set(gcf,'WindowButtonDownFcn',{}); set(gcf,'KeyPressFcn',{@KeyPressFcn0,hObject,eventdata,handles});
            end
        else
            set(handles.auto_cut,'Value',0);
        end

end
%updates the handles structure
guidata(hObject, handles);




% --- Executes on button press in seg_fill.
function seg_fill_Callback(hObject, eventdata, handles)
% hObject    handle to seg_fill (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of seg_fill
global sub_list sub_num in_folder input_name out_folder mask_name mat_tail resolution0 image_folder DAPI_channel DAPI_th0 new_DAPI new_DAPI3D bw_DAPI bw_applied bw_applied3D sel_DAPI sel_applied L_ratio DAPI_on0 applied_on0 foci_on0 new_DAPI3D_save bw_applied3D_save v_save mask_stack im_foci im_foci3D im_foci3D_save layer_I layer_num bw_cut
if get(handles.current_seg,'Value')
    sel_DAPI = imfill(sel_DAPI,'holes');
    bw_DAPI = bw_DAPI | sel_DAPI;
else
    sel_applied = imfill(sel_applied,'holes');
    bw_applied = bw_applied | sel_applied;
end
%%% Image display
overlay = im_overlay(new_DAPI,sel_DAPI,sel_applied,bw_cut,bw_DAPI,bw_applied,im_foci,handles);
axes(handles.im_show);
v = axis; 
imshow(overlay)
axis(v);




% --- Executes on button press in im_lock.
function im_lock_Callback(hObject, eventdata, handles)
% hObject    handle to im_lock (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of im_lock
global sub_list sub_num in_folder input_name out_folder mask_name mat_tail resolution0 image_folder DAPI_channel DAPI_th0 new_DAPI new_DAPI3D bw_DAPI bw_applied bw_applied3D sel_DAPI sel_applied L_ratio DAPI_on0 applied_on0 foci_on0 new_DAPI3D_save bw_applied3D_save v_save mask_stack im_foci im_foci3D im_foci3D_save layer_I layer_num area_thresh bw_cut
set(handles.im_lock,'String','Wait');
%updates the handles structure
guidata(hObject, handles);

if get(handles.im_lock,'Value')   %%% Lock
    bw_applied3D(:,:,layer_I) = bw_applied;
    bw_applied3D_save(v_save(3):v_save(4),v_save(1):v_save(2),:) = bw_applied3D;
    
    v_save = floor(axis-0.5)+[1,0,1,0];
    v_save = [max(v_save(1),1),min(v_save(2),size(bw_applied3D_save,2)),max(v_save(3),1),min(v_save(4),size(bw_applied3D_save,1))];
%    set(hObject,'toolbar','none');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    im_foci3D = im_foci3D_save(v_save(3):v_save(4),v_save(1):v_save(2),:);
    im_foci = im_foci3D(:,:,layer_I);
    
    new_DAPI3D = new_DAPI3D_save(v_save(3):v_save(4),v_save(1):v_save(2),:);
    new_DAPI = new_DAPI3D(:,:,layer_I);
    
    bw_DAPI = im2bw(new_DAPI,DAPI_th0);
    bw_DAPI = imfill(bw_DAPI,'holes');
    bw_DAPI = bwareaopen(bw_DAPI, area_thresh);

    sel_DAPI = false(size(new_DAPI));
    sel_applied = false(size(new_DAPI));
    
    bw_applied3D = bw_applied3D_save(v_save(3):v_save(4),v_save(1):v_save(2),:);
    bw_applied = bw_applied3D(:,:,layer_I);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    set(handles.im_lock,'String','Unlock');
    bw_cut = false(size(new_DAPI));
    %%% Image display
    overlay = im_overlay(new_DAPI,sel_DAPI,sel_applied,bw_cut,bw_DAPI,bw_applied,im_foci,handles);
    imshow(overlay)
    
else                                  %%% Unlock
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    im_foci3D = im_foci3D_save;
    im_foci = im_foci3D(:,:,layer_I);
    
    new_DAPI3D = new_DAPI3D_save;
    new_DAPI = new_DAPI3D(:,:,layer_I);
    
    bw_applied3D(:,:,layer_I) = bw_applied;
    bw_applied3D_save(v_save(3):v_save(4),v_save(1):v_save(2),:) = bw_applied3D;
    bw_applied3D = bw_applied3D_save;
    bw_applied = bw_applied3D(:,:,layer_I);

    bw_DAPI = im2bw(new_DAPI,DAPI_th0);
    bw_DAPI = imfill(bw_DAPI,'holes');
    bw_DAPI = bwareaopen(bw_DAPI, area_thresh);

    sel_DAPI = false(size(new_DAPI));
    sel_applied = false(size(new_DAPI));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    set(hObject,'toolbar','figure');
    set(handles.im_lock,'String','Lock');
    bw_cut = false(size(new_DAPI));
    
    overlay = im_overlay(new_DAPI,sel_DAPI,sel_applied,bw_cut,bw_DAPI,bw_applied,im_foci,handles);
    imshow(overlay)
    axis(v_save+[-0.5,0.5,-0.5,0.5]);
    v_save = [1,size(new_DAPI3D_save,2),1,size(new_DAPI3D_save,1)];
end

%updates the handles structure
guidata(hObject, handles);

    


% --- Executes on button press in apply_seg.
function apply_seg_Callback(hObject, eventdata, handles)
% hObject    handle to apply_seg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of apply_seg
global sub_list sub_num in_folder input_name out_folder mask_name mat_tail resolution0 image_folder DAPI_channel DAPI_th0 new_DAPI new_DAPI3D bw_DAPI bw_applied bw_applied3D sel_DAPI sel_applied L_ratio DAPI_on0 applied_on0 foci_on0 new_DAPI3D_save bw_applied3D_save v_save mask_stack im_foci im_foci3D im_foci3D_save layer_I layer_num bw_cut
if get(handles.current_seg,'Value')
    bw_applied = bw_applied | sel_DAPI;
    sel_applied = sel_applied | sel_DAPI;
    
    label_bw = bwlabel(bw_applied);
    label_sel = unique(label_bw(sel_applied));
    sel_applied = sel_applied | ismember(label_bw,label_sel);
    
    %%% Image display
    overlay = im_overlay(new_DAPI,sel_DAPI,sel_applied,bw_cut,bw_DAPI,bw_applied,im_foci,handles);
    axes(handles.im_show);
    v = axis; 
    imshow(overlay)
    axis(v);
end
set(gcf,'KeyPressFcn',{@KeyPressFcn0,hObject,eventdata,handles});




% --- Executes on button press in all_sel.
function all_sel_Callback(hObject, eventdata, handles)
% hObject    handle to all_sel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global sub_list sub_num in_folder input_name out_folder mask_name mat_tail resolution0 image_folder DAPI_channel DAPI_th0 new_DAPI new_DAPI3D bw_DAPI bw_applied bw_applied3D sel_DAPI sel_applied L_ratio DAPI_on0 applied_on0 foci_on0 new_DAPI3D_save bw_applied3D_save v_save mask_stack im_foci im_foci3D im_foci3D_save layer_I layer_num bw_cut
if get(handles.current_seg,'Value')
    if get(handles.sel_empty,'Value')
        bw_DAPI_prop = regionprops(bw_DAPI,bw_applied,'MaxIntensity','Area');
        bw_DAPI_empty = find(~[bw_DAPI_prop.MaxIntensity] & [bw_DAPI_prop.Area] > 500);
        bw_DAPI0 = imclearborder(ismember(bwlabel(bw_DAPI),bw_DAPI_empty));
    elseif get(handles.sel_single,'Value')
        bw_DAPI_prop = regionprops(bw_DAPI,bwlabel(bw_applied),'MaxIntensity','MinIntensity');
        bw_DAPI_prop2 = regionprops(bw_DAPI,bwlabel(bw_applied)+1000000*(~bw_applied),'MinIntensity');
        bw_DAPI_prop3 = regionprops(bw_DAPI,bw_applied,'MeanIntensity','Area');
%         bw_DAPI_single = find(([bw_DAPI_prop3.MeanIntensity] < 0.999) & ([bw_DAPI_prop.MaxIntensity] == [bw_DAPI_prop2.MinIntensity]) & ([bw_DAPI_prop.MaxIntensity] > [bw_DAPI_prop.MinIntensity]));
        bw_DAPI_single = find(([bw_DAPI_prop3.Area] >= 200) & ([bw_DAPI_prop3.MeanIntensity] < 0.85) & ([bw_DAPI_prop.MaxIntensity] == [bw_DAPI_prop2.MinIntensity]) & ([bw_DAPI_prop.MaxIntensity] > [bw_DAPI_prop.MinIntensity]));
        bw_DAPI0 = imclearborder(ismember(bwlabel(bw_DAPI),bw_DAPI_single));
        
% %         bw_DAPI_prop4 = regionprops(bw_applied,bw_DAPI,'MeanIntensity');
% %         ind_applied_sel = find([bw_DAPI_prop4.MeanIntensity] < 0.90);
% %         label_applied = bwlabel(bw_applied);
% %         true_DAPI_sel = ismember(label_applied,ind_applied_sel);
% %         label_DAPI = bwlabel(bw_DAPI);
% %         bw_DAPI_single = setdiff(unique(label_DAPI(true_DAPI_sel)),0);
% %         bw_DAPI0 = imclearborder(ismember(bwlabel(bw_DAPI),bw_DAPI_single));
        
    else
        bw_DAPI_prop = regionprops(bw_DAPI,bwlabel(bw_applied),'MaxIntensity');
        bw_DAPI_prop2 = regionprops(bw_DAPI,bwlabel(bw_applied)+1000000*(~bw_applied),'MinIntensity');
        bw_DAPI_multi = find([bw_DAPI_prop.MaxIntensity] > [bw_DAPI_prop2.MinIntensity]);
        bw_DAPI0 = imclearborder(ismember(bwlabel(bw_DAPI),bw_DAPI_multi));
    end
    sel_DAPI = bw_DAPI0 | sel_DAPI;
%     sel_DAPI = bw_DAPI | sel_DAPI;
else
    if get(handles.auto_cut,'Value')
        bw_applied_prop = regionprops(bw_applied & sel_applied,'EquivDiameter');
        r0 = max(round(mean([bw_applied_prop.EquivDiameter])/2*0.4),18);
        disk0 = getnhood(strel('disk',r0));
        disk0 = disk0/nnz(disk0);
        bw_cut = imdilate(conv2(double(sel_applied),disk0,'same') >= 0.99,strel('disk',2));
        set(gcf,'KeyPressFcn',{@KeyPressFcn2,handles});
    else
        if get(handles.sel_single,'Value')
% %     % %         bw_applied_prop3 = regionprops(bw_applied,bw_DAPI,'MeanIntensity');
% %     % %         bw_applied_single = find(([bw_applied_prop3.MeanIntensity] == 1);

% %             bw_applied_prop = regionprops(bw_applied,bwlabel(bw_DAPI),'MaxIntensity');
% %             bw_DAPI_prop3 = regionprops(bw_DAPI,bw_applied,'MeanIntensity');
% %             bw_applied_prop2 = regionprops(bw_applied,bw_DAPI,'MeanIntensity');
% %             
% %             bw_DAPI_single = find([bw_DAPI_prop3.MeanIntensity] < 0.92);
% %             bw_applied_single = find(ismember([bw_applied_prop.MaxIntensity],bw_DAPI_single));
            
% % % %             bw_applied_single = find([bw_applied_prop2.MeanIntensity] < 0.92);
% % 
% %             if any(sel_applied(:))
% %                 sel_applied = imclearborder(ismember(bwlabel(bw_applied),bw_applied_single) & sel_applied);
% %             else
% %                 sel_applied = imclearborder(ismember(bwlabel(bw_applied),bw_applied_single));
% %             end
            
            label_applied = bwlabel(bw_applied);
            bw_applied_single = setdiff(unique(label_applied(sel_DAPI)),0);
            sel_applied = imclearborder(ismember(bwlabel(bw_applied),bw_applied_single));
            
        elseif get(handles.sel_multi,'Value')
            bw_applied_prop3 = regionprops(bw_applied,bw_DAPI,'MeanIntensity');
            bw_applied_multi = find([bw_applied_prop3.MeanIntensity] < 0.90);
            sel_applied = imclearborder(ismember(bwlabel(bw_applied),bw_applied_multi));
        else
% %             sel_applied = imclearborder(bw_applied) | sel_applied;
% %     %     sel_applied = bw_applied | sel_applied;
            bw_applied_prop3 = regionprops(bw_applied,bw_DAPI,'MaxIntensity');
            bw_applied_multi = find([bw_applied_prop3.MaxIntensity] == 0);
            sel_applied = imclearborder(ismember(bwlabel(bw_applied),bw_applied_multi));
        end
    end
end
%%% Image display
overlay = im_overlay(new_DAPI,sel_DAPI,sel_applied,bw_cut,bw_DAPI,bw_applied,im_foci,handles);
axes(handles.im_show);
v = axis; 
imshow(overlay)
axis(v);



% --- Executes on button press in cls_sel.
function cls_sel_Callback(hObject, eventdata, handles)
% hObject    handle to cls_sel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global sub_list sub_num in_folder input_name out_folder mask_name mat_tail resolution0 image_folder DAPI_channel DAPI_th0 new_DAPI new_DAPI3D bw_DAPI bw_applied bw_applied3D sel_DAPI sel_applied L_ratio DAPI_on0 applied_on0 foci_on0 new_DAPI3D_save bw_applied3D_save v_save mask_stack im_foci im_foci3D im_foci3D_save layer_I layer_num bw_cut
if get(handles.current_seg,'Value')
    sel_DAPI = false(size(sel_DAPI));
else
    sel_applied = false(size(sel_applied));
end
%%% Image display
overlay = im_overlay(new_DAPI,sel_DAPI,sel_applied,bw_cut,bw_DAPI,bw_applied,im_foci,handles);
axes(handles.im_show);
v = axis; 
imshow(overlay)
axis(v);




function ButttonDownFcn1(src,event,handles)
global sub_list sub_num in_folder input_name out_folder mask_name mat_tail resolution0 image_folder DAPI_channel DAPI_th0 new_DAPI new_DAPI3D bw_DAPI bw_applied bw_applied3D sel_DAPI sel_applied L_ratio DAPI_on0 applied_on0 foci_on0 new_DAPI3D_save bw_applied3D_save v_save mask_stack im_foci im_foci3D im_foci3D_save layer_I layer_num bw_cut
pt = get(gca,'CurrentPoint');
x = floor(pt(1,1)+0.5);
y = floor(pt(1,2)+0.5);
if get(handles.current_seg,'Value')
    label_DAPI = bwlabel(bw_DAPI);
    label_sel = bwlabel(sel_DAPI);
    if label_sel(y,x)
        sel_DAPI(label_sel == label_sel(y,x)) = false;
    elseif label_DAPI(y,x)
        sel_DAPI(label_DAPI == label_DAPI(y,x)) = true;
    end
else
    label_DAPI = bwlabel(bw_applied);
    label_sel = bwlabel(sel_applied);
    if label_sel(y,x)
        sel_applied(label_sel == label_sel(y,x)) = false;
    elseif label_DAPI(y,x)
        sel_applied(label_DAPI == label_DAPI(y,x)) = true;
    end
end
%%% Image display
overlay = im_overlay(new_DAPI,sel_DAPI,sel_applied,bw_cut,bw_DAPI,bw_applied,im_foci,handles);
axes(handles.im_show);
v = axis; 
imshow(overlay)
axis(v);




function ButttonDownFcn2(src,event,handles)
global sub_list sub_num in_folder input_name out_folder mask_name mat_tail resolution0 image_folder DAPI_channel DAPI_th0 new_DAPI new_DAPI3D bw_DAPI bw_applied bw_applied3D sel_DAPI sel_applied L_ratio DAPI_on0 applied_on0 foci_on0 new_DAPI3D_save bw_applied3D_save v_save mask_stack im_foci im_foci3D im_foci3D_save layer_I layer_num bw_cut
pt = get(gca,'CurrentPoint');
x = floor(pt(1,1)+0.5);
y = floor(pt(1,2)+0.5);
rxy = 5;
% if get(handles.current_seg,'Value')
%     label_DAPI = bwlabel(bw_DAPI);
%     label_sel = bwlabel(sel_DAPI);
%     if label_sel(y,x)
%         sel_DAPI(label_sel == label_sel(y,x)) = false;
%     elseif label_DAPI(y,x)
%         sel_DAPI(label_DAPI == label_DAPI(y,x)) = true;
%     end
% else
%     label_DAPI = bwlabel(bw_applied);
%     label_sel = bwlabel(sel_applied);
%     if label_sel(y,x)
%         sel_applied(label_sel == label_sel(y,x)) = false;
%     elseif label_DAPI(y,x)
%         sel_applied(label_DAPI == label_DAPI(y,x)) = true;
%     end
% end
xmin = max(1,x-rxy);
xmax = min(size(bw_cut,2),x+rxy);
ymin = max(1,y-rxy);
ymax = min(size(bw_cut,1),y+rxy);

if any(any(bw_cut(ymin:ymax,xmin:xmax)))
    bw_cut_label = bwlabel(bw_cut);
    I00 = max(max(bw_cut_label(ymin:ymax,xmin:xmax)));
    bw_cut = bw_cut & ~ismember(bw_cut_label,I00);
else
    bw_cut(y,x) = true;
end

%%% Image display
overlay = im_overlay(new_DAPI,sel_DAPI,sel_applied,bw_cut,bw_DAPI,bw_applied,im_foci,handles);
axes(handles.im_show);
v = axis; 
imshow(overlay)
axis(v);





function KeyPressFcn0(obj,evt,hObject,eventdata,handles)
if strcmp(evt.Key,'uparrow')
    temp = str2double(get(handles.I_layer,'String'));
    set(handles.I_layer,'String',num2str(temp+1));
    I_layer_Callback(hObject, eventdata, handles);
elseif strcmp(evt.Key,'downarrow')
    temp = str2double(get(handles.I_layer,'String'));
    set(handles.I_layer,'String',num2str(temp-1));
    I_layer_Callback(hObject, eventdata, handles);
% elseif strcmp(evt.Key,'control')
%     all_sel_Callback(hObject, eventdata, handles);
% elseif strcmp(evt.Key,'shift')
%     conv_all_Callback(hObject, eventdata, handles);
% elseif strcmp(evt.Key,'space')
%     apply_seg_Callback(hObject, eventdata, handles);
end






function KeyPressFcn2(obj,evt,handles)
global sub_list sub_num in_folder input_name out_folder mask_name mat_tail resolution0 image_folder DAPI_channel DAPI_th0 new_DAPI new_DAPI3D bw_DAPI bw_applied bw_applied3D sel_DAPI sel_applied L_ratio DAPI_on0 applied_on0 foci_on0 new_DAPI3D_save bw_applied3D_save v_save mask_stack im_foci im_foci3D im_foci3D_save layer_I layer_num  bw_cut
if double(evt.Character) == 13
    quit_reply = questdlg('Show cutting?');
    if strcmp(quit_reply,'Yes')
        temp_label = bwlabel(bw_applied);
        temp_sel = ismember(temp_label,temp_label(bw_cut));
        D= bwdist(~bw_applied);
        g2 = imimposemin(-D,bw_cut);
        bw_cut = imdilate(watershed(g2) == 0,strel('disk',1)) & (temp_sel);
        %%% Image display
        overlay = im_overlay(new_DAPI,sel_DAPI,sel_applied,bw_cut,bw_DAPI,bw_applied,im_foci,handles);
        axes(handles.im_show);
        v = axis; 
        imshow(overlay)
        axis(v);

        quit_reply = questdlg('Apply cutting?');
        if strcmp(quit_reply,'Yes')
            bw_applied = bw_applied & (~ bw_cut);
            sel_applied = sel_applied & (~ bw_cut);
        end
        bw_cut = false(size(new_DAPI));
%         set(gcf,'WindowButtonDownFcn',{}); set(gcf,'KeyPressFcn',{});

        %%% Image display
        overlay = im_overlay(new_DAPI,sel_DAPI,sel_applied,bw_cut,bw_DAPI,bw_applied,im_foci,handles);
        axes(handles.im_show);
        v = axis; 
        imshow(overlay)
        axis(v);
    end
    
    
elseif double(evt.Character) == 27
    bw_cut = false(size(new_DAPI));
%     set(gcf,'WindowButtonDownFcn',{}); set(gcf,'KeyPressFcn',{});
    %%% Image display
    overlay = im_overlay(new_DAPI,sel_DAPI,sel_applied,bw_cut,bw_DAPI,bw_applied,im_foci,handles);
    axes(handles.im_show);
    v = axis; 
    imshow(overlay)
    axis(v);

end





function th_closereq(hObject, eventdata, handles)
global sub_list sub_num in_folder input_name out_folder mask_name mat_tail resolution0 image_folder DAPI_channel DAPI_th0 new_DAPI new_DAPI3D bw_DAPI bw_applied bw_applied3D sel_DAPI sel_applied L_ratio DAPI_on0 applied_on0 foci_on0 new_DAPI3D_save bw_applied3D_save v_save mask_stack im_foci im_foci3D im_foci3D_save layer_I layer_num  bw_cut

quit_reply = questdlg('Quit the program (Make sure you have saved the result)?');
if strcmp(quit_reply,'Yes')
    clear global
    delete(gcf)
else
    return
end




% --- Executes on button press in conv_all.
function conv_all_Callback(hObject, eventdata, handles)
% hObject    handle to conv_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global sub_list sub_num in_folder input_name out_folder mask_name mat_tail resolution0 image_folder DAPI_channel DAPI_th0 new_DAPI new_DAPI3D bw_DAPI bw_applied bw_applied3D sel_DAPI sel_applied L_ratio DAPI_on0 applied_on0 foci_on0 new_DAPI3D_save bw_applied3D_save v_save mask_stack im_foci im_foci3D im_foci3D_save layer_I layer_num area_thresh bw_cut

quit_reply = questdlg('Convert all nuclei to convex objects?');
    
if strcmp(quit_reply,'Yes')
    bw_applied3D(:,:,layer_I) = bw_applied;
    bw_applied3D_save(v_save(3):v_save(4),v_save(1):v_save(2),:) = bw_applied3D;


    bw_applied3D_save = conv3D(bw_applied3D_save);

    bw_applied3D = bw_applied3D_save(v_save(3):v_save(4),v_save(1):v_save(2),:);
    bw_applied = bw_applied3D(:,:,layer_I);

    if any(any(sel_applied))
        sel_prop = regionprops(sel_applied, 'Centroid');
        centroid_xy = round(cell2mat({sel_prop.Centroid}'));
        applied_label = bwlabel(bw_applied);
        sel_list = applied_label(sub2ind(size(sel_applied),centroid_xy(:,2),centroid_xy(:,1)));
        sel_applied = ismember(applied_label,sel_list);
    end

    %%% Image display
    overlay = im_overlay(new_DAPI,sel_DAPI,sel_applied,bw_cut,bw_DAPI,bw_applied,im_foci,handles);
    axes(handles.im_show);
    v = axis; 
    imshow(overlay)
    axis(v);
end






% --- Executes on button press in label_layer.
function label_layer_Callback(hObject, eventdata, handles)
% hObject    handle to label_layer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global sub_list sub_num in_folder input_name out_folder mask_name mat_tail resolution0 image_folder  DAPI_th0 new_DAPI new_DAPI3D bw_DAPI bw_applied bw_applied3D sel_DAPI sel_applied L_ratio DAPI_on0 applied_on0 foci_on0 new_DAPI3D_save bw_applied3D_save v_save mask_stack im_foci im_foci3D im_foci3D_save layer_I layer_num  bw_cut
quit_reply = questdlg('Label the segmentation result?');
if strcmp(quit_reply,'Yes')
    set(handles.label_layer,'String','Wait');
    % Update handles structure
    guidata(hObject, handles);
    
    bw_applied3D(:,:,layer_I) = bw_applied;
    bw_applied3D_save(v_save(3):v_save(4),v_save(1):v_save(2),:) = bw_applied3D;

    mask_stack = label3D(bw_applied3D_save);
    
    set(handles.label_layer,'String','Label');
    % Update handles structure
    guidata(hObject, handles);
end







% --- Executes on button press in thresh_apply.
function thresh_apply_Callback(hObject, eventdata, handles)
% hObject    handle to thresh_apply (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global sub_list sub_num in_folder input_name out_folder mask_name mat_tail resolution0 image_folder DAPI_channel DAPI_th0 new_DAPI new_DAPI3D bw_DAPI bw_applied bw_applied3D sel_DAPI sel_applied L_ratio DAPI_on0 applied_on0 foci_on0 new_DAPI3D_save bw_applied3D_save v_save mask_stack im_foci im_foci3D im_foci3D_save layer_I layer_num  bw_cut area_thresh

apply_reply = questdlg('Apply threshold to the whole stack?');
if strcmp(apply_reply,'Yes')
    set(handles.thresh_apply,'String','Wait');
    
    %get the string for the editText component
    DAPI_th0 = str2num(get(handles.DAPI_th,'String'));
    sel_DAPI = false(size(sel_DAPI));
    sel_applied = false(size(sel_applied));
    for I_layer = 1:size(bw_applied3D,3)
        bw_applied3D(:,:,I_layer) = im2bw(new_DAPI3D(:,:,I_layer),DAPI_th0);
        bw_applied3D(:,:,I_layer) = imfill(bw_applied3D(:,:,I_layer),'holes');
        bw_applied3D(:,:,I_layer) = bwareaopen(bw_applied3D(:,:,I_layer), area_thresh);
        bw_applied3D(:,:,I_layer) = imclearborder(bw_applied3D(:,:,I_layer));
    end
    bw_applied3D_save(v_save(3):v_save(4),v_save(1):v_save(2),:) = bw_applied3D;
    bw_applied = bw_applied3D(:,:,layer_I);
    
    %%% Image display
    overlay = im_overlay(new_DAPI,sel_DAPI,sel_applied,bw_cut,bw_DAPI,bw_applied,im_foci,handles);
    axes(handles.im_show);
    v = axis; 
    imshow(overlay)
    axis(v);
    
    set(handles.thresh_apply,'String','Apply to all');
end
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in foci_on.
function foci_on_Callback(hObject, eventdata, handles)
% hObject    handle to foci_on (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of foci_on
global sub_list sub_num in_folder input_name out_folder mask_name mat_tail resolution0 image_folder DAPI_channel DAPI_th0 new_DAPI new_DAPI3D bw_DAPI bw_applied bw_applied3D sel_DAPI sel_applied L_ratio DAPI_on0 applied_on0 foci_on0 new_DAPI3D_save bw_applied3D_save v_save mask_stack im_foci im_foci3D im_foci3D_save layer_I layer_num  bw_cut
foci_on0 = get(handles.foci_on,'Value');
%%% Image display
overlay = im_overlay(new_DAPI,sel_DAPI,sel_applied,bw_cut,bw_DAPI,bw_applied,im_foci,handles);
axes(handles.im_show);
v = axis; 
imshow(overlay)
axis(v);


function overlay = im_overlay(new_DAPI,sel_DAPI,sel_applied,bw_cut,bw_DAPI,bw_applied,im_foci,handles)
DAPI_on0 = get(handles.DAPI_on,'Value');
applied_on0 = get(handles.applied_on,'Value');
foci_on0 = get(handles.foci_on,'Value');
r0 = 0.3;

overlay(:,:,1) = r0*new_DAPI+0.3*DAPI_on0*sel_DAPI+applied_on0*bwperim(sel_applied)-imdilate(bw_cut,strel('disk',5))+foci_on0*(im_foci & applied_on0*(~bwperim(bw_applied)));
overlay(:,:,2) = r0*new_DAPI+0.3*DAPI_on0*(bw_DAPI & (~sel_DAPI))+applied_on0*bwperim(bw_applied & (~sel_applied))-imdilate(bw_cut,strel('disk',5));
overlay(:,:,3) = r0*new_DAPI+imdilate(bw_cut,strel('disk',5))+foci_on0*(im_foci & applied_on0*(~bwperim(bw_applied)));



function v_input_Callback(hObject, eventdata, handles)
% hObject    handle to v_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of v_input as text
%        str2double(get(hObject,'String')) returns contents of v_input as a double

global sub_list sub_num in_folder input_name out_folder mask_name mat_tail resolution0 image_folder DAPI_channel DAPI_th0 new_DAPI new_DAPI3D bw_DAPI bw_applied bw_applied3D sel_DAPI sel_applied L_ratio DAPI_on0 applied_on0 foci_on0 new_DAPI3D_save bw_applied3D_save v_save mask_stack im_foci im_foci3D im_foci3D_save layer_I layer_num  bw_cut area_thresh

v0 = size(bw_applied);

try
    xy_range0 = eval(get(handles.v_input,'String'));
    if xy_range0(2) > xy_range0(1) && xy_range0(4) > xy_range0(3) && xy_range0(1) >= 0 && xy_range0(2) <= v0(2)  && xy_range0(3) >= 0 && xy_range0(4) <= v0(1)
        is_vinput = true;
        v_update = xy_range0;
    else
        is_vinput = false;
    end
catch
    is_vinput = false;
end

if ~get(handles.im_lock,'Value') && is_vinput
    axes(handles.im_show);
    v_current = axis; 
    axis(v_update);
    
    quit_reply = questdlg('Move to this range?');
    if ~strcmp(quit_reply,'Yes')
        axes(handles.im_show);
        axis(v_current);
    end
end

set(handles.v_input,'String','')






% --- Executes during object creation, after setting all properties.
function v_input_CreateFcn(hObject, eventdata, handles)
% hObject    handle to v_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
