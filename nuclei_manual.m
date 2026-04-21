function varargout = nuclei_manual(varargin)
% NUCLEI_MANUAL M-file for nuclei_manual.fig
%      NUCLEI_MANUAL, by itself, creates a new NUCLEI_MANUAL or raises the existing
%      singleton*.
%
%      H = NUCLEI_MANUAL returns the handle to a new NUCLEI_MANUAL or the handle to
%      the existing singleton*.
%
%      NUCLEI_MANUAL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in NUCLEI_MANUAL.M with the given input arguments.
%
%      NUCLEI_MANUAL('Property','Value',...) creates a new NUCLEI_MANUAL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before nuclei_manual_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to nuclei_manual_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help nuclei_manual

% Last Modified by GUIDE v2.5 10-May-2011 11:35:52

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @nuclei_manual_OpeningFcn, ...
                   'gui_OutputFcn',  @nuclei_manual_OutputFcn, ...
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


% --- Executes just before nuclei_manual is made visible.
function nuclei_manual_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to nuclei_manual (see VARARGIN)

% Parameter setting:
set(gcf,'WindowButtonDownFcn',{});
global sub_list sub_num in_folder input_name out_folder mat_tail resolution0 image_folder nu_add seg_add figure_tail WGA_channel DAPI_channel WGA_th0 DAPI_th0 new_DAPI new_WGA0 fWGA bw_DAPI bw_WGA bw_WGA_applied sel_WGA sel_WGA_applied N_cycle L_ratio bound_threshold bound_threshold2 WGA_on0 DAPI_on0 applied_on0 new_WGA0_save sel_WGA_save sel_WGA_applied_save bw_WGA_save bw_WGA_applied_save bw_DAPI_save v_save fWGA_save new_DAPI_save max_image
in_folder = 'stacks/';
input_name = 'matchlist.xls';
out_folder = 'Results/';
mat_tail = '.mat';
nu_add = '_nucleus';
seg_add = '_seg';
figure_tail = '.fig';
resolution0 = 0.091;
J1 = 1;
WGA_th0 = zeros(0);
DAPI_th0 = zeros(0);
recalculate_WGA = false;

% File loading:
open_folder = uigetdir;
open_folder = [open_folder,'\'];
[sub_num, sub_list] = xlsread([open_folder,in_folder,input_name]);
[M1,~] = size(sub_list);
set(handles.list_J,'String',num2str(J1));
set(handles.all_J,'String',num2str(M1));
set(handles.folder_name,'String',open_folder);
set(handles.file_name,'String',sub_list{J1,3}(1:(length(sub_list{J1,3})-1)));
if exist([open_folder,out_folder,sub_list{J1,3}(1:(length(sub_list{J1,3})-1)),mat_tail],'file')
    load([open_folder,out_folder,sub_list{J1,3}(1:(length(sub_list{J1,3})-1)),mat_tail]);
else
    if ~exist([open_folder,out_folder],'dir')
        mkdir([open_folder,out_folder]);
    end
    image_folder = [open_folder,in_folder,sub_list{J1,3}];
    image_type = '*.tif';
    WGA_channel = sub_num(J1,6);
    DAPI_channel = sub_num(J1,7);
    protein_channel = sub_num(J1,8);
    RNA_channel = sub_num(J1,10);
    resolution = sub_num(J1,9);
    [seg_bw,cyto_bw,max_image,N_cycle,WGA_th0,DAPI_th0] = nuclei_segDAPI(image_folder,image_type,DAPI_channel,resolution);
    save([open_folder,out_folder,sub_list{J1,3}(1:(length(sub_list{J1,3})-1)),mat_tail],'seg_bw','N_cycle','cyto_bw','WGA_th0','DAPI_th0','max_image','image_folder','WGA_channel','DAPI_channel','protein_channel','RNA_channel','resolution');
end
resolution = sub_num(J1,9);
L_ratio = (resolution/resolution0);
set(handles.list_J,'UserData', L_ratio);

% Choose default command line output for nuclei_manual
handles.output = hObject;
%set(hObject,'toolbar','figure');
set(hObject,'Resize','On');
WGA_channel = DAPI_channel;

% Initialization of GUI:
if isempty(WGA_th0) || isempty(DAPI_th0)% || WGA_th0 <= 0 || DAPI_th0 <= 0 % || isempty(fWGA) || recalculate_WGA
    WGA_th0 = 0.1;
    DAPI_th0 = 0.1;
    %[~,~,max_image,N_cycle,WGA_th0,DAPI_th0] = nuclei_seg4(image_folder,image_type,WGA_channel,DAPI_channel,resolution);
    %[~,~,max_image,N_cycle,WGA_th0,DAPI_th0] = nuclei_segDAPI(image_folder,image_type,DAPI_channel,resolution);
    %close(51)
end
set(handles.WGA_default,'UserData',WGA_th0);
set(handles.DAPI_default,'UserData',DAPI_th0);
set(handles.WGA_th,'String', num2str(WGA_th0));
set(handles.DAPI_th,'String', num2str(DAPI_th0));
set(handles.WGA_on,'Value',true);
set(handles.applied_on,'Value',true);
set(handles.DAPI_on,'Value',true);
set(handles.cancel_seg,'UserData',seg_bw);
set(handles.current_seg,'Value',true);
set(handles.im_lock,'Value',false);
set(handles.im_lock,'String','Lock');
WGA_on0 = 1;
DAPI_on0 = 1;
applied_on0 = 1;

% Show image:
H = -fspecial('log',15,5);
bound_threshold = 10000;
bound_threshold2 = 50;
se = strel('disk',floor(20/L_ratio));

new_WGA0 = double(max_image(:,:,WGA_channel));
new_WGA0 = new_WGA0/max(max(new_WGA0));
new_WGA = imfilter(max_image(:,:,WGA_channel),fspecial('gaussian',10,3/L_ratio),'symmetric','conv'); %%% tophat filtering for WGA channel
new_DAPI = imtophat(max_image(:,:,DAPI_channel), se); %%% tophat filtering for DAPI channel
new_DAPI = double(new_DAPI);
new_DAPI = new_DAPI/max(max(new_DAPI));
bw_DAPI = im2bw(new_DAPI,DAPI_th0);
bw_DAPI = imfill(bw_DAPI,'holes');
bw_DAPI = imopen(bw_DAPI, strel('disk',floor(5/L_ratio)));
bw_DAPI = imerode(bw_DAPI, strel('disk',floor(3/L_ratio)));

new_WGA = double(new_WGA);
fWGA = new_WGA;
fWGA = fWGA/max(max(fWGA));
maskall = im2bw(fWGA,WGA_th0);
maskall = imopen(maskall,strel('disk',round(3/L_ratio)));
maskall = imclose(maskall,strel('disk',round(3/L_ratio)));
maskall = imopen(maskall,strel('disk',round(20/L_ratio)));
maskall = bwareaopen(maskall,round(1000/L_ratio/L_ratio));
se_mask = regionprops(maskall,'Area');
se_area = [se_mask.Area];
area1 = geomean(se_area(se_area > 1000/L_ratio/L_ratio));   %%% calculate the mean area for a single nucleus
bw_WGA = reseg(maskall,area1);

bw_WGA_applied = seg_bw;
sel_WGA = false(size(bw_WGA));
sel_WGA_applied = sel_WGA;

fWGA_save = fWGA;
new_DAPI_save = new_DAPI;
new_WGA0_save = new_WGA0;
sel_WGA_save = sel_WGA;
sel_WGA_applied_save = sel_WGA_applied;
bw_WGA_save = bw_WGA;
bw_WGA_applied_save = bw_WGA_applied;
bw_DAPI_save = bw_DAPI;

%%% Image display
overlay(:,:,1) = new_WGA0*0.8+0.3*WGA_on0*sel_WGA+applied_on0*bwperim(sel_WGA_applied);
overlay(:,:,2) = new_WGA0*0.8+0.3*WGA_on0*(bw_WGA & (~sel_WGA))+applied_on0*bwperim(bw_WGA_applied & (~sel_WGA_applied));
overlay(:,:,3) = new_WGA0*0.8+DAPI_on0*bwperim(bw_DAPI);
axes(handles.im_show);
imshow(overlay)
v_save = [1,size(fWGA_save,2),1,size(fWGA_save,1)];

%%% Add button groups:
set(handles.seg_group,'SelectionChangeFcn',@seg_group_SelectionChangeFcn);
set(gcf,'CloseRequestFcn',@th_closereq);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes nuclei_manual wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = nuclei_manual_OutputFcn(hObject, eventdata, handles) 
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
set(gcf,'WindowButtonDownFcn',{});
global sub_list sub_num in_folder input_name out_folder mat_tail resolution0 image_folder nu_add seg_add figure_tail WGA_channel DAPI_channel WGA_th0 DAPI_th0 new_DAPI new_WGA0 fWGA bw_DAPI bw_WGA bw_WGA_applied sel_WGA sel_WGA_applied N_cycle L_ratio bound_threshold bound_threshold2 WGA_on0 DAPI_on0 applied_on0 new_WGA0_save sel_WGA_save sel_WGA_applied_save bw_WGA_save bw_WGA_applied_save bw_DAPI_save v_save fWGA_save new_DAPI_save max_image
fig = gcf;
JJ = str2num(get(handles.list_J,'String'));
quit_reply = questdlg(['Move to file #',sub_list{JJ,3}(1:(length(sub_list{JJ,3})-1)),' (Make sure you have saved the result)?']);
if strcmp(quit_reply,'Yes')
    set(handles.next_file,'String','Wait');
    % Update handles structure
    guidata(hObject, handles);
    % File loading:
    M1 = str2num(get(handles.all_J,'String'));
    
    if JJ <= M1
        open_folder = get(handles.folder_name,'String');
        set(handles.list_J,'String',num2str(JJ));
        set(handles.file_name,'String',sub_list{JJ,3}(1:(length(sub_list{JJ,3})-1)));
        WGA_th0 = zeros(0);
        DAPI_th0 = zeros(0);
        if exist([open_folder,out_folder,sub_list{JJ,3}(1:(length(sub_list{JJ,3})-1)),mat_tail],'file')
            load([open_folder,out_folder,sub_list{JJ,3}(1:(length(sub_list{JJ,3})-1)),mat_tail]);
        else
            if ~exist([open_folder,out_folder],'dir')
                mkdir([open_folder,out_folder]);
            end
            image_folder = [open_folder,in_folder,sub_list{JJ,3}];
            image_type = '*.tif';
            WGA_channel = sub_num(JJ,6);
            DAPI_channel = sub_num(JJ,7);
            protein_channel = sub_num(JJ,8);
            RNA_channel = sub_num(JJ,10);
            resolution = sub_num(JJ,9);
            [seg_bw,cyto_bw,max_image,N_cycle,WGA_th0,DAPI_th0] = nuclei_segDAPI(image_folder,image_type,DAPI_channel,resolution);
            save([open_folder,out_folder,sub_list{JJ,3}(1:(length(sub_list{JJ,3})-1)),mat_tail],'seg_bw','N_cycle','cyto_bw','WGA_th0','DAPI_th0','max_image','image_folder','WGA_channel','DAPI_channel','protein_channel','RNA_channel','resolution');
        end
        resolution = sub_num(JJ,9);
        L_ratio = (resolution/resolution0);
        WGA_channel = DAPI_channel;
        
%         set(handles.ps_value,'UserData', L_ratio);

        % Initialization of GUI:
        if isempty(WGA_th0) || isempty(DAPI_th0)% || WGA_th0 <= 0 || DAPI_th0 <= 0 
            WGA_th0 = 0.1;
            DAPI_th0 = 0.1;
            %[~,~,max_image,N_cycle,WGA_th0,DAPI_th0] = nuclei_seg4(image_folder,image_type,WGA_channel,DAPI_channel,resolution);
            %[~,~,max_image,N_cycle,WGA_th0,DAPI_th0] = nuclei_segDAPI(image_folder,image_type,DAPI_channel,resolution);
            %close(51)
        end
        set(handles.WGA_default,'UserData',WGA_th0);
        set(handles.DAPI_default,'UserData',DAPI_th0);
        set(handles.WGA_th,'String', num2str(WGA_th0));
        set(handles.DAPI_th,'String', num2str(DAPI_th0));
        set(handles.WGA_on,'Value',true);
        set(handles.applied_on,'Value',true);
        set(handles.DAPI_on,'Value',true);
        set(handles.cancel_seg,'UserData',seg_bw);
        set(handles.im_lock,'Value',false);
        set(handles.im_lock,'String','Lock');
        WGA_on0 = 1;
        DAPI_on0 = 1;
        applied_on0 = 1;

        % Show image:
        H = -fspecial('log',15,5);
        
        se = strel('disk',floor(20/L_ratio));

        new_WGA0 = double(max_image(:,:,WGA_channel));
        new_WGA0 = new_WGA0/max(max(new_WGA0));
        new_WGA = imfilter(max_image(:,:,WGA_channel),fspecial('gaussian',10,3/L_ratio),'symmetric','conv'); %%% tophat filtering for WGA channel
        new_DAPI = imtophat(max_image(:,:,DAPI_channel), se); %%% tophat filtering for DAPI channel
        new_DAPI = double(new_DAPI);
        new_DAPI = new_DAPI/max(max(new_DAPI));
        bw_DAPI = im2bw(new_DAPI,DAPI_th0);
        bw_DAPI = imfill(bw_DAPI,'holes');
        bw_DAPI = imopen(bw_DAPI, strel('disk',floor(5/L_ratio)));
        bw_DAPI = imerode(bw_DAPI, strel('disk',floor(3/L_ratio)));

        new_WGA = double(new_WGA);
        fWGA = new_WGA;
        fWGA = fWGA/max(max(fWGA));
        maskall = im2bw(fWGA,WGA_th0);
        maskall = imopen(maskall,strel('disk',round(3/L_ratio)));
        maskall = imclose(maskall,strel('disk',round(3/L_ratio)));
        maskall = imopen(maskall,strel('disk',round(20/L_ratio)));
        maskall = bwareaopen(maskall,round(1000/L_ratio/L_ratio));
        se_mask = regionprops(maskall,'Area');
        se_area = [se_mask.Area];
        area1 = geomean(se_area(se_area > 1000/L_ratio/L_ratio));   %%% calculate the mean area for a single nucleus
        bw_WGA = reseg(maskall,area1);

        bw_WGA_applied = seg_bw;
        sel_WGA = false(size(bw_WGA));
        sel_WGA_applied = sel_WGA;

        fWGA_save = fWGA;
        new_DAPI_save = new_DAPI;
        new_WGA0_save = new_WGA0;
        sel_WGA_save = sel_WGA;
        sel_WGA_applied_save = sel_WGA_applied;
        bw_WGA_save = bw_WGA;
        bw_WGA_applied_save = bw_WGA_applied;
        bw_DAPI_save = bw_DAPI;

        %%% Image display
        overlay(:,:,1) = new_WGA0*0.8+0.3*WGA_on0*sel_WGA+applied_on0*bwperim(sel_WGA_applied);
        overlay(:,:,2) = new_WGA0*0.8+0.3*WGA_on0*(bw_WGA & (~sel_WGA))+applied_on0*bwperim(bw_WGA_applied & (~sel_WGA_applied));
        overlay(:,:,3) = new_WGA0*0.8+DAPI_on0*bwperim(bw_DAPI);
        axes(handles.im_show);
        imshow(overlay)
        v_save = [1,size(fWGA_save,2),1,size(fWGA_save,1)];
        %set(hObject,'toolbar','figure');

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
set(gcf,'WindowButtonDownFcn',{});
global sub_list sub_num in_folder input_name out_folder mat_tail resolution0 image_folder nu_add seg_add figure_tail WGA_channel DAPI_channel WGA_th0 DAPI_th0 new_DAPI new_WGA0 fWGA bw_DAPI bw_WGA bw_WGA_applied sel_WGA sel_WGA_applied N_cycle L_ratio bound_threshold bound_threshold2 WGA_on0 DAPI_on0 applied_on0 new_WGA0_save sel_WGA_save sel_WGA_applied_save bw_WGA_save bw_WGA_applied_save bw_DAPI_save v_save fWGA_save new_DAPI_save max_image
fig = gcf;
quit_reply = questdlg('Move to the next file (Make sure you have saved the result)?');
if strcmp(quit_reply,'Yes')
    set(handles.next_file,'String','Wait');
    % Update handles structure
    guidata(hObject, handles);
    % File loading:
    JJ = str2num(get(handles.list_J,'String'))+1;
    M1 = str2num(get(handles.all_J,'String'));
    
    if JJ <= M1
        open_folder = get(handles.folder_name,'String');
        set(handles.list_J,'String',num2str(JJ));
        set(handles.file_name,'String',sub_list{JJ,3}(1:(length(sub_list{JJ,3})-1)));
        WGA_th0 = zeros(0);
        DAPI_th0 = zeros(0);
        if exist([open_folder,out_folder,sub_list{JJ,3}(1:(length(sub_list{JJ,3})-1)),mat_tail],'file')
            load([open_folder,out_folder,sub_list{JJ,3}(1:(length(sub_list{JJ,3})-1)),mat_tail]);
        else
            if ~exist([open_folder,out_folder],'dir')
                mkdir([open_folder,out_folder]);
            end
            image_folder = [open_folder,in_folder,sub_list{JJ,3}];
            image_type = '*.tif';
            WGA_channel = sub_num(JJ,6);
            DAPI_channel = sub_num(JJ,7);
            protein_channel = sub_num(JJ,8);
            RNA_channel = sub_num(JJ,10);
            resolution = sub_num(JJ,9);
            [seg_bw,cyto_bw,max_image,N_cycle,WGA_th0,DAPI_th0] = nuclei_segDAPI(image_folder,image_type,DAPI_channel,resolution);
            save([open_folder,out_folder,sub_list{JJ,3}(1:(length(sub_list{JJ,3})-1)),mat_tail],'seg_bw','N_cycle','cyto_bw','WGA_th0','DAPI_th0','max_image','image_folder','WGA_channel','DAPI_channel','protein_channel','RNA_channel','resolution');
        end
        resolution = sub_num(JJ,9);
        L_ratio = (resolution/resolution0);
        WGA_channel = DAPI_channel;
        
        % Initialization of GUI:
        if isempty(WGA_th0) || isempty(DAPI_th0)% || WGA_th0 <= 0 || DAPI_th0 <= 0 
            WGA_th0 = 0.1;
            DAPI_th0 = 0.1;
            %[~,~,max_image,N_cycle,WGA_th0,DAPI_th0] = nuclei_seg4(image_folder,image_type,WGA_channel,DAPI_channel,resolution);
            %[~,~,max_image,N_cycle,WGA_th0,DAPI_th0] = nuclei_segDAPI(image_folder,image_type,DAPI_channel,resolution);
            %close(51)
        end
        set(handles.WGA_default,'UserData',WGA_th0);
        set(handles.DAPI_default,'UserData',DAPI_th0);
        set(handles.WGA_th,'String', num2str(WGA_th0));
        set(handles.DAPI_th,'String', num2str(DAPI_th0));
        set(handles.WGA_on,'Value',true);
        set(handles.applied_on,'Value',true);
        set(handles.DAPI_on,'Value',true);
        set(handles.cancel_seg,'UserData',seg_bw);
        set(handles.im_lock,'Value',false);
        set(handles.im_lock,'String','Lock');
        WGA_on0 = 1;
        DAPI_on0 = 1;
        applied_on0 = 1;

        % Show image:
        H = -fspecial('log',15,5);
        se = strel('disk',floor(20/L_ratio));

        new_WGA0 = double(max_image(:,:,WGA_channel));
        new_WGA0 = new_WGA0/max(max(new_WGA0));
        new_WGA = imfilter(max_image(:,:,WGA_channel),fspecial('gaussian',10,3/L_ratio),'symmetric','conv'); %%% tophat filtering for WGA channel
        new_DAPI = imtophat(max_image(:,:,DAPI_channel), se); %%% tophat filtering for DAPI channel
        new_DAPI = double(new_DAPI);
        new_DAPI = new_DAPI/max(max(new_DAPI));
        bw_DAPI = im2bw(new_DAPI,DAPI_th0);
        bw_DAPI = imfill(bw_DAPI,'holes');
        bw_DAPI = imopen(bw_DAPI, strel('disk',floor(5/L_ratio)));
        bw_DAPI = imerode(bw_DAPI, strel('disk',floor(3/L_ratio)));

        new_WGA = double(new_WGA);
        fWGA = new_WGA;
        fWGA = fWGA/max(max(fWGA));
        maskall = im2bw(fWGA,WGA_th0);
        maskall = imopen(maskall,strel('disk',round(3/L_ratio)));
        maskall = imclose(maskall,strel('disk',round(3/L_ratio)));
        maskall = imopen(maskall,strel('disk',round(20/L_ratio)));
        maskall = bwareaopen(maskall,round(1000/L_ratio/L_ratio));
        se_mask = regionprops(maskall,'Area');
        se_area = [se_mask.Area];
        area1 = geomean(se_area(se_area > 1000/L_ratio/L_ratio));   %%% calculate the mean area for a single nucleus
        bw_WGA = reseg(maskall,area1);

        bw_WGA_applied = seg_bw;
        sel_WGA = false(size(bw_WGA));
        sel_WGA_applied = sel_WGA;

        fWGA_save = fWGA;
        new_DAPI_save = new_DAPI;
        new_WGA0_save = new_WGA0;
        sel_WGA_save = sel_WGA;
        sel_WGA_applied_save = sel_WGA_applied;
        bw_WGA_save = bw_WGA;
        bw_WGA_applied_save = bw_WGA_applied;
        bw_DAPI_save = bw_DAPI;

        %%% Image display
        overlay(:,:,1) = new_WGA0*0.8+0.3*WGA_on0*sel_WGA+applied_on0*bwperim(sel_WGA_applied);
        overlay(:,:,2) = new_WGA0*0.8+0.3*WGA_on0*(bw_WGA & (~sel_WGA))+applied_on0*bwperim(bw_WGA_applied & (~sel_WGA_applied));
        overlay(:,:,3) = new_WGA0*0.8+DAPI_on0*bwperim(bw_DAPI);
        axes(handles.im_show);
        imshow(overlay)
        v_save = [1,size(fWGA_save,2),1,size(fWGA_save,1)];
        %set(hObject,'toolbar','figure');

    else
        close(gcf)
    end
    set(handles.next_file,'String','Next');
end
% Update handles structure
guidata(hObject, handles);



% --- Executes on button press in save_seg.
function save_seg_Callback(hObject, eventdata, handles)
% hObject    handle to save_seg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global sub_list sub_num in_folder input_name out_folder mat_tail resolution0 image_folder nu_add seg_add figure_tail WGA_channel DAPI_channel WGA_th0 DAPI_th0 new_DAPI new_WGA0 fWGA bw_DAPI bw_WGA bw_WGA_applied sel_WGA sel_WGA_applied N_cycle L_ratio bound_threshold bound_threshold2 WGA_on0 DAPI_on0 applied_on0 new_WGA0_save sel_WGA_save sel_WGA_applied_save bw_WGA_save bw_WGA_applied_save bw_DAPI_save v_save fWGA_save new_DAPI_save max_image
quit_reply = questdlg('Save the segmentation result?');
if strcmp(quit_reply,'Yes')
    set(handles.save_seg,'String','Wait');
    % Update handles structure
    guidata(hObject, handles);
    
    %%% Data loading
    area_ratio2 = 2;
    bw_WGA_save(v_save(3):v_save(4),v_save(1):v_save(2)) = bw_WGA;
    bw_WGA_applied_save(v_save(3):v_save(4),v_save(1):v_save(2)) = bw_WGA_applied;
    sel_WGA_save(v_save(3):v_save(4),v_save(1):v_save(2)) = sel_WGA;
    sel_WGA_applied_save(v_save(3):v_save(4),v_save(1):v_save(2)) = sel_WGA_applied;
    bw_DAPI_save(v_save(3):v_save(4),v_save(1):v_save(2)) = bw_DAPI;
    
    seg_bw = bw_WGA_applied_save;
    all_prop = regionprops(seg_bw,'Area');
    N_cycle = round(log2(size(all_prop(:),1)))+2;   %%% Calculate the nuclear cycle number
    
    %%% Cytoplasm region recognition: %%%======================================
    embryo_region = bwconvhull(seg_bw);
    cyto_bw = embryo_region & (~seg_bw);
    cyto_bw = imerode(cyto_bw, strel('disk',floor(5/L_ratio)));
%     embryo_region = false(size(seg_bw));
%     if any(any(seg_bw))
%         temp_label = 2*seg_bw;
%         temp_label(1,1) = 1;
%         embryo_prop = regionprops(temp_label,'ConvexImage','BoundingBox');
%         x1 = uint16(embryo_prop(2).BoundingBox(2));
%         x2 = uint16(x1+embryo_prop(2).BoundingBox(4)-1);
%         y1 = uint16(embryo_prop(2).BoundingBox(1));
%         y2 = uint16(y1+embryo_prop(2).BoundingBox(3)-1);
%         embryo_region(x1:x2,y1:y2) = embryo_prop(2).ConvexImage;
%     end
% 
%     DAPI_prop = regionprops(logical(bw_DAPI_save),seg_bw,'MaxIntensity','Centroid');
%     r_more = ceil(sqrt(area_ratio2*median([all_prop.Area])/pi));
%     more_bw = zeros(size(seg_bw));
%     if size(DAPI_prop(:),1)>0
%         keepIdx = find([DAPI_prop.MaxIntensity] == 0);
%         for I_more = 1:length(keepIdx)
%             more_bw(ceil(DAPI_prop(keepIdx(I_more)).Centroid(2)),ceil(DAPI_prop(keepIdx(I_more)).Centroid(1))) = 1;
%         end
%     end
%     more_bw = logical(conv2(more_bw,double(getnhood(strel('disk',r_more))),'same'));
%     cyto_bw = embryo_region & (~seg_bw) & (~more_bw);
%     cyto_bw = imerode(cyto_bw, strel('disk',floor(5/L_ratio)));
%     %bw_DAPI = imopen(bw_DAPI, strel('disk',floor(10/L_ratio)));
    %%% =======================================================================
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(51)
    bw_perim_WGA = bwperim(seg_bw);
    bw_perim_DAPI = bwperim(bw_DAPI_save);
    overlay = imoverlay(adapthisteq(new_WGA0_save), bw_perim_WGA, [0,1,0]);
    overlay = imoverlay(overlay, bw_perim_DAPI, [1,0,0]);
    imshow(overlay)
    title([image_folder,' (green: nucleus recognition, red: DAPI signal), cycle = ',num2str(N_cycle)],'Interpreter','none');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    open_folder = get(handles.folder_name,'String');
    sub_name = get(handles.file_name,'String');
    if ~exist([open_folder,out_folder],'dir')
        mkdir([open_folder,out_folder]);
    end
    saveas(51,[open_folder,out_folder,sub_name,nu_add,seg_add,figure_tail]);
    close(51)
%    if exist([open_folder,out_folder,sub_name])
        save([open_folder,out_folder,sub_name],'seg_bw','N_cycle','cyto_bw','WGA_th0','DAPI_th0','max_image','-append');
%    else
%        save([open_folder,out_folder,sub_name],'seg_bw','N_cycle','cyto_bw','WGA_th0','DAPI_th0','max_image');
%    end
    set(handles.save_seg,'String','Save');
end
% Update handles structure
guidata(hObject, handles);



% --- Executes on button press in cancel_seg.
function cancel_seg_Callback(hObject, eventdata, handles)
% hObject    handle to cancel_seg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(gcf,'WindowButtonDownFcn',{});
global sub_list sub_num in_folder input_name out_folder mat_tail resolution0 image_folder nu_add seg_add figure_tail WGA_channel DAPI_channel WGA_th0 DAPI_th0 new_DAPI new_WGA0 fWGA bw_DAPI bw_WGA bw_WGA_applied sel_WGA sel_WGA_applied N_cycle L_ratio bound_threshold bound_threshold2 WGA_on0 DAPI_on0 applied_on0 new_WGA0_save sel_WGA_save sel_WGA_applied_save bw_WGA_save bw_WGA_applied_save bw_DAPI_save v_save fWGA_save new_DAPI_save
quit_reply = questdlg('All parameters return to default values?');
if strcmp(quit_reply,'Yes')
    set(handles.cancel_seg,'String','Wait');
    % Update handles structure
    guidata(hObject, handles);
    %%% Data loading
    WGA_th0 = get(handles.WGA_default,'UserData');
    DAPI_th0 = get(handles.DAPI_default,'UserData');
    seg_bw = get(handles.cancel_seg,'UserData');

    %%% Data update
    set(handles.WGA_th,'String', num2str(WGA_th0));
    set(handles.DAPI_th,'String', num2str(DAPI_th0));
    set(handles.WGA_on,'Value',true);
    set(handles.applied_on,'Value',true);
    set(handles.DAPI_on,'Value',true);
    set(handles.current_seg,'Value',true);
    set(handles.im_lock,'Value',false);
    set(handles.im_lock,'String','Lock');
    WGA_on0 = 1;
    DAPI_on0 = 1;
    applied_on0 = 1;

    %%% Segmentation update
    bw_DAPI = im2bw(new_DAPI_save,DAPI_th0);
    bw_DAPI = imfill(bw_DAPI,'holes');
    bw_DAPI = imopen(bw_DAPI, strel('disk',floor(5/L_ratio)));
    bw_DAPI = imerode(bw_DAPI, strel('disk',floor(3/L_ratio)));
    
    maskall = im2bw(fWGA_save,WGA_th0);
    maskall = imopen(maskall,strel('disk',round(3/L_ratio)));
    maskall = imclose(maskall,strel('disk',round(3/L_ratio)));
    maskall = imopen(maskall,strel('disk',round(20/L_ratio)));
    maskall = bwareaopen(maskall,round(1000/L_ratio/L_ratio));
    se_mask = regionprops(maskall,'Area');
    se_area = [se_mask.Area];
    area1 = geomean(se_area(se_area > 1000/L_ratio/L_ratio));   %%% calculate the mean area for a single nucleus
    bw_WGA = reseg(maskall,area1);

    bw_WGA_applied = seg_bw;
    sel_WGA = false(size(bw_WGA));
    sel_WGA_applied = sel_WGA;

    fWGA = fWGA_save;
    new_DAPI = new_DAPI_save;
    new_WGA0 = new_WGA0_save;
    sel_WGA_save = sel_WGA;
    sel_WGA_applied_save = sel_WGA_applied;
    bw_WGA_save = bw_WGA;
    bw_WGA_applied_save = bw_WGA_applied;
    bw_DAPI_save = bw_DAPI;

    %%% Image display
    overlay(:,:,1) = new_WGA0*0.8+0.3*WGA_on0*sel_WGA+applied_on0*bwperim(sel_WGA_applied);
    overlay(:,:,2) = new_WGA0*0.8+0.3*WGA_on0*(bw_WGA & (~sel_WGA))+applied_on0*bwperim(bw_WGA_applied & (~sel_WGA_applied));
    overlay(:,:,3) = new_WGA0*0.8+DAPI_on0*bwperim(bw_DAPI);
    axes(handles.im_show);
    imshow(overlay)
    v_save = [1,size(fWGA_save,2),1,size(fWGA_save,1)];
    %set(hObject,'toolbar','figure');
    set(handles.cancel_seg,'String','Cancel');
end
%%% Update handles structure
guidata(hObject, handles);




function WGA_th_Callback(hObject, eventdata, handles)
% hObject    handle to WGA_th (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of WGA_th as text
%        str2double(get(hObject,'String')) returns contents of WGA_th as a double
global sub_list sub_num in_folder input_name out_folder mat_tail resolution0 image_folder nu_add seg_add figure_tail WGA_channel DAPI_channel WGA_th0 DAPI_th0 new_DAPI new_WGA0 fWGA bw_DAPI bw_WGA bw_WGA_applied sel_WGA sel_WGA_applied N_cycle L_ratio bound_threshold bound_threshold2 WGA_on0 DAPI_on0 applied_on0 new_WGA0_save sel_WGA_save sel_WGA_applied_save bw_WGA_save bw_WGA_applied_save bw_DAPI_save v_save fWGA_save new_DAPI_save

minT = get(handles.WGA_th,'Min');
maxT = get(handles.WGA_th,'Max');

%get the string for the editText component
WGA_th0 = str2num(get(handles.WGA_th,'String'));
 
%if user inputs something is not a number, or if the input is less than 0
%or greater than 100, then the slider value defaults to 0
if (isempty(WGA_th0) || WGA_th0 < minT)
    WGA_th0 = minT;
    set(handles.WGA_th,'String',num2str(WGA_th0));
elseif WGA_th0 > maxT
    WGA_th0 = maxT;
    set(handles.WGA_th,'String',num2str(WGA_th0));
end

maskall = im2bw(fWGA,WGA_th0*max(max(fWGA)));
maskall = imopen(maskall,strel('disk',round(3/L_ratio)));
maskall = imclose(maskall,strel('disk',round(3/L_ratio)));
maskall = imopen(maskall,strel('disk',round(20/L_ratio)));
maskall = bwareaopen(maskall,round(1000/L_ratio/L_ratio));
se_mask = regionprops(maskall,'Area');
se_area = [se_mask.Area];
area1 = geomean(se_area(se_area > 1000/L_ratio/L_ratio));   %%% calculate the mean area for a single nucleus
bw_WGA = reseg(maskall,area1);

overlay(:,:,1) = new_WGA0*0.8+0.3*WGA_on0*sel_WGA+applied_on0*bwperim(sel_WGA_applied);
overlay(:,:,2) = new_WGA0*0.8+0.3*WGA_on0*(bw_WGA & (~sel_WGA))+applied_on0*bwperim(bw_WGA_applied & (~sel_WGA_applied));
overlay(:,:,3) = new_WGA0*0.8+DAPI_on0*bwperim(bw_DAPI);
axes(handles.im_show);
v = axis; 
imshow(overlay)
axis(v);
% Update handles structure
%guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function WGA_th_CreateFcn(hObject, eventdata, handles)
% hObject    handle to WGA_th (see GCBO)
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
global sub_list sub_num in_folder input_name out_folder mat_tail resolution0 image_folder nu_add seg_add figure_tail WGA_channel DAPI_channel WGA_th0 DAPI_th0 new_DAPI new_WGA0 fWGA bw_DAPI bw_WGA bw_WGA_applied sel_WGA sel_WGA_applied N_cycle L_ratio bound_threshold bound_threshold2 WGA_on0 DAPI_on0 applied_on0 new_WGA0_save sel_WGA_save sel_WGA_applied_save bw_WGA_save bw_WGA_applied_save bw_DAPI_save v_save fWGA_save new_DAPI_save

minT = get(handles.DAPI_th,'Min');
maxT = get(handles.DAPI_th,'Max');

%get the string for the editText component
DAPI_th0 = str2num(get(handles.DAPI_th,'String'));
 
%if user inputs something is not a number, or if the input is less than 0
%or greater than 100, then the slider value defaults to 0
if (isempty(DAPI_th0) || DAPI_th0 < minT)
    DAPI_th0 = minT;
    set(handles.DAPI_th,'String',num2str(DAPI_th0));
elseif WGA_th0 > maxT
    DAPI_th0 = maxT;
    set(handles.DAPI_th,'String',num2str(DAPI_th0));
end

bw_DAPI = im2bw(new_DAPI,DAPI_th0*max(max(new_DAPI)));
bw_DAPI = imfill(bw_DAPI,'holes');
bw_DAPI = imopen(bw_DAPI, strel('disk',floor(5/L_ratio)));
bw_DAPI = imerode(bw_DAPI, strel('disk',floor(3/L_ratio)));

overlay(:,:,1) = new_WGA0*0.8+0.3*WGA_on0*sel_WGA+applied_on0*bwperim(sel_WGA_applied);
overlay(:,:,2) = new_WGA0*0.8+0.3*WGA_on0*(bw_WGA & (~sel_WGA))+applied_on0*bwperim(bw_WGA_applied & (~sel_WGA_applied));
overlay(:,:,3) = new_WGA0*0.8+DAPI_on0*bwperim(bw_DAPI);
axes(handles.im_show);
v = axis; 
imshow(overlay)
axis(v);
% Update handles structure
%guidata(hObject, handles);



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



% --- Executes on button press in WGA_default.
function WGA_default_Callback(hObject, eventdata, handles)
% hObject    handle to WGA_default (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global sub_list sub_num in_folder input_name out_folder mat_tail resolution0 image_folder nu_add seg_add figure_tail WGA_channel DAPI_channel WGA_th0 DAPI_th0 new_DAPI new_WGA0 fWGA bw_DAPI bw_WGA bw_WGA_applied sel_WGA sel_WGA_applied N_cycle L_ratio bound_threshold bound_threshold2 WGA_on0 DAPI_on0 applied_on0 new_WGA0_save sel_WGA_save sel_WGA_applied_save bw_WGA_save bw_WGA_applied_save bw_DAPI_save v_save fWGA_save new_DAPI_save
WGA_th0 = get(handles.WGA_default,'UserData');

% data update
set(handles.WGA_th,'String',num2str(WGA_th0));
set(handles.WGA_on,'Value',true);
set(handles.applied_on,'Value',true);
set(handles.DAPI_on,'Value',true);
set(handles.current_seg,'Value',true);
WGA_on0 = 1;
DAPI_on0 = 1;
applied_on0 = 1;

maskall = im2bw(fWGA,WGA_th0*max(max(fWGA)));
maskall = imopen(maskall,strel('disk',round(3/L_ratio)));
maskall = imclose(maskall,strel('disk',round(3/L_ratio)));
maskall = imopen(maskall,strel('disk',round(20/L_ratio)));
maskall = bwareaopen(maskall,round(1000/L_ratio/L_ratio));
se_mask = regionprops(maskall,'Area');
se_area = [se_mask.Area];
area1 = geomean(se_area(se_area > 1000/L_ratio/L_ratio));   %%% calculate the mean area for a single nucleus
bw_WGA = reseg(maskall,area1);

%update the foci recognition
overlay(:,:,1) = new_WGA0*0.8+0.3*WGA_on0*sel_WGA+applied_on0*bwperim(sel_WGA_applied);
overlay(:,:,2) = new_WGA0*0.8+0.3*WGA_on0*(bw_WGA & (~sel_WGA))+applied_on0*bwperim(bw_WGA_applied & (~sel_WGA_applied));
overlay(:,:,3) = new_WGA0*0.8+DAPI_on0*bwperim(bw_DAPI);
axes(handles.im_show);
v = axis; 
imshow(overlay)
axis(v);

% Update handles structure
guidata(hObject, handles);




% --- Executes on button press in DAPI_default.
function DAPI_default_Callback(hObject, eventdata, handles)
% hObject    handle to DAPI_default (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global sub_list sub_num in_folder input_name out_folder mat_tail resolution0 image_folder nu_add seg_add figure_tail WGA_channel DAPI_channel WGA_th0 DAPI_th0 new_DAPI new_WGA0 fWGA bw_DAPI bw_WGA bw_WGA_applied sel_WGA sel_WGA_applied N_cycle L_ratio bound_threshold bound_threshold2 WGA_on0 DAPI_on0 applied_on0 new_WGA0_save sel_WGA_save sel_WGA_applied_save bw_WGA_save bw_WGA_applied_save bw_DAPI_save v_save fWGA_save new_DAPI_save
DAPI_th0 = get(handles.DAPI_default,'UserData');

% data update
set(handles.DAPI_th,'String',num2str(DAPI_th0));
set(handles.WGA_on,'Value',true);
set(handles.applied_on,'Value',true);
set(handles.DAPI_on,'Value',true);
set(handles.current_seg,'Value',true);
WGA_on0 = 1;
DAPI_on0 = 1;
applied_on0 = 1;

bw_DAPI = im2bw(new_DAPI,DAPI_th0*max(max(new_DAPI)));
bw_DAPI = imfill(bw_DAPI,'holes');
bw_DAPI = imopen(bw_DAPI, strel('disk',floor(5/L_ratio)));
bw_DAPI = imerode(bw_DAPI, strel('disk',floor(3/L_ratio)));

%update the foci recognition
overlay(:,:,1) = new_WGA0*0.8+0.3*WGA_on0*sel_WGA+applied_on0*bwperim(sel_WGA_applied);
overlay(:,:,2) = new_WGA0*0.8+0.3*WGA_on0*(bw_WGA & (~sel_WGA))+applied_on0*bwperim(bw_WGA_applied & (~sel_WGA_applied));
overlay(:,:,3) = new_WGA0*0.8+DAPI_on0*bwperim(bw_DAPI);
axes(handles.im_show);
v = axis; 
imshow(overlay)
axis(v);

% Update handles structure
guidata(hObject, handles);



% --- Executes on button press in WGA_on.
function WGA_on_Callback(hObject, eventdata, handles)
% hObject    handle to WGA_on (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of WGA_on
global sub_list sub_num in_folder input_name out_folder mat_tail resolution0 image_folder nu_add seg_add figure_tail WGA_channel DAPI_channel WGA_th0 DAPI_th0 new_DAPI new_WGA0 fWGA bw_DAPI bw_WGA bw_WGA_applied sel_WGA sel_WGA_applied N_cycle L_ratio bound_threshold bound_threshold2 WGA_on0 DAPI_on0 applied_on0 new_WGA0_save sel_WGA_save sel_WGA_applied_save bw_WGA_save bw_WGA_applied_save bw_DAPI_save v_save fWGA_save new_DAPI_save
WGA_on0 = get(handles.WGA_on,'Value');
%update the foci recognition
overlay(:,:,1) = new_WGA0*0.8+0.3*WGA_on0*sel_WGA+applied_on0*bwperim(sel_WGA_applied);
overlay(:,:,2) = new_WGA0*0.8+0.3*WGA_on0*(bw_WGA & (~sel_WGA))+applied_on0*bwperim(bw_WGA_applied & (~sel_WGA_applied));
overlay(:,:,3) = new_WGA0*0.8+DAPI_on0*bwperim(bw_DAPI);
axes(handles.im_show);
v = axis; 
imshow(overlay)
axis(v);



% --- Executes on button press in DAPI_on.
function DAPI_on_Callback(hObject, eventdata, handles)
% hObject    handle to DAPI_on (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of DAPI_on
global sub_list sub_num in_folder input_name out_folder mat_tail resolution0 image_folder nu_add seg_add figure_tail WGA_channel DAPI_channel WGA_th0 DAPI_th0 new_DAPI new_WGA0 fWGA bw_DAPI bw_WGA bw_WGA_applied sel_WGA sel_WGA_applied N_cycle L_ratio bound_threshold bound_threshold2 WGA_on0 DAPI_on0 applied_on0 new_WGA0_save sel_WGA_save sel_WGA_applied_save bw_WGA_save bw_WGA_applied_save bw_DAPI_save v_save fWGA_save new_DAPI_save
DAPI_on0 = get(handles.DAPI_on,'Value');
%update the foci recognition
overlay(:,:,1) = new_WGA0*0.8+0.3*WGA_on0*sel_WGA+applied_on0*bwperim(sel_WGA_applied);
overlay(:,:,2) = new_WGA0*0.8+0.3*WGA_on0*(bw_WGA & (~sel_WGA))+applied_on0*bwperim(bw_WGA_applied & (~sel_WGA_applied));
overlay(:,:,3) = new_WGA0*0.8+DAPI_on0*bwperim(bw_DAPI);
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
global sub_list sub_num in_folder input_name out_folder mat_tail resolution0 image_folder nu_add seg_add figure_tail WGA_channel DAPI_channel WGA_th0 DAPI_th0 new_DAPI new_WGA0 fWGA bw_DAPI bw_WGA bw_WGA_applied sel_WGA sel_WGA_applied N_cycle L_ratio bound_threshold bound_threshold2 WGA_on0 DAPI_on0 applied_on0 new_WGA0_save sel_WGA_save sel_WGA_applied_save bw_WGA_save bw_WGA_applied_save bw_DAPI_save v_save fWGA_save new_DAPI_save
if get(handles.current_seg,'Value')
    sel_WGA = convex_seg(sel_WGA);
    bw_WGA = bw_WGA | sel_WGA;
else
    sel_WGA_applied = convex_seg(sel_WGA_applied);
    bw_WGA_applied = bw_WGA_applied | sel_WGA_applied;
end
%update the foci recognition
overlay(:,:,1) = new_WGA0*0.8+0.3*WGA_on0*sel_WGA+applied_on0*bwperim(sel_WGA_applied);
overlay(:,:,2) = new_WGA0*0.8+0.3*WGA_on0*(bw_WGA & (~sel_WGA))+applied_on0*bwperim(bw_WGA_applied & (~sel_WGA_applied));
overlay(:,:,3) = new_WGA0*0.8+DAPI_on0*bwperim(bw_DAPI);
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
global sub_list sub_num in_folder input_name out_folder mat_tail resolution0 image_folder nu_add seg_add figure_tail WGA_channel DAPI_channel WGA_th0 DAPI_th0 new_DAPI new_WGA0 fWGA bw_DAPI bw_WGA bw_WGA_applied sel_WGA sel_WGA_applied N_cycle L_ratio bound_threshold bound_threshold2 WGA_on0 DAPI_on0 applied_on0 new_WGA0_save sel_WGA_save sel_WGA_applied_save bw_WGA_save bw_WGA_applied_save bw_DAPI_save v_save fWGA_save new_DAPI_save
if get(handles.current_seg,'Value')
    %sel_WGA = imdilate(sel_WGA, strel('disk',2));
    label_WGA = bwlabel(bwmorph(bw_WGA,'thicken',1));
    sel_WGA = ismember(label_WGA,unique(label_WGA(sel_WGA)));
    bw_WGA = bw_WGA | sel_WGA;
else
    %sel_WGA_applied = imdilate(sel_WGA_applied, strel('disk',2));
    %sel_WGA_applied = bwmorph(sel_WGA_applied,'thicken',1);
    label_WGA_applied = bwlabel(bwmorph(bw_WGA_applied,'thicken',1));
    sel_WGA_applied = ismember(label_WGA_applied,unique(label_WGA_applied(sel_WGA_applied)));
    bw_WGA_applied = bw_WGA_applied | sel_WGA_applied;
end
%update the foci recognition
overlay(:,:,1) = new_WGA0*0.8+0.3*WGA_on0*sel_WGA+applied_on0*bwperim(sel_WGA_applied);
overlay(:,:,2) = new_WGA0*0.8+0.3*WGA_on0*(bw_WGA & (~sel_WGA))+applied_on0*bwperim(bw_WGA_applied & (~sel_WGA_applied));
overlay(:,:,3) = new_WGA0*0.8+DAPI_on0*bwperim(bw_DAPI);
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
global sub_list sub_num in_folder input_name out_folder mat_tail resolution0 image_folder nu_add seg_add figure_tail WGA_channel DAPI_channel WGA_th0 DAPI_th0 new_DAPI new_WGA0 fWGA bw_DAPI bw_WGA bw_WGA_applied sel_WGA sel_WGA_applied N_cycle L_ratio bound_threshold bound_threshold2 WGA_on0 DAPI_on0 applied_on0 new_WGA0_save sel_WGA_save sel_WGA_applied_save bw_WGA_save bw_WGA_applied_save bw_DAPI_save v_save fWGA_save new_DAPI_save
if get(handles.current_seg,'Value')
    bw_WGA = bw_WGA & (~sel_WGA);
    sel_WGA = imerode(sel_WGA, strel('disk',2));
    bw_WGA = bw_WGA | sel_WGA;
else
    bw_WGA_applied = bw_WGA_applied & (~sel_WGA_applied);
    sel_WGA_applied = imerode(sel_WGA_applied, strel('disk',2));
    bw_WGA_applied = bw_WGA_applied | sel_WGA_applied;
end
%update the foci recognition
overlay(:,:,1) = new_WGA0*0.8+0.3*WGA_on0*sel_WGA+applied_on0*bwperim(sel_WGA_applied);
overlay(:,:,2) = new_WGA0*0.8+0.3*WGA_on0*(bw_WGA & (~sel_WGA))+applied_on0*bwperim(bw_WGA_applied & (~sel_WGA_applied));
overlay(:,:,3) = new_WGA0*0.8+DAPI_on0*bwperim(bw_DAPI);
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
global sub_list sub_num in_folder input_name out_folder mat_tail resolution0 image_folder nu_add seg_add figure_tail WGA_channel DAPI_channel WGA_th0 DAPI_th0 new_DAPI new_WGA0 fWGA bw_DAPI bw_WGA bw_WGA_applied sel_WGA sel_WGA_applied N_cycle L_ratio bound_threshold bound_threshold2 WGA_on0 DAPI_on0 applied_on0 new_WGA0_save sel_WGA_save sel_WGA_applied_save bw_WGA_save bw_WGA_applied_save bw_DAPI_save v_save fWGA_save new_DAPI_save
applied_on0 = get(handles.applied_on,'Value');
%update the foci recognition
overlay(:,:,1) = new_WGA0*0.8+0.3*WGA_on0*sel_WGA+applied_on0*bwperim(sel_WGA_applied);
overlay(:,:,2) = new_WGA0*0.8+0.3*WGA_on0*(bw_WGA & (~sel_WGA))+applied_on0*bwperim(bw_WGA_applied & (~sel_WGA_applied));
overlay(:,:,3) = new_WGA0*0.8+DAPI_on0*bwperim(bw_DAPI);
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
global sub_list sub_num in_folder input_name out_folder mat_tail resolution0 image_folder nu_add seg_add figure_tail WGA_channel DAPI_channel WGA_th0 DAPI_th0 new_DAPI new_WGA0 fWGA bw_DAPI bw_WGA bw_WGA_applied sel_WGA sel_WGA_applied N_cycle L_ratio bound_threshold bound_threshold2 WGA_on0 DAPI_on0 applied_on0 new_WGA0_save sel_WGA_save sel_WGA_applied_save bw_WGA_save bw_WGA_applied_save bw_DAPI_save v_save fWGA_save new_DAPI_save
if get(handles.current_seg,'Value')
    %label_WGA = bwlabel(bw_WGA);
    %bw_WGA = logical(ismember(label_WGA,setdiff(unique(label_WGA),unique(label_WGA(sel_WGA)))));
    %sel_WGA = false(size(sel_WGA));
    bw_WGA(sel_WGA) = false;
    sel_WGA(sel_WGA) = false;
else
    %label_WGA_applied = bwlabel(bw_WGA_applied);
    %bw_WGA_applied = logical(ismember(label_WGA_applied,unique(label_WGA_applied(~sel_WGA_applied))));
    %sel_WGA_applied = false(size(sel_WGA_applied));
    bw_WGA_applied(sel_WGA_applied) = false;
    sel_WGA_applied(sel_WGA_applied) = false;
end
%update the foci recognition
overlay(:,:,1) = new_WGA0*0.8+0.3*WGA_on0*sel_WGA+applied_on0*bwperim(sel_WGA_applied);
overlay(:,:,2) = new_WGA0*0.8+0.3*WGA_on0*(bw_WGA & (~sel_WGA))+applied_on0*bwperim(bw_WGA_applied & (~sel_WGA_applied));
overlay(:,:,3) = new_WGA0*0.8+DAPI_on0*bwperim(bw_DAPI);
axes(handles.im_show);
v = axis; 
imshow(overlay)
axis(v);




function seg_group_SelectionChangeFcn(hObject, eventdata)
 
%retrieve GUI data, i.e. the handles structure
handles = guidata(hObject); 
global sub_list sub_num in_folder input_name out_folder mat_tail resolution0 image_folder nu_add seg_add figure_tail WGA_channel DAPI_channel WGA_th0 DAPI_th0 new_DAPI new_WGA0 fWGA bw_DAPI bw_WGA bw_WGA_applied sel_WGA sel_WGA_applied N_cycle L_ratio bound_threshold bound_threshold2 WGA_on0 DAPI_on0 applied_on0 new_WGA0_save sel_WGA_save sel_WGA_applied_save bw_WGA_save bw_WGA_applied_save bw_DAPI_save v_save fWGA_save new_DAPI_save
 
switch get(eventdata.NewValue,'Tag')   % Get Tag of selected object
    case 'sel_cls'
      %execute this code when "select/clear" button is selected
        if ~get(handles.manual_draw,'Value')
            if get(handles.sel_cls,'Value')
                axes(handles.im_show);
                zoom off
                pan off
                set(gcf,'WindowButtonDownFcn',{@ButttonDownFcn1,handles});
            else
                axes(handles.im_show);
                set(gcf,'WindowButtonDownFcn',{});
            end
        else
            set(handles.sel_cls,'Value',0);
        end
        
    case 'manual_draw'
      %execute this code when "manual draw" button is selected        
        if get(handles.manual_draw,'Value')
            set(handles.manual_draw,'ForegroundColor',[1,0,0]);
            axes(handles.im_show);
            zoom off
            pan off
            set(gcf,'WindowButtonDownFcn',{});
            add_switch = true;
            while add_switch
                overlay(:,:,1) = new_WGA0*0.8+0.3*WGA_on0*sel_WGA+applied_on0*bwperim(sel_WGA_applied);
                overlay(:,:,2) = new_WGA0*0.8+0.3*WGA_on0*(bw_WGA & (~sel_WGA))+applied_on0*bwperim(bw_WGA_applied & (~sel_WGA_applied));
                overlay(:,:,3) = new_WGA0*0.8+DAPI_on0*bwperim(bw_DAPI);
                axes(handles.im_show);
                v = axis; 
                imshow(overlay)
                axis(v);
                bw_add = roipoly;
                add_switch = any(any(bw_add));
                if add_switch
                    if get(handles.current_seg,'Value')
                        bw_WGA = bw_WGA | bw_add;
                        label_WGA = bwlabel(bw_WGA);
                        sel_WGA = sel_WGA | ismember(label_WGA,unique(label_WGA(logical(bw_add))));
                    else
                        bw_WGA_applied = bw_WGA_applied | bw_add;
                        label_WGA_applied = bwlabel(bw_WGA_applied);
                        sel_WGA_applied = sel_WGA_applied | ismember(label_WGA_applied,unique(label_WGA_applied(logical(bw_add))));
                    end
                end
                add_switch = false;
            end
            overlay(:,:,1) = new_WGA0*0.8+0.3*WGA_on0*sel_WGA+applied_on0*bwperim(sel_WGA_applied);
            overlay(:,:,2) = new_WGA0*0.8+0.3*WGA_on0*(bw_WGA & (~sel_WGA))+applied_on0*bwperim(bw_WGA_applied & (~sel_WGA_applied));
            overlay(:,:,3) = new_WGA0*0.8+DAPI_on0*bwperim(bw_DAPI);
            axes(handles.im_show);
            v = axis; 
            imshow(overlay)
            axis(v);
            set(handles.manual_draw,'Value',0);
            set(handles.manual_draw,'ForegroundColor',[0,0,0]);
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
global sub_list sub_num in_folder input_name out_folder mat_tail resolution0 image_folder nu_add seg_add figure_tail WGA_channel DAPI_channel WGA_th0 DAPI_th0 new_DAPI new_WGA0 fWGA bw_DAPI bw_WGA bw_WGA_applied sel_WGA sel_WGA_applied N_cycle L_ratio bound_threshold bound_threshold2 WGA_on0 DAPI_on0 applied_on0 new_WGA0_save sel_WGA_save sel_WGA_applied_save bw_WGA_save bw_WGA_applied_save bw_DAPI_save v_save fWGA_save new_DAPI_save
if get(handles.current_seg,'Value')
    sel_WGA = imfill(sel_WGA,'holes');
    bw_WGA = bw_WGA | sel_WGA;
else
    sel_WGA_applied = imfill(sel_WGA_applied,'holes');
    bw_WGA_applied = bw_WGA_applied | sel_WGA_applied;
end
%update the foci recognition
overlay(:,:,1) = new_WGA0*0.8+0.3*WGA_on0*sel_WGA+applied_on0*bwperim(sel_WGA_applied);
overlay(:,:,2) = new_WGA0*0.8+0.3*WGA_on0*(bw_WGA & (~sel_WGA))+applied_on0*bwperim(bw_WGA_applied & (~sel_WGA_applied));
overlay(:,:,3) = new_WGA0*0.8+DAPI_on0*bwperim(bw_DAPI);
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
global sub_list sub_num in_folder input_name out_folder mat_tail resolution0 image_folder nu_add seg_add figure_tail WGA_channel DAPI_channel WGA_th0 DAPI_th0 new_DAPI new_WGA0 fWGA bw_DAPI bw_WGA bw_WGA_applied sel_WGA sel_WGA_applied N_cycle L_ratio bound_threshold bound_threshold2 WGA_on0 DAPI_on0 applied_on0 new_WGA0_save sel_WGA_save sel_WGA_applied_save bw_WGA_save bw_WGA_applied_save bw_DAPI_save v_save fWGA_save new_DAPI_save
set(handles.im_lock,'String','Wait');
%updates the handles structure
guidata(hObject, handles);

if get(handles.im_lock,'Value')   %%% Unlock
    v_save = floor(axis-0.5)+[1,0,1,0];
%    set(hObject,'toolbar','none');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fWGA = fWGA_save(v_save(3):v_save(4),v_save(1):v_save(2));
    new_DAPI = new_DAPI_save(v_save(3):v_save(4),v_save(1):v_save(2));
    new_WGA0 = new_WGA0_save(v_save(3):v_save(4),v_save(1):v_save(2));
    
    sel_WGA_save = sel_WGA;
    sel_WGA_applied_save = sel_WGA_applied;
    bw_WGA_save = bw_WGA;
    bw_WGA_applied_save = bw_WGA_applied;
    bw_DAPI_save = bw_DAPI;
    
    sel_WGA = sel_WGA(v_save(3):v_save(4),v_save(1):v_save(2));
    sel_WGA_applied = sel_WGA_applied(v_save(3):v_save(4),v_save(1):v_save(2));
    bw_WGA = bw_WGA(v_save(3):v_save(4),v_save(1):v_save(2));
    bw_WGA_applied = bw_WGA_applied(v_save(3):v_save(4),v_save(1):v_save(2));
    bw_DAPI = bw_DAPI(v_save(3):v_save(4),v_save(1):v_save(2));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    set(handles.im_lock,'String','Unlock');
    
    overlay(:,:,1) = new_WGA0*0.8+0.3*WGA_on0*sel_WGA+applied_on0*bwperim(sel_WGA_applied);
    overlay(:,:,2) = new_WGA0*0.8+0.3*WGA_on0*(bw_WGA & (~sel_WGA))+applied_on0*bwperim(bw_WGA_applied & (~sel_WGA_applied));
    overlay(:,:,3) = new_WGA0*0.8+DAPI_on0*bwperim(bw_DAPI);
    imshow(overlay)
    
else                                  %%% Lock
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fWGA = fWGA_save;
    new_DAPI = new_DAPI_save;
    new_WGA0 = new_WGA0_save;
    
    bw_WGA_save(v_save(3):v_save(4),v_save(1):v_save(2)) = bw_WGA;
    bw_WGA_applied_save(v_save(3):v_save(4),v_save(1):v_save(2)) = bw_WGA_applied;
    sel_WGA_save(v_save(3):v_save(4),v_save(1):v_save(2)) = sel_WGA;
    sel_WGA_applied_save(v_save(3):v_save(4),v_save(1):v_save(2)) = sel_WGA_applied;
    bw_DAPI_save(v_save(3):v_save(4),v_save(1):v_save(2)) = bw_DAPI;
    
    sel_WGA = sel_WGA_save;
    sel_WGA_applied = sel_WGA_applied_save;
    bw_WGA = bw_WGA_save;
    bw_WGA_applied = bw_WGA_applied_save;
    bw_DAPI = bw_DAPI_save;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    set(hObject,'toolbar','figure');
    set(handles.im_lock,'String','Lock');
    
    overlay(:,:,1) = new_WGA0*0.8+0.3*WGA_on0*sel_WGA+applied_on0*bwperim(sel_WGA_applied);
    overlay(:,:,2) = new_WGA0*0.8+0.3*WGA_on0*(bw_WGA & (~sel_WGA))+applied_on0*bwperim(bw_WGA_applied & (~sel_WGA_applied));
    overlay(:,:,3) = new_WGA0*0.8+DAPI_on0*bwperim(bw_DAPI);
    imshow(overlay)
    axis(v_save+[-0.5,0.5,-0.5,0.5]);
    v_save = [1,size(fWGA_save,2),1,size(fWGA_save,1)];
end

%updates the handles structure
guidata(hObject, handles);

    

% --- Executes on button press in apply_seg.
function apply_seg_Callback(hObject, eventdata, handles)
% hObject    handle to apply_seg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of apply_seg
global sub_list sub_num in_folder input_name out_folder mat_tail resolution0 image_folder nu_add seg_add figure_tail WGA_channel DAPI_channel WGA_th0 DAPI_th0 new_DAPI new_WGA0 fWGA bw_DAPI bw_WGA bw_WGA_applied sel_WGA sel_WGA_applied N_cycle L_ratio bound_threshold bound_threshold2 WGA_on0 DAPI_on0 applied_on0 new_WGA0_save sel_WGA_save sel_WGA_applied_save bw_WGA_save bw_WGA_applied_save bw_DAPI_save v_save fWGA_save new_DAPI_save
if get(handles.current_seg,'Value')
    bw_WGA_applied = bw_WGA_applied | sel_WGA;
    sel_WGA_applied = sel_WGA_applied | sel_WGA;
    %update the foci recognition
    overlay(:,:,1) = new_WGA0*0.8+0.3*WGA_on0*sel_WGA+applied_on0*bwperim(sel_WGA_applied);
    overlay(:,:,2) = new_WGA0*0.8+0.3*WGA_on0*(bw_WGA & (~sel_WGA))+applied_on0*bwperim(bw_WGA_applied & (~sel_WGA_applied));
    overlay(:,:,3) = new_WGA0*0.8+DAPI_on0*bwperim(bw_DAPI);
    axes(handles.im_show);
    v = axis; 
    imshow(overlay)
    axis(v);
end



% --- Executes on button press in auto_sel.
function auto_sel_Callback(hObject, eventdata, handles)
% hObject    handle to auto_sel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
convex_threshold = 0.3;
global sub_list sub_num in_folder input_name out_folder mat_tail resolution0 image_folder nu_add seg_add figure_tail WGA_channel DAPI_channel WGA_th0 DAPI_th0 new_DAPI new_WGA0 fWGA bw_DAPI bw_WGA bw_WGA_applied sel_WGA sel_WGA_applied N_cycle L_ratio bound_threshold bound_threshold2 WGA_on0 DAPI_on0 applied_on0 new_WGA0_save sel_WGA_save sel_WGA_applied_save bw_WGA_save bw_WGA_applied_save bw_DAPI_save v_save fWGA_save new_DAPI_save
if get(handles.current_seg,'Value')
    WGA_prop = regionprops(logical(bw_WGA),bw_DAPI,'MaxIntensity','Area','ConvexArea');
    if size(WGA_prop(:),1)>0
        keepIdx = find(([WGA_prop.ConvexArea]-[WGA_prop.Area])./[WGA_prop.ConvexArea] < convex_threshold);
        keepIdx2 = find([WGA_prop.MaxIntensity] > 0);
        keepIdx3 = find(([WGA_prop.Area] > 0.25*median([WGA_prop.Area]))&([WGA_prop.Area] < 2*median([WGA_prop.Area])));
        keepIdx = intersect(keepIdx,keepIdx2);
        keepIdx = intersect(keepIdx,keepIdx3);
        WGA_all = logical(ismember(bwlabel(bw_WGA),keepIdx)); %stores objects that fit in the keepIdx
    end
    sel_WGA = WGA_all | sel_WGA;
    
    %update the foci recognition
    overlay(:,:,1) = new_WGA0*0.8+0.3*WGA_on0*sel_WGA+applied_on0*bwperim(sel_WGA_applied);
    overlay(:,:,2) = new_WGA0*0.8+0.3*WGA_on0*(bw_WGA & (~sel_WGA))+applied_on0*bwperim(bw_WGA_applied & (~sel_WGA_applied));
    overlay(:,:,3) = new_WGA0*0.8+DAPI_on0*bwperim(bw_DAPI);
    axes(handles.im_show);
    v = axis; 
    imshow(overlay)
    axis(v);
end



% --- Executes on button press in all_sel.
function all_sel_Callback(hObject, eventdata, handles)
% hObject    handle to all_sel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global sub_list sub_num in_folder input_name out_folder mat_tail resolution0 image_folder nu_add seg_add figure_tail WGA_channel DAPI_channel WGA_th0 DAPI_th0 new_DAPI new_WGA0 fWGA bw_DAPI bw_WGA bw_WGA_applied sel_WGA sel_WGA_applied N_cycle L_ratio bound_threshold bound_threshold2 WGA_on0 DAPI_on0 applied_on0 new_WGA0_save sel_WGA_save sel_WGA_applied_save bw_WGA_save bw_WGA_applied_save bw_DAPI_save v_save fWGA_save new_DAPI_save
if get(handles.current_seg,'Value')
    sel_WGA = bw_WGA | sel_WGA;
else
    sel_WGA_applied = bw_WGA_applied | sel_WGA_applied;
    
end
%update the foci recognition
overlay(:,:,1) = new_WGA0*0.8+0.3*WGA_on0*sel_WGA+applied_on0*bwperim(sel_WGA_applied);
overlay(:,:,2) = new_WGA0*0.8+0.3*WGA_on0*(bw_WGA & (~sel_WGA))+applied_on0*bwperim(bw_WGA_applied & (~sel_WGA_applied));
overlay(:,:,3) = new_WGA0*0.8+DAPI_on0*bwperim(bw_DAPI);
axes(handles.im_show);
v = axis; 
imshow(overlay)
axis(v);



% --- Executes on button press in cls_sel.
function cls_sel_Callback(hObject, eventdata, handles)
% hObject    handle to cls_sel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global sub_list sub_num in_folder input_name out_folder mat_tail resolution0 image_folder nu_add seg_add figure_tail WGA_channel DAPI_channel WGA_th0 DAPI_th0 new_DAPI new_WGA0 fWGA bw_DAPI bw_WGA bw_WGA_applied sel_WGA sel_WGA_applied N_cycle L_ratio bound_threshold bound_threshold2 WGA_on0 DAPI_on0 applied_on0 new_WGA0_save sel_WGA_save sel_WGA_applied_save bw_WGA_save bw_WGA_applied_save bw_DAPI_save v_save fWGA_save new_DAPI_save
if get(handles.current_seg,'Value')
    sel_WGA = false(size(sel_WGA));
else
    sel_WGA_applied = false(size(sel_WGA_applied));
end
%update the foci recognition
overlay(:,:,1) = new_WGA0*0.8+0.3*WGA_on0*sel_WGA+applied_on0*bwperim(sel_WGA_applied);
overlay(:,:,2) = new_WGA0*0.8+0.3*WGA_on0*(bw_WGA & (~sel_WGA))+applied_on0*bwperim(bw_WGA_applied & (~sel_WGA_applied));
overlay(:,:,3) = new_WGA0*0.8+DAPI_on0*bwperim(bw_DAPI);
axes(handles.im_show);
v = axis; 
imshow(overlay)
axis(v);



function ButttonDownFcn1(src,event,handles)
global sub_list sub_num in_folder input_name out_folder mat_tail resolution0 image_folder nu_add seg_add figure_tail WGA_channel DAPI_channel WGA_th0 DAPI_th0 new_DAPI new_WGA0 fWGA bw_DAPI bw_WGA bw_WGA_applied sel_WGA sel_WGA_applied N_cycle L_ratio bound_threshold bound_threshold2 WGA_on0 DAPI_on0 applied_on0 new_WGA0_save sel_WGA_save sel_WGA_applied_save bw_WGA_save bw_WGA_applied_save bw_DAPI_save v_save fWGA_save new_DAPI_save
pt = get(gca,'CurrentPoint');
x = floor(pt(1,1)+0.5);
y = floor(pt(1,2)+0.5);
if get(handles.current_seg,'Value')
    label_WGA = bwlabel(bw_WGA);
    label_sel = bwlabel(sel_WGA);
    if label_sel(y,x)
        sel_WGA(label_sel == label_sel(y,x)) = false;
    elseif label_WGA(y,x)
        sel_WGA(label_WGA == label_WGA(y,x)) = true;
    end
else
    label_WGA = bwlabel(bw_WGA_applied);
    label_sel = bwlabel(sel_WGA_applied);
    if label_sel(y,x)
        sel_WGA_applied(label_sel == label_sel(y,x)) = false;
    elseif label_WGA(y,x)
        sel_WGA_applied(label_WGA == label_WGA(y,x)) = true;
    end
end
%update the foci recognition
overlay(:,:,1) = new_WGA0*0.8+0.3*WGA_on0*sel_WGA+applied_on0*bwperim(sel_WGA_applied);
overlay(:,:,2) = new_WGA0*0.8+0.3*WGA_on0*(bw_WGA & (~sel_WGA))+applied_on0*bwperim(bw_WGA_applied & (~sel_WGA_applied));
overlay(:,:,3) = new_WGA0*0.8+DAPI_on0*bwperim(bw_DAPI);
axes(handles.im_show);
v = axis; 
imshow(overlay)
axis(v);


function th_closereq(hObject, eventdata, handles)
global sub_list sub_num in_folder input_name out_folder mat_tail resolution0 image_folder nu_add seg_add figure_tail WGA_channel DAPI_channel WGA_th0 DAPI_th0 new_DAPI new_WGA0 fWGA bw_DAPI bw_WGA bw_WGA_applied sel_WGA sel_WGA_applied N_cycle L_ratio bound_threshold bound_threshold2 WGA_on0 DAPI_on0 applied_on0 new_WGA0_save sel_WGA_save sel_WGA_applied_save bw_WGA_save bw_WGA_applied_save bw_DAPI_save v_save fWGA_save new_DAPI_save

quit_reply = questdlg('Quit the program (Make sure you have saved the result)?');
if strcmp(quit_reply,'Yes')
    clear global
    delete(gcf)
else
    return
end
