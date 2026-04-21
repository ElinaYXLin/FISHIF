function varargout = emmask_manual(varargin)
% EMMASK_MANUAL M-file for emmask_manual.fig
%      EMMASK_MANUAL, by itself, creates a new EMMASK_MANUAL or raises the existing
%      singleton*.
%
%      H = EMMASK_MANUAL returns the handle to a new EMMASK_MANUAL or the handle to
%      the existing singleton*.
%
%      EMMASK_MANUAL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EMMASK_MANUAL.M with the given input arguments.
%
%      EMMASK_MANUAL('Property','Value',...) creates a new EMMASK_MANUAL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before emmask_manual_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to emmask_manual_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help emmask_manual

% Last Modified by GUIDE v2.5 08-Sep-2013 15:03:11

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @emmask_manual_OpeningFcn, ...
                   'gui_OutputFcn',  @emmask_manual_OutputFcn, ...
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


% --- Executes just before emmask_manual is made visible.
function emmask_manual_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to emmask_manual (see VARARGIN)

% Parameter setting:
global in_folder input_name out_folder out_folder2 mat_tail DAPI_th0
in_folder = 'stacks/';
input_name = 'matchlist.xls';
out_folder = 'Results/';
out_folder2 = 'Results_protein2/';
mat_tail = '.mat';
DAPI_th0 = 0.6;

%%% Add button groups:
set(handles.seg_group,'SelectionChangeFcn',@seg_group_SelectionChangeFcn);
set(gcf,'CloseRequestFcn',@th_closereq);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes emmask_manual wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = emmask_manual_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
% varargout{1} = handles.output;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% --- Executes on button press in load_file.
function load_file_Callback(hObject, eventdata, handles)
% hObject    handle to load_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% File loading:
global in_folder input_name out_folder mat_tail sub_list
set(handles.load_file,'String','Wait');

open_folder = uigetdir;
open_folder = [open_folder,'\'];
[~, sub_list] = xlsread([open_folder,in_folder,input_name]);
[M1,~] = size(sub_list);
J1 = 1;
set(handles.folder_name,'String',open_folder);
set(handles.list_J,'String',num2str(J1));
set(handles.all_J,'String',num2str(M1));
list_J_Callback(hObject, eventdata, handles)
set(handles.load_file,'String','Load');

% Update handles structure
guidata(hObject, handles);







function list_J_Callback(hObject, eventdata, handles)
% hObject    handle to list_J (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of list_J as text
%        str2double(get(hObject,'String')) returns contents of list_J as a double
global in_folder input_name out_folder mat_tail sub_list em_mask DAPI_all DAPI_end DAPI_th0 DAPI_th DAPI_min0 DAPI_max0 ch0

fig = gcf;
JJ = str2num(get(handles.list_J,'String'));
M1 = str2num(get(handles.all_J,'String'));

if JJ > M1
    JJ = M1;
elseif JJ < 1
    JJ = 1;
end
quit_reply = questdlg(['Move to file #',sub_list{JJ,3}(1:(length(sub_list{JJ,3})-1)),' (Make sure you have saved the result)?']);
    
if strcmp(quit_reply,'Yes')
    set(handles.next_file,'String','Wait');
    % Update handles structure
    guidata(hObject, handles);
    % File loading:

    open_folder = get(handles.folder_name,'String');
    set(handles.list_J,'String',num2str(JJ));
    set(handles.file_name,'String',sub_list{JJ,3}(1:(length(sub_list{JJ,3})-1)));

    %%% load embryo mask/image: %%% ============================================
    load([open_folder,out_folder,sub_list{JJ,3}(1:end-1),mat_tail],'em_mask','DAPI_channel','max_image');
    ch0 = DAPI_channel;
    DAPI_all = max_image(:,:,DAPI_channel);
%     ch0 = [2,5];
%     DAPI_all = mean(max_image(:,:,ch0),3);
    image_folder = [open_folder,in_folder,sub_list{JJ,3}];
    fname = dir([image_folder,'*.tif']);
    temp0 = imread([image_folder,fname(end).name]);
    DAPI_end = temp0(:,:,DAPI_channel);
    DAPI_th = DAPI_th0;
    DAPI_min0 = 0;
    DAPI_max0 = 1;
    %%% =======================================================================


    % Choose default command line output for emmask_manual
    handles.output = hObject;
    %set(hObject,'toolbar','figure');
%        set(hObject,'Resize','On');

    % Initialization of GUI:
    set(handles.DAPI_min,'String', '0');
    set(handles.DAPI_max,'String', '1');
    set(handles.DAPI_thresh,'String', num2str(DAPI_th));


    % Show image:
    overlay = im_out(DAPI_all,em_mask,DAPI_min0,DAPI_max0);
    axes(handles.im_show);
    imshow(overlay)
    
    EL_info = get_EL(em_mask);
    hold on
    plot(EL_info([1,3]),EL_info([2,4]),'g')
    hold off

    set(handles.next_file,'String','Next');
    set(handles.list_J,'String',num2str(JJ));
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
global in_folder input_name out_folder out_folder2 mat_tail sub_list em_mask DAPI_all DAPI_end DAPI_th0 DAPI_th
emmask_add = '_emmask';
figure_tail = '.fig';

quit_reply = questdlg('Save the segmentation result?');
if strcmp(quit_reply,'Yes')
    set(handles.save_seg,'String','Wait');
    % Update handles structure
    guidata(hObject, handles);
    
    %%% Data storage
    open_folder = get(handles.folder_name,'String');
    JJ = str2num(get(handles.list_J,'String'));
    save([open_folder,out_folder,sub_list{JJ,3}(1:end-1),mat_tail],'em_mask','-append','-v7.3');
    h1 = figure;
        overlay = im_out(DAPI_all,em_mask,0,1);
        imshow(overlay)
    
        EL_info = get_EL(em_mask);
        hold on
        plot(EL_info([1,3]),EL_info([2,4]),'g')
        hold off
    
        title([open_folder,in_folder,sub_list{JJ,3}],'Interpreter','none')
% %     saveas(h1,[open_folder,out_folder,sub_list{JJ,3}(1:(end-1)),emmask_add,figure_tail]);
    hgsave(h1,[open_folder,out_folder,sub_list{JJ,3}(1:(end-1)),emmask_add,figure_tail],'-v7.3');
    
    if exist([open_folder,out_folder2])
        save([open_folder,out_folder2,sub_list{JJ,3}(1:end-1),mat_tail],'em_mask','-append','-v7.3');
        saveas(h1,[open_folder,out_folder2,sub_list{JJ,3}(1:(end-1)),emmask_add,figure_tail]);
    end
    
    close(h1)
    
    set(handles.save_seg,'String','Save');
    guidata(hObject, handles);
end
% Update handles structure
guidata(hObject, handles);








function DAPI_min_Callback(hObject, eventdata, handles)
% hObject    handle to DAPI_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DAPI_min as text
%        str2double(get(hObject,'String')) returns contents of DAPI_min as a double
global in_folder input_name out_folder mat_tail sub_list em_mask DAPI_all DAPI_end DAPI_th0 DAPI_th DAPI_min0 DAPI_max0

minT = get(handles.DAPI_min,'Min');
maxT = get(handles.DAPI_min,'Max');

%get the string for the editText component
DAPI_temp = str2num(get(handles.DAPI_min,'String'));
 
%if user inputs something is not a number, or if the input is less than 0
%or greater than 100, then the slider value defaults to 0
if (isempty(DAPI_temp) || DAPI_temp < minT)
    DAPI_temp = minT;
elseif DAPI_temp > DAPI_max0
    DAPI_temp = DAPI_max0;
end
DAPI_min0 = DAPI_temp;
set(handles.DAPI_min,'String',num2str(DAPI_min0));

%%% Image display
overlay = im_out(DAPI_all,em_mask,DAPI_min0,DAPI_max0);
axes(handles.im_show);
v = axis; 
imshow(overlay)
    
    EL_info = get_EL(em_mask);
    hold on
    plot(EL_info([1,3]),EL_info([2,4]),'g')
    hold off
    
axis(v);
% Update handles structure
guidata(hObject, handles);







% --- Executes during object creation, after setting all properties.
function DAPI_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DAPI_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end







function DAPI_max_Callback(hObject, eventdata, handles)
% hObject    handle to DAPI_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DAPI_max as text
%        str2double(get(hObject,'String')) returns contents of DAPI_max as a double
global in_folder input_name out_folder mat_tail sub_list em_mask DAPI_all DAPI_end DAPI_th0 DAPI_th DAPI_min0 DAPI_max0

minT = get(handles.DAPI_max,'Min');
maxT = get(handles.DAPI_max,'Max');

%get the string for the editText component
DAPI_temp = str2num(get(handles.DAPI_max,'String'));
 
%if user inputs something is not a number, or if the input is less than 0
%or greater than 100, then the slider value defaults to 0
if (isempty(DAPI_temp) || DAPI_temp > maxT)
    DAPI_temp = maxT;
elseif DAPI_temp < DAPI_min0
    DAPI_temp = DAPI_min0;
end
DAPI_max0 = DAPI_temp;
set(handles.DAPI_max,'String',num2str(DAPI_max0));

%%% Image display
overlay = im_out(DAPI_all,em_mask,DAPI_min0,DAPI_max0);
axes(handles.im_show);
v = axis; 
imshow(overlay)
    
    EL_info = get_EL(em_mask);
    hold on
    plot(EL_info([1,3]),EL_info([2,4]),'g')
    hold off
    
axis(v);
% Update handles structure
guidata(hObject, handles);






% --- Executes during object creation, after setting all properties.
function DAPI_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DAPI_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end








function DAPI_thresh_Callback(hObject, eventdata, handles)
% hObject    handle to DAPI_thresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DAPI_thresh as text
%        str2double(get(hObject,'String')) returns contents of DAPI_thresh as a double
global in_folder input_name out_folder mat_tail sub_list em_mask DAPI_all DAPI_end DAPI_th0 DAPI_th DAPI_min0 DAPI_max0 ch0

th0 = str2num(get(handles.DAPI_thresh,'String'));
if isempty(th0) || th0 < 0
    th0 = DAPI_th0;
end
DAPI_th = th0;

JJ = str2num(get(handles.list_J,'String'));
open_folder = get(handles.folder_name,'String');

em_mask = get_emmask([open_folder,in_folder,sub_list{JJ,3}],th0,ch0);
% % em_mask = imopen(em_mask,strel('Disk',30));

%%% Image display
overlay = im_out(DAPI_all,em_mask,DAPI_min0,DAPI_max0);
axes(handles.im_show);
v = axis; 
imshow(overlay)
    
    EL_info = get_EL(em_mask);
    hold on
    plot(EL_info([1,3]),EL_info([2,4]),'g')
    hold off
    
axis(v);
% Update handles structure
guidata(hObject, handles);









% --- Executes during object creation, after setting all properties.
function DAPI_thresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DAPI_thresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end









function seg_group_SelectionChangeFcn(hObject, eventdata)
 
%retrieve GUI data, i.e. the handles structure
handles = guidata(hObject); 
global in_folder input_name out_folder mat_tail sub_list em_mask DAPI_all DAPI_end DAPI_th0 DAPI_th DAPI_min0 DAPI_max0
 
switch get(eventdata.NewValue,'Tag')   % Get Tag of selected object
    case 'manual_draw'
      %execute this code when "manual draw" button is selected        
        if get(handles.manual_draw,'Value')
            set(handles.manual_draw,'ForegroundColor',[1,0,0]);
            axes(handles.im_show);
            zoom off
            pan off
                
            bw_add = roipoly;
            if nnz(bw_add)
                em_mask = em_mask | bw_add;

                %%% Image display
                overlay = im_out(DAPI_all,em_mask,DAPI_min0,DAPI_max0);
                axes(handles.im_show);
                v = axis; 
                imshow(overlay)
    
                EL_info = get_EL(em_mask);
                hold on
                plot(EL_info([1,3]),EL_info([2,4]),'g')
                hold off
    
                axis(v);
            end
            set(handles.manual_draw,'Value',0);
            set(handles.manual_draw,'ForegroundColor',[0,0,0]);
        end
 
    case 'manual_erase'
      %execute this code when "manual erase" button is selected        
        if get(handles.manual_erase,'Value')
            set(handles.manual_erase,'ForegroundColor',[1,0,0]);
            axes(handles.im_show);
            zoom off
            pan off

            bw_cut = roipoly;
            if nnz(bw_cut)
                em_mask = em_mask & (~bw_cut);

                %%% Image display
                overlay = im_out(DAPI_all,em_mask,DAPI_min0,DAPI_max0);
                axes(handles.im_show);
                v = axis; 
                imshow(overlay)
    
                EL_info = get_EL(em_mask);
                hold on
                plot(EL_info([1,3]),EL_info([2,4]),'g')
                hold off
    
                axis(v);
            end
            set(handles.manual_erase,'Value',0);
            set(handles.manual_erase,'ForegroundColor',[0,0,0]);
        end
end
%updates the handles structure
guidata(hObject, handles);









function th_closereq(hObject, eventdata, handles)

quit_reply = questdlg('Quit the program (Make sure you have saved the result)?');
if strcmp(quit_reply,'Yes')
    clear global
    delete(gcf)
else
    return
end







function overlay = im_out(DAPI_all,em_mask,DAPI_min0,DAPI_max0)
%% A function to create output image:
DAPI_out = double(DAPI_all)/65535;
DAPI_out = (DAPI_out-DAPI_min0)/(DAPI_max0-DAPI_min0);
em_perim = imdilate(bwperim(em_mask),strel('disk',1));
overlay = zeros([size(em_mask),3]);
overlay(:,:,1) = em_perim;
overlay(:,:,2) = em_perim;
overlay(:,:,3) = em_perim+DAPI_out;





