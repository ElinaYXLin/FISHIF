function varargout = stack_th(varargin)
% STACK_TH M-file for stack_th.fig
%      STACK_TH, by itself, creates a new STACK_TH or raises the existing
%      singleton*.
%
%      H = STACK_TH returns the handle to a new STACK_TH or the handle to
%      the existing singleton*.
%
%      STACK_TH('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in STACK_TH.M with the given input arguments.
%
%      STACK_TH('Property','Value',...) creates a new STACK_TH or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before stack_th_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to stack_th_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help stack_th

% Last Modified by GUIDE v2.5 20-May-2011 17:16:22

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @stack_th_OpeningFcn, ...
                   'gui_OutputFcn',  @stack_th_OutputFcn, ...
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


% --- Executes just before stack_th is made visible.
function stack_th_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to stack_th (see VARARGIN)
global imstack immask dlayer dcontrast dcontrmin dfsize dimthresh dmask dpeak mstatus resolution0 L_ratio g vlayer vcontrast vcontrmin vfsize vimthresh vmask vpeak limsize limthresh limcontrast bk_label sub_folder bk_name sg_name probe_name

%% Default value setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dlayer = 1;
dcontrast = 1;
dcontrmin = 0;
dfsize = 1.5;
dimthresh = 0.01;
dmask = true;
dpeak = true;
bk_label = 'OreR';
mstatus = {'Mask on','Mask off'};
resolution0 = 0.13;
limsize = [0,10];
limthresh = [0,1];
limcontrast = [0,1];
sub_folder = 'Histogram/';
bk_name = 'bk_record.xls';
sg_name = 'signal_record.xls';
probe_name = '_hb';
set(handles.show_bk,'Value',false);
set(handles.presult_on,'Value',true);
set(handles.fit_on,'Value',true);
% set(handles.sigma_show,'Color','none');
% set(handles.sigma_zoom,'Color','none');
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Choose default command line output for stack_th
handles.output = hObject;
set(hObject,'Resize','On');
%load_file_Callback(hObject, eventdata, handles)

set(gcf,'WindowButtonDownFcn',{@ButttonDownFcn2});
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes stack_th wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = stack_th_OutputFcn(hObject, eventdata, handles) 
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

overlay = cf_show1(imstack,g,vlayer,vimthresh,vpeak,vcontrast,vcontrmin,vmask);
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

overlay = cf_show1(imstack,g,vlayer,vimthresh,vpeak,vcontrast,vcontrmin,vmask);
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

[overlay,g] = cf_show(imstack,vlayer,vfsize,L_ratio,vimthresh,vpeak,vcontrast,vcontrmin,vmask);
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

overlay = cf_show1(imstack,g,vlayer,vimthresh,vpeak,vcontrast,vcontrmin,vmask);
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
global imstack imstack0 immask dlayer dcontrast dcontrmin dfsize dimthresh dmask dpeak mstatus resolution0 L_ratio g vlayer vcontrast vcontrmin vfsize vimthresh vmask vpeak limsize limthresh limcontrast max_spot n0 xout0 n xout imname imfolder bk_label bk_use mask_out
imstack = zeros(0);
imstack0 = zeros(0);
immask = false(0);
max_spot = zeros(0);
n0 = zeros(0);
xout0 = zeros(0);
n = zeros(0);
xout = zeros(0);
bk_use = cell(0);
Idim = 2;
standard_record = 'Calibration/Results/standard.mat';
global standard_data
standard_data = load(standard_record);

%% Image stack loading: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(varargin)
    [imname,imfolder,~] = uigetfile('*.*','Select a file to load');   %%% Pick up a lsm file
end

set(handles.load_file,'String','Busy');
% Update handles structure
guidata(hObject, handles);

lsm_stack = tiffread([imfolder,imname]); %%% Lsm file loading
resolution = lsm_stack(1).lsm.VoxelSizeX/(1e-6);
L_ratio = (resolution/resolution0);

immax = length(lsm_stack);
tempmax = zeros(1,immax);
if isa(lsm_stack(1).data,'cell')
    %%%%% Intensity correction mask generation: ---------------------------
    minsize = min(size(lsm_stack(1).data{1}));
    corr_size = [minsize*[1,1],1];
    imcorr = corr_mask(corr_size,{'TMR'},resolution);
    %%%%% -----------------------------------------------------------------
    for I_layer = 1:immax
        temp0 = lsm_stack(I_layer).data;
        temp = temp0{Idim};
        temp = uint16(double(temp)./repmat(imcorr,size(temp)./minsize));
        imstack0(:,:,I_layer) = temp;
        imstack(:,:,I_layer) = imfilter(temp,fspecial('gaussian',5,1),'symmetric','conv');
        tempmax(I_layer) = max(max(imstack(:,:,I_layer)));
    end
else
    %%%%% Intensity correction mask generation: ---------------------------
    minsize = min(size(lsm_stack(1).data));
    corr_size = [minsize*[1,1],1];
    imcorr = corr_mask(corr_size,{'A488'},resolution);
    %%%%% -----------------------------------------------------------------
    for I_layer = 1:immax
        temp = lsm_stack(I_layer).data;
        temp = uint16(double(temp)./repmat(imcorr,size(temp)./minsize));
        imstack0(:,:,I_layer) = temp;
        imstack(:,:,I_layer) = imfilter(temp,fspecial('gaussian',5,1),'symmetric','conv');
        tempmax(I_layer) = max(max(imstack(:,:,I_layer)));
    end
end
temp0 = fspecial('gaussian',5,0.3);
kernel0 = zeros(1,1,5);
kernel0(1,1,:) = temp0(3,:);
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
mask_out = spmask(imstack,1);   %%% 3D local maximal mask generation (Primary filter)
[overlay,g] = cf_show(imstack,vlayer,vfsize,L_ratio,vimthresh,vpeak,vcontrast,vcontrmin,vmask);
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

overlay = cf_show1(imstack,g,vlayer,vimthresh,vpeak,vcontrast,vcontrmin,vmask);
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
global imstack0 imstack immask dlayer dcontrast dcontrmin dfsize dimthresh dmask dpeak mstatus resolution0 L_ratio g vlayer vcontrast vcontrmin vfsize vimthresh vmask vpeak limsize limthresh limcontrast
minT = limcontrast(1);
maxT = limcontrast(2);
vcontrast = str2num(get(handles.contrast_value,'String'));
if (isempty(vcontrast) || vcontrast < minT)
    vcontrast = minT;
elseif vcontrast > maxT
    vcontrast = maxT;
end
set(handles.contrast_value,'String',num2str(vcontrast));


%imstack1 = imstack0./max(max(max(imstack0)));

overlay = cf_show1(imstack,g,vlayer,vimthresh,vpeak,vcontrast,vcontrmin,vmask);
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

overlay = cf_show1(imstack,g,vlayer,vimthresh,vpeak,vcontrast,vcontrmin,vmask);
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
global imstack imstack0 immask dlayer dcontrast dcontrmin dfsize dimthresh dmask dpeak mstatus resolution0 L_ratio g vlayer vcontrast vcontrmin vfsize vimthresh vmask vpeak limsize limthresh limcontrast max_spot n0 xout0 n xout sg_name bk_use imname imfolder sub_folder mask_out m_data probe_name
n = zeros(0);
xout = zeros(0);
set(handles.hist_on,'String','Wait')
guidata(hObject, handles);

xls_tail = '.xls';
xls_name2 = [imname(1:find(imname == '.',1,'last')-1),'_raw',xls_tail];

if get(handles.presult_on,'Value') && exist([imfolder,sub_folder,xls_name2],'file')
    [max_spot,~,~] = xlsread([imfolder,sub_folder,xls_name2]);
else
    % for I_layer = 1:get(handles.layer_slide,'Max')
    %     immask(:,:,I_layer) = im2bw(g(:,:,I_layer),vimthresh);
    % end
    % STATS = regionprops(immask,imstack*65535,'MaxIntensity');
    % max_spot = imstack(imregionalmax(g)&immask)*65535;
    mask_out = spmask(imstack,1);   %%% 3D local maximal mask generation (Primary filter)
    mask_out = mask_out;% & (imstack*65535 >= 5000);   %%% Prefilter using peak intensity threshold
    % max_spot = imstack(mask_out)*65535;
    peak_parameter = ptrack(imstack0,imstack,mask_out);   %%% Track/fit mRNA spots
    max_spot = spfilter(peak_parameter);   %%% Fine filter to get rid of noise
end
M_spot = max_spot(:,1).*max_spot(:,2).*max_spot(:,3)*2*pi;

%% Background calculation: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if get(handles.show_bk,'Value') && ~get(handles.background_on,'Value')
    if isempty(bk_use)
        if exist([imfolder,sg_name],'file')
            [~,sg_record,~] = xlsread([imfolder,sg_name]);
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
            [bk_temp,~,~] = xlsread([imfolder,sub_folder,bk_use{I_bk}]);
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
%set(gca,'YScale','log')
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
low_th = 102000;
low_bin = 500;
low_th2 = 100000;
low_bin2 = 5000;

%%%%% Figure plot:
[n,xout] = hist(M_spot,[0:low_bin:low_th]);
%[n,xout] = hist(M_spot(:,1),[0:1000:60000]);
%n = n/sum(n)*100;
n = n(1:(end-1));
xout = xout(1:(end-1));
bar(xout,n)
legend_text = {'Raw distribution'};

%%%%% Statistics:
N_start = strfind(imname,probe_name)+length(probe_name);
N_end = strfind(imname,'_');
if isempty(N_start)
    N_site = 48;
else
    N_site = str2num(imname(N_start(1):(N_end(find(N_end > N_start(1),1))-1)));
end
m_temp = M_spot(M_spot(:,1)<low_th,1);
m_mean = mean(m_temp);
m_std = std(m_temp);
m_median = median(m_temp);
m_peak = xout(find(n == max(n),1));
m_p =  1./(N_site/(m_mean/m_std)^2+1);
m_data = [m_mean,m_std,m_median,m_peak,m_p];


if get(handles.show_bk,'Value') && ~isempty(n_bk)
    [~,Imax] = max(n_bk);
    %[~,Imin] = min(abs(xout-xout_bk(Imax)));
    [~,Imin] = max(n);
    nlevel = n(Imin);
    hold on
    plot(xout_bk,n_bk/max(n_bk)*nlevel,'r')
    hold off
    legend_text = cat(2,legend_text,{'Background'});
end

if get(handles.fit_on,'Value')
    Nx = [0:N_site];
    py = binopdf(Nx,N_site,m_p);
    Aprobe = m_mean./m_p./N_site;
    Ay = max(n)./max(py);
    hold on
    plot(Aprobe*Nx,Ay*py,'g')
    hold off
    legend_text = cat(2,legend_text,{'Binomial fit'});
end

legend(legend_text)

xlim([0,low_th2]);
%xlim([0,10000]);
set(gca,'XTick',[0:low_bin2:low_th2])%,'YScale','log')
%set(gca,'XTick',[0:5000:50000],'YScale','log')

ylabel('#')
ylim0 = ylim;
ylim0(1) = 1;
ylim(ylim0);
xlim0 = xlim;

r_bin = 0.1;
r_max = 3+r_bin;
[nrx,rxout] = hist(max_spot(:,2),[0:r_bin:r_max]);
nrx = nrx(1:(end-1));
rxout = rxout(1:(end-1));
rx_peak = rxout(find(nrx == max(nrx),1));

[nry,ryout] = hist(max_spot(:,3),[0:r_bin:r_max]);
nry = nry(1:(end-1));
ryout = ryout(1:(end-1));
ry_peak = ryout(find(nry == max(nry),1));


set(handles.text_input,'String',['Low intensity zoom in (mean = ',num2str(m_mean),', std = ',num2str(m_std),', median = ',num2str(m_median),', peak = ',num2str(m_peak),', p0 = ',num2str(m_p),', r_x = ',num2str(rx_peak),', r_y = ',num2str(ry_peak)])

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

overlay = cf_show1(imstack,g,vlayer,vimthresh,vpeak,vcontrast,vcontrmin,vmask);
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
global imstack immask dlayer dcontrast dcontrmin dfsize dimthresh dmask dpeak mstatus resolution0 L_ratio g vlayer vcontrast vcontrmin vfsize vimthresh vmask vpeak limsize limthresh limcontrast imname imfolder n0 xout0 n xout bk_name sg_name bk_use sub_folder max_spot m_data
hist_tail = '.fig';
tif_tail = '.tif';
xls_tail = '.xls';
hist_name = [imname(1:find(imname == '.',1,'last')-1),hist_tail];
tif_name = [imname(1:find(imname == '.',1,'last')-1),tif_tail];
xls_name0 = [imname(1:find(imname == '.',1,'last')-1),'_all',xls_tail];
xls_name1 = [imname(1:find(imname == '.',1,'last')-1),'_low',xls_tail];
xls_name2 = [imname(1:find(imname == '.',1,'last')-1),'_raw',xls_tail];

%% Save images and data: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist([imfolder,sub_folder],'dir')
    mkdir([imfolder,sub_folder]);
end
hgsave([imfolder,sub_folder,hist_name]);
saveas(gcf,[imfolder,sub_folder,tif_name]);
if ~isempty(n0)
    xlswrite([imfolder,sub_folder,xls_name0],[n0',xout0']);
    xlswrite([imfolder,sub_folder,xls_name1],[n',xout']);
    xlswrite([imfolder,sub_folder,xls_name2],max_spot);
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Save the record: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if get(handles.background_on,'Value')
        %%% Save background record: %%% =======================================
        if exist([imfolder,bk_name],'file')
            [bk_data,bk_record,~] = xlsread([imfolder,bk_name]);
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
        xlswrite([imfolder,bk_name],cat(2,bk_record,num2cell(bk_data)));
        %%% ===================================================================
    else
        %%% Save data record: %%% =============================================
        if exist([imfolder,sg_name],'file')
            [sg_data,sg_record,~] = xlsread([imfolder,sg_name]);
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
        xlswrite([imfolder,sg_name],cat(2,sg_record,num2cell(sg_data)));
        %%% ===================================================================
    end
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

% Update handles structure
guidata(hObject, handles);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% --- Executes on button press in hist_auto.
function hist_auto_Callback(hObject, eventdata, handles)
% hObject    handle to hist_auto (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global imstack immask dlayer dcontrast dcontrmin dfsize dimthresh dmask dpeak mstatus resolution0 L_ratio g vlayer vcontrast vcontrmin vfsize vimthresh vmask vpeak limsize limthresh limcontrast imname imfolder
lsm_type = '*.lsm';

set(handles.hist_auto,'String','Busy');
% Update handles structure
guidata(hObject, handles);

imfolder = uigetdir;
imfolder = [imfolder,'/'];
lsm_name = dir([imfolder,lsm_type]);

if length(lsm_name) > 0
    for I_file = 1:length(lsm_name)
        imname = lsm_name(I_file).name;
        load_file_Callback(hObject, eventdata, handles, true);
        hist_on_Callback(hObject, eventdata, handles);
        hist_save_Callback(hObject, eventdata, handles);
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
if exist([imfolder,bk_name],'file')
    [~,bk_record,~] = xlsread([imfolder,bk_name]);
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







function ButttonDownFcn2(src,event)
global imstack immask dlayer dcontrast dcontrmin dfsize dimthresh dmask dpeak mstatus resolution0 L_ratio g vlayer vcontrast vcontrmin vfsize vimthresh vmask vpeak limsize limthresh limcontrast
pt = get(gca,'CurrentPoint');
x = floor(pt(1,1)+0.5);
y = floor(pt(1,2)+0.5);
xy_radius = 10;
foci_sub = imstack((y-xy_radius):(y+xy_radius),(x-xy_radius):(x+xy_radius),vlayer)*65535;
%[X,Y] = meshgrid([-xy_radius:xy_radius],[-xy_radius:xy_radius]);
figure
subplot(1,2,1)
    surf(foci_sub);
subplot(1,2,2)
    plot(reshape(imstack(y,x,:),1,size(imstack,3))*65535,'ob')
    hold on
    plot(vlayer*[1,1],imstack(y,x,vlayer)*65535*[0,1],'--r')
    hold off
    


% --- Executes on button press in fit_on.
function fit_on_Callback(hObject, eventdata, handles)
% hObject    handle to fit_on (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of fit_on
