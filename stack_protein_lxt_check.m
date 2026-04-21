function stack_protein_lxt_check(file_list,single_switch)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to stack_RNA_lx (see VARARGIN)
global imstack immask dlayer dcontrast dcontrmin dfsize dimthresh dmask dpeak mstatus resolution0 L_ratio g vlayer vcontrast vcontrmin vfsize vimthresh vmask vpeak limsize limthresh limcontrast bk_label sub_folder bk_name sg_name stack_th stack_th0 single_switch0

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
bk_name = 'bk_record.xls';
sg_name = 'signal_record.xls';

single_switch0 = single_switch;
if single_switch0
    sub_folder = 'Histogram_protein_A/';
else
    sub_folder = 'Histogram_protein_P/';
end

if ischar(file_list)
    [~, folder_list] = xlsread(file_list);
    folder_list = folder_list(strcmpi('b',folder_list(:,6)),1);
    folder_list = folder_list([3:end]);
    hist_auto(folder_list);
else
    hist_auto(file_list);
end    

end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function hist_auto(imfolder_list)
global imstack immask dlayer dcontrast dcontrmin dfsize dimthresh dmask dpeak mstatus resolution0 L_ratio g vlayer vcontrast vcontrmin vfsize vimthresh vmask vpeak limsize limthresh limcontrast imname imfolder

list_name = 'matchlist.xls';


for I_folder = 1:length(imfolder_list)
    imfolder = imfolder_list{I_folder};
    if ~(imfolder(end) == '/' || imfolder(end) == '\')
        imfolder = [imfolder,'/'];
    end
    
    I_cut = strfind(imfolder,'stacks');
    if I_cut
        imname = imfolder(I_cut+7:end-1);
        imfolder = imfolder(1:I_cut+6);
        load_file;
        hist_on;
        hist_save;
    else           
        imfolder = [imfolder,'stacks/'];
        [~,file_list] = xlsread([imfolder,list_name]);
        if size(file_list,1) > 0
            for I_file = 1:size(file_list,1)
                imname = file_list{I_file,3};
                imname = imname(1:(end-1));
                load_file;
% %                 hist_on;
% %                 hist_save;
            end
        end
    end
    
end

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function load_file
global imstack imstack0 immask dlayer dcontrast dcontrmin dfsize dimthresh dmask dpeak mstatus resolution0 L_ratio g vlayer vcontrast vcontrmin vfsize vimthresh vmask vpeak limsize limthresh limcontrast max_spot n0 xout0 n xout imname imfolder bk_label bk_use mask_out stack_th single_switch0 mask_stack
imstack = zeros(0);
imstack0 = zeros(0);
immask = false(0);
max_spot = zeros(0);
n0 = zeros(0);
xout0 = zeros(0);
n = zeros(0);
xout = zeros(0);
bk_use = cell(0);
standard_record = 'Calibration/Results/standard.mat';
global standard_data
standard_data = load(standard_record);
list_name = 'matchlist.xls';
image_tail = '*.tif';

%% Image stack loading: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[num_list, file_list] = xlsread([imfolder,list_name]);
I_file = strcmp([imname,'/'],file_list(:,3));
I_file = I_file | strcmp([imname,'/'],file_list(:,3));

lsm_stack = dir([imfolder,imname,'/',image_tail]); %%% image file list loading
resolution = num_list(I_file,9);
RNA_channel = num_list(I_file,8);
%RNA_channel = num_list(I_file,12);
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

load([imfolder(1:find((imfolder(1:end-1) == '/') | (imfolder(1:end-1) == '\'),1,'last')),'masks/',imname,'/mask.mat']);
    mask1D = max(max(logical(mask_stack),[],3),[],1);
    ELmin = find(mask1D,1);
    ELmax = find(mask1D,1,'last');

    rxmin1 = 2/16;
    rxmax1 = 4/16;
    rymin1 = 2/8;
    rymax1 = 6/8;

    xmin1 = round(ELmin+(ELmax-ELmin)*rxmin1);
    xmax1 = round(ELmin+(ELmax-ELmin)*rxmax1);
    ymin1 = round(rymin1*size(imstack,1));
    ymax1 = round(rymax1*size(imstack,1));


    rxmin2 = 1-rxmax1;
    rxmax2 = 1-rxmin1;
    rymin2 = rymin1;
    rymax2 = rymax1;

    xmin2 = round(ELmin+(ELmax-ELmin)*rxmin2);
    xmax2 = round(ELmin+(ELmax-ELmin)*rxmax2);
    ymin2 = round(rymin2*size(imstack,1));
    ymax2 = round(rymax2*size(imstack,1));
    
    [imfolder,imname]
    mean(mean(mean(imstack0(ymin1:ymax1,xmin1:xmax1,:)))) >= mean(mean(mean(imstack0(ymin2:ymax2,xmin2:xmax2,:))))


end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function hist_on
global imstack imstack0 immask dlayer dcontrast dcontrmin dfsize dimthresh dmask dpeak mstatus resolution0 L_ratio g vlayer vcontrast vcontrmin vfsize vimthresh vmask vpeak limsize limthresh limcontrast max_spot n0 xout0 n xout sg_name bk_use imname imfolder sub_folder mask_out m_data stack_th stack_th0 RNA_mask lim t0 dt0 single_switch0 mask_stack
n = zeros(0);
xout = zeros(0);
low_th = 500;
% lim = [0:10:500];
rxy = 10;
xls_tail = '.xls';
xls_name2 = [imname,'_raw',xls_tail];
k0 = strfind(imfolder,'stacks');
imfolder0 = imfolder(1:(k0(end)-1));

    lim = [0:20:2000];
    decrease_th = 0.9^((lim(2)-lim(1))/20);
    [stack_th,t0,dt0] = th_find_low(spmask1(imstack,1),lim,decrease_th,low_th);

% % % if (stack_th > lim(end)) ||(stack_th < lim(1)) || isempty(stack_th) || isnan(stack_th)
% % %     stack_th = stack_th0;
% % % end

stack_th = 200;

% if get(handles.presult_on,'Value') && exist([imfolder0,sub_folder,xls_name2],'file')
%     [max_spot,~,~] = xlsread([imfolder0,sub_folder,xls_name2]);
%     max_spot = spfilter(max_spot,rxy,rxy,1);   %%% Fine filter to get rid of noise
% else
    [mask_out,mask_out2D] = spmask1(imstack,1,stack_th);   %%% 3D local maximal mask generation (Primary filter)
    
    embryo_mask = repmat(imdilate(bwconvhull(max(logical(mask_stack),[],3)),strel('disk',20)),[1,1,size(mask_out,3)]);
    embryo_mask = embryo_mask & (~imdilate(mask_stack,strel('disk',5)));
    mask_out = mask_out & embryo_mask;
    mask_out2D = mask_out2D & embryo_mask;
    
    peak_parameter = ptrack(imstack0,imstack,mask_out,mask_out2D);   %%% Track/fit mRNA spots
    max_spot = spfilter(peak_parameter,rxy,rxy,1);   %%% Fine filter to get rid of noise
% end
M_spot = max_spot(:,1).*max_spot(:,2).*max_spot(:,3)*2*pi;

RNA_mask = false(size(imstack,1),size(imstack,2),size(imstack,3));
for ii = 1:size(max_spot,1)
    x_spot = min(max(round(max_spot(ii,6)),1),size(imstack,1));
    y_spot = min(max(round(max_spot(ii,7)),1),size(imstack,2));
    RNA_mask(x_spot,y_spot,round(max_spot(ii,8))) = true;
end

[n0,xout0] = hist(M_spot,60);

low_th = 1.01e5;
low_bin = 1e3;
low_th2 = 1e5;
low_bin2 = 5e3;
[n,xout] = hist(M_spot,[0:low_bin:low_th]);

end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function hist_save
global imstack0 imstack imname imfolder n0 xout0 n xout sub_folder max_spot stack_th lim t0 dt0
hist_tail = '.fig';
tif_tail = '.tif';
xls_tail = '.xls';
hist_name = [imname,hist_tail];
tif_name = [imname,tif_tail];
xls_name0 = [imname,'_all',xls_tail];
xls_name1 = [imname,'_low',xls_tail];
xls_name2 = [imname,'_raw',xls_tail];
mat_name = [imname,'_value.mat'];
k0 = strfind(imfolder,'stacks');
imfolder0 = imfolder(1:(k0(end)-1));

%% Save images and data: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist([imfolder0,sub_folder],'dir')
    mkdir([imfolder0,sub_folder]);
end
% hgsave([imfolder0,sub_folder,hist_name]);
% saveas(gcf,[imfolder0,sub_folder,tif_name]);
stack_th0 = stack_th;
save([imfolder0,sub_folder,mat_name],'n0','xout0','n','xout','max_spot','stack_th0','lim','t0','dt0')
% % % % % xlswrite([imfolder0,sub_folder,xls_name0],[n0',xout0']);
% % % % % xlswrite([imfolder0,sub_folder,xls_name1],[n',xout']);
% % % % % xlswrite([imfolder0,sub_folder,xls_name2],max_spot);
% % % % % 
% % % % % %%% Foci recognition output: %%% ==========================================
% % % % % foci_add = '_foci';
% % % % % seg_add = '_seg';
% % % % % figure_tail = '.fig';
% % % % % out_folder = 'Results/';
% % % % % result_folder = [imfolder(1:find((imfolder(1:(end-1)) == '/')|(imfolder(1:(end-1)) == '\'),1,'last')),out_folder];
% % % % % foci_bw = false(size(imstack0,1),size(imstack0,2));
% % % % % 
% % % % % for ii = 1:size(max_spot,1)
% % % % %     x_spot = min(max(round(max_spot(ii,6)),1),size(imstack,1));
% % % % %     y_spot = min(max(round(max_spot(ii,7)),1),size(imstack,2));
% % % % %     foci_bw(x_spot,y_spot) = true;
% % % % % end
% % % % % 
% % % % % foci_bw0 = imdilate(foci_bw,strel('disk',4));
% % % % % bw_perim_g = bwperim(foci_bw0);
% % % % % %bw_perim_WGA = bwperim(seg_bw);
% % % % % new_image(:,:,1) = max(imstack0,[],3)./max(max(max(imstack0,[],3)));
% % % % % new_image(:,:,2) = 0;
% % % % % new_image(:,:,3) = 0;
% % % % % overlay = imoverlay(new_image, bw_perim_g, [1,1,1]);
% % % % % %overlay = imoverlay(overlay,bw_perim_WGA,[0,0,1]);
% % % % % 
% % % % % % figure(4)
% % % % % % imshow(overlay)
% % % % % % title([result_folder,imname,', white: transcription foci recognition'],'Interpreter','none')
% % % % % % saveas(4,[imfolder0,sub_folder,imname,foci_add,seg_add,figure_tail]);
% % % % % % close(4)
% % % % % % imwrite(overlay,[imfolder0,sub_folder,imname,foci_add,seg_add,figure_tail]);
% % % % % %%% =======================================================================
% % % % % 
% % % % % 
% % % % % %%% Foci recognition mask selection: %%% ==================================
% % % % % th_add = '_th';
% % % % % figure(3)
% % % % % [AX,H1,H2] = plotyy(lim,t0,lim,dt0);
% % % % % hold on
% % % % % plot(stack_th*[1,1],ylim,'--')
% % % % % title(['Threshold - # profile and the selected threshold value (',imfolder0,') threshold = ',num2str(stack_th)],'Interpreter','none')
% % % % % xlabel('Threshold value (A.U.)')
% % % % % set(get(AX(1),'Ylabel'),'String','#') 
% % % % % set(get(AX(2),'Ylabel'),'String','#/#') 
% % % % % %legend('# of recognized regions','#(n)/#(n-1)','Threshold value')
% % % % % saveas(3,[imfolder0,sub_folder,imname,foci_add,th_add,figure_tail]);
% % % % % close(3)
% % % % % %%% =======================================================================
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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

end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



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
try
    lim2 = lim(I2(find(lim(I2) >= lim1,1))+2);
catch
    lim2 = [];
end
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

end
