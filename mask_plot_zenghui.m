clear all
close all

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A program to plot the nuclear mask on top of the image %%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic
%% Parameter setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
list_name = 'data_zenghui.xls';
in_folder = 'stacks/';
input_name = 'matchlist.xls';
mismatch_name = 'mismatch.xls';
mismatch_name2 = 'mismatch_60X.xls';
xymismatch_name = 'xymismatch.mat';
xymismatch_name2 = 'xymismatch_60X.mat';
flip_list = {'Cad'};
image_type = '*.tif';
mask_folder = 'masks/';
mask_name = 'mask.mat';
out_folder = 'Results/';
% out_folder0 = 'Results_3Dz1/';
out_folder0 = '';
hist_folder = 'Histogram/';
fit_folder = 'Histogram_A/';
hist_folder2 = 'Histogram_RNA2/';
fit_folder2 = 'Histogram_A_RNA2/';
fit_add = '_spot_fit';
hist_tail = '_raw.xls';
N_thresh = 3;
output_tail = '.xls';
figure_tail = '.fig';
mat_tail = '.mat';
seg_add = '_seg';
th_add = '_th';
foci_add = '_foci';
fish_add = '_fish';
num_add = '_num';
int_add = '_int';
cmp_add = '_cmp';
sel_add = '_sel';
protein_add = '_protein';
RNA_add = '_RNA';
signal2_add = '_RNA2';
nu_add = '_nucleus';
cyto_add = '_cytoplasm';
fate_add = '_fate';
reg_add = '_regulation';
fluc_add = '_fluc';
bino_add = '_bino';
local_add = '_local';
hist_add = '_hist';
fake_add = '_fake';
abs_add = '_abs';
rad_add = '_rad';
D3_add = '_3D';
compare_add = '_compare';
DAPI_add = '_DAPI';
noise_add = '_noise';
sub_pos = [3,3];
cycle_pos0 = zeros(1,sub_pos(1));
cycle_range = [11,12,13];

sub_i = 2;
out_name = 'mask_im';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Data loading: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[num_list, folder_list] = xlsread(list_name);
folder_list = folder_list(strcmpi('T',folder_list(:,6)),:);
[N1,N2] = size(folder_list);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Image analysis: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for list_I = 1:N1
    [sub_num, sub_list] = xlsread([folder_list{list_I,1},in_folder,input_name]);
    [M1,M2] = size(sub_list);
    
    figure(list_I)
    maximize(list_I)
    sub_j = ceil(M1/sub_i);
    
    for list_J = 1:M1
        %%% Load image stack:
        image_folder = [folder_list{list_I,1},in_folder,sub_list{list_J,3}];
        temp0 = zeros(0);
        file0 = dir([image_folder,image_type]);
        for list_k = 1:length(file0)
            temp1 = imread([image_folder,file0(list_k).name]);
            temp0 = cat(3,temp0,temp1);
        end
%         max_image = max(temp0,[],3);
        
        %%% Load 3D mask:
        load([folder_list{list_I,1},mask_folder,sub_list{list_J,3},mask_name]);   
        
        %%% Plot the overlay image:
        subaxis(sub_i,sub_j,list_J, 'Spacing', 0, 'PaddingRight',0.03, 'PaddingLeft',0.05, 'PaddingTop',0.03, 'PaddingBottom',0.03, 'MarginLeft', 0.01, 'MarginRight', 0.001, 'MarginTop', 0.015, 'MarginBottom', 0.025);
        im0 = max(temp0,[],3);
        im0 = double(im0);
        im0 = im0/max(im0(:));
        im0 = imadjust(im0);
        out0 = repmat(im0,[1,1,3]);
        mask0 = logical(max(mask_stack,[],3));
        out0(:,:,2) = out0(:,:,2)+double(bwperim(mask0));
        N_nu = max(mask_stack(:));
        A_nu = nnz(mask0(:));
        A_all = numel(mask0);
        nu_cyto = A_nu/(A_all-A_nu);
        
        imshow(out0)
%         title([folder_list{list_I,1},in_folder,sub_list{list_J,3},': N_n_u = ',num2str(N_nu),', nu/cyto = ',num2str(nu_cyto)])
        title([sub_list{list_J,3},': N_n_u = ',num2str(N_nu),', nu/cyto = ',num2str(nu_cyto)])
        
    end
    
    result_folder = [folder_list{list_I,1},out_folder];
    if ~exist(result_folder)
        mkdir(result_folder)
    end
    saveas(list_I,[result_folder,out_name,figure_tail])
end
toc