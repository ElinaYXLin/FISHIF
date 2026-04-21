clear all
close all

tic
xy_mismatch_check = true;
use_linear = true; %%% whether to use linear representation for max_image00 and SS

%% Parameter setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
list_name = 'Duallist_bhh.xls';
in_folder = 'stacks/';
% old_add = '_old';
old_add = '';
input_name = 'matchlist.xls';

mismatch_name = 'mismatch.xls';
mismatch_name2 = 'mismatch_60X.xls';
% mismatch_name3 = 'mismatch_60X_FA_02142020.xls';
mismatch_name3 = 'mismatch_60X_FA.xls';
xymismatch_name = 'xymismatch.mat';
xymismatch_name2 = 'xymismatch_60X.mat';
xymismatch_name3 = 'xymismatch_60X_FA.mat';
fast_add = 'F';
airy_add = 'A';
double_add = 'D';
stitch_name = 'stitch_parameter.mat';

flip_list = {'Cad'};
% flip0 = [];
flip0 = false;
image_type = '*.tif';
mask_folder = 'masks/';
mask_name = 'mask';
out_folder = 'Results/';
% out_folder0 = 'Results_3Dz1/';
% out_folder0 = 'Results_decross/';
% out_folder0 = 'Results_original/';
% out_folder0 = 'Results_noalignment/';
out_folder0 = '';
hist_folder = 'Histogram_alignment/';
fit_folder = 'Histogram_A/';
hist_folder2 = 'Histogram_alignment_RNA2/';
fit_folder2 = 'Histogram_A_RNA2/';
fit_folder_protein = 'Histogram_protein_A/';
fit_add = '_spot_fit';
hist_tail = '_raw.xls';
N_thresh = 3;
N_thresh2 = 0;
EL_check_name = 'RNA_stack';
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
original_add = '_original';

bin0 = 1;
r_edge00 = 20;
dr = 1;
r00 = -10;
r_edge0 = [(-r_edge00-1)*dr:dr:r_edge00*dr]+r00;
expand_th = 0.3;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Data loading: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[num_list, folder_list] = xlsread(list_name);
% folder_list = folder_list(strcmpi('b',folder_list(:,6)) | strcmpi('ba',folder_list(:,6)),:);
% folder_list = folder_list(strcmpi('b3',folder_list(:,6)),:);
% folder_list = folder_list(strcmpi('h',folder_list(:,6)),:);
% folder_list = folder_list(strcmpi('b3n',folder_list(:,6)),:);
folder_list = folder_list(strcmpi('T',folder_list(:,6)),:);
% % folder_list = folder_list(strcmpi('h3k27acgal4new',folder_list(:,6)),:);
% folder_list = folder_list(strcmpi('cgal4new',folder_list(:,6)),:);
% folder_list = folder_list(strcmpi('h4k5ac',folder_list(:,6)),:);
% folder_list = folder_list(strcmpi('cbgal4',folder_list(:,6)),:);
% folder_list = folder_list(strcmpi('bgal4',folder_list(:,6)),:);
[N1,N2] = size(folder_list);
% % run_list = {2,2,1};
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Image analysis: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for list_I = 1:N1
    [sub_num, sub_list] = xlsread([folder_list{list_I,1},in_folder,input_name]);
    [M1,M2] = size(sub_list);
    
    for list_J = 2%1:M1%run_list{list_I}

        result_folder = [folder_list{list_I,1},out_folder];
        load([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mat_tail]);
        load([folder_list{list_I,1},mask_folder,sub_list{list_J,3},mask_name,mat_tail]);   %%% load 3D mask
        
        image_folder = [folder_list{list_I,1},in_folder,sub_list{list_J,3}];
        
        
        mask_stack = double(mask_stack);
        size0 = size(mask_stack);
        
        DAPI_stack = stack3D_single(image_folder,DAPI_channel); toc   %%% load 3D image stacks
% %         [nucleus_DAPI_profile,DNA_mask,nucleus_z_profile] = DAPI_profile3D(mask_stack,DAPI_stack,image_folder,N_cycle,resolution,5);

        if bin0 > 1
            DAPI_stack = imresize(DAPI_stack,1/bin0);
            mask_stack = imresize(mask_stack,1/bin0,'nearest');
        end; toc
        
        [nucleus_expand_profile,r_edge] = DAPI_profile3D_expand(mask_stack,DAPI_stack,r_edge0); toc
        clear DAPI_stack
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        


%% Adjust the size of nucleus mask: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        mask_stack = uint16(mask_adjust(mask_stack,nucleus_expand_profile,r_edge,expand_th)); toc
        
        if bin0 > 1
            mask_stack = imresize(mask_stack,bin0,'nearest');
            mask_stack = mask_stack(1:size0(1),1:size0(2),:);
%             [size0;size(mask_stack)]
        end
        
        mask_name_old = [folder_list{list_I,1},mask_folder,sub_list{list_J,3},mask_name,original_add,mat_tail];
        mask_name_new = [folder_list{list_I,1},mask_folder,sub_list{list_J,3},mask_name,mat_tail];
        if ~exist(mask_name_old)
            copyfile(mask_name_new,mask_name_old);
        end
        save(mask_name_new,'mask_stack','-v7.3')
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    end
end
toc

