clear all
close all

tic
xy_mismatch_check = true;
use_linear = true; %%% whether to use linear representation for max_image00 and SS

%% Parameter setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
list_name = 'Duallist.xls';
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
flip0 = true;
image_type = '*.tif';
mask_folder = 'masks/';
mask_name = 'mask.mat';
out_folder = 'Results2/';
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
sub_pos = [3,3];
cycle_pos0 = zeros(1,sub_pos(1));
cycle_range = [11,12,13];
% EL_range = [];
EL_range = [0,0.35,0.8,1,0.7,0.8]; tbk0 = false;  %%% Bcd
% % EL_range = [0.2,0.55,0.7,0.9,0.7,0.9]; tbk0 = false;  %%% Hb
% % EL_range = [0.3,0.7,0,0.1,0,0.1]; tbk0 = false;  %%% Cad
% % EL_range = [0.1,0.9,1,1,1,1]; tbk0 = true;  %%% Histone
sigmaxz = [1.35,1.3];
I_qp = 5; %%% type of protein quantification: 1. fluctuation method; 5. spot method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Data loading: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[num_list, folder_list] = xlsread(list_name);
% folder_list = folder_list(strcmpi('b',folder_list(:,6)) | strcmpi('ba',folder_list(:,6)),:);
% folder_list = folder_list(strcmpi('b3',folder_list(:,6)),:);
% folder_list = folder_list(strcmpi('h',folder_list(:,6)),:);
% folder_list = folder_list(strcmpi('b3n',folder_list(:,6)),:);
folder_list = folder_list(strcmpi('bcgal4new',folder_list(:,6)),:);
% folder_list = folder_list(strcmpi('cgal4new',folder_list(:,6)),:);
% folder_list = folder_list(strcmpi('h4k5ac',folder_list(:,6)),:);
% folder_list = folder_list(strcmpi('cbgal4',folder_list(:,6)),:);
% folder_list = folder_list(strcmpi('bgal4',folder_list(:,6)),:);
[N1,N2] = size(folder_list);
% % run_list = {[1,2,4,5]};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Image analysis: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for list_I = 3%1:N1
    [sub_num, sub_list] = xlsread([folder_list{list_I,1},in_folder,input_name]);
    [M1,M2] = size(sub_list);
    
    for list_J = 1%:M1%run_list{list_I}
        if isempty(strfind(sub_list{list_J,3},'_60X')) && isempty(strfind(sub_list{list_J,2},'60X')) && isempty(strfind(sub_list{list_J,2},airy_add)) && isempty(strfind(sub_list{list_J,2},fast_add))
            [~,~,mismatch_matrix] = xlsread(mismatch_name);
            load(xymismatch_name)
            mismatch_name
        elseif (~isempty(strfind(sub_list{list_J,3},'_60X')) || ~isempty(strfind(sub_list{list_J,2},'60X'))) && isempty(strfind(sub_list{list_J,2},airy_add)) && isempty(strfind(sub_list{list_J,2},fast_add))
            [~,~,mismatch_matrix] = xlsread(mismatch_name2);
            load(xymismatch_name2)
            mismatch_name2
        elseif (~isempty(strfind(sub_list{list_J,3},'_60X')) || ~isempty(strfind(sub_list{list_J,2},'60X'))) && ~isempty(strfind(sub_list{list_J,2},fast_add))
            [~,~,mismatch_matrix] = xlsread(mismatch_name3);
            load(xymismatch_name3)
            mismatch_name3
        else
            error('No mismatch information')
        end
        
%         image_folder = [folder_list{list_I,1},in_folder,sub_list{list_J,3}];
        image_folder = [folder_list{list_I,1},in_folder,sub_list{list_J,3}(1:end-1),old_add,'/'];
        
        if isempty(flip0)
            flip_axis = any(cellfun(@(x) ~isempty(strfind(image_folder,x)),flip_list));   %%% Check whether the protein profile has an opposite direction (to Bcd)
        else
            flip_axis = flip0;
        end
        
        result_folder = [folder_list{list_I,1},out_folder];
        
        load([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mat_tail]);
        image_folder = [folder_list{list_I,1},in_folder,sub_list{list_J,3}(1:end-1),old_add,'/'];
        
        if ~isempty(strfind(sub_list{list_J,2},double_add))
            Nbin = sub_num(list_J,2:3);
            Mdim = 1:2;
        else
            Nbin = ones(1,2);
            Mdim = sub_num(list_J,3);
            Nbin(Mdim) = sub_num(list_J,2);
            Mdim = 1:2;
        end
        
        %%% Stitch parameter loading: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [x_start_im,y_start_im,x_center_im,y_center_im] = tile_info_extract(image_folder,Nbin,max_image);
        %%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        load([folder_list{list_I,1},mask_folder,sub_list{list_J,3},mask_name]);   %%% load 3D mask
        mask_stack = double(mask_stack);
        
        all_color = eval(folder_list{list_I,5});
        protein_color = all_color{protein_channel};
        RNA_color = all_color{RNA_channel};
        protein_RNA_mismatch = mismatch_matrix{strcmp(protein_color,mismatch_matrix(1,:)),strcmp(RNA_color,mismatch_matrix(:,1))};
        if xy_mismatch_check
            protein_RNA_xymismatch = {eval([protein_color,'_',RNA_color]),eval([protein_color,'_',RNA_color,'_con']),eval([protein_color,'_',RNA_color,'_x0']),Nbin,Mdim,resolution,resolution_mismatch};
        else
            protein_RNA_xymismatch = [];
        end
        
% %         [signal_stack,RNA_stack,DAPI_stack] = stack3D(imclearborder(seg_bw),protein_channel,DAPI_channel,RNA_channel,image_folder,protein_RNA_mismatch);   %%% load 3D image stacks
% %         [nucleus_DAPI_profile,DNA_mask] = DAPI_profile3D(mask_stack,DAPI_stack,image_folder,N_cycle,resolution);
% %         clear DAPI_stack
        
        EL_info = get_EL(em_mask);   %%% get EL information (extreme points, EL length) from the embryo mask
% %         flip_EL = EL_orientation(nucleus_DAPI_profile,EL_info,mask_stack,eval(EL_check_name),flip_axis,[],resolution);
        flip_EL = EL_orientation2(EL_info,mask_stack,[folder_list{list_I,1},hist_folder,sub_list{list_J,3}(1:(end-1)),hist_tail],protein_RNA_mismatch,flip_axis,[]);

% %         if exist([folder_list{list_I,1},hist_folder,sub_list{list_J,3}(1:(end-1)),hist_tail],'file')
% %             [foci_list,~,~] = xlsread([folder_list{list_I,1},hist_folder,sub_list{list_J,3}(1:(end-1)),hist_tail]);
% %         else
% %             [foci_list,~,~] = xlsread([folder_list{list_I,1},hist_folder,sub_list{list_J,3}(1:(end-1)),hist_tail2]);
% %         end
% %         I_foci = prod(foci_list(:,1:3),2);
% %         
% %         x0 = EL_info(1);
% %         y0 = EL_info(2);
% %         x1 = EL_info(3);
% %         y1 = EL_info(4);
% %         L2_extreme = EL_info(5);
% %         foci_xy = foci_list(:,[7,6]);
% %         foci_EL = 1-dot((foci_xy-repmat([x0,y0],size(foci_xy,1),1)),repmat([(x1-x0),(y1-y0)],size(foci_xy,1),1),2)/L2_extreme;
% %         
% %         if flip_EL
% %             foci_EL = 1-foci_EL;
% %         end
% %         
% %         figure
% %             plot(foci_EL,I_foci,'.')
        
        save([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mat_tail],'flip_EL','-append');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        
    end
end
toc