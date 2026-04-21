clear all
close all

tic
xy_mismatch_check = true;
epsilonxy = 0.55;
epsilonxyz = 0.55;
epsilonz = 1;

%% Parameter setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
list_name = 'Duallist.xls';
in_folder = 'stacks/';
% old_add = '_old';
old_add = '';
input_name = 'matchlist.xls';
mismatch_name = 'mismatch.xls';
mismatch_name2 = 'mismatch_60X.xls';
xymismatch_name = 'xymismatch.mat';
xymismatch_name2 = 'xymismatch_60X.mat';
flip_list = {'Cad'};
flip0 = false;
image_type = '*.tif';
mask_folder = 'masks/';
mask_name = 'mask.mat';
out_folder = 'Results/';
% out_folder0 = 'Results_3Dz1/';
% out_folder0 = 'Results_decross/';
% out_folder0 = 'Results_original/';
out_folder0 = 'Results0/';
% hist_folder = 'Histogram/';
hist_folder = 'Histogram_alignment/';
fit_folder = 'Histogram_A/';
hist_folder2 = 'Histogram_en/';
fit_folder2 = 'Histogram_protein_A/';
% hist_folder_couple = 'Histogram_alignment/';
% hist_folder2_couple = 'Histogram_alignment_RNA2/';
fit_add = '_spot_fit';
hist_tail = '_raw.xls';
hist_tail2 = '_raw.xlsx';
hist_link_tail = '_link.xls';
N_thresh = 3;
output_tail = '.xls';
figure_tail = '.fig';
figure_tail2 = '.png';
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
enrich_add = '_enrich';
spot_add = '_spot';
hist_add = '_hist';

sub_pos = [4,5];
ana_folder = 'enrichspot_result/';
embryo_name = '11012014_10652-1-3M_FISHIF_lacZ_gal4_Bcd_par2';
embryo_type = {'b3n'};
% % embryo_name = '09302015_16249-1-5M_FISHIF_hb_gal4_Bcd_par2';
% % embryo_type = {'bgal4'};
var_list = {'nucleus_protein_profile','nucleus_protein_profile_ab','quanti_p','max_image','em_mask'};
dzmax = 1;
dz = -dzmax:dzmax;
inf0 = 1e8;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Data loading: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[num_list, folder_list] = xlsread(list_name);
folder_list = folder_list(strcmpi(embryo_type,folder_list(:,6)),:);
[N1,N2] = size(folder_list);
% run_list = {'[1]','[2]'};
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Index loading: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
list_ind = zeros(0);
em_index_all = 0;
for list_I = 1:N1
    [sub_num, sub_list] = xlsread([folder_list{list_I,1},in_folder,input_name]);
    [M1,M2] = size(sub_list);
    
    for list_J = 1:M1%eval(run_list{list_I})
        em_index_all = em_index_all+1;
        list_ind = cat(1,list_ind,[list_I,list_J]);
    end
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Output initialization: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xyz_all0_cell = cell(em_index_all,1);
dxy0_cell = cell(em_index_all,1);
Imin_all0_cell = cell(em_index_all,1);
I_RNA0_cell = cell(em_index_all,1);
I_en0_cell = cell(em_index_all,1);
N_en0_cell = cell(em_index_all,1);
EL_RNA0_cell = cell(em_index_all,1);
C_RNA0_cell = cell(em_index_all,1);
foci_list0_cell = cell(em_index_all,1);
Nu_en_cell = cell(em_index_all,1);

xyz_all1_cell = cell(em_index_all,1);
dxy1_cell = cell(em_index_all,1);
Imin_all1_cell = cell(em_index_all,1);
I_RNA1_cell = cell(em_index_all,1);
I_en1_cell = cell(em_index_all,1);
I_en2_cell = cell(em_index_all,1);
N_en1_cell = cell(em_index_all,1);
EL_RNA1_cell = cell(em_index_all,1);
C_RNA1_cell = cell(em_index_all,1);
foci_list1_cell = cell(em_index_all,1);
foci_circle_area_cell = cell(em_index_all,2);
Nu_RNA_cell = cell(em_index_all,1);

ind_em =zeros(em_index_all,1);
b0_em = zeros(em_index_all,1);
b2_em = zeros(em_index_all,1);
cycle_em = zeros(em_index_all,1);
em_name = cell(em_index_all,1);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Data analysis: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pool_name = parpool(10);
% matlabpool(9)
parfor em_index = 1:size(list_ind,1)
    list_I = list_ind(em_index,1);
    list_J = list_ind(em_index,2);
    [sub_num, sub_list] = xlsread([folder_list{list_I,1},in_folder,input_name]);
        
% %     if isempty(strfind(sub_list{list_J,3},'_60X'))
% %         [~,~,mismatch_matrix] = xlsread(mismatch_name);
% %         xymismatch_name0 = xymismatch_name;
% %     else
% %         [~,~,mismatch_matrix] = xlsread(mismatch_name2);
% %         xymismatch_name0 = xymismatch_name2;
% %     end

%         image_folder = [folder_list{list_I,1},in_folder,sub_list{list_J,3}];
        image_folder = [folder_list{list_I,1},in_folder,sub_list{list_J,3}(1:end-1),old_add,'/'];
        

        
% %         if exist([folder_list{list_I,1},hist_folder,sub_list{list_J,3}(1:(end-1)),hist_tail],'file')
% %             [foci_list,~,~] = xlsread([folder_list{list_I,1},hist_folder,sub_list{list_J,3}(1:(end-1)),hist_tail]);
% %         else
% %             [foci_list,~,~] = xlsread([folder_list{list_I,1},hist_folder,sub_list{list_J,3}(1:(end-1)),hist_tail2]);
% %         end
        if exist([folder_list{list_I,1},hist_folder2,sub_list{list_J,3}(1:(end-1)),hist_tail],'file')
            [foci_list2,~,~] = xlsread([folder_list{list_I,1},hist_folder2,sub_list{list_J,3}(1:(end-1)),hist_tail]);
        else
            [foci_list2,~,~] = xlsread([folder_list{list_I,1},hist_folder2,sub_list{list_J,3}(1:(end-1)),hist_tail2]);
        end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
end
delete(pool_name)
% matlabpool close
toc

