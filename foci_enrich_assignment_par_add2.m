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

sub_pos = [5,7];
ana_folder = 'enrichspot_result/';
% % embryo_name = '11012014_10652-1-3M_FISHIF_lacZ_gal4_Bcd';
% % embryo_type = {'b3n'};
% % % embryo_name = '09302015_16249-1-5M_FISHIF_hb_gal4_Bcd_par2';
% % % embryo_type = {'bgal4'};
% embryo_name = '12072016b_1_16249-1-5M_FISHIF_hb_gal4_H3_par3';
% embryo_type = {'h3'};
% % embryo_name = '12072016b_2_16249-1-5M_FISHIF_hb_gal4_H4K5ac_par3';
% % embryo_type = {'h4k5ac'};
embryo_name = '12072016b_3_16249-1-5M_FISHIF_hb_gal4_H3K27ac_par3';
embryo_type = {'h3k27ac'};
var_list = {'nucleus_protein_profile','nucleus_protein_profile_ab','quanti_p','max_image','em_mask'};
dzmax = 1;
dz = -dzmax:dzmax;
inf0 = 1e8;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Data loading: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([ana_folder,embryo_name,'/',embryo_name,mat_tail]);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Data analysis: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for em_index = 1:size(foci_circle_area_cell,1)
    dx = 0.02;
    xl0 = 0:dx:(3+dx);
    foci_circle_area_cell{em_index,2} = xl0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
end
toc


%% Output: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if exist([ana_folder,embryo_name,'/']) ~= 7
    mkdir([ana_folder,embryo_name,'/']);
end
save([ana_folder,embryo_name,'/',embryo_name,mat_tail],'foci_circle_area_cell','-append');

