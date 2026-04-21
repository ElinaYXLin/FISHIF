clear all
close all

tic
%% Parameter setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
list_name = 'Duallist.xls';
in_folder = 'stacks/';
% old_add = '_old';
old_add = '';
input_name = 'matchlist.xls';

image_type = '*.tif';
mask_folder = 'masks/';
mask_name = 'mask.mat';
out_folder = 'Results/';
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
mask_image_name = 'mask_image';
mask_overlap_name = 'mask_overlap';
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

sub_ij = [3,5];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Data loading: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[num_list, folder_list] = xlsread(list_name);
% folder_list = folder_list(strcmpi('b',folder_list(:,6)) | strcmpi('ba',folder_list(:,6)),:);
% folder_list = folder_list(strcmpi('b3',folder_list(:,6)),:);
% folder_list = folder_list(strcmpi('h',folder_list(:,6)),:);
% folder_list = folder_list(strcmpi('b3n',folder_list(:,6)),:);
folder_list = folder_list(strcmpi('T',folder_list(:,6)),:);
% folder_list = folder_list(strcmpi('cgal4new',folder_list(:,6)),:);
% folder_list = folder_list(strcmpi('h4k5ac',folder_list(:,6)),:);
% folder_list = folder_list(strcmpi('cbgal4',folder_list(:,6)),:);
% folder_list = folder_list(strcmpi('bgal4',folder_list(:,6)),:);
[N1,N2] = size(folder_list);
% run_list = {5};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Image analysis: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for list_I = 1%:N1
    [sub_num, sub_list] = xlsread([folder_list{list_I,1},in_folder,input_name]);
    [M1,M2] = size(sub_list);
    
    for list_J = 1:M1%run_list{list_I}
        image_folder = [folder_list{list_I,1},in_folder,sub_list{list_J,3}(1:end-1),'/'];
        result_folder = [folder_list{list_I,1},out_folder];
        
        open([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,D3_add,figure_tail]);
        open([folder_list{list_I,1},fit_folder_protein,sub_list{list_J,3}(1:end-1),fit_add,figure_tail])
        
% %         open([folder_list{list_I,1},fit_folder,sub_list{list_J,3}(1:end-1),fit_add,figure_tail])
% %         open([folder_list{list_I,1},fit_folder2,sub_list{list_J,3}(1:end-1),fit_add,figure_tail])
        
        
% %         open([folder_list{list_I,1},mask_folder,sub_list{list_J,3},mask_image_name,figure_tail])
% %         open([folder_list{list_I,1},mask_folder,sub_list{list_J,3},mask_overlap_name,figure_tail])
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        
    end
end
toc