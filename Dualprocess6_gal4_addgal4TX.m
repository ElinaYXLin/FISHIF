clear all
close all

tic
%% Parameter setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
list_name = 'Duallist.xls';
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
hist_folder2 = 'Histogram_gal4/';
fit_folder2 = 'Histogram_A_gal4/';
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
signal2_add = '_gal4';
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Data loading: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[num_list, folder_list] = xlsread(list_name);
% folder_list = folder_list(strcmpi('b',folder_list(:,6)) | strcmpi('ba',folder_list(:,6)),:);
% folder_list = folder_list(strcmpi('b3',folder_list(:,6)),:);
% folder_list = folder_list(strcmpi('h',folder_list(:,6)),:);
% folder_list = folder_list(strcmpi('b3h',folder_list(:,6)),:);
folder_list = folder_list(strcmpi('T',folder_list(:,6)),:);
[N1,N2] = size(folder_list);
% run_list = {'[1]','[2]'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Image analysis: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for list_I = 1:N1
    [sub_num, sub_list] = xlsread([folder_list{list_I,1},in_folder,input_name]);
    [M1,M2] = size(sub_list);
    
    for list_J = 1:M1%eval(run_list{list_I})
        result_folder = [folder_list{list_I,1},out_folder];

        load([folder_list{list_I,1},fit_folder,sub_list{list_J,3}(1:end-1),fit_add,mat_tail]);   %%% load spot intensity fitting result
        single_Inten = b;
        load([folder_list{list_I,1},fit_folder2,sub_list{list_J,3}(1:end-1),fit_add,mat_tail]);   %%% load spot intensity fitting result
        single_Inten2 = b;
        
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        nucleus_signal2_profile = xlsread([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),signal2_add,nu_add,output_tail]);
        nucleus_signal2_profile(:,4) = ceil(nucleus_signal2_profile(:,4)*single_Inten/single_Inten2);
        
        cytoplasmic_signal2_profile = xlsread([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),signal2_add,cyto_add,output_tail]);
        
        foci_signal2_profile = xlsread([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),signal2_add,foci_add,output_tail]);
        foci_signal2_profile(:,3:4) = ceil(foci_signal2_profile(:,3:4)*single_Inten/single_Inten2);
        
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
%%% Output: %%%============================================================
        save([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mat_tail],'nucleus_signal2_profile','foci_signal2_profile','cytoplasmic_signal2_profile','-append');
        
        clear ('nucleus_signal2_profile','foci_signal2_profile','cytoplasmic_signal2_profile');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        
    end
end
toc




