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
nf_folder = 'nuclei_false/';
% embryo_name = '11012014_10652-1-3M_FISHIF_lacZ_gal4_Bcd';
% embryo_type = {'b3n'};
% % embryo_name = '09302015_16249-1-5M_FISHIF_hb_gal4_Bcd';
% % embryo_type = {'bgal4'};
% % % embryo_name = '12072016b_1_16249-1-5M_FISHIF_hb_gal4_H3_par3';
% % % embryo_type = {'h3'};
% % % % embryo_name = '12072016b_2_16249-1-5M_FISHIF_hb_gal4_H4K5ac_par3';
% % % % embryo_type = {'h4k5ac'};
embryo_name = '12072016b_3_16249-1-5M_FISHIF_hb_gal4_H3K27ac_par3';
embryo_type = {'h3k27ac'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Output initialization: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Data loading: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[num_list, folder_list] = xlsread(list_name);
folder_list = folder_list(strcmpi(embryo_type,folder_list(:,6)),:);
[N1,N2] = size(folder_list);
% run_list = {'[1]','[2]'};
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nuclei_false_cell = cell(0);
%% Data analysis: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for list_I = 1:N1
    [sub_num, sub_list] = xlsread([folder_list{list_I,1},in_folder,input_name]);
    [M1,M2] = size(sub_list);
    
    for list_J = 1:M1%eval(run_list{list_I})
        load([folder_list{list_I,1},mask_folder,sub_list{list_J,3},mask_name]);   %%% load 3D mask
        z_size = size(mask_stack,3);
        nuclei_false = zeros(0);
        false_im = zeros(size(mask_stack));
        
        clear temp0
        
        for zz = 1:z_size
            stats = regionprops(mask_stack(:,:,zz),bwlabel(logical(mask_stack(:,:,zz))),'MaxIntensity','MinIntensity');
            N_list = sort(unique(mask_stack(:,:,zz)));N_list = N_list(2:end);
            nuclei_false = union(nuclei_false,N_list([stats.MaxIntensity] ~= [stats.MinIntensity]));
        end
        for zz = 1:z_size
            false_im(:,:,zz) = ismember(mask_stack(:,:,zz),nuclei_false);
        end
        
        mask2D = any(mask_stack,3);
        false2D = any(false_im,3);
        temp0(:,:,1) = mask2D;
        temp0(:,:,2) = mask2D-false2D;
        temp0(:,:,3) = mask2D-false2D;
        if ~isempty(nuclei_false)
            figure
                imshow(double(temp0))
                title([folder_list{list_I,1},sub_list{list_J,3}],'Interpreter','none')
            if ~exist([ana_folder,embryo_name,'/',nf_folder])
                mkdir([ana_folder,embryo_name,'/',nf_folder]);
            end
            saveas(gcf,[ana_folder,embryo_name,'/',nf_folder,sub_list{list_J,3}(1:end-1),figure_tail])
        end
        
        if ~isempty(nuclei_false)
            nuclei_false_cell = cat(1,nuclei_false_cell,{folder_list{list_I,1},sub_list{list_J,3},nuclei_false});
        end
    end
end
xlswrite([ana_folder,embryo_name,'/',nf_folder,nf_folder(1:end-1),output_tail],nuclei_false_cell)

toc


%% Output: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        