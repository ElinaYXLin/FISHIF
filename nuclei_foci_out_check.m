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
nf_folder = 'foci_out/';
embryo_name = '11012014_10652-1-3M_FISHIF_lacZ_gal4_Bcd_par2';
embryo_type = {'b3n'};
% % embryo_name = '09302015_16249-1-5M_FISHIF_hb_gal4_Bcd_par2';
% % embryo_type = {'bgal4'};
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
        false_im = false(size(mask_stack));
        true_im = false(size(mask_stack));
        
        clear temp0
        
        if exist([folder_list{list_I,1},hist_folder,sub_list{list_J,3}(1:(end-1)),hist_tail])
            [foci_list,~,~] = xlsread([folder_list{list_I,1},hist_folder,sub_list{list_J,3}(1:(end-1)),hist_tail]);
        else
            [foci_list,~,~] = xlsread([folder_list{list_I,1},hist_folder,sub_list{list_J,3}(1:(end-1)),hist_tail2]);
        end
        
        mask_expand = dilate3D(mask_stack,8);
%         mask_expand = mask_stack;
        Nuclei_foci = mask_expand(sub2ind(size(mask_stack),round(foci_list(:,6)),round(foci_list(:,7)),round(foci_list(:,8))));
        false_id = Nuclei_foci == 0;
        false_im(sub2ind(size(mask_stack),round(foci_list(false_id,6)),round(foci_list(false_id,7)),round(foci_list(false_id,8)))) = true;
        true_im(sub2ind(size(mask_stack),round(foci_list(~false_id,6)),round(foci_list(~false_id,7)),round(foci_list(~false_id,8)))) = true;
        
        mask2D = imdilate(bwperim(any(mask_expand,3)),strel('disk',3));
        false2D = any(false_im,3);
        false2D_show = bwthicken(false2D,8);
        true2D = any(true_im,3);
        true2D_show = bwthicken(true2D,8);
        temp0(:,:,1) = mask2D+false2D_show;
        temp0(:,:,2) = mask2D+true2D_show;
        temp0(:,:,3) = mask2D;
        if any(false_id)
            figure
                imshow(double(temp0))
                title([folder_list{list_I,1},sub_list{list_J,3}],'Interpreter','none')
            if ~exist([ana_folder,embryo_name,'/',nf_folder])
                mkdir([ana_folder,embryo_name,'/',nf_folder]);
            end
            saveas(gcf,[ana_folder,embryo_name,'/',nf_folder,sub_list{list_J,3}(1:end-1),figure_tail])
            xlswrite([ana_folder,embryo_name,'/',nf_folder,sub_list{list_J,3}(1:end-1),output_tail],foci_list(false_id,:))
        end
    end
end

toc


%% Output: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        