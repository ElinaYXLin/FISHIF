clear all
close all

tic

%% Parameter setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
list_name = 'Duallist.xls';
in_folder = 'stacks/';
input_name = 'matchlist.xls';

image_type = '*.tif';
mask_folder = 'masks/';
mask_name = 'mask';
out_folder = 'Results/';
fit_folder_protein = {'Histogram_protein_A/','Histogram_protein_M/','Histogram_protein_P/'};
mask2_add = '_mask_area';
fit_add = '_spot_fit';
hist_tail = '_raw.xls';
hist_tail2 = '_raw.xlsx';
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
purify_add = '_purify';
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Data loading: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[num_list, folder_list] = xlsread(list_name);
% folder_list = folder_list(strcmpi('b',folder_list(:,6)) | strcmpi('ba',folder_list(:,6)),:);
% folder_list = folder_list(strcmpi('b3',folder_list(:,6)),:);
% folder_list = folder_list(strcmpi('h',folder_list(:,6)),:);
% folder_list = folder_list(strcmpi('b3n',folder_list(:,6)),:);
folder_list = folder_list(strcmpi('S',folder_list(:,6)),:);
% folder_list = folder_list(strcmpi('cgal4new',folder_list(:,6)),:);
% folder_list = folder_list(strcmpi('h4k5ac',folder_list(:,6)),:);
% folder_list = folder_list(strcmpi('cbgal4',folder_list(:,6)),:);
% folder_list = folder_list(strcmpi('bgal4',folder_list(:,6)),:);
[N1,N2] = size(folder_list);
% run_list = {2:4,[],[]};
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for list_fit = 1:length(fit_folder_protein)
%% Data loaing: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for list_I = 1:N1
        [sub_num, sub_list] = xlsread([folder_list{list_I,1},in_folder,input_name]);
        [M1,M2] = size(sub_list);

        for list_J = 1:M1%run_list{list_I}
    % %         load([folder_list{list_I,1},mask_folder,sub_list{list_J,3},mask_name,mat_tail]);   %%% load 3D mask
    % %         mask1 = any(mask_stack,3);

            load([folder_list{list_I,1},fit_folder_protein{list_fit},sub_list{list_J,3}(1:end-1),mask2_add,mat_tail]);   %%% load 3D mask from the spot fitting result
% %             mask_new = imdilate(~embryo_mask,strel('disk',400)) & ~imdilate(~embryo_mask,strel('disk',10));
            
            load([folder_list{list_I,1},mask_folder,sub_list{list_J,3},mask_name,mat_tail]);   %%% load 3D mask
            embryo_mask = embryo_mask & ~imdilate(mask_stack(xylim(1):xylim(2),xylim(3):xylim(4),:),strel('disk',2));
            mask_new = imdilate(~embryo_mask,strel('disk',400)) & ~imdilate(~embryo_mask,strel('disk',10));
            clear mask_stack

            if exist([folder_list{list_I,1},fit_folder_protein{list_fit},sub_list{list_J,3}(1:(end-1)),hist_tail],'file')
                Inten0 = xlsread([folder_list{list_I,1},fit_folder_protein{list_fit},sub_list{list_J,3}(1:(end-1)),hist_tail]);
            else
                Inten0 = xlsread([folder_list{list_I,1},fit_folder_protein{list_fit},sub_list{list_J,3}(1:(end-1)),hist_tail2]);
            end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        


%% Purify the identified spots: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            size0 = size(mask_new);
            x0 = min(max(round(Inten0(:,6)),1),size0(1));
            y0 = min(max(round(Inten0(:,7)),1),size0(2));
            z0 = min(max(round(Inten0(:,8)),1),size0(3));

            ind0 = sub2ind(size0,x0,y0,z0);
            spot_true = mask_new(ind0);

            save([folder_list{list_I,1},fit_folder_protein{list_fit},sub_list{list_J,3}(1:end-1),purify_add,mat_tail],'spot_true')
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        end
    end
end
toc

