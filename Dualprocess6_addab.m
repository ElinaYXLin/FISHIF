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
hist_folder2 = 'Histogram_nullo/';
fit_folder2 = 'Histogram_default/';
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
signal2_add = '_nullo';
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
sub_pos = [3,3];
cycle_pos0 = zeros(1,sub_pos(1));
cycle_range = [11,12,13];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Data loading: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[num_list, folder_list] = xlsread(list_name);
folder_list = folder_list(strcmpi('b3h',folder_list(:,6)),:);
[N1,N2] = size(folder_list);
% run_list = {'[1,5]','1:M1'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Image analysis: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for list_I = 1:N1
    [sub_num, sub_list] = xlsread([folder_list{list_I,1},in_folder,input_name]);
    [M1,M2] = size(sub_list);
    if ~isempty(out_folder0)
        copyfile([folder_list{list_I,1},out_folder],[folder_list{list_I,1},out_folder0]);
    end
    
    for list_J = 1:M1%eval(run_list{list_I})
        if isempty(strfind(sub_list{list_J,3},'_60X'))
            [~,~,mismatch_matrix] = xlsread(mismatch_name);
            load(xymismatch_name);
        else
            [~,~,mismatch_matrix] = xlsread(mismatch_name2);
            load(xymismatch_name2);
        end
        
        flip_axis = any(cellfun(@(x) ~isempty(strfind(folder_list{list_I,1},x)),flip_list));   %%% Check whether the protein profile has an opposite direction (to Bcd)
        
        image_folder = [folder_list{list_I,1},in_folder,sub_list{list_J,3}];
        result_folder = [folder_list{list_I,1},out_folder];
        resolutionz = sub_num(list_J,11);
        signal2_channel = sub_num(list_J,12);
        Nbin = sub_num(list_J,2);
        Mdim = sub_num(list_J,3);
        load([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mat_tail]);
        load([folder_list{list_I,1},mask_folder,sub_list{list_J,3},mask_name]);   %%% load 3D mask
        z_size = size(mask_stack,3);
        load([folder_list{list_I,1},fit_folder,sub_list{list_J,3}(1:end-1),fit_add,mat_tail]);   %%% load spot intensity fitting result

        Inten_thresh = b*N_thresh;   %%% set foci intensity threshold
%         single_Inten = b;
%         load([folder_list{list_I,1},fit_folder2,sub_list{list_J,3}(1:end-1),fit_add,mat_tail]);   %%% load spot intensity fitting result
%         Inten_thresh2 = b*N_thresh;   %%% set foci intensity threshold

        all_color = eval(folder_list{list_I,5});
        protein_color = all_color{protein_channel};
        RNA_color = all_color{RNA_channel};
        protein_RNA_mismatch = mismatch_matrix{strcmp(protein_color,mismatch_matrix(1,:)),strcmp(RNA_color,mismatch_matrix(:,1))};
%         protein_RNA_xymismatch = {eval([protein_color,'_',RNA_color]),eval([protein_color,'_',RNA_color,'_con']),eval([protein_color,'_',RNA_color,'_x0']),Nbin,Mdim,resolution,resolution_mismatch};
%         signal2_color = all_color{signal2_channel};
%         protein_signal2_mismatch = mismatch_matrix{strcmp(protein_color,mismatch_matrix(1,:)),strcmp(signal2_color,mismatch_matrix(:,1))};
%         protein_signal2_xymismatch = {eval([protein_color,'_',signal2_color]),eval([protein_color,'_',signal2_color,'_con']),eval([protein_color,'_',signal2_color,'_x0']),Nbin,Mdim,resolution,resolution_mismatch};

%         [signal_stack,RNA_stack,DAPI_stack] = stack3D(imclearborder(seg_bw),protein_channel,DAPI_channel,RNA_channel,image_folder,protein_RNA_mismatch);   %%% load 3D image stacks
%         [nucleus_DAPI_profile,DNA_mask] = DAPI_profile3D(mask_stack,DAPI_stack,image_folder,N_cycle);
%         clear DAPI_stack
        
        
        [foci_bw3D,~,~] = modified_foci3D([folder_list{list_I,1},hist_folder,sub_list{list_J,3}(1:(end-1)),hist_tail],max_image,RNA_channel,mask_stack,N_cycle,image_folder,protein_RNA_mismatch,Inten_thresh);
%         [foci_bw3D,max_image00,SS] = modified_foci3D2([folder_list{list_I,1},hist_folder,sub_list{list_J,3}(1:(end-1)),hist_tail],max_image,RNA_channel,mask_stack,N_cycle,image_folder,protein_RNA_mismatch,protein_RNA_xymismatch,Inten_thresh);
        foci_bw2D = max(foci_bw3D,[],3);
        
%         index_true = (nucleus_protein_profile(:,4) > 1) & (nucleus_protein_profile(:,4) < z_size);   %%% check whether the nuclei are on the first or last slice
% 
%         if ~flip_axis
%             mean_post = mean(nucleus_protein_profile((nucleus_protein_profile(:,1) > 0.8) & index_true,2));
%         else
%             mean_post = mean(nucleus_protein_profile((nucleus_protein_profile(:,1) < 0.2) & index_true,2));
%         end
% 
% %         nucleus_protein_ab = (nucleus_protein_profile(:,2)-mean_post)/quanti_p(1);
% %         nucleus_protein_ab(nucleus_protein_ab < 0) = 0;
% %         nucleus_protein_profile_ab = [nucleus_protein_profile(:,1),nucleus_protein_ab,nucleus_protein_profile(:,4)];
% %         save([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mat_tail],'nucleus_protein_profile_ab','-append');
% 
%         for ir = 1:length(foci_data)
%             foci_data{ir}(:,13) = (foci_data{ir}(:,13)*quanti_p(1)+quanti_p(2)+mean_post)/quanti_p(1);
%             fake_data{ir}(:,13) = (fake_data{ir}(:,13)*quanti_p(1)+quanti_p(2)+mean_post)/quanti_p(1);
%             foci_data2{ir}(:,13) = (foci_data2{ir}(:,13)*quanti_p(1)+quanti_p(2)+mean_post)/quanti_p(1);
%             fake_data2{ir}(:,13) = (fake_data2{ir}(:,13)*quanti_p(1)+quanti_p(2)+mean_post)/quanti_p(1);
%         end
%         save([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mat_tail],'foci_data','fake_data','foci_data2','fake_data2','-append');
% 
% %         if any(nucleus_protein_profile(:,1) ~= nucleus_RNA_profile(:,1))
% %             nucleus_RNA_profile(:,1) = 1-nucleus_RNA_profile(:,1);
% % %             nucleus_signal2_profile(:,1) = 1-nucleus_signal2_profile(:,1);
% %             foci_RNA_profile(:,1) = 1-foci_RNA_profile(:,1);
% % %             foci_signal2_profile(:,1) = 1-foci_signal2_profile(:,1);
% %             cytoplasmic_RNA_profile(:,1) = 1-cytoplasmic_RNA_profile(:,1);
% % %             cytoplasmic_signal2_profile(:,1) = 1-cytoplasmic_signal2_profile(:,1);
% %         end
% % %         save([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mat_tail],'nucleus_RNA_profile','nucleus_signal2_profile','foci_RNA_profile','foci_signal2_profile','cytoplasmic_RNA_profile','cytoplasmic_signal2_profile','-append');
% %         save([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mat_tail],'nucleus_RNA_profile','foci_RNA_profile','cytoplasmic_RNA_profile','-append');
        save([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mat_tail],'foci_bw2D','-append');
%         rmvar([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mat_tail],'foci_bw3D')
%         
    end
end
toc