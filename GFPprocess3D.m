clear all
close all
tic
%% Parameter setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
list_name = 'GFPlist.xls';
in_folder = 'stacks/';
input_name = 'matchlist.xls';
mismatch_name = 'mismatch.xls';
image_type = '*.tif';
out_folder = 'Results/';
mask_folder = 'masks/';
mask_name = 'mask.mat';
hist_folder = 'Histogram/';
hist_tail = '_raw.xls';
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
pro_add = '_pro';
RNA_add = '_RNA';
nu_add = '_nucleus';
cyto_add = '_cytoplasm';
GFPim_add = '_GFP_im';
immu_add = '_immu';
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
GFP_add = '_GFP';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Data loading: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[num_list, folder_list] = xlsread(list_name);
folder_list = folder_list(strcmpi('T',folder_list(:,6)),:);
% folder_list = folder_list(strcmpi('ggo',folder_list(:,6)) | strcmpi('gbn',folder_list(:,6)) | strcmpi('gg',folder_list(:,6)),:);
[N1,N2] = size(folder_list);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Image analysis: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for list_I = 1:N1
    [sub_num, sub_list] = xlsread([folder_list{list_I,1},in_folder,input_name]);
    [M1,M2] = size(sub_list);
    [~,~,mismatch_matrix] = xlsread(mismatch_name);
    
    for list_J = 1:M1
        image_folder = [folder_list{list_I,1},in_folder,sub_list{list_J,3}];
        result_folder = [folder_list{list_I,1},out_folder];
        load([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mat_tail]);

        load([folder_list{list_I,1},mask_folder,sub_list{list_J,3},mask_name]);   %%% load 3D mask
        N_cycle = round(log2(max(mask_stack(:))))+2;   %%% Calculate the nuclear cycle number
        resolutionz = sub_num(list_J,11);

        all_color = eval(folder_list{list_I,5});
        protein_color = all_color{protein_channel};
        GFP_color = all_color{GFP_channel};
        protein_GFP_mismatch = mismatch_matrix{strcmp(protein_color,mismatch_matrix(1,:)),strcmp(GFP_color,mismatch_matrix(:,1))};
        
        [signal_stack,GFP_stack,DAPI_stack] = stack3D(imclearborder(seg_bw),protein_channel,DAPI_channel,GFP_channel,image_folder,protein_GFP_mismatch);   %%% load 3D image stacks
        [nucleus_DAPI_profile,~] = DAPI_profile3D(mask_stack,DAPI_stack,image_folder,N_cycle);
        clear DAPI_stack
        
        EL_info = get_EL(em_mask);   %%% get EL information (extreme points, EL length) from the embryo mask
        [nucleus_protein_profile,cytoplasmic_protein_profile,quanti_p,nucleus_protein_ab,~] = protein_profile3D(nucleus_DAPI_profile,max_image,EL_info,protein_channel,mask_stack,signal_stack,image_folder,N_cycle,resolution,resolutionz,'',[],[],[],[],[1.8,1]);
        nucleus_protein_profile_ab = [nucleus_protein_profile(:,1),nucleus_protein_ab,nucleus_protein_profile(:,4),nucleus_protein_profile(:,4)];
        [nucleus_GFP_profile,cytoplasmic_GFP_profile,quanti_GFP_p,nucleus_GFP_ab,~] = protein_profile3D(nucleus_DAPI_profile,max_image,EL_info,GFP_channel,mask_stack,GFP_stack,image_folder,N_cycle,resolution,resolutionz,'GFP',22,24,52,53);
        nucleus_GFP_profile_ab = [nucleus_GFP_profile(:,1),nucleus_GFP_ab,nucleus_GFP_profile(:,4),nucleus_GFP_profile(:,4)];

        clear mask_stack signal_stack GFP_stack

        GFPIF_profile(nucleus_protein_profile,nucleus_GFP_profile,image_folder,N_cycle);
        
        figure(16)
            temp_image = zeros(size(max_image,1),size(max_image,2),3);
            temp_image(:,:,1) = double(max_image(:,:,protein_channel))./max(max(double(max_image(:,:,protein_channel))));
            try
                imshow(temp_image)
            catch
                imshow_lxt(temp_image)
            end
            title(['Immunofluorescence signal: ',image_folder,', cycle = ',num2str(N_cycle)]);
        figure(17)
            temp_image = zeros(size(max_image,1),size(max_image,2),3);
            temp_image(:,:,2) = double(max_image(:,:,GFP_channel))./max(max(double(max_image(:,:,GFP_channel))));
            try
                imshow(temp_image)
            catch
                imshow_lxt(temp_image)
            end
            title(['GFP signal: ',image_folder,', cycle = ',num2str(N_cycle)]);
            
%%% Output: %%%============================================================
        result_folder = [folder_list{list_I,1},out_folder];
        if exist(result_folder) ~= 7
            mkdir(result_folder);
        end
        saveas(2,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,D3_add,figure_tail]);
        saveas(14,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,fate_add,D3_add,figure_tail]);
        saveas(22,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),GFP_add,D3_add,figure_tail]);
        saveas(24,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),GFP_add,fate_add,D3_add,figure_tail]);
        saveas(25,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),GFP_add,protein_add,D3_add,figure_tail]);
        saveas(16,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),immu_add,D3_add,figure_tail]);
        saveas(17,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),GFPim_add,D3_add,figure_tail]);
        saveas(62,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,fluc_add,D3_add,figure_tail]);
        saveas(63,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,bino_add,D3_add,figure_tail]);
        saveas(52,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),GFP_add,fluc_add,D3_add,figure_tail]);
        saveas(53,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),GFP_add,bino_add,D3_add,figure_tail]);
        
        %save([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mat_tail],'image_folder','image_type','WGA_channel','DAPI_channel','RNA_channel','protein_channel','resolution','seg_bw','cyto_bw','max_image','foci_bw','nucleus_RNA_profile','foci_RNA_profile','cytoplasmic_RNA_profile','nucleus_protein_profile','cytoplasmic_protein_profile','N_cycle','foci_data','fake_data','h','t_absolute','-append');
        save([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mat_tail],'image_folder','image_type','WGA_channel','DAPI_channel','protein_channel','GFP_channel','resolution','resolutionz','seg_bw','cyto_bw','max_image','nucleus_protein_profile','cytoplasmic_protein_profile','nucleus_GFP_profile','cytoplasmic_GFP_profile','N_cycle','resolutionz','nucleus_protein_profile_ab','nucleus_GFP_profile_ab','quanti_p','-append');
        
        if ~isempty(nucleus_protein_profile)
           xlswrite([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,nu_add,output_tail],nucleus_protein_profile);
        end
        if ~isempty(cytoplasmic_protein_profile)
           xlswrite([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,cyto_add,output_tail],cytoplasmic_protein_profile);
        end
        if ~isempty(nucleus_GFP_profile)
           xlswrite([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),GFP_add,nu_add,output_tail],nucleus_GFP_profile);
        end
        if ~isempty(cytoplasmic_GFP_profile)
           xlswrite([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),GFP_add,cyto_add,output_tail],cytoplasmic_GFP_profile);
        end
       
        sub_num(list_J,13) = N_cycle;
        
        clear ('seg_bw','cyto_bw','max_image','foci_bw','nucleus_RNA_profile','foci_RNA_profile','cytoplasmic_RNA_profile','nucleus_protein_profile','cytoplasmic_protein_profile','N_cycle','foci_data','fake_data','nucleus_protein_profile_ab','nucleus_GFP_profile_ab','quanti_p');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        
    end
    xlswrite([folder_list{list_I,1},in_folder,input_name],cat(2,sub_list,num2cell(sub_num)));
end
toc