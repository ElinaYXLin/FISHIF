clear

%% Parameter setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
list_name = 'proteinlist.xls';
in_folder = 'stacks/';
input_name = 'matchlist.xls';
image_type = '*.tif';
out_folder = 'Results/';
output_tail = '.xls';
figure_tail = '.fig';
mat_tail = '.mat';
seg_add = '_seg';
pro_add = '_pro';
nu_add = '_nucleus';
cyto_add = '_cytoplasm';
protein_add = '_protein';
fate_add = '_fate';
fluc_add = '_fluc';
bino_add = '_bino';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Data loading: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[num_list, folder_list] = xlsread(list_name);
folder_list = folder_list(strcmpi('T',folder_list(:,6)),:);
[N1,N2] = size(folder_list);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Image analysis: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for list_I = 1:N1
    [sub_num, sub_list] = xlsread([folder_list{list_I,1},in_folder,input_name]);
    [M1,M2] = size(sub_list);
    for list_J = 1:M1
        %image_folder = [folder_list{list_I,1},in_folder,sub_list{list_J,3}];
        %WGA_channel = sub_num(list_J,6);
        %DAPI_channel = sub_num(list_J,7);
        %protein_channel = sub_num(list_J,8);
        %resolution = sub_num(list_J,9);
        result_folder = [folder_list{list_I,1},out_folder];
        load([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mat_tail]);
        
        %[seg_bw,cyto_bw,max_image,N_cycle] = nuclei_seg4(image_folder,image_type,WGA_channel,DAPI_channel,resolution);
        [mask_stack,signal_stack] = mask3D(seg_bw,protein_channel,DAPI_channel,image_folder);
        %[seg_bw,cyto_bw,max_image,N_cycle] = nuclei_seg(image_folder,image_type,WGA_channel,DAPI_channel,resolution);
%         [nucleus_protein_profile,cytoplasmic_protein_profile] = protein_profile(seg_bw,cyto_bw,max_image,signal_channel,image_folder,N_cycle)
        [nucleus_protein_profile,cytoplasmic_protein_profile] = protein_profile3(seg_bw,cyto_bw,max_image,protein_channel,mask_stack,signal_stack,image_folder,N_cycle,resolution);
        clear mask_stack signal_stack
        
%%% Output: %%%============================================================
        result_folder = [folder_list{list_I,1},out_folder];
        if exist(result_folder) ~= 7
            mkdir(result_folder);
        end
        %saveas(1,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),seg_add,figure_tail]);
        saveas(2,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),pro_add,figure_tail]);
        saveas(14,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,fate_add,figure_tail]);
        saveas(62,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),pro_add,fluc_add,figure_tail]);
        saveas(63,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),pro_add,bino_add,figure_tail]);
        save([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mat_tail],'image_folder','image_type','WGA_channel','DAPI_channel','protein_channel','resolution','seg_bw','cyto_bw','max_image','nucleus_protein_profile','cytoplasmic_protein_profile','N_cycle');
        xlswrite([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),nu_add,output_tail],nucleus_protein_profile);
        xlswrite([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),cyto_add,output_tail],cytoplasmic_protein_profile);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        sub_list{list_J,13} = N_cycle;
    end
    %xlswrite([folder_list{list_I,1},in_folder,input_name],sub_list);
end