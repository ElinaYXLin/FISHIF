clear

%% Parameter setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
list_name = 'GFPlist.xls';
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
GFP_add = '_GFP';
GFPim_add = '_GFP_im';
fate_add = '_fate';
immu_add = '_immu';
fluc_add = '_fluc';
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
    for list_J = M1:M1
        image_folder = [folder_list{list_I,1},in_folder,sub_list{list_J,3}];
        WGA_channel = sub_num(list_J,6);
        DAPI_channel = sub_num(list_J,7);
        protein_channel = sub_num(list_J,8);
        GFP_channel = sub_num(list_J,10);
        resolution = sub_num(list_J,9);
        result_folder = [folder_list{list_I,1},out_folder];
        load([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mat_tail]);
        N_cycle = round(log2(max(max(bwlabel(seg_bw)))))+2;   %%% Calculate the nuclear cycle number
        [mask_stack,signal_stack] = mask3D(seg_bw,protein_channel,DAPI_channel,image_folder);
        [nucleus_protein_profile,cytoplasmic_protein_profile] = protein_profile3(seg_bw,cyto_bw,max_image,protein_channel,mask_stack,signal_stack,image_folder,N_cycle,resolution);
        clear mask_stack signal_stack
        [mask_stack,signal_stack2] = mask3D(seg_bw,GFP_channel,DAPI_channel,image_folder);
        [nucleus_GFP_profile,cytoplasmic_GFP_profile] = protein_profile3(seg_bw,cyto_bw,max_image,GFP_channel,mask_stack,signal_stack2,image_folder,N_cycle,resolution,'GFP',22,24);
        clear mask_stack signal_stack2
        %[seg_bw,cyto_bw,max_image,N_cycle] = nuclei_seg(image_folder,image_type,WGA_channel,DAPI_channel,resolution);
        %[nucleus_protein_profile,cytoplasmic_protein_profile] = protein_profile(seg_bw,cyto_bw,max_image,protein_channel,image_folder,N_cycle);
        %[nucleus_GFP_profile,cytoplasmic_GFP_profile] = GFP_profile(seg_bw,cyto_bw,max_image,GFP_channel,image_folder,N_cycle);
        GFPIF_profile(nucleus_protein_profile,nucleus_GFP_profile,image_folder,N_cycle);
        
        figure(16)
            temp_image = zeros(size(max_image,1),size(max_image,2),3);
            temp_image(:,:,1) = double(max_image(:,:,protein_channel))./max(max(double(max_image(:,:,protein_channel))));
            imshow(temp_image);
            title(['Immunofluorescence signal: ',image_folder,', cycle = ',num2str(N_cycle)]);
        figure(17)
            temp_image = zeros(size(max_image,1),size(max_image,2),3);
            temp_image(:,:,2) = double(max_image(:,:,GFP_channel))./max(max(double(max_image(:,:,GFP_channel))));
            imshow(temp_image);
            title(['GFP signal: ',image_folder,', cycle = ',num2str(N_cycle)]);
            
%%% Output: %%%============================================================
        result_folder = [folder_list{list_I,1},out_folder];
        if exist(result_folder) ~= 7
            mkdir(result_folder);
        end
        %saveas(51,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),seg_add,figure_tail]);
        saveas(2,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),pro_add,figure_tail]);
        saveas(14,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,fate_add,figure_tail]);
        saveas(22,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),GFP_add,figure_tail]);
        saveas(24,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),GFP_add,fate_add,figure_tail]);
        saveas(25,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),GFP_add,protein_add,figure_tail]);
        saveas(16,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),immu_add,figure_tail]);
        saveas(17,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),GFPim_add,figure_tail]);
        saveas(62,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),pro_add,fluc_add,figure_tail]);

        
        
        save([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mat_tail],'image_folder','image_type','WGA_channel','DAPI_channel','protein_channel','GFP_channel','resolution','seg_bw','cyto_bw','max_image','nucleus_protein_profile','cytoplasmic_protein_profile','nucleus_GFP_profile','cytoplasmic_GFP_profile','N_cycle','WGA_th0','DAPI_th0');
        xlswrite([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),nu_add,protein_add,output_tail],nucleus_protein_profile);
        xlswrite([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),cyto_add,protein_add,output_tail],cytoplasmic_protein_profile);
        xlswrite([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),nu_add,GFP_add,output_tail],nucleus_GFP_profile);
        xlswrite([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),cyto_add,GFP_add,output_tail],cytoplasmic_GFP_profile);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        sub_list{list_J,13} = N_cycle;
    end
    %xlswrite([folder_list{list_I,1},in_folder,input_name],sub_list);
end