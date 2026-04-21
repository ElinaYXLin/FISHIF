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
emmask_add = '_emmask';
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
        image_folder = [folder_list{list_I,1},in_folder,sub_list{list_J,3}];
        WGA_channel = sub_num(list_J,6);
        DAPI_channel = sub_num(list_J,7);
        protein_channel = sub_num(list_J,8);
        RNA_channel = protein_channel;
        resolution = sub_num(list_J,9);
        [seg_bw,cyto_bw,max_image,N_cycle,WGA_th0,DAPI_th0] = nuclei_segDAPI(image_folder,image_type,DAPI_channel,resolution);
%         [seg_bw,cyto_bw,max_image,N_cycle,WGA_th0,DAPI_th0] = nuclei_seg4(image_folder,image_type,WGA_channel,DAPI_channel,resolution);
%         [seg_bw,cyto_bw,max_image,N_cycle] = nuclei_seg(image_folder,image_type,WGA_channel,DAPI_channel,resolution);
%         [nucleus_protein_profile,cytoplasmic_protein_profile] = protein_profile(seg_bw,cyto_bw,max_image,signal_channel,image_folder,N_cycle)
        em_mask = get_emmask(image_folder);
        figure(1)
            out_image = double(repmat(bwperim(em_mask),[1,1,3]));
            out_image(:,:,3) = out_image(:,:,3)+double(max_image(:,:,1))/max(max(double(max_image(:,:,1))));
            imshow(out_image)
            title(image_folder,'Interpreter','none')
            
        if max(max(bwlabel(seg_bw))) >= 20
            [mask_stack,signal_stack,~] = mask3D(seg_bw,protein_channel,DAPI_channel,RNA_channel,image_folder);
            [nucleus_protein_profile,cytoplasmic_protein_profile] = protein_profile3(seg_bw,cyto_bw,max_image,protein_channel,mask_stack,signal_stack,image_folder,N_cycle,resolution);
        else
            nucleus_protein_profile = 0;
            cytoplasmic_protein_profile = 0;
            figure(2)
            figure(14)
            figure(62)
            figure(63)
        end
        clear mask_stack signal_stack
        
%%% Output: %%%============================================================
        result_folder = [folder_list{list_I,1},out_folder];
        if exist(result_folder) ~= 7
            mkdir(result_folder);
        end
        saveas(1,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),emmask_add,figure_tail]);
        saveas(51,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),seg_add,figure_tail]);
        saveas(2,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),pro_add,figure_tail]);
        saveas(14,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,fate_add,figure_tail]);
        saveas(62,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),pro_add,fluc_add,figure_tail]);
        saveas(63,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),pro_add,bino_add,figure_tail]);
        save([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mat_tail],'em_mask','image_folder','image_type','WGA_channel','DAPI_channel','protein_channel','resolution','seg_bw','cyto_bw','max_image','nucleus_protein_profile','cytoplasmic_protein_profile','N_cycle','WGA_th0','DAPI_th0');
        xlswrite([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),nu_add,output_tail],nucleus_protein_profile);
        xlswrite([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),cyto_add,output_tail],cytoplasmic_protein_profile);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        sub_list{list_J,13} = N_cycle;
    end
    %xlswrite([folder_list{list_I,1},in_folder,input_name],sub_list);
end