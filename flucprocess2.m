clear

%% Parameter setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
list_name = 'fluclist.xls';
in_folder = 'stacks/';
input_name = 'matchlist.xls';
image_type = '*.tif';
out_folder = 'Results/';
image_source = 'stacks/';
output_tail = '.xls';
figure_tail = '.fig';
mat_tail = '.mat';
mark_list = 'rb';
seg_add = '_seg';
pro_add = '_pro';
nu_add = '_nucleus';
cyto_add = '_cytoplasm';
protein_add = '_protein';
fluc_add = '_fluc';
bino_add = '_bino';
slide_add = '_slide';
sample_add = '_sample';
down_add = '_down';
corr_add = '_corr';
cluster_add = '_cluster';
fate_add = '_fate';
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
    close all
    for list_J = 1:M1
        %image_folder = [folder_list{list_I,1},in_folder,sub_list{list_J,3}];
        %WGA_channel = sub_num(list_J,6);
        %DAPI_channel = sub_num(list_J,7);
        %protein_channel = sub_num(list_J,8);
        %resolution = sub_num(list_J,9);
        result_folder = [folder_list{list_I,1},out_folder];
        load([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mat_tail]);
        ind1 = find((image_folder(1:(end-1))=='/')|(image_folder(1:(end-1))=='\'),1,'last');
        image_folder = [folder_list{list_I,1},image_source,sub_list{list_J,3}];
        %[seg_bw,cyto_bw,max_image,N_cycle,WGA_th0,DAPI_th0] = nuclei_segDAPI(image_folder,image_type,DAPI_channel,resolution);
        %[seg_bw,cyto_bw,max_image,N_cycle] = nuclei_seg4(image_folder,image_type,WGA_channel,DAPI_channel,resolution);
        [mask_stack,signal_stack,~] = mask3D(seg_bw,protein_channel,DAPI_channel,RNA_channel,image_folder);
        %[seg_bw,cyto_bw,max_image,N_cycle] = nuclei_seg(image_folder,image_type,WGA_channel,DAPI_channel,resolution);
%         [nucleus_protein_profile,cytoplasmic_protein_profile] = protein_profile(seg_bw,cyto_bw,max_image,signal_channel,image_folder,N_cycle)

%         figure(62)
%             clf
%             h1 = axes;
%         figure(63)
%             clf
%             h2 = axes;
%         figure(64)
%             h3 = subplot(2,ceil(M1./4),ceil(list_J./2));
%         figure(65)
%             window_show = 11:10:41;
%             mh = length(window_show);
%             h4 = zeros(1,mh);
%             for i_h = 1:mh
%                 h4(i_h) = subplot(mh,ceil(M1./2),ceil(list_J./2)+(i_h-1)*ceil(M1./2));
%             end
% %         figure(66)
% %             h5 = subplot(2,ceil(M1./4),ceil(list_J./2));
% 
%         figure(66)
%             h5 = subplot(2,ceil(M1./2),list_J);
%             
%         figure(67)
%             k_show = 2:5:20;
%             mh = length(k_show);
%             h6 = zeros(1,mh);
%             for i_h = 1:mh
%                 h6(i_h) = subplot(mh,ceil(M1./2),ceil(list_J./2)+(i_h-1)*ceil(M1./2));
%             end
%             
%         figure(68)
%             h7 = subplot(2,ceil(M1./4),ceil(list_J./2));
            
        %[nucleus_protein_profile,fluc_protein_profile] = protein_fluc4(seg_bw,cyto_bw,max_image,protein_channel,mask_stack,signal_stack,image_folder,resolution,h1,h2,[],[],h3,mark_list(mod(list_J-1,2)+1),h4,mark_list(mod(list_J-1,2)+1),window_show,h5,h6,mark_list(mod(list_J-1,2)+1),k_show,h7);
        [nucleus_protein_profile,fluc_protein_profile] = protein_fluc4(seg_bw,cyto_bw,max_image,protein_channel,mask_stack,signal_stack,image_folder,resolution);
        %[nucleus_protein_profile,fluc_protein_profile] = both_fluc3(seg_bw,cyto_bw,max_image,protein_channel,mask_stack,signal_stack,image_folder,resolution);
        clear mask_stack signal_stack
        
%%% Output: %%%============================================================
        result_folder = [folder_list{list_I,1},out_folder];
        if exist(result_folder) ~= 7
            mkdir(result_folder);
        end
        %saveas(1,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),seg_add,figure_tail]);
        saveas(62,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),pro_add,fluc_add,figure_tail]);
        saveas(63,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),pro_add,bino_add,figure_tail]);
        saveas(64,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),pro_add,slide_add,figure_tail]);
        saveas(65,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),pro_add,slide_add,sample_add,figure_tail]);
        saveas(66,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),pro_add,corr_add,figure_tail]);
        saveas(67,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),pro_add,cluster_add,sample_add,figure_tail]);
        saveas(68,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),pro_add,cluster_add,figure_tail]);
        
        save([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mat_tail],'image_folder','image_type','WGA_channel','DAPI_channel','protein_channel','resolution','seg_bw','cyto_bw','max_image','nucleus_protein_profile','fluc_protein_profile','N_cycle','-append');
        xlswrite([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),nu_add,output_tail],nucleus_protein_profile);
        xlswrite([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),fluc_add,output_tail],fluc_protein_profile);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        sub_list{list_J,13} = N_cycle;
    end
    %xlswrite([folder_list{list_I,1},in_folder,input_name],sub_list);
end