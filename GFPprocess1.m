clear all
close all
tic
%% Parameter setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
list_name = 'GFPlist.xls';
in_folder = 'stacks/';
input_name = 'matchlist.xls';
image_type = '*.tif';
out_folder = 'Results/';
hist_folder = 'Histogram/';
mask_folder = 'masks/';
mask_name = 'mask.mat';
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
a3D_add = '_3D';
RNA_add = '_RNA';
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
        GFP_channel = sub_num(list_J,10);
        protein_channel = sub_num(list_J,8);
        resolution = sub_num(list_J,9);
        if size(sub_num,2) > 10
            resolutionz = sub_num(list_J,11);
        else
            resolutionz = 4*resolution;
        end
        %[seg_bw,cyto_bw,max_image] = nuclei_seg4(image_folder,image_type,WGA_channel,DAPI_channel,resolution);

        %foci_bw = RNA_seg(seg_bw,max_image,RNA_channel,WGA_channel,resolution,image_folder);
        result_folder = [folder_list{list_I,1},out_folder];
        load([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mat_tail]);
        %[seg_bw,cyto_bw,max_image,N_cycle,WGA_th0,DAPI_th0] = nuclei_segDAPI(image_folder,image_type,DAPI_channel,resolution);
        [signal_stack,GFP_stack] = stack3D(imclearborder(seg_bw),protein_channel,DAPI_channel,GFP_channel,image_folder);
        load([folder_list{list_I,1},mask_folder,sub_list{list_J,3},mask_name]);
        %N_cycle = sub_num(list_J,11);
%        foci_bw = RNA_seg_manual(seg_bw,max_image,RNA_channel,WGA_channel,resolution,image_folder,N_cycle);
        N_cycle = round(log2(max(mask_stack(:))))+2;   %%% Calculate the nuclear cycle number

        %nucleus_RNA_profile = xlsread([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),RNA_add,nu_add,output_tail]);
        %cytoplasmic_RNA_profile = xlsread([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),RNA_add,cyto_add,output_tail]);
        %foci_RNA_profile = xlsread([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),RNA_add,foci_add,output_tail]);
        %nucleus_protein_profile = xlsread([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,nu_add,output_tail]);
        %cytoplasmic_protein_profile = xlsread([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,cyto_add,output_tail]);
        
        %%% Re-imaging of the foci region: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %foci_bw0 = imdilate(foci_bw,strel('disk',4));
        %bw_perim_g = bwperim(foci_bw0);
        %bw_perim_WGA = bwperim(seg_bw);
        %s_max = size(max_image);
        %new_image = zeros(s_max(1),s_max(2),3);
        %new_image(:,:,1) = double(max_image(:,:,RNA_channel))/double(max(max(max_image(:,:,RNA_channel))));
        %new_image(:,:,2) = 0;
        %new_image(:,:,3) = 0;
        %overlay = imoverlay(new_image, bw_perim_g, [1,1,1]);
        %overlay = imoverlay(overlay,bw_perim_WGA,[0,0,1]);
        %figure(4)
        %imshow(overlay)
        %title([image_folder,' (cycle = ',num2str(N_cycle),', white: transcription foci recognition, blue: nucleus recognition)'])
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %[foci_bw00,max_image00,SS,foci_layer] = modified_foci([folder_list{list_I,1},hist_folder,sub_list{list_J,3}(1:(end-1)),hist_tail],max_image,RNA_channel,imclearborder(seg_bw),mask_stack,N_cycle,image_folder);
        %[nucleus_RNA_profile,foci_RNA_profile,cytoplasmic_RNA_profile] = RNA_profile(imclearborder(seg_bw),cyto_bw,foci_bw00,max_image00,RNA_channel,resolution,image_folder,N_cycle,SS);   %%% Input the foci area matrix from modified_foci.m if gaussian fitting method is used
        
        %[nucleus_protein_profile,cytoplasmic_protein_profile] = protein_profile3(imclearborder(seg_bw),cyto_bw,max_image,protein_channel,mask_stack,signal_stack,image_folder,N_cycle);
        [nucleus_protein_profile,cytoplasmic_protein_profile,quanti_p] = protein_profile3D(imclearborder(seg_bw),cyto_bw,max_image,protein_channel,mask_stack,signal_stack,image_folder,N_cycle,resolution);
        [nucleus_GFP_profile,cytoplasmic_GFP_profile,quanti_GFP_p] = protein_profile3D(imclearborder(seg_bw),cyto_bw,max_image,GFP_channel,mask_stack,GFP_stack,image_folder,N_cycle,resolution,'GFP',22,24);
        %signal_stack = (signal_stack+quanti_p(2))/quanti_p(1);
        %nucleus_protein_profile_ab = [nucleus_protein_profile(:,1),(nucleus_protein_profile(:,2)+quanti_p(2))/quanti_p(1)];
        %[foci_data,fake_data,h,t_absolute,r_size] = dual_local_local(foci_bw00,max_image00(:,:,RNA_channel),SS,foci_layer,mask_stack,signal_stack,RNA_stack,nucleus_protein_profile_ab,image_folder,N_cycle,resolution,resolutionz);
        %clear foci_bw00 max_image00 SS foci_layer
        clear mask_stack signal_stack GFP_stack 
        
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

%         [nucleus_protein_profile,cytoplasmic_protein_profile] = protein_profile(imclearborder(seg_bw),cyto_bw,max_image,signal_channel,image_folder,N_cycle)
        %dual_profile(nucleus_protein_profile,nucleus_RNA_profile,foci_RNA_profile,image_folder,N_cycle);
        
            
%%% Output: %%%============================================================
        result_folder = [folder_list{list_I,1},out_folder];
        if exist(result_folder) ~= 7
            mkdir(result_folder);
        end
        %saveas(51,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),nu_add,seg_add,figure_tail]);
        saveas(2,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,a3D_add,figure_tail]);
        %saveas(3,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,th_add,figure_tail]);
%         saveas(4,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,seg_add,figure_tail]);
%         saveas(5,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),fish_add,int_add,figure_tail]);
%         saveas(6,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),fish_add,num_add,figure_tail]);
%         saveas(7,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,num_add,figure_tail]);
%         saveas(8,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,int_add,figure_tail]);
%         saveas(9,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,int_add,cmp_add,figure_tail]);
%         %saveas(10,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,th_add,sel_add,figure_tail]);
%         saveas(11,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),nu_add,fate_add,figure_tail]);
%         saveas(12,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),nu_add,reg_add,figure_tail]);
%         saveas(13,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,reg_add,figure_tail]);
        saveas(14,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,fate_add,a3D_add,figure_tail]);
%         saveas(15,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),nu_add,foci_add,reg_add,figure_tail]);
        saveas(22,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),GFP_add,figure_tail]);
        saveas(24,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),GFP_add,fate_add,figure_tail]);
        saveas(25,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),GFP_add,protein_add,figure_tail]);
        saveas(16,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),immu_add,figure_tail]);
        saveas(17,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),GFPim_add,figure_tail]);
        saveas(62,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,fluc_add,a3D_add,figure_tail]);
        saveas(63,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,bino_add,a3D_add,figure_tail]);
%         saveas(73,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,local_add,foci_add,figure_tail]);
%         saveas(74,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,local_add,hist_add,figure_tail]);
%         saveas(75,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,local_add,RNA_add,figure_tail]);
%         saveas(76,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,local_add,abs_add,foci_add,figure_tail]);
%         saveas(77,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,local_add,abs_add,hist_add,figure_tail]);
%         saveas(80,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,local_add,rad_add,figure_tail]);
%         saveas(81,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,local_add,abs_add,rad_add,figure_tail]);
        
        save([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mat_tail],'image_folder','image_type','WGA_channel','DAPI_channel','protein_channel','GFP_channel','resolution','resolutionz','seg_bw','cyto_bw','max_image','nucleus_protein_profile','cytoplasmic_protein_profile','nucleus_GFP_profile','cytoplasmic_GFP_profile','N_cycle','-append');
%         xlswrite([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),nu_add,protein_add,output_tail],nucleus_protein_profile);
%         xlswrite([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),cyto_add,protein_add,output_tail],cytoplasmic_protein_profile);
%         xlswrite([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),nu_add,GFP_add,output_tail],nucleus_GFP_profile);
%         xlswrite([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),cyto_add,GFP_add,output_tail],cytoplasmic_GFP_profile);
%         save([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mat_tail],'image_folder','image_type','WGA_channel','DAPI_channel','RNA_channel','protein_channel','resolution','seg_bw','cyto_bw','max_image','nucleus_RNA_profile','foci_RNA_profile','cytoplasmic_RNA_profile','nucleus_protein_profile','cytoplasmic_protein_profile','N_cycle','foci_data','fake_data','h','t_absolute','r_size','-append');
        
%         if ~isempty(nucleus_RNA_profile)
%            xlswrite([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),RNA_add,nu_add,output_tail],nucleus_RNA_profile);
%         end
%         if ~isempty(cytoplasmic_RNA_profile)
%            xlswrite([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),RNA_add,cyto_add,output_tail],cytoplasmic_RNA_profile);
%         end
%         if ~isempty(foci_RNA_profile)
%            xlswrite([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),RNA_add,foci_add,output_tail],foci_RNA_profile);
%         end
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
%         for I_data = 1:length(foci_data)
%             if ~isempty(foci_data{I_data})
%                 xlswrite([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,local_add,foci_add,output_tail],foci_data{I_data},I_data);
%             end
%         end
%         for I_data = 1:length(fake_data)
%             if ~isempty(fake_data{I_data})
%                 xlswrite([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,local_add,fake_add,foci_add,output_tail],fake_data{I_data},I_data);
%             end
%         end
       
        sub_num(list_J,13) = N_cycle;
        
        clear ('seg_bw','cyto_bw','max_image','foci_bw','nucleus_RNA_profile','foci_RNA_profile','cytoplasmic_RNA_profile','nucleus_protein_profile','cytoplasmic_protein_profile','N_cycle','foci_data','fake_data');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        
    end
    xlswrite([folder_list{list_I,1},in_folder,input_name],cat(2,sub_list,num2cell(sub_num)));
end
toc