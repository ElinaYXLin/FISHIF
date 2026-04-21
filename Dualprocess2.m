clear

%% Parameter setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
list_name = 'Duallist.xls';
in_folder = 'stacks/';
input_name = 'matchlist.xls';
image_type = '*.tif';
out_folder = 'Results/';
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
nu_add = '_nucleus';
cyto_add = '_cytoplasm';
fate_add = '_fate';
reg_add = '_regulation';
immu_add = '_immu';
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
        psize0 = zeros(0);
        g_threshold0 = zeros(0);
        WGA_th0 = zeros(0);
        DAPI_th0 = zeros(0);
        result_folder = [folder_list{list_I,1},out_folder];
        load([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mat_tail]);
        N_cycle = round(log2(max(max(bwlabel(seg_bw)))))+2;   %%% Calculate the nuclear cycle number
        %image_folder = [folder_list{list_I,1},in_folder,sub_list{list_J,3}];
        %WGA_channel = sub_num(list_J,6);
        %DAPI_channel = sub_num(list_J,7);
        %RNA_channel = sub_num(list_J,10);
        %protein_channel = sub_num(list_J,8);
        %resolution = sub_num(list_J,9);
        %[seg_bw,cyto_bw,max_image,N_cycle,WGA_th0,DAPI_th0] = nuclei_seg4(image_folder,image_type,WGA_channel,DAPI_channel,resolution)
        %[foci_bw,psize0,g_threshold0] = RNA_seg(seg_bw,max_image,RNA_channel,WGA_channel,resolution,image_folder,N_cycle);
        %max_image = imagepack(image_folder,image_type);
        
        [nucleus_RNA_profile,foci_RNA_profile,cytoplasmic_RNA_profile] = RNA_profile(seg_bw,cyto_bw,foci_bw,max_image,RNA_channel,resolution,image_folder,N_cycle);
        [nucleus_protein_profile,cytoplasmic_protein_profile] = protein_profile(seg_bw,cyto_bw,max_image,protein_channel,image_folder,N_cycle);
        dual_profile(nucleus_protein_profile,nucleus_RNA_profile,foci_RNA_profile,image_folder,N_cycle);
        
        figure(16)
            temp_image = zeros(size(max_image,1),size(max_image,2),3);
            temp_image(:,:,2) = double(max_image(:,:,protein_channel))./max(max(double(max_image(:,:,protein_channel))));
            imshow(temp_image);
            title(['Immunofluorescence signal: ',image_folder,', cycle = ',num2str(N_cycle)]);
        figure(17)
            imshow(double(max_image(:,:,RNA_channel))./max(max(double(max_image(:,:,RNA_channel)))));
            title(['FISH signal: ',image_folder,', cycle = ',num2str(N_cycle)]);
        
            
%%% Output: %%%============================================================
        result_folder = [folder_list{list_I,1},out_folder];
        if exist(result_folder) ~= 7
            mkdir(result_folder);
        end
%        saveas(51,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),nu_add,seg_add,figure_tail]);
        saveas(2,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,figure_tail]);
%        saveas(3,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,th_add,figure_tail]);
%        saveas(4,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,seg_add,figure_tail]);
        saveas(5,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),fish_add,int_add,figure_tail]);
        saveas(6,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),fish_add,num_add,figure_tail]);
        saveas(7,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,num_add,figure_tail]);
        saveas(8,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,int_add,figure_tail]);
        saveas(9,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,int_add,cmp_add,figure_tail]);
%        saveas(10,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,th_add,sel_add,figure_tail]);
        saveas(11,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),nu_add,fate_add,figure_tail]);
        saveas(12,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),nu_add,reg_add,figure_tail]);
        saveas(13,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,reg_add,figure_tail]);
        saveas(14,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,fate_add,figure_tail]);
        saveas(15,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),nu_add,foci_add,reg_add,figure_tail]);
        saveas(16,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),immu_add,figure_tail]);
        saveas(17,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),fish_add,figure_tail]);
        save([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mat_tail],'image_folder','image_type','WGA_channel','DAPI_channel','RNA_channel','protein_channel','resolution','seg_bw','cyto_bw','max_image','foci_bw','nucleus_RNA_profile','foci_RNA_profile','cytoplasmic_RNA_profile','nucleus_protein_profile','cytoplasmic_protein_profile','N_cycle','psize0','g_threshold0','WGA_th0','DAPI_th0');
        
         if ~isempty(nucleus_RNA_profile)
             xlswrite([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),RNA_add,nu_add,output_tail],nucleus_RNA_profile);
         end
         if ~isempty(cytoplasmic_RNA_profile)
             xlswrite([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),RNA_add,cyto_add,output_tail],cytoplasmic_RNA_profile);
         end
         if ~isempty(foci_RNA_profile)
             xlswrite([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),RNA_add,foci_add,output_tail],foci_RNA_profile);
         end
         if ~isempty(nucleus_protein_profile)
             xlswrite([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,nu_add,output_tail],nucleus_protein_profile);
         end
         if ~isempty(cytoplasmic_protein_profile)
             xlswrite([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,cyto_add,output_tail],cytoplasmic_protein_profile);
         end
         sub_num(list_J,11) = N_cycle;
         
         clear ('seg_bw','cyto_bw','max_image','foci_bw','nucleus_RNA_profile','foci_RNA_profile','cytoplasmic_RNA_profile','nucleus_protein_profile','cytoplasmic_protein_profile','N_cycle','WGA_th0','DAPI_th0');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        
    end
    xlswrite([folder_list{list_I,1},in_folder,input_name],cat(2,sub_list,num2cell(sub_num)));
end