clear all
close all

%% Parameter setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
list_name = 'Duallist.xls';
in_folder = 'stacks/';
input_name = 'matchlist.xls';
standard_record = 'Calibration/Results/standard.mat';
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
fluc_add = '_fluc';
bino_add = '_bino';
emmask_add = '_emmask';
double_add = 'D';
% size0 = 5000;
bin0 = 3;
Smin = 100;
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
    channel_name = eval(folder_list{list_I,5});
    for list_J = 1:M1
        image_folder = [folder_list{list_I,1},in_folder,sub_list{list_J,3}];
        
        if ~isempty(strfind(sub_list{list_J,2},double_add))
            Nbin = sub_num(list_J,2:3);
            Mdim = 1:2;
        else
            Mdim = sub_num(list_J,3);
            Nbin = sub_num(list_J,2);
        end
        
        WGA_channel = sub_num(list_J,6);
        DAPI_channel = sub_num(list_J,7);
        RNA_channel = sub_num(list_J,10);
        protein_channel = sub_num(list_J,8);
        signal2_channel = sub_num(list_J,12);
        resolution = sub_num(list_J,9);
        resolutionz = sub_num(list_J,11);
        all_color = eval(folder_list{list_I,5});


        temp_name = dir([image_folder,image_type]);
        temp0 = imread([image_folder,temp_name(1).name]);
% %         if size0 > 0
% %             scale0 = min(size0/max(size(temp0)),1);
% %         else
% %             scale0 = 1;
% %         end
        if bin0 > 1
            scale0 = 1/bin0;
        else
            scale0 = 1;
        end
        clear temp_name temp0
        
        %[seg_bw,cyto_bw,max_image,N_cycle,WGA_th0,DAPI_th0] = nuclei_seg4(image_folder,image_type,WGA_channel,DAPI_channel,resolution);
        [seg_bw,cyto_bw,max_image,N_cycle,WGA_th0,DAPI_th0] = nuclei_segDAPI(image_folder,image_type,DAPI_channel,resolution,scale0,Smin);
        
        result_folder = [folder_list{list_I,1},out_folder];
        save([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mat_tail],'max_image','-append');
        clear ('em_mask','seg_bw','cyto_bw','max_image','foci_bw','nucleus_RNA_profile','foci_RNA_profile','cytoplasmic_RNA_profile','nucleus_protein_profile','cytoplasmic_protein_profile','N_cycle','WGA_th0','DAPI_th0');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        
    end
end

close all

% clear global standard_data