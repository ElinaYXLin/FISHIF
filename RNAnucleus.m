clear

%% Parameter setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
list_name = 'RNAlist.xls';
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
nu_add = '_nucleus';
cyto_add = '_cytoplasm';
fate_add = '_fate';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Data loading: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[num_list, folder_list] = xlsread(list_name);
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
        signal_channel = sub_num(list_J,8);
        resolution = sub_num(list_J,9);
        [seg_bw,cyto_bw,max_image] = nuclei_seg4(image_folder,image_type,WGA_channel,DAPI_channel,resolution);
        %foci_bw = RNA_seg(seg_bw,max_image,signal_channel,WGA_channel,resolution,image_folder);
        %[nucleus_profile,foci_profile,cytoplasmic_profile] = RNA_profile(seg_bw,cyto_bw,foci_bw,max_image,signal_channel,resolution,image_folder);
        
%%% Output: %%%============================================================
        result_folder = [folder_list{list_I,1},out_folder];
        if exist(result_folder) ~= 7
            mkdir(result_folder);
        end
        saveas(1,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),nu_add,seg_add,figure_tail]);
        %saveas(3,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,th_add,figure_tail]);
        %saveas(4,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,seg_add,figure_tail]);
        %saveas(5,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),fish_add,int_add,figure_tail]);
        %saveas(6,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),fish_add,num_add,figure_tail]);
        %saveas(7,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,num_add,figure_tail]);
        %saveas(8,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,int_add,figure_tail]);
        %saveas(9,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,int_add,cmp_add,figure_tail]);
        %saveas(10,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,th_add,sel_add,figure_tail]);
        %saveas(11,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),nu_add,fate_add,figure_tail]);
        save([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mat_tail],'image_folder','image_type','WGA_channel','DAPI_channel','signal_channel','resolution','seg_bw','cyto_bw','max_image');
        
        %if ~isempty(nucleus_profile)
        %    xlswrite([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),nu_add,output_tail],nucleus_profile);
        %end
        %if ~isempty(cytoplasmic_profile)
        %    xlswrite([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),cyto_add,output_tail],cytoplasmic_profile);
        %end
        %if ~isempty(foci_profile)
        %    xlswrite([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,output_tail],foci_profile);
        %end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        
    end
end