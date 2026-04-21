clear

%% Parameter setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
list_name = 'RNAlist.xls';
in_folder = 'stacks/';
input_name = 'matchlist.xls';
gene_table = 'mRNA_intensity.xls';
image_type = '*.tif';
out_folder = 'Results/';
output_tail = '.xls';
figure_tail = '.fig';
mat_tail = '.mat';
hist_folder = 'Histogram/';
hist_tail = '_raw.xls';
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
probe_name = '_hb';
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
        RNA_channel = signal_channel;
        N_cycle = round(log2(max(max(bwlabel(seg_bw)))))+2;   %%% Calculate the nuclear cycle number
        %[seg_bw,cyto_bw,max_image,N_cycle,WGA_th0,DAPI_th0] = nuclei_seg4(image_folder,image_type,WGA_channel,DAPI_channel,resolution);
        %[foci_bw,psize0,g_threshold0] = RNA_seg(seg_bw,max_image,signal_channel,WGA_channel,resolution,image_folder,N_cycle);
        %[nucleus_profile,foci_profile,cytoplasmic_profile] = RNA_profile(seg_bw,cyto_bw,foci_bw,max_image,signal_channel,resolution,image_folder,N_cycle);
        [foci_bw00,max_image00,SS] = modified_foci([folder_list{list_I,1},hist_folder,sub_list{list_J,3}(1:(end-1)),hist_tail],max_image,RNA_channel,seg_bw,N_cycle,image_folder);
        
        N_start = strfind(sub_list{list_J,3},probe_name)+1;
        N_end = strfind(sub_list{list_J,3},'_');
        N_site = sub_list{list_J,3}(N_start(1):(N_end(find(N_end > N_start(1),1))-1));
        [foci_list,gene_list] = xlsread(gene_table);
        I_gene = strcmp(N_site, gene_list);
        Ifoci0 = foci_list(I_gene,:);
        
        [nucleus_RNA_profile,foci_RNA_profile,cytoplasmic_RNA_profile] = RNA_profile(seg_bw,cyto_bw,foci_bw00,max_image00,RNA_channel,resolution,image_folder,N_cycle,SS,Ifoci0);   %%% Input the foci area matrix from modified_foci.m if gaussian fitting method is used
        clear foci_bw00 max_image00 SS
 
        
        figure(17)
            imshow(double(max_image(:,:,RNA_channel))./max(max(double(max_image(:,:,RNA_channel)))));
            title(['FISH signal: ',image_folder,', cycle = ',num2str(N_cycle)]);
        
%%% Output: %%%============================================================
        result_folder = [folder_list{list_I,1},out_folder];
        if exist(result_folder) ~= 7
            mkdir(result_folder);
        end
        %saveas(51,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),nu_add,seg_add,figure_tail]);
        %saveas(3,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,th_add,figure_tail]);
        saveas(4,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,seg_add,figure_tail]);
        saveas(5,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),fish_add,int_add,figure_tail]);
        saveas(6,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),fish_add,num_add,figure_tail]);
        saveas(7,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,num_add,figure_tail]);
        saveas(8,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,int_add,figure_tail]);
        saveas(9,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,int_add,cmp_add,figure_tail]);
        %saveas(10,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,th_add,sel_add,figure_tail]);
        saveas(11,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),nu_add,fate_add,figure_tail]);
        saveas(17,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),fish_add,figure_tail]);
        save([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mat_tail],'image_folder','image_type','WGA_channel','DAPI_channel','signal_channel','resolution','seg_bw','cyto_bw','max_image','foci_bw','nucleus_profile','foci_profile','cytoplasmic_profile','N_cycle','psize0','g_threshold0','WGA_th0','DAPI_th0');
        
        if ~isempty(nucleus_profile)
            xlswrite([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),nu_add,output_tail],nucleus_profile);
        end
        if ~isempty(cytoplasmic_profile)
            xlswrite([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),cyto_add,output_tail],cytoplasmic_profile);
        end
        if ~isempty(foci_profile)
            xlswrite([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,output_tail],foci_profile);
        end
        sub_num(list_J,11) = N_cycle;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        
    end
    xlswrite([folder_list{list_I,1},in_folder,input_name],cat(2,sub_list,num2cell(sub_num)));
end