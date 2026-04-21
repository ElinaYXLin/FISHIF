clear all
close all

tic
%% Parameter setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
list_name = 'Duallist.xls';
in_folder = 'stacks/';
input_name = 'matchlist.xls';
mismatch_name = 'mismatch.xls';
image_type = '*.tif';
mask_folder = 'masks/';
mask_name = 'mask.mat';
out_folder = 'Results/';
hist_folder = 'Histogram/';
fit_folder = 'Histogram_A/';
hist_folder2 = 'Histogram_nullo/';
fit_folder2 = 'Histogram_default/';
fit_add = '_spot_fit';
hist_tail = '_raw.xls';
N_thresh = 3;
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
signal2_add = '_nullo';
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
EL_add = '_EL';
D3_add = '_3D';
compare_add = '_compare';
DAPI_add = '_DAPI';
mean_add = '_mean';
slope_add = '_slope';
sub_pos = [3,3];
cycle_pos0 = zeros(1,sub_pos(1));
cycle_range = [11,12,13];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Data loading: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[num_list, folder_list] = xlsread(list_name);
folder_list = folder_list(strcmpi('T',folder_list(:,6)),:);
[N1,N2] = size(folder_list);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Image analysis: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for list_I = [2,10]%N1
    [sub_num, sub_list] = xlsread([folder_list{list_I,1},in_folder,input_name]);
    [M1,M2] = size(sub_list);
    [~,~,mismatch_matrix] = xlsread(mismatch_name);

    for list_J = 1:M1
        image_folder = [folder_list{list_I,1},in_folder,sub_list{list_J,3}];
        result_folder = [folder_list{list_I,1},out_folder];
        resolutionz = sub_num(list_J,11);
        signal2_channel = sub_num(list_J,12);
        load([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mat_tail]);

        load([folder_list{list_I,1},mask_folder,sub_list{list_J,3},mask_name]);   %%% load 3D mask
        load([folder_list{list_I,1},fit_folder,sub_list{list_J,3}(1:end-1),fit_add,mat_tail]);   %%% load spot intensity fitting result
        Inten_thresh = b*N_thresh;   %%% set foci intensity threshold
%         load([folder_list{list_I,1},fit_folder2,sub_list{list_J,3}(1:end-1),fit_add,mat_tail]);   %%% load spot intensity fitting result
%         Inten_thresh2 = b*N_thresh;   %%% set foci intensity threshold
%         N_cycle = round(log2(max(mask_stack(:))))+2;   %%% Calculate the nuclear cycle number
%         N_cycle = round(log2(max(mask_stack(:))))+3;   %%% Calculate the nuclear cycle number
        N_cycle = sub_num(list_J,13);

        all_color = eval(folder_list{list_I,5});
        protein_color = all_color{protein_channel};
        RNA_color = all_color{RNA_channel};
        protein_RNA_mismatch = mismatch_matrix{strcmp(protein_color,mismatch_matrix(1,:)),strcmp(RNA_color,mismatch_matrix(:,1))};
        signal2_color = all_color{signal2_channel};
        protein_signal2_mismatch = mismatch_matrix{strcmp(protein_color,mismatch_matrix(1,:)),strcmp(signal2_color,mismatch_matrix(:,1))};

        p = polyfit(nucleus_protein_profile(:,4),sqrt(nucleus_protein_profile(:,5)),1);
        k_DAPI = p(1);
        
        figure(111)
            imshow(double(max_image(:,:,DAPI_channel))/double(max(max(max_image(:,:,DAPI_channel)))))
            title([image_folder,': DAPI image'],'Interpreter','none')
        

%        [signal_stack,RNA_stack] = stack3D(imclearborder(seg_bw),protein_channel,DAPI_channel,RNA_channel,image_folder,protein_RNA_mismatch);   %%% load 3D image stacks
%        [nucleus_protein_profile,cytoplasmic_protein_profile,quanti_p] = protein_profile3D(max_image,protein_channel,mask_stack,signal_stack,image_folder,N_cycle,resolution);
        I_cycle = find(N_cycle == cycle_range);
        if N_cycle > max(cycle_range)
            I_cycle = find(max(cycle_range) == cycle_range);
        elseif N_cycle < min(cycle_range)
            I_cycle = find(min(cycle_range) == cycle_range);
        end
        if ~isempty(I_cycle)
            cycle_pos0(I_cycle) = cycle_pos0(I_cycle)+1;
            sub_pos(3) = sub2ind(sub_pos([1,2]),cycle_pos0(I_cycle),I_cycle);
        else
            sub_pos(3) = 0;
        end
        
        figure(1011)
            subplot(sub_pos(1),sub_pos(2),sub_pos(3))
            imshow(double(max_image(:,:,DAPI_channel))/double(max(max(max_image(:,:,DAPI_channel)))))
            title([image_folder,char(10),'DAPI image, cycle = ',num2str(N_cycle),', mean TX = ',num2str(mean(nucleus_RNA_profile(:,4))),', DAPI slope = ',num2str(k_DAPI)],'Interpreter','none')
       
        
        enrichment_compare(nucleus_protein_profile,nucleus_RNA_profile/b,foci_data,fake_data,r_size,foci_data2,fake_data2,signal2_add(2:end),image_folder,N_cycle,sub_pos,k_DAPI,true)
            
%%% Output: %%%============================================================
        result_folder = [folder_list{list_I,1},out_folder];
        if exist(result_folder) ~= 7
            mkdir(result_folder);
        end
        saveas(111,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),DAPI_add,figure_tail]);
%        saveas(14,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,fate_add,D3_add,figure_tail]);
        saveas(576,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,local_add,abs_add,foci_add,D3_add,compare_add,figure_tail]);
        saveas(581,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,local_add,abs_add,rad_add,D3_add,compare_add,figure_tail]);
        saveas(582,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,local_add,abs_add,EL_add,D3_add,compare_add,figure_tail]);
        
        saveas(5076,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,local_add,abs_add,foci_add,D3_add,compare_add,figure_tail]);
        saveas(5077,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,local_add,abs_add,foci_add,mean_add,D3_add,compare_add,figure_tail]);
        saveas(5078,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,local_add,abs_add,foci_add,slope_add,D3_add,compare_add,figure_tail]);
        saveas(5081,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,local_add,abs_add,rad_add,D3_add,compare_add,figure_tail]);
        saveas(5082,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,local_add,abs_add,EL_add,D3_add,compare_add,figure_tail]);

        
        clear ('seg_bw','cyto_bw','max_image','foci_bw','nucleus_RNA_profile','foci_RNA_profile','cytoplasmic_RNA_profile','nucleus_protein_profile','cytoplasmic_protein_profile','N_cycle','foci_data','fake_data','foci_signal2_profile','cytoplasmic_signal2_profile','foci_data2','fake_data2');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        
    end
    xlswrite([folder_list{list_I,1},in_folder,input_name],cat(2,sub_list,num2cell(sub_num)));
end
toc