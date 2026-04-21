clear all
close all

tic
%% Parameter setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
list_name = 'Duallist.xls';
in_folder = 'stacks/';
input_name = 'matchlist.xls';
mismatch_name = 'mismatch.xls';
image_type = '*.tif';
out_folder = 'Results/';
mask_folder = 'masks/';
mask_name = 'mask.mat';
hist_folder = 'Histogram/';
fit_folder = 'Histogram_A/';
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
D3_add = '_3D';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Data loading: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[num_list, folder_list] = xlsread(list_name);
folder_list = folder_list(strcmpi('T',folder_list(:,6)),:);
[N1,N2] = size(folder_list);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Image analysis: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for list_I = 1:1%N1
    [sub_num, sub_list] = xlsread([folder_list{list_I,1},in_folder,input_name]);
    [M1,M2] = size(sub_list);
    [~,~,mismatch_matrix] = xlsread(mismatch_name);

    for list_J = 1:1%1:M1
        image_folder = [folder_list{list_I,1},in_folder,sub_list{list_J,3}];
        result_folder = [folder_list{list_I,1},out_folder];
        resolutionz = sub_num(list_J,11);
        signal2_channel = sub_num(list_J,12);
        load([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mat_tail]);

        load([folder_list{list_I,1},mask_folder,sub_list{list_J,3},mask_name]);   %%% load 3D mask
        load([folder_list{list_I,1},fit_folder,sub_list{list_J,3}(1:end-1),fit_add,mat_tail]);   %%% load spot intensity fitting result
        Inten_thresh = b*N_thresh;   %%% set foci intensity threshold
        N_cycle = round(log2(max(mask_stack(:))))+2;   %%% Calculate the nuclear cycle number
        
        all_color = eval(folder_list{list_I,5});
        protein_color = all_color{protein_channel};
        RNA_color = all_color{RNA_channel};
        signal2_color = all_color{signal2_channel};
        protein_RNA_mismatch = mismatch_matrix{strcmp(protein_color,mismatch_matrix(1,:)),strcmp(RNA_color,mismatch_matrix(:,1))};

        [signal_stack,RNA_stack] = stack3D(imclearborder(seg_bw),protein_channel,DAPI_channel,RNA_channel,image_folder,protein_RNA_mismatch);   %%% load 3D image stacks
        
        [foci_bw3D,max_image00,SS] = modified_foci3D([folder_list{list_I,1},hist_folder,sub_list{list_J,3}(1:(end-1)),hist_tail],max_image,RNA_channel,mask_stack,N_cycle,image_folder,protein_RNA_mismatch,Inten_thresh);
        [nucleus_RNA_profile,foci_RNA_profile,cytoplasmic_RNA_profile] = RNA_profile3D(mask_stack,cyto_bw,foci_bw3D,max_image00,SS,RNA_stack,resolution,image_folder,N_cycle);   %%% Input the foci area matrix from modified_foci.m if gaussian fitting method is used

        [nucleus_protein_profile,cytoplasmic_protein_profile,quanti_p] = protein_profile3D(max_image,protein_channel,mask_stack,signal_stack,image_folder,N_cycle,resolution);
        signal_stack = (signal_stack+quanti_p(2))/quanti_p(1);
        nucleus_protein_profile_ab = [nucleus_protein_profile(:,1),(nucleus_protein_profile(:,2)+quanti_p(2))/quanti_p(1)];

        foci_bw00 = find(foci_bw3D);
        max_image00 = max_image00.*SS;
        max_image_list = max_image00(foci_bw00);
        clear foci_bw3D max_image00 SS 
        
        [foci_data,fake_data,h,t_absolute,r_size] = dual_local_local3D2(foci_bw00,max_image_list,mask_stack,signal_stack,RNA_stack,nucleus_protein_profile_ab,image_folder,N_cycle,resolution,resolutionz);
        clear mask_stack signal_stack RNA_stack nucleus_protein_profile_ab

        dual_profile(nucleus_protein_profile,nucleus_RNA_profile,foci_RNA_profile,image_folder,N_cycle);
        
            
%%% Output: %%%============================================================
        result_folder = [folder_list{list_I,1},out_folder];
        if exist(result_folder) ~= 7
            mkdir(result_folder);
        end
        %saveas(51,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),nu_add,seg_add,figure_tail]);
        saveas(2,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,D3_add,figure_tail]);
        %saveas(3,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,th_add,figure_tail]);
        saveas(4,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,seg_add,D3_add,figure_tail]);
        saveas(5,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),fish_add,int_add,D3_add,figure_tail]);
        saveas(6,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),fish_add,num_add,D3_add,figure_tail]);
        saveas(7,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,num_add,D3_add,figure_tail]);
        saveas(8,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,int_add,D3_add,figure_tail]);
        saveas(9,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,int_add,D3_add,cmp_add,figure_tail]);
        %saveas(10,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,th_add,sel_add,figure_tail]);
        saveas(11,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),nu_add,fate_add,D3_add,figure_tail]);
        saveas(12,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),nu_add,reg_add,D3_add,figure_tail]);
        saveas(13,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,reg_add,D3_add,figure_tail]);
        saveas(14,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,fate_add,D3_add,figure_tail]);
        saveas(15,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),nu_add,foci_add,reg_add,D3_add,figure_tail]);
        saveas(62,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,fluc_add,D3_add,figure_tail]);
        saveas(63,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,bino_add,D3_add,figure_tail]);
        saveas(73,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,local_add,foci_add,D3_add,figure_tail]);
        saveas(74,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,local_add,hist_add,D3_add,figure_tail]);
        saveas(75,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,local_add,RNA_add,D3_add,figure_tail]);
        saveas(76,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,local_add,abs_add,foci_add,D3_add,figure_tail]);
        saveas(77,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,local_add,abs_add,hist_add,D3_add,figure_tail]);
        saveas(80,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,local_add,rad_add,D3_add,figure_tail]);
        saveas(81,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,local_add,abs_add,rad_add,D3_add,figure_tail]);
        
        %save([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mat_tail],'image_folder','image_type','WGA_channel','DAPI_channel','RNA_channel','protein_channel','resolution','seg_bw','cyto_bw','max_image','foci_bw','nucleus_RNA_profile','foci_RNA_profile','cytoplasmic_RNA_profile','nucleus_protein_profile','cytoplasmic_protein_profile','N_cycle','foci_data','fake_data','h','t_absolute','-append');
        save([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mat_tail],'image_folder','image_type','WGA_channel','DAPI_channel','RNA_channel','protein_channel','resolution','seg_bw','cyto_bw','max_image','nucleus_RNA_profile','foci_RNA_profile','cytoplasmic_RNA_profile','nucleus_protein_profile','cytoplasmic_protein_profile','N_cycle','foci_data','fake_data','h','t_absolute','r_size','resolutionz','quanti_p','-append');
        
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
        for I_data = 1:length(foci_data)
            if ~isempty(foci_data{I_data})
                xlswrite([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,local_add,foci_add,output_tail],foci_data{I_data},I_data);
            end
        end
        for I_data = 1:length(fake_data)
            if ~isempty(fake_data{I_data})
                xlswrite([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,local_add,fake_add,foci_add,output_tail],fake_data{I_data},I_data);
            end
        end
       
        sub_num(list_J,13) = N_cycle;
        
        clear ('seg_bw','cyto_bw','max_image','foci_bw','nucleus_RNA_profile','foci_RNA_profile','cytoplasmic_RNA_profile','nucleus_protein_profile','cytoplasmic_protein_profile','N_cycle','foci_data','fake_data');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        
    end
    xlswrite([folder_list{list_I,1},in_folder,input_name],cat(2,sub_list,num2cell(sub_num)));
end
toc