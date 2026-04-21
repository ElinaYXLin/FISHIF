clear all
close all

tic
%% Parameter setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
list_name = 'RNAlist.xls';
in_folder = 'stacks/';
input_name = 'matchlist.xls';
mismatch_name = 'mismatch.xls';
mismatch_name2 = 'mismatch_60X.xls';
xymismatch_name = 'xymismatch.mat';
xymismatch_name2 = 'xymismatch_60X.mat';
flip_list = {'Cad'};
image_type = '*.tif';
mask_folder = 'masks/';
mask_name = 'mask.mat';
out_folder = 'Results/';
% out_folder0 = 'Results_3Dz1/';
out_folder0 = '';
hist_folder = 'Histogram/';
fit_folder = 'Histogram_A/';
hist_folder2 = 'Histogram_nullo/';
fit_folder2 = 'Histogram_default/';
fit_add = '_spot_fit';
hist_tail = '_raw.xls';
N_thresh = 1;
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
D3_add = '_3D';
compare_add = '_compare';
DAPI_add = '_DAPI';
sub_pos = [3,3];
cycle_pos0 = zeros(1,sub_pos(1));
cycle_range = [11,12,13];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Data loading: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[num_list, folder_list] = xlsread(list_name);
folder_list = folder_list(strcmpi('T',folder_list(:,6)),:);
[N1,N2] = size(folder_list);
% run_list = {'[1]','[2]'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Image analysis: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for list_I = 1:N1
    [sub_num, sub_list] = xlsread([folder_list{list_I,1},in_folder,input_name]);
    [M1,M2] = size(sub_list);
%     if ~isempty(out_folder0)
%         copyfile([folder_list{list_I,1},out_folder],[folder_list{list_I,1},out_folder0]);
%     end
    
    for list_J = 1:M1%eval(run_list{list_I})
%         if isempty(strfind(sub_list{list_J,3},'_60X'))
%             [~,~,mismatch_matrix] = xlsread(mismatch_name);
%             load(xymismatch_name);
%         else
%             [~,~,mismatch_matrix] = xlsread(mismatch_name2);
%             load(xymismatch_name2);
%         end
        
        image_folder = [folder_list{list_I,1},in_folder,sub_list{list_J,3}];
        result_folder = [folder_list{list_I,1},out_folder];
        resolutionz = sub_num(list_J,11);
        Nbin = sub_num(list_J,2);
        Mdim = sub_num(list_J,3);
        load([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mat_tail]);

        load([folder_list{list_I,1},mask_folder,sub_list{list_J,3},mask_name]);   %%% load 3D mask
        RNA_channel = signal_channel;
        z_size = size(mask_stack,3);
        load([folder_list{list_I,1},fit_folder,sub_list{list_J,3}(1:end-1),fit_add,mat_tail]);   %%% load spot intensity fitting result
        Inten_thresh = b*N_thresh;   %%% set foci intensity threshold
        single_Inten = b;
        N_cycle = round(log2(max(mask_stack(:))))+2;   %%% Calculate the nuclear cycle number
        

        [~,RNA_stack,DAPI_stack] = stack3D(imclearborder(seg_bw),RNA_channel,DAPI_channel,RNA_channel,image_folder,0);   %%% load 3D image stacks
        [nucleus_DAPI_profile,DNA_mask] = DAPI_profile3D(mask_stack,DAPI_stack,image_folder,N_cycle);
        clear DAPI_stack
        
        
        [foci_bw3D,max_image00,SS] = modified_foci3D([folder_list{list_I,1},hist_folder,sub_list{list_J,3}(1:(end-1)),hist_tail],max_image,RNA_channel,mask_stack,N_cycle,image_folder,0,Inten_thresh);
        foci_bw2D = max(foci_bw3D,[],3);
        max_image00 = max_image00/single_Inten;

        [nucleus_RNA_profile,foci_RNA_profile,cytoplasmic_RNA_profile] = RNA_profile3D(mask_stack,nucleus_DAPI_profile,cyto_bw,foci_bw3D,max_image00,SS,RNA_stack,resolution,image_folder,N_cycle,[]);   %%% Input the foci area matrix from modified_foci.m if gaussian fitting method is used

        clear foci_bw3D max_image00 SS 
        clear RNA_stack
        
            
%%% Output: %%%============================================================
        result_folder = [folder_list{list_I,1},out_folder];
        if exist(result_folder) ~= 7
            mkdir(result_folder);
        end
        %saveas(51,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),nu_add,seg_add,figure_tail]);
%         saveas(2,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,D3_add,figure_tail]);
        %saveas(3,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,th_add,figure_tail]);
        saveas(4,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,seg_add,D3_add,figure_tail]);
        saveas(5,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),fish_add,int_add,D3_add,figure_tail]);
        saveas(6,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),fish_add,num_add,D3_add,figure_tail]);
        saveas(7,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,num_add,D3_add,figure_tail]);
        saveas(8,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,int_add,D3_add,figure_tail]);
        saveas(9,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,int_add,D3_add,cmp_add,figure_tail]);
        %saveas(10,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,th_add,sel_add,figure_tail]);
        saveas(11,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),nu_add,fate_add,D3_add,figure_tail]);

%         saveas(44,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,signal2_add,seg_add,D3_add,figure_tail]);
%         saveas(45,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),fish_add,signal2_add,int_add,D3_add,figure_tail]);
%         saveas(46,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),fish_add,signal2_add,num_add,D3_add,figure_tail]);
%         saveas(47,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,signal2_add,num_add,D3_add,figure_tail]);
%         saveas(48,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,signal2_add,int_add,D3_add,figure_tail]);
%         saveas(49,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,signal2_add,int_add,D3_add,cmp_add,figure_tail]);
%         saveas(411,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),nu_add,signal2_add,fate_add,D3_add,figure_tail]);
%         
%         saveas(12,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),nu_add,reg_add,D3_add,figure_tail]);
%         saveas(13,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,reg_add,D3_add,figure_tail]);
%         saveas(14,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,fate_add,D3_add,figure_tail]);
%         saveas(15,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),nu_add,foci_add,reg_add,D3_add,figure_tail]);
        saveas(16,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),nu_add,DAPI_add,D3_add,figure_tail]);
        
%         saveas(62,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,fluc_add,D3_add,figure_tail]);
%         saveas(63,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,bino_add,D3_add,figure_tail]);
%         
%         saveas(73,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,local_add,foci_add,D3_add,figure_tail]);
%         saveas(74,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,local_add,hist_add,D3_add,figure_tail]);
%         saveas(75,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,local_add,RNA_add,D3_add,figure_tail]);
%         saveas(76,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,local_add,abs_add,foci_add,D3_add,figure_tail]);
%         saveas(77,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,local_add,abs_add,hist_add,D3_add,figure_tail]);
%         saveas(80,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,local_add,rad_add,D3_add,figure_tail]);
%         saveas(81,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,local_add,abs_add,rad_add,D3_add,figure_tail]);
%         
%         saveas(473,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,local_add,foci_add,D3_add,signal2_add,figure_tail]);
%         saveas(474,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,local_add,hist_add,D3_add,signal2_add,figure_tail]);
%         saveas(475,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,local_add,RNA_add,D3_add,signal2_add,figure_tail]);
%         saveas(476,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,local_add,abs_add,foci_add,D3_add,signal2_add,figure_tail]);
%         saveas(477,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,local_add,abs_add,hist_add,D3_add,signal2_add,figure_tail]);
%         saveas(480,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,local_add,rad_add,D3_add,signal2_add,figure_tail]);
%         saveas(481,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,local_add,abs_add,rad_add,D3_add,signal2_add,figure_tail]);

        %saveas(576,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,local_add,abs_add,foci_add,D3_add,compare_add,figure_tail]);
        %saveas(581,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,local_add,abs_add,rad_add,D3_add,compare_add,figure_tail]);

        %save([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mat_tail],'image_folder','image_type','WGA_channel','DAPI_channel','RNA_channel','protein_channel','resolution','seg_bw','cyto_bw','max_image','foci_bw','nucleus_RNA_profile','foci_RNA_profile','cytoplasmic_RNA_profile','nucleus_protein_profile','cytoplasmic_protein_profile','N_cycle','foci_data','fake_data','h','t_absolute','-append');
        save([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mat_tail],'image_folder','image_type','DAPI_channel','RNA_channel','resolution','seg_bw','cyto_bw','max_image','foci_bw2D','nucleus_RNA_profile','foci_RNA_profile','cytoplasmic_RNA_profile','N_cycle','resolutionz','-append');
        
        if ~isempty(nucleus_RNA_profile)
           xlswrite([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),RNA_add,nu_add,output_tail],nucleus_RNA_profile);
        end
        if ~isempty(cytoplasmic_RNA_profile)
           xlswrite([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),RNA_add,cyto_add,output_tail],cytoplasmic_RNA_profile);
        end
        if ~isempty(foci_RNA_profile)
           xlswrite([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),RNA_add,foci_add,output_tail],foci_RNA_profile);
        end
%         if ~isempty(nucleus_protein_profile)
%            xlswrite([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,nu_add,output_tail],nucleus_protein_profile);
%         end
%         if ~isempty(cytoplasmic_protein_profile)
%            xlswrite([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,cyto_add,output_tail],cytoplasmic_protein_profile);
%         end
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
%         
%         if ~isempty(nucleus_signal2_profile)
%            xlswrite([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),signal2_add,nu_add,output_tail],nucleus_signal2_profile);
%         end
%         if ~isempty(cytoplasmic_signal2_profile)
%            xlswrite([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),signal2_add,cyto_add,output_tail],cytoplasmic_signal2_profile);
%         end
%         if ~isempty(foci_signal2_profile)
%            xlswrite([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),signal2_add,foci_add,output_tail],foci_signal2_profile);
%         end
%         for I_data = 1:length(foci_data2)
%             if ~isempty(foci_data2{I_data})
%                 xlswrite([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,local_add,signal2_add,foci_add,output_tail],foci_data2{I_data},I_data);
%             end
%         end
%         for I_data = 1:length(fake_data2)
%             if ~isempty(fake_data2{I_data})
%                 xlswrite([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,local_add,signal2_add,fake_add,foci_add,output_tail],fake_data2{I_data},I_data);
%             end
%         end
% 
        
        sub_num(list_J,13) = N_cycle;
        
        clear ('seg_bw','cyto_bw','max_image','foci_bw','nucleus_RNA_profile','foci_RNA_profile','cytoplasmic_RNA_profile','nucleus_protein_profile','cytoplasmic_protein_profile','N_cycle','foci_data','fake_data','foci_signal2_profile','cytoplasmic_signal2_profile','foci_data2','fake_data2');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        
    end
%     try
%         xlswrite([folder_list{list_I,1},in_folder,input_name],cat(2,sub_list,num2cell(sub_num)));
%     catch
% %         xlwrite([folder_list{list_I,1},in_folder,input_name],cat(2,sub_list,num2cell(sub_num)));
%     end
end
toc