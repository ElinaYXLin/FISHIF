clear all
close all

tic
xy_mismatch_check = true;
use_linear = true; %%% whether to use linear representation for max_image00 and SS

%% Parameter setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
list_name = 'Duallist.xls';
in_folder = 'stacks/';
% old_add = '_old';
old_add = '';
input_name = 'matchlist.xls';

mismatch_name = 'mismatch.xls';
mismatch_name2 = 'mismatch_60X.xls';
% mismatch_name3 = 'mismatch_60X_FA_02142020.xls';
mismatch_name3 = 'mismatch_60X_FA.xls';
xymismatch_name = 'xymismatch.mat';
xymismatch_name2 = 'xymismatch_60X.mat';
xymismatch_name3 = 'xymismatch_60X_FA.mat';
fast_add = 'F';
airy_add = 'A';
double_add = 'D';
stitch_name = 'stitch_parameter.mat';

flip_list = {'Cad'};
% flip0 = [];
flip0 = false;
image_type = '*.tif';
mask_folder = 'masks/';
mask_name = 'mask.mat';
out_folder = 'Results2/';
% out_folder0 = 'Results_3Dz1/';
% out_folder0 = 'Results_decross/';
% out_folder0 = 'Results_original/';
% out_folder0 = 'Results_noalignment/';
out_folder0 = '';
hist_folder = 'Histogram_alignment/';
fit_folder = 'Histogram_A/';
hist_folder2 = 'Histogram_alignment_RNA2/';
fit_folder2 = 'Histogram_A_RNA2/';
fit_folder_protein = 'Histogram_protein2_A/';
fit_add = '_spot_fit';
hist_tail = '_raw.xls';
N_thresh = 0;
N_thresh2 = 0;
EL_check_name = 'RNA_stack';
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
signal2_add = '_RNA2';
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
noise_add = '_noise';
sub_pos = [3,3];
cycle_pos0 = zeros(1,sub_pos(1));
cycle_range = [11,12,13];
% EL_range = [];
% % EL_range = [0,0.35,0.8,1,0.7,0.8];  %%% Bcd
EL_range = [0.2,0.55,0.7,0.9,0.7,0.9];  %%% Hb
% % EL_range = [0.3,0.7,0,0.1,0,0.1];  %%% Cad
% EL_range = [0.1,0.9,1,1,1,1];  %%% Histone
% % sigmaxz = [1.35,1.3];
sigmaxz = [1.6,1.5];
I_qp = 5; %%% type of protein quantification: 1. fluctuation method; 5. spot method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Data loading: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[num_list, folder_list] = xlsread(list_name);
% folder_list = folder_list(strcmpi('b',folder_list(:,6)) | strcmpi('ba',folder_list(:,6)),:);
% folder_list = folder_list(strcmpi('b3',folder_list(:,6)),:);
% folder_list = folder_list(strcmpi('h',folder_list(:,6)),:);
% folder_list = folder_list(strcmpi('b3n',folder_list(:,6)),:);
folder_list = folder_list(strcmpi('T',folder_list(:,6)),:);
% folder_list = folder_list(strcmpi('cgal4new',folder_list(:,6)),:);
% folder_list = folder_list(strcmpi('h4k5ac',folder_list(:,6)),:);
% folder_list = folder_list(strcmpi('cbgal4',folder_list(:,6)),:);
% folder_list = folder_list(strcmpi('bgal4',folder_list(:,6)),:);
[N1,N2] = size(folder_list);
run_list = {[],5};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Image analysis: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for list_I = 1:N1
    [sub_num, sub_list] = xlsread([folder_list{list_I,1},in_folder,input_name]);
    [M1,M2] = size(sub_list);
    if ~isempty(out_folder0)
        copyfile([folder_list{list_I,1},out_folder],[folder_list{list_I,1},out_folder0]);
    end
    
    for list_J = run_list{list_I}
        if isempty(strfind(sub_list{list_J,3},'_60X')) && isempty(strfind(sub_list{list_J,2},'60X')) && isempty(strfind(sub_list{list_J,2},airy_add)) && isempty(strfind(sub_list{list_J,2},fast_add))
            [~,~,mismatch_matrix] = xlsread(mismatch_name);
            load(xymismatch_name)
            mismatch_name
        elseif (~isempty(strfind(sub_list{list_J,3},'_60X')) || ~isempty(strfind(sub_list{list_J,2},'60X'))) && isempty(strfind(sub_list{list_J,2},airy_add)) && isempty(strfind(sub_list{list_J,2},fast_add))
            [~,~,mismatch_matrix] = xlsread(mismatch_name2);
            load(xymismatch_name2)
            mismatch_name2
        elseif (~isempty(strfind(sub_list{list_J,3},'_60X')) || ~isempty(strfind(sub_list{list_J,2},'60X'))) && ~isempty(strfind(sub_list{list_J,2},fast_add))
            [~,~,mismatch_matrix] = xlsread(mismatch_name3);
            load(xymismatch_name3)
            mismatch_name3
        else
            error('No mismatch information')
        end
        
%         image_folder = [folder_list{list_I,1},in_folder,sub_list{list_J,3}];
        image_folder = [folder_list{list_I,1},in_folder,sub_list{list_J,3}(1:end-1),old_add,'/'];
        
        if isempty(flip0)
            flip_axis = any(cellfun(@(x) ~isempty(strfind(image_folder,x)),flip_list));   %%% Check whether the protein profile has an opposite direction (to Bcd)
        else
            flip_axis = flip0;
        end
        
        result_folder = [folder_list{list_I,1},out_folder];
        
        load([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mat_tail]);
        image_folder = [folder_list{list_I,1},in_folder,sub_list{list_J,3}(1:end-1),old_add,'/'];
        
        resolutionz = sub_num(list_J,11);
        signal2_channel = sub_num(list_J,12);
        protein_channel = sub_num(list_J,14);

        if ~isempty(strfind(sub_list{list_J,2},double_add))
            Nbin = sub_num(list_J,2:3);
            Mdim = 1:2;
        else
            Nbin = ones(1,2);
            Mdim = sub_num(list_J,3);
            Nbin(Mdim) = sub_num(list_J,2);
            Mdim = 1:2;
        end
        
        %%% Stitch parameter loading: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [x_start_im,y_start_im,x_center_im,y_center_im] = tile_info_extract(image_folder,Nbin,max_image);
        %%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        load([folder_list{list_I,1},mask_folder,sub_list{list_J,3},mask_name]);   %%% load 3D mask
        mask_stack = double(mask_stack);
        z_size = size(mask_stack,3);
        load([folder_list{list_I,1},fit_folder,sub_list{list_J,3}(1:end-1),fit_add,mat_tail]);   %%% load spot intensity fitting result
% %         b = 1;
        Inten_thresh = b*N_thresh;   %%% set foci intensity threshold
        single_Inten = b;
        load([folder_list{list_I,1},fit_folder2,sub_list{list_J,3}(1:end-1),fit_add,mat_tail]);   %%% load spot intensity fitting result
        Inten_thresh2 = b*N_thresh2;   %%% set foci intensity threshold
        single_Inten2 = b;
        
        load([folder_list{list_I,1},fit_folder_protein,sub_list{list_J,3}(1:end-1),fit_add,mat_tail]);   %%% load spot intensity fitting result
        Inten_protein = b;
        
        N_cycle = round(log2(max(mask_stack(:))))+2;   %%% Calculate the nuclear cycle number
%         N_cycle = sub_num(list_J,13);

        all_color = eval(folder_list{list_I,5});
        protein_color = all_color{protein_channel};
        RNA_color = all_color{RNA_channel};
        protein_RNA_mismatch = mismatch_matrix{strcmp(protein_color,mismatch_matrix(1,:)),strcmp(RNA_color,mismatch_matrix(:,1))};
        if xy_mismatch_check
            protein_RNA_xymismatch = {eval([protein_color,'_',RNA_color]),eval([protein_color,'_',RNA_color,'_con']),eval([protein_color,'_',RNA_color,'_x0']),Nbin,Mdim,resolution,resolution_mismatch};
        else
            protein_RNA_xymismatch = [];
        end
        signal2_color = all_color{signal2_channel};
        protein_signal2_mismatch = mismatch_matrix{strcmp(protein_color,mismatch_matrix(1,:)),strcmp(signal2_color,mismatch_matrix(:,1))};
        if xy_mismatch_check
            protein_signal2_xymismatch = {eval([protein_color,'_',signal2_color]),eval([protein_color,'_',signal2_color,'_con']),eval([protein_color,'_',signal2_color,'_x0']),Nbin,Mdim,resolution,resolution_mismatch};
        else
            protein_signal2_xymismatch = [];
        end
        
        [signal_stack,RNA_stack,DAPI_stack] = stack3D(imclearborder(seg_bw),protein_channel,DAPI_channel,RNA_channel,image_folder,protein_RNA_mismatch);   %%% load 3D image stacks
        [nucleus_DAPI_profile,DNA_mask] = DAPI_profile3D(mask_stack,DAPI_stack,image_folder,N_cycle,resolution);
        clear DAPI_stack
        
        
        [foci_bw3D,max_image00,SS] = modified_foci3D([folder_list{list_I,1},hist_folder,sub_list{list_J,3}(1:(end-1)),hist_tail],max_image,RNA_channel,mask_stack,N_cycle,image_folder,protein_RNA_mismatch,Inten_thresh,[],resolution,use_linear);
        foci_bw2D = max(foci_bw3D,[],3);
        max_image00 = max_image00/single_Inten;
        EL_info = get_EL(em_mask);   %%% get EL information (extreme points, EL length) from the embryo mask
% %         flip_EL = EL_orientation(nucleus_DAPI_profile,EL_info,mask_stack,eval(EL_check_name),flip_axis,[],resolution);
        flip_EL = EL_orientation2(EL_info,mask_stack,[folder_list{list_I,1},hist_folder,sub_list{list_J,3}(1:(end-1)),hist_tail],protein_RNA_mismatch,flip_axis,[]);
        [nucleus_protein_profile,cytoplasmic_protein_profile,quanti_p,nucleus_protein_ab] = protein_profile3D2(nucleus_DAPI_profile,max_image,EL_info,protein_channel,mask_stack,signal_stack,image_folder,N_cycle,resolution,resolutionz,flip_EL,[],[],[],[],[],sigmaxz,flip_axis,EL_range,Inten_protein);
        
%         [nucleus_RNA_profile,foci_RNA_profile,cytoplasmic_RNA_profile] = RNA_profile3D(mask_stack,nucleus_DAPI_profile,cyto_bw,foci_bw3D,max_image00,EL_info,SS,RNA_stack,resolution,image_folder,N_cycle,[],flip_EL);   %%% Input the foci area matrix from modified_foci.m if gaussian fitting method is used
        [nucleus_RNA_profile,foci_RNA_profile,cytoplasmic_RNA_profile] = RNA_profile3D(mask_stack,nucleus_DAPI_profile,em_mask,foci_bw3D,max_image00,EL_info,SS,RNA_stack,resolution,image_folder,N_cycle,[],flip_EL,use_linear);   %%% Input the foci area matrix from modified_foci.m if gaussian fitting method is used

%         signal_stack = (signal_stack+quanti_p(2))/quanti_p(1);
%         nucleus_protein_profile_ab = [nucleus_protein_profile(:,1),(nucleus_protein_profile(:,2)+quanti_p(2))/quanti_p(1),nucleus_DAPI_profile];
%        nuclei_variance = [0;nucleus_protein_profile(:,3)];
%        variance_stack = nuclei_variance(mask_stack+1);
        nucleus_protein_profile_ab = [nucleus_protein_profile(:,1),nucleus_protein_ab,nucleus_protein_profile(:,4),nucleus_protein_profile(:,4),nucleus_protein_ab*quanti_p(1)/quanti_p(5)];
%         pro_temp = [0;nucleus_protein_ab];
%         signal_stack = pro_temp(mask_stack+1);
% % %         signal_stack = (double(signal_stack)+quanti_p(2))/quanti_p(1);
        signal_stack = (double(signal_stack)+quanti_p(2))/quanti_p(I_qp);
%         signal_stack = signal_stack/quanti_p(1);
        nucleus_protein_profile = [nucleus_protein_profile,nucleus_DAPI_profile];
        RNA_DNA(nucleus_protein_profile(:,5),nucleus_RNA_profile(:,[3,4]),image_folder,N_cycle,100);
        TXnoise_plot(nucleus_RNA_profile,N_cycle,image_folder,[101,102]);
        
        [foci_bw3D,max_image00,SS] = modified_foci3D2([folder_list{list_I,1},hist_folder,sub_list{list_J,3}(1:(end-1)),hist_tail],max_image,RNA_channel,mask_stack,N_cycle,image_folder,protein_RNA_mismatch,protein_RNA_xymismatch,Inten_thresh,[],resolution,{x_start_im,y_start_im,x_center_im,y_center_im},use_linear);
        max_image00 = max_image00/single_Inten;
        foci_bw00 = find(foci_bw3D);
        max_image00 = max_image00.*SS;
        if ~use_linear
            max_image_list = max_image00(foci_bw00);
        else
            max_image_list = max_image00;
        end
        foci_mask = logical(max(foci_bw3D,[],3));
        clear foci_bw3D max_image00 SS 
        [foci_data,fake_data,h,t_absolute,r_size] = dual_local_local3D3(foci_bw00,max_image_list,mask_stack,DNA_mask,signal_stack,RNA_stack,nucleus_protein_profile_ab,image_folder,N_cycle,resolution,resolutionz,[],[],sigmaxz);
%         clear RNA_stack
        
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        [~,signal2_stack] = stack3D(imclearborder(seg_bw),protein_channel,DAPI_channel,signal2_channel,image_folder,protein_signal2_mismatch);   %%% load 3D image stacks

        [foci_bw3D2,max_image002,SS2] = modified_foci3D([folder_list{list_I,1},hist_folder2,sub_list{list_J,3}(1:(end-1)),hist_tail],max_image,signal2_channel,mask_stack,N_cycle,image_folder,protein_signal2_mismatch,Inten_thresh2,44,resolution,use_linear);
        max_image002 = max_image002/single_Inten2;
        [nucleus_signal2_profile,foci_signal2_profile,cytoplasmic_signal2_profile] = RNA_profile3D(mask_stack,nucleus_DAPI_profile,em_mask,foci_bw3D2,max_image002,EL_info,SS2,signal2_stack,resolution,image_folder,N_cycle,[45,46,47,48,49,411],flip_EL,use_linear);   %%% Input the foci area matrix from modified_foci.m if gaussian fitting method is used
        nucleus_signal2_profile(:,1) = nucleus_protein_profile(:,1);
        
        [foci_bw3D2,max_image002,SS2] = modified_foci3D2([folder_list{list_I,1},hist_folder2,sub_list{list_J,3}(1:(end-1)),hist_tail],max_image,signal2_channel,mask_stack,N_cycle,image_folder,protein_signal2_mismatch,protein_signal2_xymismatch,Inten_thresh2,44,resolution,{x_start_im,y_start_im,x_center_im,y_center_im},use_linear);
        max_image002 = max_image002/single_Inten2;
        foci_bw002 = find(foci_bw3D2);
        max_image002 = max_image002.*SS2;
        if ~use_linear
            max_image_list2 = max_image002(foci_bw002);
        else
            max_image_list2 = max_image002;
        end
        clear foci_bw3D2 max_image002 SS2 
%         [foci_data2,fake_data2,h2,t_absolute2,r_size2] = dual_local_local3D3(foci_bw002,max_image_list2,mask_stack,DNA_mask,signal_stack,signal2_stack,nucleus_protein_profile_ab,image_folder,N_cycle,resolution,resolutionz,foci_mask,[473:477,480,481]);
%         [foci_data2,fake_data2,h2,t_absolute2,r_size2] = dual_local_local3D3(foci_bw002,max_image_list2,mask_stack,DNA_mask,signal_stack,signal2_stack,nucleus_protein_profile_ab,image_folder,N_cycle,resolution,resolutionz,[],[473:477,480,481]);
        [foci_data2,fake_data2,h2,t_absolute2,r_size2] = dual_local_local3D3(foci_bw002,max_image_list2,mask_stack,DNA_mask,signal_stack,signal2_stack,nucleus_protein_profile_ab,image_folder,N_cycle,resolution,resolutionz,[],[473:477,480,481],sigmaxz,[],RNA_stack);
        clear mask_stack signal_stack signal2_stack RNA_stack% nucleus_protein_profile_ab

        dual_profile(nucleus_protein_profile_ab,nucleus_RNA_profile,foci_RNA_profile,image_folder,N_cycle,z_size,{'M','#'});
        %enrichment_compare(foci_data,fake_data,r_size,foci_data2,fake_data2,signal2_add(2:end),image_folder,N_cycle,true)
            
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
% %         saveas(11,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),nu_add,fate_add,D3_add,figure_tail]);
        hgsave(11,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),nu_add,fate_add,D3_add,figure_tail],'-v7.3');

        saveas(44,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,signal2_add,seg_add,D3_add,figure_tail]);
        saveas(45,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),fish_add,signal2_add,int_add,D3_add,figure_tail]);
        saveas(46,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),fish_add,signal2_add,num_add,D3_add,figure_tail]);
        saveas(47,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,signal2_add,num_add,D3_add,figure_tail]);
        saveas(48,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,signal2_add,int_add,D3_add,figure_tail]);
        saveas(49,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,signal2_add,int_add,D3_add,cmp_add,figure_tail]);
% %         saveas(411,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),nu_add,signal2_add,fate_add,D3_add,figure_tail]);
        hgsave(411,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),nu_add,signal2_add,fate_add,D3_add,figure_tail],'-v7.3');
        
        saveas(12,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),nu_add,reg_add,D3_add,figure_tail]);
        saveas(13,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,reg_add,D3_add,figure_tail]);
% %         saveas(14,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,fate_add,D3_add,figure_tail]);
        hgsave(14,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,fate_add,D3_add,figure_tail],'-v7.3');
        saveas(15,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),nu_add,foci_add,reg_add,D3_add,figure_tail]);
        saveas(16,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),nu_add,DAPI_add,D3_add,figure_tail]);
        
        saveas(62,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,fluc_add,D3_add,figure_tail]);
        saveas(63,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,bino_add,D3_add,figure_tail]);
        
        saveas(100,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),DAPI_add,RNA_add,figure_tail]);
        saveas(101,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),RNA_add,noise_add,figure_tail]);

        saveas(73,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,local_add,foci_add,D3_add,figure_tail]);
        saveas(74,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,local_add,hist_add,D3_add,figure_tail]);
        saveas(75,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,local_add,RNA_add,D3_add,figure_tail]);
        saveas(76,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,local_add,abs_add,foci_add,D3_add,figure_tail]);
        saveas(77,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,local_add,abs_add,hist_add,D3_add,figure_tail]);
        saveas(80,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,local_add,rad_add,D3_add,figure_tail]);
        saveas(81,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,local_add,abs_add,rad_add,D3_add,figure_tail]);
        
        saveas(473,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,local_add,foci_add,D3_add,signal2_add,figure_tail]);
        saveas(474,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,local_add,hist_add,D3_add,signal2_add,figure_tail]);
        saveas(475,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,local_add,RNA_add,D3_add,signal2_add,figure_tail]);
        saveas(476,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,local_add,abs_add,foci_add,D3_add,signal2_add,figure_tail]);
        saveas(477,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,local_add,abs_add,hist_add,D3_add,signal2_add,figure_tail]);
        saveas(480,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,local_add,rad_add,D3_add,signal2_add,figure_tail]);
        saveas(481,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,local_add,abs_add,rad_add,D3_add,signal2_add,figure_tail]);

        %saveas(576,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,local_add,abs_add,foci_add,D3_add,compare_add,figure_tail]);
        %saveas(581,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,local_add,abs_add,rad_add,D3_add,compare_add,figure_tail]);

        %save([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mat_tail],'image_folder','image_type','WGA_channel','DAPI_channel','RNA_channel','protein_channel','resolution','seg_bw','cyto_bw','max_image','foci_bw','nucleus_RNA_profile','foci_RNA_profile','cytoplasmic_RNA_profile','nucleus_protein_profile','cytoplasmic_protein_profile','N_cycle','foci_data','fake_data','h','t_absolute','-append');
        save([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mat_tail],'image_folder','image_type','WGA_channel','DAPI_channel','RNA_channel','protein_channel','resolution','seg_bw','cyto_bw','max_image','foci_bw2D','nucleus_RNA_profile','foci_RNA_profile','cytoplasmic_RNA_profile','nucleus_signal2_profile','foci_signal2_profile','cytoplasmic_signal2_profile','nucleus_protein_profile','cytoplasmic_protein_profile','N_cycle','foci_data','fake_data','h','t_absolute','r_size','foci_data2','fake_data2','h2','t_absolute2','r_size2','resolutionz','quanti_p','nucleus_protein_profile_ab','all_color','Nbin','Mdim','flip_EL','-append');
        
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
        
        if ~isempty(nucleus_signal2_profile)
           xlswrite([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),signal2_add,nu_add,output_tail],nucleus_signal2_profile);
        end
        if ~isempty(cytoplasmic_signal2_profile)
           xlswrite([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),signal2_add,cyto_add,output_tail],cytoplasmic_signal2_profile);
        end
        if ~isempty(foci_signal2_profile)
           xlswrite([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),signal2_add,foci_add,output_tail],foci_signal2_profile);
        end
        for I_data = 1:length(foci_data2)
            if ~isempty(foci_data2{I_data})
                xlswrite([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,local_add,signal2_add,foci_add,output_tail],foci_data2{I_data},I_data);
            end
        end
        for I_data = 1:length(fake_data2)
            if ~isempty(fake_data2{I_data})
                xlswrite([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,local_add,signal2_add,fake_add,foci_add,output_tail],fake_data2{I_data},I_data);
            end
        end

        
        sub_num(list_J,13) = N_cycle;
        
        clear ('seg_bw','cyto_bw','max_image','foci_bw','nucleus_RNA_profile','foci_RNA_profile','cytoplasmic_RNA_profile','nucleus_signal2_profile','foci_signal2_profile','cytoplasmic_signal2_profile','nucleus_protein_profile','cytoplasmic_protein_profile','N_cycle','foci_data','fake_data','foci_signal2_profile','cytoplasmic_signal2_profile','foci_data2','fake_data2');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        
    end
    
    sub_list(:,4:4+size(sub_num,2)-1) = num2cell(sub_num);
    try
        xlswrite([folder_list{list_I,1},in_folder,input_name],sub_list);
    catch
        xlwrite([folder_list{list_I,1},in_folder,input_name],sub_list);
    end
end
toc