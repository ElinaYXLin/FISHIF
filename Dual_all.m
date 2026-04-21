function Dual_all(total_name,varargin)

% clear all
close all

if ~isempty(varargin)
    simbol_name = varargin{1};
    real_gene = varargin{2};
    control_gene = varargin{3};
else
    simbol_name = {'b','ba'};
    real_gene = 'hb';
    control_gene = 'act5C';
end

tic
%% Parameter setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
real_name = real_gene;
control_name = control_gene;
total_folder = 'Regulation_result/';
sub_folder = 'Individual/';
list_name = 'Duallist.xls';
in_folder = 'stacks/';
input_name = 'matchlist.xls';
mismatch_name = 'mismatch.xls';
mismatch_name2 = 'mismatch_60X.xls';
image_type = '*.tif';
mask_folder = 'masks/';
mask_name = 'mask.mat';
out_folder = 'Results/';
hist_folder = 'Histogram/';
fit_folder = 'Histogram_A/';
fit_folder_protein = 'Histogram_protein_A/';
hist_folder2 = 'Histogram_nullo/';
fit_folder2 = 'Histogram_default/';
fit_add = '_spot_fit';
hist_tail = '_raw.xls';
N_thresh = 3;
output_tail = '.xls';
figure_tail = '.fig';
figure_tail2 = '.jpg';
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
average_add = '_average';
spot_add = '_spot';
TX_add = '_TX';
var_add = '_var';
abs_add = '_abs';
para_add = '_para';
p_add = '_scale';
null_add = '_null';
all_add = '_all';
con_add = '_con';
vm_add = '_vm';
low_add = '_low';
k2s_add = '_k2s';
% sub_pos = [3,3];
sub_pos = [4,11];
cycle_pos0 = zeros(1,sub_pos(1));
cycle_range = [11,12,13,14];
Nthresh00 = 80;   %%% foci number threshold for enrichment averaging
th_good = 0;%0.5;
cycle_choose = [11,12];
% cycle_choose = [1:20];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Data loading: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[num_list, folder_list] = xlsread(list_name);
% folder_list = folder_list(strcmpi('b',folder_list(:,6)) | strcmpi('ba',folder_list(:,6)),:);
% folder_list = folder_list(strcmpi('b3',folder_list(:,6)),:);
Ifolder = false(size(folder_list),1);
for  Icheck = 1:length(simbol_name)
    Ifolder = Ifolder | strcmpi(simbol_name{Icheck},folder_list(:,6));
end
folder_list = folder_list(Ifolder,:);

[N1,N2] = size(folder_list);

%%%%%%%%%%%%%%%%%%%
figure(30);clf

figure(2);clf
figure(22);clf
figure(14);clf
figure(62);clf
figure(20);clf
figure(220);clf
figure(140);clf
figure(620);clf

figure(5);clf
figure(31);clf
figure(33);clf
figure(35);clf
figure(34);clf
figure(36);clf
figure(7);clf
figure(6);clf
figure(8);clf
figure(11);clf
figure(10);clf
figure(13);clf
figure(50);clf
figure(310);clf
figure(330);clf
figure(350);clf
figure(340);clf
figure(360);clf
figure(70);clf
figure(110);clf
figure(100);clf
figure(130);clf

figure(4);clf
figure(40);clf

figure(12);clf
figure(91);clf
figure(15);clf
figure(17);clf
figure(18);clf
figure(51);clf
figure(120);clf
figure(910);clf
figure(150);clf
figure(170);clf
figure(510);clf
%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%
maximize(30)

maximize(2)
maximize(22)
maximize(14)
maximize(62)
maximize(20)
maximize(220)
maximize(140)
maximize(620)

maximize(5)
maximize(31)
maximize(33)
maximize(35)
maximize(34)
maximize(36)
maximize(7)
maximize(6)
maximize(8)
maximize(11)
maximize(10)
maximize(13)
maximize(50)
maximize(310)
maximize(330)
maximize(350)
maximize(340)
maximize(360)
maximize(70)
maximize(110)
maximize(100)
maximize(130)

maximize(4)
maximize(40)

maximize(12)
maximize(91)
maximize(15)
maximize(17)
maximize(18)
maximize(51)
maximize(120)
maximize(910)
maximize(150)
maximize(510)
%%%%%%%%%%%%%%%%%%%

out_para = zeros(0);
out_p = zeros(0);
out_name = cell(0);
out_title = {'File name','Cycle','Imax','lambda_I','Imin','NImax','Cmax','lambda_C','Cmin','NCmax','NCmax_par','H_I','C0_I','Amax_I','Amin_I','H_N','C0_N','Amax_N','Amin_N','H_A','C0_A','Pmax_A','Amin_A'};
out_titlep = {'#','File name','Cycle','k1','mean(post)','DAPI k1','DAPI mean(post)','# nuclei','nuclei/cyto','I0RNA'};
out_titles = {'#','File name','Cycle','mean(S)','mean(C)','mean(S*C)','max(S)','max(C)','max(S*C)','max(Sm)','max(Cm)','max(Sm*Cm)','RNA_bgA','RNA_bgP','Sembryo'};
out_titlef = {'#','File name','Cycle','max(FISHI)','std(FISHI)','err(FISHI)','EL(FISHI)','EL(1/2I)','max(FISHN)','std(FISHN)','err(FISHN)','EL(FISHN)','EL(1/2N)','max(FISHN0)','std(FISHN0)','err(FISHN0)','EL(FISHN0)','EL(1/2N0)','a_EL','b_EL','kon*T_EL','kTX*T_EL','mean corr_EL','mean corr control_EL','a_Bcd','b_Bcd','kon*T_Bcd','kTX*T_Bcd','mean corr_Bcd','mean corr control_Bcd','h_k','c_k'};
N_plus = 0;

pro_out = cell(0);
S_out = zeros(0);
FISHI_out = cell(0);
FISHN_out = cell(0);
fociN_out = cell(0);
out_fi = zeros(0);
RNAI_out = zeros(0);
RNAN_out = zeros(0);
null_out = zeros(0);
embryo_info = zeros(0);
Nccode = 6;
color_code = mycolors(Nccode);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if total_name(end) == '/' || total_name(end) == '\'
    total_name = total_name(1:end-1);
end
result_folder0 = [total_folder,total_name,'/'];
total_name((total_name == '/') | (total_name == '\')) = '_';

if exist([result_folder0,sub_folder]) ~= 7
    mkdir([result_folder0,sub_folder]);
end



%% Image analysis: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for list_I = 1:N1
    [sub_num, sub_list, sub_all] = xlsread([folder_list{list_I,1},in_folder,input_name]);
    if size(sub_all,2) >= 17
        I_se = ~strcmp('F',sub_all(:,17));
        sub_num = sub_num(I_se,:);
        sub_list = sub_list(I_se,:);
    end
    [M1,M2] = size(sub_list);

    for list_J = 1:M1        
        %%% Data loading: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        image_folder = [folder_list{list_I,1},in_folder,sub_list{list_J,3}];
        result_folder = [folder_list{list_I,1},out_folder];
        load([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mat_tail]);   %%% load analysis result
        load([folder_list{list_I,1},mask_folder,sub_list{list_J,3},mask_name]);   %%% load 3D mask
        z_size = size(mask_stack,3);
        nucleus_protein_profile0 = nucleus_protein_profile;
        nucleus_protein_profile(:,4) = (nucleus_protein_profile(:,4) > 1) & (nucleus_protein_profile(:,4) < z_size);   %%% check whether the nuclei are on the first or last slice
        nucleus_protein_profile_ab(:,3) = nucleus_protein_profile(:,4);
        
        load([folder_list{list_I,1},fit_folder,sub_list{list_J,3}(1:end-1),fit_add,mat_tail]);   %%% load spot intensity fitting result
        nucleus_RNA_profile = foci_refine(nucleus_RNA_profile,foci_RNA_profile,N_thresh);
        nucleus_RNA_profile(:,[2,5,6]) = nucleus_RNA_profile(:,[2,5,6])/b;
        nucleus_RNA_profile = nucleus_RNA_profile(:,[1:4,6]);
        cytoplasmic_RNA_profile(:,2) = cytoplasmic_RNA_profile(:,2)/b;
        
        start_limit = 0.4;
        end_limit = 1-start_limit;
        index_true = logical(nucleus_protein_profile(:,4));
        flip_EL =  mean(nucleus_protein_profile(index_true & (nucleus_protein_profile(:,1) <= start_limit),2)) < mean(nucleus_protein_profile(index_true & (nucleus_protein_profile(:,1) >= end_limit),2));
        
        out_name = cat(1,out_name,{image_folder});
        N_plus = N_plus+1;
        %%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%% subplot organization: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        N_cycle = sub_num(list_J,13);
        I_cycle = find(N_cycle == cycle_range);
        if N_cycle > max(cycle_range)
            I_cycle = find(max(cycle_range) == cycle_range);
        elseif N_cycle < min(cycle_range)
            I_cycle = find(min(cycle_range) == cycle_range);
        end
        if ~isempty(I_cycle)
            cycle_pos0(I_cycle) = cycle_pos0(I_cycle)+1;
            sub_pos(3) = sub2ind(sub_pos([2,1]),cycle_pos0(I_cycle),I_cycle);
        else
            sub_pos(3) = 0;
        end
        %%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%% Pre-modifications: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %         load([folder_list{list_I,1},fit_folder_protein,sub_list{list_J,3}(1:end-1),fit_add,mat_tail]);   %%% load spot intensity fitting result
%         new_p = quanti_p(1)*(1.35/1.2)^2;
%         quanti_p(1) = new_p;
% %         par_p = b_p/sqrt(2).*resolution.*resolution.*6.02e8;
%         nucleus_protein_profile_ab(:,2) = nucleus_protein_profile_ab(:,2)*(1.2/1.35)^2;
% %         nucleus_protein_profile_ab_p = nucleus_protein_profile_ab(:,2)/par_p*new_p;
% %         par_C = max(nucleus_protein_profile_ab_p);
%         par_C = max(nucleus_protein_profile_ab(:,2))*(1.35/1.2)^2;
        [nucleus_protein_profile_ab,p_temp] = protein_rescale(nucleus_protein_profile,nucleus_protein_profile_ab,quanti_p);
        [~,pD_temp] = protein_rescale(nucleus_protein_profile(:,[1,5:end]),nucleus_protein_profile_ab,quanti_p);
        mean_ratio = protein_profile3D_short2(mask_stack);
        out_p = cat(1,out_p,[N_cycle,p_temp,pD_temp,size(nucleus_protein_profile,1),mean_ratio,b]);
        par_C = max(nucleus_protein_profile_ab(:,2));
        %%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%% Nuclear image plot: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        figure(30)
        im_out = max_image(:,:,DAPI_channel);
        imshow(double(im_out)/max(double(im_out(:))));
        title(image_folder);
        %%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%% Data re-plotting/analysis: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [lambda_I,lambda_C,pro_temp] = protein_profile3D_show(nucleus_protein_profile,cytoplasmic_protein_profile,quanti_p,nucleus_protein_profile_ab,mask_stack,['Sample ',num2str(N_plus)],N_cycle,sub_pos,flip_EL);   %%% protein profile plot (2,22,14,62)
        S_temp = protein_profile3D_short(nucleus_protein_profile0,nucleus_protein_profile_ab,mask_stack,resolution,lambda_C,index_true);   %%% Nuclear area extraction
        [FISHI_temp,FISHN_temp,fociN_temp,FISHI_temp2,FISHN_temp2,fociN_temp2,phat_temp,S2_temp,S3_temp] = RNA_profile3D_show(nucleus_RNA_profile,foci_RNA_profile,cytoplasmic_RNA_profile,mask_stack,['Sample ',num2str(N_plus)],N_cycle,sub_pos,flip_EL,nucleus_protein_profile_ab);   %%% RNA profile plot (5,31,32,6,8,11,10,13)
        RNA_profile3D_show(nucleus_RNA_profile,foci_RNA_profile,cytoplasmic_RNA_profile,mask_stack,['Sample ',num2str(N_plus)],N_cycle,sub_pos,flip_EL,nucleus_protein_profile_ab,[0,0,35,0,0,350,0,0,0,0,0,0,0,0,0,0],[0.5,0.6]);   %%% RNA profile plot (5,31,32,6,8,11,10,13)
        modified_foci3D_show(foci_bw2D,max_image,RNA_channel,mask_stack,N_cycle,image_folder,sub_pos);   %%% foci location plot (4)
        [H_I,H_N,pro_center0,RNAI_rate,RNAN_rate,H_A,pro_center,null_rate,phat2_temp] = dual_profile_show(nucleus_protein_profile_ab,nucleus_RNA_profile,foci_RNA_profile,['Sample ',num2str(N_plus)],N_cycle,sub_pos,color_code(N_plus-Nccode*(ceil(N_plus/Nccode)-1),:));   %%% Protein-RNA regulation plot (12,91,15,17)
        dual_profile_show(nucleus_protein_profile_ab,nucleus_RNA_profile,foci_RNA_profile,['Sample ',num2str(N_plus)],N_cycle,sub_pos,color_code(N_plus-Nccode*(ceil(N_plus/Nccode)-1),:),'Protein',[0,0,0,0,0,0,0,0,0,36,360,0,0],[2e-9,5e-9]);   %%% Protein-RNA regulation plot (12,91,15,17)
        
        out_para = cat(1,out_para,[N_cycle,lambda_I,lambda_C,par_C,H_I,H_N,H_A]);
        pro_out = cat(1,pro_out,pro_temp);
        S_out = cat(1,S_out,[N_cycle,S_temp,S2_temp/resolution/resolution/resolutionz,S3_temp]);
        FISHI_out = cat(1,FISHI_out,FISHI_temp);
        FISHN_out = cat(1,FISHN_out,FISHN_temp);
        fociN_out = cat(1,fociN_out,fociN_temp);
        out_fi = cat(1,out_fi,[N_cycle,FISHI_temp2,FISHN_temp2,fociN_temp2,phat_temp,phat2_temp]);
        RNAI_out = cat(1,RNAI_out,RNAI_rate);
        RNAN_out = cat(1,RNAN_out,RNAN_rate);
        null_out = cat(1,null_out,null_rate);
        embryo_info = cat(1,embryo_info,[N_cycle,size(foci_RNA_profile,1)]);
        %%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%% Image output for individual embryos: %%%%%%%%%%%%%%%%%%%%%%%%%%
        image_name = image_folder(1:end-1);
        image_name((image_name == '/') | (image_name == '\')) = '_';
        
        saveas(30,[result_folder0,sub_folder,image_name,DAPI_add,figure_tail]);

        saveas(20,[result_folder0,sub_folder,image_name,protein_add,D3_add,figure_tail]);
        saveas(220,[result_folder0,sub_folder,image_name,protein_add,abs_add,D3_add,figure_tail]);
        saveas(140,[result_folder0,sub_folder,image_name,protein_add,fate_add,D3_add,figure_tail]);
        saveas(620,[result_folder0,sub_folder,image_name,protein_add,fluc_add,D3_add,figure_tail]);

        saveas(110,[result_folder0,sub_folder,image_name,nu_add,fate_add,D3_add,figure_tail]);
        saveas(100,[result_folder0,sub_folder,image_name,nu_add,int_add,fate_add,D3_add,figure_tail]);
        saveas(50,[result_folder0,sub_folder,image_name,fish_add,int_add,D3_add,figure_tail]);
        saveas(310,[result_folder0,sub_folder,image_name,fish_add,int_add,vm_add,EL_add,D3_add,figure_tail]);
        saveas(330,[result_folder0,sub_folder,image_name,fish_add,int_add,hist_add,D3_add,figure_tail]);
        saveas(350,[result_folder0,sub_folder,image_name,fish_add,int_add,hist_add,low_add,D3_add,figure_tail]);
        saveas(340,[result_folder0,sub_folder,image_name,fish_add,int_add,hist_add,protein_add,D3_add,figure_tail]);
        saveas(360,[result_folder0,sub_folder,image_name,fish_add,int_add,hist_add,protein_add,low_add,D3_add,figure_tail]);
        saveas(70,[result_folder0,sub_folder,image_name,fish_add,num_add,D3_add,figure_tail]);
        saveas(130,[result_folder0,sub_folder,image_name,fish_add,int_add,cmp_add,D3_add,figure_tail]);

        saveas(40,[result_folder0,sub_folder,image_name,foci_add,seg_add,D3_add,figure_tail]);

        saveas(120,[result_folder0,sub_folder,image_name,nu_add,reg_add,D3_add,figure_tail]);
        saveas(910,[result_folder0,sub_folder,image_name,fish_add,int_add,vm_add,protein_add,D3_add,figure_tail]);
        saveas(150,[result_folder0,sub_folder,image_name,nu_add,foci_add,reg_add,D3_add,figure_tail]);
        saveas(170,[result_folder0,sub_folder,image_name,foci_add,null_add,D3_add,figure_tail]);
        saveas(510,[result_folder0,sub_folder,image_name,foci_add,k2s_add,D3_add,figure_tail]);
        %%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    close(30)
    close(20)
    close(220)
    close(140)
    close(620)
    close(110)
    close(100)
    close(50)
    close(310)
    close(330)
    close(340)
    close(70)
    close(130)
    close(40)
    close(120)
    close(910)
    close(150)
    close(170)
    close(510)


%% Mean plot of I regulation: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cycle_sel = cat(2,num2cell(cycle_range),{cycle_choose});
name_sel = cat(2,cellstr([num2str(cycle_range'),repmat('/',length(cycle_range),1)])',{''});
out_namef = out_name;
out_name0 = out_name;
HH_I = out_para(:,end-11:end-8);
% cycle_sel = {cycle_choose};
% name_sel = {''};

for I_sel = 1:length(cycle_sel)
    good_data = (embryo_info(:,2) >= Nthresh00) & ismember(embryo_info(:,1),cycle_sel{I_sel}) & (max(null_out,[],2) >= th_good);
    if nnz(good_data)
        sub_folder2 = name_sel{I_sel};
        % % RNAI_out0 = RNAI_out(good_data,:);
        % % RNAI_value = RNAI_out0;
        % % RNAI_value(isnan(RNAI_out0)) = 0;
        % % RNAI_mean = sum(RNAI_value)./sum(~isnan(RNAI_out0));
        % % RNAI_std0 = sqrt((sum(RNAI_value.^2)./sum(~isnan(RNAI_out0))-RNAI_mean.^2)./sum(~isnan(RNAI_out0)));
        % % good_bin = (~isnan(RNAI_std0)) & (RNAI_std0 > 0);

        figure(122);clf
        maximize(122)
        figure(105);clf
        maximize(105)
        figure(131);clf
        maximize(131)
        figure(107);clf
        maximize(107)
        figure(112);clf
        maximize(112)
        figure(191);clf
        maximize(191)
        figure(115);clf
        maximize(115)
        figure(117);clf
        maximize(117)
        figure(101);clf
        maximize(101)

        [H_I,H_N,H_A,para_EL,para_Bcd] = dual_all_show([total_name,'/',sub_folder2],out_name,[good_data,embryo_info(:,1)],HH_I,pro_out,FISHI_out,FISHN_out,fociN_out,RNAI_out,RNAN_out,null_out,pro_center0,pro_center,protein_add(2:end),4,[],color_code);
        if length(cycle_sel{I_sel}) == 1
            out_fi = cat(1,out_fi,[cycle_sel{I_sel}*ones(2,1),nan(2,15),para_EL,para_Bcd]);
            out_namef = cat(1,out_namef,{'Fit';'Fit error'});
        end

        out_para = cat(1,out_para,[nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,H_I,H_N,H_A]);
        out_para = cat(1,out_para,[nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,zeros(size(H_I)),zeros(size(H_N)),zeros(size(H_A))]);
        out_name = cat(1,out_name,['Cycle ',num2str(cycle_sel{I_sel}),' mean'],['Cycle ',num2str(cycle_sel{I_sel}),' err']);
        
        if exist([result_folder0,sub_folder2]) ~= 7
            mkdir([result_folder0,sub_folder2]);
        end
        saveas(122,[result_folder0,sub_folder2,total_name,protein_add,abs_add,mean_add,D3_add,figure_tail]);
        saveas(105,[result_folder0,sub_folder2,total_name,fish_add,int_add,mean_add,D3_add,figure_tail]);
        saveas(131,[result_folder0,sub_folder2,total_name,fish_add,int_add,vm_add,EL_add,mean_add,D3_add,figure_tail]);
        saveas(107,[result_folder0,sub_folder2,total_name,fish_add,num_add,mean_add,D3_add,figure_tail]);
        saveas(112,[result_folder0,sub_folder2,total_name,nu_add,reg_add,mean_add,D3_add,figure_tail]);
        saveas(191,[result_folder0,sub_folder2,total_name,fish_add,int_add,vm_add,protein_add,mean_add,D3_add,figure_tail]);
        saveas(115,[result_folder0,sub_folder2,total_name,nu_add,foci_add,reg_add,mean_add,D3_add,figure_tail]);
        saveas(117,[result_folder0,sub_folder2,total_name,foci_add,null_add,mean_add,D3_add,figure_tail]);
        saveas(101,[result_folder0,sub_folder2,total_name,nu_add,fate_add,mean_add,D3_add,figure_tail]);
    end
end
close(122)
close(105)
close(131)
close(107)
close(112)
close(191)
close(115)
close(117)
close(101)

% % %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % out_para = cat(1,out_para,[nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,H_I,H_N,H_A]);
% % out_name = cat(1,out_name,{'Total'});
% % 
% % for I_cycle = 1:length(cycle_range)
% %     if cycle_range(I_cycle) == min(cycle_range)
% %         ind_cycle = out_para(:,1) <= cycle_range(I_cycle);
% %     elseif cycle_range(I_cycle) == max(cycle_range)
% %         ind_cycle = out_para(:,1) >= cycle_range(I_cycle);
% %     else
% %         ind_cycle = out_para(:,1) == cycle_range(I_cycle);
% %     end
% %     if nnz(ind_cycle)
% %         out_para = cat(1,out_para,mean(out_para(ind_cycle,:)),std0(out_para(ind_cycle,:)));
% %         out_name = cat(1,out_name,['Cycle ',num2str(cycle_range(I_cycle)),' mean'],['Cycle ',num2str(cycle_range(I_cycle)),' err']);
% %     end
% % end




%% Output: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
saveas(2,[result_folder0,total_name,protein_add,D3_add,figure_tail]);
saveas(22,[result_folder0,total_name,protein_add,abs_add,D3_add,figure_tail]);
% hgsave(14,[result_folder0,total_name,protein_add,fate_add,D3_add,figure_tail],'-v7.3');
saveas(14,[result_folder0,total_name,protein_add,fate_add,D3_add,figure_tail2]);
saveas(62,[result_folder0,total_name,protein_add,fluc_add,D3_add,figure_tail]);

saveas(5,[result_folder0,total_name,fish_add,int_add,D3_add,figure_tail]);
saveas(31,[result_folder0,total_name,fish_add,int_add,vm_add,EL_add,D3_add,figure_tail]);
saveas(33,[result_folder0,total_name,fish_add,int_add,hist_add,D3_add,figure_tail]);
saveas(35,[result_folder0,total_name,fish_add,int_add,hist_add,low_add,D3_add,figure_tail]);
saveas(34,[result_folder0,total_name,fish_add,int_add,hist_add,protein_add,D3_add,figure_tail]);
saveas(36,[result_folder0,total_name,fish_add,int_add,hist_add,protein_add,low_add,D3_add,figure_tail]);
saveas(7,[result_folder0,total_name,fish_add,num_add,mean_add,D3_add,figure_tail]);
saveas(6,[result_folder0,total_name,fish_add,num_add,D3_add,figure_tail]);
saveas(8,[result_folder0,total_name,foci_add,int_add,D3_add,figure_tail]);
% hgsave(11,[result_folder0,total_name,nu_add,fate_add,D3_add,figure_tail],'-v7.3');
saveas(11,[result_folder0,total_name,nu_add,fate_add,D3_add,figure_tail2]);
saveas(10,[result_folder0,total_name,nu_add,int_add,fate_add,D3_add,figure_tail2]);
saveas(13,[result_folder0,total_name,fish_add,int_add,cmp_add,D3_add,figure_tail]);

saveas(4,[result_folder0,total_name,foci_add,seg_add,D3_add,figure_tail2]);

saveas(12,[result_folder0,total_name,nu_add,reg_add,D3_add,figure_tail]);
saveas(91,[result_folder0,total_name,fish_add,int_add,vm_add,protein_add,D3_add,figure_tail]);
saveas(15,[result_folder0,total_name,nu_add,foci_add,reg_add,D3_add,figure_tail]);
saveas(17,[result_folder0,total_name,foci_add,null_add,D3_add,figure_tail]);
saveas(18,[result_folder0,total_name,foci_add,null_add,all_add,D3_add,figure_tail]);
saveas(51,[result_folder0,total_name,foci_add,k2s_add,D3_add,figure_tail]);


out_data = cat(1,out_title,cat(2,out_name,num2cell(out_para)));
out_pdata = cat(1,out_titlep,cat(2,num2cell([1:N_plus]'),out_name0,num2cell(out_p)));
out_sdata = cat(1,out_titles,cat(2,num2cell([1:N_plus]'),out_name0,num2cell(S_out)));
out_fdata = cat(1,out_titlef,cat(2,num2cell([[1:N_plus],nan(1,length(out_namef)-N_plus)]'),out_namef,num2cell(out_fi)));

xlswrite([result_folder0,total_name,para_add,output_tail],out_data);
xlswrite([result_folder0,total_name,p_add,output_tail],out_pdata);
xlswrite([result_folder0,total_name,con_add,output_tail],out_sdata);
xlswrite([result_folder0,total_name,fish_add,output_tail],out_fdata);



toc
close all
end


function y = Hill(beta,x)
 y = beta(3)*x.^beta(1)./(x.^beta(1)+beta(2)^beta(1))+beta(4);
end