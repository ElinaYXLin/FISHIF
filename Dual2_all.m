function Dual2_all(total_name,varargin)

% clear all
close all

if ~isempty(varargin)
    simbol_name = varargin{1};
    protein_add = varargin{2};
    signal2_add = varargin{3};
else
    simbol_name = 'T';
    protein_add = '_Hb';
    signal2_add = '_Bcd';
end
tic
%% Parameter setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
real_name = 'hb';
control_name = 'act5C';
total_folder = 'Regulation_result/';
sub_folder = 'Individual/';
list_name = 'Dual2list.xls';
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
% protein_add = '_Hb';
RNA_add = '_RNA';
% signal2_add = '_Bcd';
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
null_add = '_null';
all_add = '_all';
% sub_pos = [3,3];
sub_pos = [4,11];
cycle_pos0 = zeros(1,sub_pos(1));
cycle_range = [11,12,13,14];
Nthresh00 = 80;   %%% foci number threshold for enrichment averaging
th_good = 0;%0.5;
cycle_choose = [12,13];
% cycle_choose = [1:20];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Data loading: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[num_list, folder_list] = xlsread(list_name);
folder_list = folder_list(strcmpi(simbol_name,folder_list(:,6)),:);
[N1,N2] = size(folder_list);

%%%%%%%%%%%%%%%%%%%
figure(30);clf

figure(2);clf
figure(22);clf
figure(14);clf
figure(16);clf
figure(62);clf
figure(20);clf
figure(220);clf
figure(140);clf
figure(160);clf
figure(620);clf

figure(5);clf
figure(7);clf
figure(6);clf
figure(8);clf
figure(11);clf
figure(10);clf
figure(50);clf
figure(70);clf
figure(110);clf
figure(100);clf

figure(4);clf
figure(40);clf

figure(12);clf
figure(15);clf
figure(17);clf
figure(18);clf
figure(120);clf
figure(150);clf
figure(170);clf

figure(32);clf
figure(35);clf
figure(37);clf
figure(38);clf
figure(320);clf
figure(350);clf
figure(370);clf

figure(42);clf
figure(45);clf
figure(47);clf
figure(420);clf
figure(450);clf
figure(470);clf
%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%
maximize(30)

maximize(2)
maximize(22)
maximize(14)
maximize(16)
maximize(62)
maximize(20)
maximize(220)
maximize(140)
maximize(160)
maximize(620)

maximize(5)
maximize(7)
maximize(6)
maximize(8)
maximize(11)
maximize(10)
maximize(50)
maximize(70)
maximize(110)
maximize(100)

maximize(4)
maximize(40)

maximize(12)
maximize(15)
maximize(17)
maximize(18)
maximize(120)
maximize(150)
maximize(170)

maximize(32)
maximize(35)
maximize(37)
maximize(38)
maximize(320)
maximize(350)
maximize(370)

maximize(42)
maximize(45)
maximize(47)
maximize(420)
maximize(450)
maximize(470)
%%%%%%%%%%%%%%%%%%%

out_para = zeros(0);
out_name = cell(0);
p_name = protein_add(2:end);
p2_name = signal2_add(2:end);

out_title = {'File name','Cycle',[p2_name,' Imax'],[p2_name,' lambda_I'],[p2_name,' Imin'],[p2_name,' NImax'],[p2_name,' Cmax'],[p2_name,' lambda_C'],[p2_name,' Cmin'],[p2_name,' NCmax'],[p2_name,' NCmax_par'],[p2_name,' H_I'],[p2_name,' C0_I'],[p2_name,' Amax_I'],[p2_name,' Amin_I'],[p2_name,' H_N'],[p2_name,' C0_N'],[p2_name,' Amax_N'],[p2_name,' Amin_N'],[p2_name,' H_A'],[p2_name,' C0_A'],[p2_name,' Pmax_A'],[p2_name,' Amin_A'],...
                         'Cycle',[p_name,' Hp_I'],[p_name,' EL0_I'],[p_name,' Imax'],[p_name,' Imin'],[p_name,' NImax'],[p_name,' Hp_C'],[p_name,' EL0_C'],[p_name,' Cmax'],[p_name,' Cmin'],[p_name,' NCmax'],[p_name,' NCmax_par'],[p_name,' H_I'],[p_name,' C0_I'],[p_name,' Amax_I'],[p_name,' Amin_I'],[p_name,' H_N'],[p_name,' C0_N'],[p_name,' Amax_N'],[p_name,' Amin_N'],[p_name,' H_A'],[p_name,' C0_A'],[p_name,' Pmax_A'],[p_name,' Amin_A']};
N_plus = 0;

pro1_out = cell(0);
pro2_out = cell(0);
FISHI_out = cell(0);
FISHN_out = cell(0);
RNAI_out = zeros(0);
RNAN_out = zeros(0);
null_out = zeros(0);
RNAI2_out = zeros(0);
RNAN2_out = zeros(0);
null2_out = zeros(0);
RNAI_out1 = zeros(0);
RNAN_out1 = zeros(0);
null_out1 = zeros(0);
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
    [sub_num, sub_list] = xlsread([folder_list{list_I,1},in_folder,input_name]);
    [M1,M2] = size(sub_list);

    for list_J = 1:M1        
        %%% Data loading: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        image_folder = [folder_list{list_I,1},in_folder,sub_list{list_J,3}];
        result_folder = [folder_list{list_I,1},out_folder];
        load([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mat_tail]);   %%% load analysis result
        load([folder_list{list_I,1},mask_folder,sub_list{list_J,3},mask_name]);   %%% load 3D mask
        z_size = size(mask_stack,3);
        nucleus_protein_profile(:,4) = (nucleus_protein_profile(:,4) > 1) & (nucleus_protein_profile(:,4) < z_size) & (nucleus_protein2_profile(:,4) > 1) & (nucleus_protein2_profile(:,4) < z_size);   %%% check whether the nuclei are on the first or last slice
        nucleus_protein2_profile(:,4) = nucleus_protein_profile(:,4);   %%% check whether the nuclei are on the first or last slice
        nucleus_protein_profile_ab(:,3) = nucleus_protein_profile(:,4);
        nucleus_protein2_profile_ab(:,3) = nucleus_protein2_profile(:,4);
        
        load([folder_list{list_I,1},fit_folder,sub_list{list_J,3}(1:end-1),fit_add,mat_tail]);   %%% load spot intensity fitting result
        nucleus_RNA_profile = foci_refine(nucleus_RNA_profile,foci_RNA_profile,N_thresh);
        nucleus_RNA_profile(:,[2,5]) = nucleus_RNA_profile(:,[2,5])/b;
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
        par_C = max(nucleus_protein_profile_ab(:,2));
        par_C2 = max(nucleus_protein2_profile_ab(:,2));
        %%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%% Nuclear image plot: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        figure(30)
        im_out = max_image(:,:,DAPI_channel);
        imshow(double(im_out)/max(double(im_out(:))));
        title(image_folder);
        %%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%% Data re-plotting/analysis: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [lambda_I,lambda_C,Hp_I,Hp_C,pro2_temp,pro1_temp] = protein2_profile3D_show(nucleus_protein2_profile,cytoplasmic_protein2_profile,quanti_p2,nucleus_protein2_profile_ab,nucleus_protein_profile,cytoplasmic_protein_profile,quanti_p,nucleus_protein_profile_ab,mask_stack,['Sample ',num2str(N_plus)],N_cycle,sub_pos,flip_EL,{signal2_add(2:end),protein_add(2:end)});   %%% protein profile plot (2,22,14,16,62)
        [FISHI_temp,FISHN_temp] = RNA_profile3D_show(nucleus_RNA_profile,foci_RNA_profile,cytoplasmic_RNA_profile,mask_stack,['Sample ',num2str(N_plus)],N_cycle,sub_pos,flip_EL);   %%% RNA profile plot (5,6,8,11)
        modified_foci3D_show(foci_bw2D,max_image,RNA_channel,mask_stack,N_cycle,image_folder,sub_pos);   %%% foci location plot (4)
        [H_I,H_N,pro_center0,RNAI_rate,RNAN_rate,H_A,pro_center,null_rate] = dual_profile_show(nucleus_protein_profile_ab,nucleus_RNA_profile,foci_RNA_profile,['Sample ',num2str(N_plus)],N_cycle,sub_pos,color_code(N_plus-Nccode*(ceil(N_plus/Nccode)-1),:),protein_add(2:end));   %%% Protein-RNA regulation plot (12,15,17,18)
        [H2_I,H2_N,pro2_center0,RNAI2_rate,RNAN2_rate,H2_A,pro2_center,null2_rate] = dual_profile_show(nucleus_protein2_profile_ab,nucleus_RNA_profile,foci_RNA_profile,['Sample ',num2str(N_plus)],N_cycle,sub_pos,color_code(N_plus-Nccode*(ceil(N_plus/Nccode)-1),:),signal2_add(2:end),[32,320,35,350,37,370,38]);   %%% Protein-RNA regulation plot (32,35,37,38)
        [pro_center01,pro2_center01,RNAI_rate1,RNAN_rate1,pro_center1,pro2_center1,null_rate1] = dual2_profile_show(nucleus_protein_profile_ab,nucleus_protein2_profile_ab,nucleus_RNA_profile,foci_RNA_profile,image_folder,N_cycle,sub_pos,{protein_add(2:end),signal2_add(2:end)});
        
        out_para = cat(1,out_para,[[N_cycle,lambda_I,lambda_C,par_C2,H2_I,H2_N,H2_A],[N_cycle,Hp_I,Hp_C,par_C,H_I,H_N,H_A]]);
        pro1_out = cat(1,pro1_out,pro1_temp);
        pro2_out = cat(1,pro2_out,pro2_temp);
        FISHI_out = cat(1,FISHI_out,FISHI_temp);
        FISHN_out = cat(1,FISHN_out,FISHN_temp);
        RNAI_out = cat(1,RNAI_out,RNAI_rate);
        RNAN_out = cat(1,RNAN_out,RNAN_rate);
        null_out = cat(1,null_out,null_rate);
        RNAI2_out = cat(1,RNAI2_out,RNAI2_rate);
        RNAN2_out = cat(1,RNAN2_out,RNAN2_rate);
        null2_out = cat(1,null2_out,null2_rate);
        RNAI_out1 = cat(3,RNAI_out1,RNAI_rate1);
        RNAN_out1 = cat(3,RNAN_out1,RNAN_rate1);
        null_out1 = cat(3,null_out1,null_rate1);
        embryo_info = cat(1,embryo_info,[N_cycle,size(foci_RNA_profile,1)]);
        %%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%% Image output for individual embryos: %%%%%%%%%%%%%%%%%%%%%%%%%%
        image_name = image_folder(1:end-1);
        image_name((image_name == '/') | (image_name == '\')) = '_';
        
        saveas(30,[result_folder0,sub_folder,image_name,DAPI_add,figure_tail]);

        saveas(20,[result_folder0,sub_folder,image_name,protein_add,signal2_add,D3_add,figure_tail]);
        saveas(220,[result_folder0,sub_folder,image_name,protein_add,signal2_add,abs_add,D3_add,figure_tail]);
        saveas(140,[result_folder0,sub_folder,image_name,signal2_add,fate_add,D3_add,figure_tail]);
        saveas(160,[result_folder0,sub_folder,image_name,protein_add,fate_add,D3_add,figure_tail]);
        saveas(620,[result_folder0,sub_folder,image_name,protein_add,signal2_add,fluc_add,D3_add,figure_tail]);

        saveas(110,[result_folder0,sub_folder,image_name,nu_add,fate_add,D3_add,figure_tail]);
        saveas(100,[result_folder0,sub_folder,image_name,nu_add,int_add,fate_add,D3_add,figure_tail]);
        saveas(50,[result_folder0,sub_folder,image_name,fish_add,int_add,D3_add,figure_tail]);
        saveas(70,[result_folder0,sub_folder,image_name,fish_add,num_add,mean_add,D3_add,figure_tail]);

        saveas(40,[result_folder0,sub_folder,image_name,foci_add,seg_add,D3_add,figure_tail]);

        saveas(120,[result_folder0,sub_folder,image_name,nu_add,reg_add,protein_add,D3_add,figure_tail]);
        saveas(150,[result_folder0,sub_folder,image_name,nu_add,foci_add,reg_add,protein_add,D3_add,figure_tail]);
        saveas(170,[result_folder0,sub_folder,image_name,foci_add,null_add,protein_add,D3_add,figure_tail]);

        saveas(320,[result_folder0,sub_folder,image_name,nu_add,reg_add,signal2_add,D3_add,figure_tail]);
        saveas(350,[result_folder0,sub_folder,image_name,nu_add,foci_add,reg_add,signal2_add,D3_add,figure_tail]);
        saveas(370,[result_folder0,sub_folder,image_name,foci_add,null_add,signal2_add,D3_add,figure_tail]);

        saveas(420,[result_folder0,sub_folder,image_name,nu_add,reg_add,protein_add,signal2_add,D3_add,figure_tail]);
        saveas(450,[result_folder0,sub_folder,image_name,nu_add,foci_add,reg_add,protein_add,signal2_add,D3_add,figure_tail]);
        saveas(470,[result_folder0,sub_folder,image_name,foci_add,null_add,protein_add,signal2_add,D3_add,figure_tail]);
        %%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Mean analysis: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cycle_sel = cat(2,num2cell(cycle_range),{cycle_choose});
name_sel = cat(2,cellstr([num2str(cycle_range'),repmat('/',length(cycle_range),1)])',{''});

for I_sel = 1:length(cycle_sel)
    good_data = (embryo_info(:,2) >= Nthresh00) & ismember(embryo_info(:,1),cycle_sel{I_sel}) & (max(null_out,[],2) >= th_good);
    sub_folder2 = name_sel{I_sel};

    figure(122);clf
    maximize(122)
    figure(105);clf
    maximize(105)
    figure(107);clf
    maximize(107)
    figure(112);clf
    maximize(112)
    figure(115);clf
    maximize(115)
    figure(117);clf
    maximize(117)
    [H_I,H_N,H_A] = dual_all_show(total_name,good_data,pro1_out,FISHI_out,FISHN_out,RNAI_out,RNAN_out,null_out,pro_center0,pro_center,p_name,1);

    figure(322);clf
    maximize(322)
    figure(305);clf
    maximize(305)
    figure(307);clf
    maximize(307)
    figure(312);clf
    maximize(312)
    figure(315);clf
    maximize(315)
    figure(317);clf
    maximize(317)
    [H2_I,H2_N,H2_A] = dual_all_show(total_name,good_data,pro2_out,FISHI_out,FISHN_out,RNAI2_out,RNAN2_out,null2_out,pro2_center0,pro2_center,p2_name,4,[322,305,307,312,315,317]);

    figure(412);clf
    maximize(412)
    figure(415);clf
    maximize(415)
    figure(417);clf
    maximize(417)
    dual2_all_show(total_name,good_data,RNAI_out1,RNAN_out1,null_out1,pro_center01,pro2_center01,pro_center1,pro2_center1,{protein_add(2:end),signal2_add(2:end)});

    if exist([result_folder0,sub_folder2]) ~= 7
        mkdir([result_folder0,sub_folder2]);
    end

    saveas(122,[result_folder0,sub_folder2,total_name,protein_add,abs_add,mean_add,D3_add,figure_tail]);
    saveas(105,[result_folder0,sub_folder2,total_name,fish_add,int_add,mean_add,D3_add,figure_tail]);
    saveas(107,[result_folder0,sub_folder2,total_name,fish_add,num_add,mean_add,D3_add,figure_tail]);
    saveas(112,[result_folder0,sub_folder2,total_name,nu_add,reg_add,protein_add,mean_add,D3_add,figure_tail]);
    saveas(115,[result_folder0,sub_folder2,total_name,nu_add,foci_add,reg_add,protein_add,mean_add,D3_add,figure_tail]);
    saveas(117,[result_folder0,sub_folder2,total_name,foci_add,null_add,protein_add,mean_add,D3_add,figure_tail]);

    saveas(322,[result_folder0,sub_folder2,total_name,signal2_add,abs_add,mean_add,D3_add,figure_tail]);
    saveas(312,[result_folder0,sub_folder2,total_name,nu_add,reg_add,signal2_add,mean_add,D3_add,figure_tail]);
    saveas(315,[result_folder0,sub_folder2,total_name,nu_add,foci_add,reg_add,signal2_add,mean_add,D3_add,figure_tail]);
    saveas(317,[result_folder0,sub_folder2,total_name,foci_add,null_add,signal2_add,mean_add,D3_add,figure_tail]);

    saveas(412,[result_folder0,sub_folder2,total_name,nu_add,reg_add,protein_add,signal2_add,mean_add,D3_add,figure_tail]);
    saveas(415,[result_folder0,sub_folder2,total_name,nu_add,foci_add,reg_add,protein_add,signal2_add,mean_add,D3_add,figure_tail]);
    saveas(417,[result_folder0,sub_folder2,total_name,foci_add,null_add,protein_add,signal2_add,mean_add,D3_add,figure_tail]);
end


out_para = cat(1,out_para,[[nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,H2_I,H2_N,H2_A],[nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,H_I,H_N,H_A]]);
out_name = cat(1,out_name,{'Total'});

for I_cycle = 1:length(cycle_range)
    if cycle_range(I_cycle) == min(cycle_range)
        ind_cycle = out_para(:,1) <= cycle_range(I_cycle);
    elseif cycle_range(I_cycle) == max(cycle_range)
        ind_cycle = out_para(:,1) >= cycle_range(I_cycle);
    else
        ind_cycle = out_para(:,1) == cycle_range(I_cycle);
    end
    out_para = cat(1,out_para,mean(out_para(ind_cycle,:),1),std0(out_para(ind_cycle,:),1));
    out_name = cat(1,out_name,['Cycle ',num2str(cycle_range(I_cycle)),' mean'],['Cycle ',num2str(cycle_range(I_cycle)),' err']);
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Output: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
saveas(2,[result_folder0,total_name,protein_add,signal2_add,D3_add,figure_tail]);
saveas(22,[result_folder0,total_name,protein_add,signal2_add,abs_add,D3_add,figure_tail]);
% hgsave(14,[result_folder0,total_name,protein_add,fate_add,D3_add,figure_tail],'-v7.3');
saveas(14,[result_folder0,total_name,signal2_add,fate_add,D3_add,figure_tail2]);
saveas(16,[result_folder0,total_name,protein_add,fate_add,D3_add,figure_tail2]);
saveas(62,[result_folder0,total_name,protein_add,signal2_add,fluc_add,D3_add,figure_tail]);

saveas(5,[result_folder0,total_name,fish_add,int_add,D3_add,figure_tail]);
saveas(7,[result_folder0,total_name,fish_add,num_add,mean_add,D3_add,figure_tail]);
saveas(6,[result_folder0,total_name,fish_add,num_add,D3_add,figure_tail]);
saveas(8,[result_folder0,total_name,foci_add,int_add,D3_add,figure_tail]);
% hgsave(11,[result_folder0,total_name,nu_add,fate_add,D3_add,figure_tail],'-v7.3');
saveas(11,[result_folder0,total_name,nu_add,fate_add,D3_add,figure_tail2]);
saveas(10,[result_folder0,total_name,nu_add,int_add,fate_add,D3_add,figure_tail2]);

saveas(4,[result_folder0,total_name,foci_add,seg_add,D3_add,figure_tail2]);

saveas(12,[result_folder0,total_name,nu_add,reg_add,protein_add,D3_add,figure_tail]);
saveas(15,[result_folder0,total_name,nu_add,foci_add,reg_add,protein_add,D3_add,figure_tail]);
saveas(17,[result_folder0,total_name,foci_add,null_add,protein_add,D3_add,figure_tail]);
saveas(18,[result_folder0,total_name,foci_add,null_add,all_add,protein_add,D3_add,figure_tail]);

saveas(32,[result_folder0,total_name,nu_add,reg_add,signal2_add,D3_add,figure_tail]);
saveas(35,[result_folder0,total_name,nu_add,foci_add,reg_add,signal2_add,D3_add,figure_tail]);
saveas(37,[result_folder0,total_name,foci_add,null_add,signal2_add,D3_add,figure_tail]);
saveas(38,[result_folder0,total_name,foci_add,null_add,all_add,signal2_add,D3_add,figure_tail]);

saveas(42,[result_folder0,total_name,nu_add,reg_add,protein_add,signal2_add,D3_add,figure_tail]);
saveas(45,[result_folder0,total_name,nu_add,foci_add,reg_add,protein_add,signal2_add,D3_add,figure_tail]);
saveas(47,[result_folder0,total_name,foci_add,null_add,protein_add,signal2_add,D3_add,figure_tail]);

% saveas(112,[result_folder0,total_name,nu_add,reg_add,protein_add,mean_add,D3_add,figure_tail]);
% saveas(115,[result_folder0,total_name,nu_add,foci_add,reg_add,protein_add,mean_add,D3_add,figure_tail]);
% saveas(117,[result_folder0,total_name,foci_add,null_add,protein_add,mean_add,D3_add,figure_tail]);
% 
% saveas(312,[result_folder0,total_name,nu_add,reg_add,signal2_add,mean_add,D3_add,figure_tail]);
% saveas(315,[result_folder0,total_name,nu_add,foci_add,reg_add,signal2_add,mean_add,D3_add,figure_tail]);
% saveas(317,[result_folder0,total_name,foci_add,null_add,signal2_add,mean_add,D3_add,figure_tail]);
% 
% saveas(412,[result_folder0,total_name,nu_add,reg_add,protein_add,signal2_add,mean_add,D3_add,figure_tail]);
% saveas(415,[result_folder0,total_name,nu_add,foci_add,reg_add,protein_add,signal2_add,mean_add,D3_add,figure_tail]);
% saveas(417,[result_folder0,total_name,foci_add,null_add,protein_add,signal2_add,mean_add,D3_add,figure_tail]);

out_data = cat(1,out_title,cat(2,out_name,num2cell(out_para)));
xlswrite([result_folder0,total_name,para_add,output_tail],out_data);



toc

end


function y = Hill(beta,x)
 y = beta(3)*x.^beta(1)./(x.^beta(1)+beta(2)^beta(1))+beta(4);
end