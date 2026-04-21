clear all
close all

tic

%% Parameter setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% path0 = '\\202.121.182.19\Heng 64\';
path0 = '';
addpath(path0);
list_name = 'Duallist.xls';
in_folder = 'stacks/';
% old_add = '_old';
old_add = '';
input_name = 'matchlist.xls';
mismatch_name = 'mismatch.xls';
mismatch_name2 = 'mismatch_60X.xls';
xymismatch_name = 'xymismatch.mat';
xymismatch_name2 = 'xymismatch_60X.mat';
flip_list = {'Cad'};
flip0 = false;
image_type = '*.tif';
mask_folder = 'masks/';
mask_name = 'mask.mat';
out_folder = 'Results/';
% out_folder0 = 'Results_3Dz1/';
% out_folder0 = 'Results_decross/';
% out_folder0 = 'Results_original/';
out_folder0 = 'Results0/';
% hist_folder = 'Histogram/';
hist_folder = 'Histogram_alignment/';
hist_folder_RNA2 = 'Histogram_alignment_RNA2/';
fit_folder = 'Histogram_A/';
hist_folder2 = {'Histogram_en/','Histogram_en2/'};
fit_folder2 = {'Histogram_protein_A/','Histogram_protein2_A/'};
% hist_folder_couple = 'Histogram_alignment/';
% hist_folder2_couple = 'Histogram_alignment_RNA2/';
fit_add = '_spot_fit';
fit_add2 = '_peak_fit';
hist_tail = '_raw.xls';
hist_link_tail = '_link.xls';
N_thresh = 10; %0;
output_tail = '.xls';
figure_tail = '.fig';
figure_tail2 = '.png';
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
spot_add = '_spot';
sel_add = '_sel';
peak_add = '_peak';
fit_add2 = '_fit';
binding_add = '_binding';
fit_tail = '.sfit';
fit_name = 'fit_result';
combine_add = '_combine';

r_th = 0.28;
rout_min = 0.28;
rout_max = 0.34;
rfit_min =0.34;
rfit_max = 0.8;
rdist_max = 0.8;
r_ad = 1;
rxymax = [2.05,5.25];
% % rxymax = [0,10];
rxmax = 5;
rymax = 5;
% % rxmax = 4;
% % rymax = 4;
resolution0 = 0.083;
relip = 0.95;
% Imin0 = 3e3;
Imin0 = 0;
% dz_new = [];
dz_new = -1:1;
sub_pos = [1,2];
ana_folder = [path0,'enrichspot_result\'];
% embryo_name = 'test';
% % embryo_name = '22445_2_hb_cds_TMR_gal4-647_Hb_488_par3';
% embryo_name = '22445_2_hb_cds_TMR_gal4-647_bcd_combine_488_par3';
% % embryo_name = '22445_2_hb_cds_TMR_gal4-647_bcd_combine_488_par3_new';
% % % embryo_name = '22445_2_hb_cds_TMR_gal4-647_bcd_combine_488_par3_new_add';
embryo_name = {'22445_2_gal4-647_bcd_488_cad_555_par3_new_par3','22445_2_gal4-647_bcd_488_cad_555_par3_new_ch2_par3'};
% % % % % embryo_name = '22445_2_hb_cds_TMR_gal4-647_bcd_combine_488_par3_new2';
% % embryo_name = '22445_2_hb_cds_TMR_gal4-647_bcd_488_par3';
% % embryo_name = '22445_2_hb_cds_TMR_gal4-647_Cad_488_par3';
% % embryo_name = '22445_2_hb_cds_TMR_gal4-647_bcd-rat_488_par3';
% % embryo_name = '08022019_2_22445_2_hb_cds_TMR_gal4-647_bcd_488_par3';
% embryo_name = '09302015_16249-1-5M_FISHIF_hb_gal4_Bcd_par2';
% % embryo_name = '11012014_10652-1-3M_FISHIF_lacZ_gal4_Bcd_par2';
% embryo_name = '02072016_16249-1-5M_FISHIF_hb_gal4_Cad_Bcd_PCad_par2';
% embryo_name = '12072016b_1_16249-1-5M_FISHIF_hb_gal4_H3_par3';
% % embryo_name = '12072016b_2_16249-1-5M_FISHIF_hb_gal4_H4K5ac_par3';
% % % embryo_name = '12072016b_3_16249-1-5M_FISHIF_hb_gal4_H3K27ac_par3';
fit_name = 'enrich_hist_fit';
I_plot = 1:7;%[6,8:16];%[1,2,5,8,9,10];%[3:5,8:11,22:25,27:34];%[3:5,8:11,14:19,22:25,27:34];
I_cycle = 10:12;
epsilon = 1e-4;
smin = 0.99;
smin2 = 0.99;
smin3 = 0.99999;
rsize0 = 4;
epsilonxyz = 0.55;
epsilonxyz2 = 0.1;
r0_cali = 1;%0.2697/3;a
r_rescale0 = [ones(1,12)*3/2.8,ones(1,5)*3/3.5];

EL_range_fluc0 = [0,0.35];
EL_range_fluc1 = [0,0.35];
dz_range = [0,0];
% z_range = [1,100];%[2,22];
zdiff_nucleus_RNA = [-6,4];

gau2_gen = @(n) ['m(',num2str(n),',1)/sqrt(2*pi)/m(',num2str(n),',3)*exp(-(x-m(',num2str(n),',2)).^2/2/m(',num2str(n),',3)^2)'];   %%% single-gaussian model function text generator
gau3_gen = @(n) ['a',num2str(n),'/sqrt(2*pi)/c*exp(-(x-',num2str(n),'*b).^2/2/c^2)'];   %%% single-gaussian model function text generator
lb = [0,0,0.4];
ub = [1,inf,10];
ib = [1,1,1];

IC = 4;

fsize = 7;
fsize2 = 7;
% nfont = 'Helvetica';
nfont = 'Arial';

result_folder = ['enrichspot_result/',embryo_name{1},'/analysis results_test_combine/'];
out_name = [embryo_name{1},combine_add,'_C',num2str(I_cycle)];
enrich_add = '_enrich';
Inten_add = '_Inten';
distri_add = '_distri';
mean_add = '_mean';
TF_add = '_TF';
EL_add = '_EL';
Pon_add = '_Pon';
fig_tail = '.fig';
fig_tail2 = '.eps';

% % C0 = [55.37,64.49,69.81,16.04,73.06];
% C0 = [51,55,46,11,38,19,20,31,21,39];
% C0 = [51,55,46,11,38,50,50,45,45,33,24,29,19,20,31,21,39];
C0 =[];

t_gate = true;
% t_gate = false;
rescale_en = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
temp_scale = zeros(0);

xyz_all0_cell2 = cell(size(embryo_name));
Imin_all0_cell2 = cell(size(embryo_name));
dxy0_cell2 = cell(size(embryo_name));
% foci_circle_area_cell2 = cell(size(embryo_name));

for pp = 1:length(embryo_name)
%% Data loading and adjusting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    load([ana_folder,embryo_name{pp},'/',embryo_name{pp},mat_tail]);
    N_outmin_cell = cell(size(dxy0_cell));
    N_outmax_cell = cell(size(dxy0_cell));
    en_circle_area_cell = cell(size(dxy0_cell));
    dz_en_RNA_cell = cell(size(dxy0_cell));
    N_en2_cell = cell(size(dxy0_cell));
    enrich_RNA_cell = cell(size(dxy0_cell));
    enrich_all_cell = cell(size(dxy0_cell));
    dxyz_all0_cell = cell(size(dxy0_cell));
    xyz_all2_cell = cell(size(dxy0_cell));

    for ii = 1:length(dxy0_cell)
% %         load([ana_folder,embryo_name{pp},'/',out_folder0,out_folder0(1:end-2),num2str(ii),mat_tail])
        load([path0,em_name{ii}(1:end-1),mat_tail],'foci_data2','fake_data2','max_image','h','r_size','resolution','resolutionz','quanti_p','nucleus_protein_profile');
        id0 = strfind(em_name{ii},out_folder(1:end-1));
        load([path0,em_name{ii}(1:id0-1),fit_folder2{pp},em_name{ii}(id0+8:end-1),fit_add,mat_tail],'b');

        Inten_temp = xlsread([path0,em_name{ii}(1:id0-1),hist_folder_RNA2,em_name{ii}(id0+8:end-1),hist_tail]);
        xyz2 = Inten_temp(:,6:8).*repmat([resolution,resolution,resolutionz],size(Inten_temp,1),1);


    % %     if ~exist([ana_folder,embryo_name{pp},'/',out_folder0])
    % %         mkdir([ana_folder,embryo_name{pp},'/',out_folder0])
    % %     end
    % %     save([ana_folder,embryo_name{pp},'/',out_folder0,out_folder0(1:end-2),num2str(ii),mat_tail],'foci_data2','fake_data2','max_image','h','r_size','resolution','resolutionz','quanti_p','nucleus_protein_profile','b','Inten_temp')


        z_range0 = length(dir([path0,em_name{ii}(1:id0-1),in_folder,em_name{ii}(id0+8:end),image_type]));
    %     r_rescale2 = quanti_p(1)/(b*sqrt(2*pi)*resolution^2*resolutionz*6.02e8);
        r_rescale2 = b2_em(ii)/b;%b2_em(1)/b2_em(ii);%/sqrt(2*pi);
        %%% Rescale protein concentration:
        if ~all(EL_range_fluc0 == EL_range_fluc1)
            Iz = nucleus_protein_profile(:,4) >= 1+dz_range(1) & nucleus_protein_profile(:,4) <= z_range0-dz_range(2);
            Itrue0 = (nucleus_protein_profile(:,1) >= EL_range_fluc0(1)) & (nucleus_protein_profile(:,1) <= EL_range_fluc0(2)) & Iz;
            Itrue = (nucleus_protein_profile(:,1) >= EL_range_fluc1(1)) & (nucleus_protein_profile(:,1) <= EL_range_fluc1(2)) & Iz;

            p0_protein = polyfit(nucleus_protein_profile(Itrue0,2),nucleus_protein_profile(Itrue0,3),1);
            p_protein = polyfit(nucleus_protein_profile(Itrue,2),nucleus_protein_profile(Itrue,3),1);
            r_rescale = p0_protein(1)/p_protein(1)*r_rescale2;
        elseif ~isempty(C0)
            r_rescale = C0(1)/C0(ii)*r_rescale2;
        else
            r_rescale = 1*r_rescale2;
        end
        temp_scale = cat(2,temp_scale,r_rescale);
    % %     r_rescale = 1;
    % %     r_rescale = r_rescale*r_rescale0(ii);
        %%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if isempty(dz_new)
            dz0 = dz;
        else
            dz0 = dz_new;
        end

        z000 = nucleus_protein_profile(Nu_RNA_cell{ii},4)-foci_list1_cell{ii}(:,8);
        I_real_RNA = z000 >= zdiff_nucleus_RNA(1) & z000 <= zdiff_nucleus_RNA(2) & I_RNA1_cell{ii} >= N_thresh;
        [~,ind00] = sort(find(I_real_RNA));
        ind_RNA_new = zeros(size(I_real_RNA));
        ind_RNA_new(I_real_RNA) = ind00;

        [~,I_real1] = spfilter(foci_list0_cell{ii},rxymax,[],relip,Imin0);
    % %     [~,I_real1] = spfilter(foci_list0_cell{ii},ceil(rxmax*resolution0/resolution),ceil(rymax*resolution0/resolution),relip,Imin0);
    % % %     I_real1 = true(size(foci_list0_cell{ii},1),1);
        I_real2 = Nu_en_cell{ii} == Nu_RNA_cell{ii}(Imin_all0_cell{ii});
        I_real3 = ismember(foci_list0_cell{ii}(:,8)-foci_list1_cell{ii}(Imin_all0_cell{ii},8),dz0);
        if t_gate
            load([path0,em_name{ii}(1:id0-1),hist_folder2{pp},em_name{ii}(id0+8:end-1),fit_add,mat_tail],'p0');
            I_real4 = polyval(p0,sqrt(foci_list0_cell{ii}(:,2).*foci_list0_cell{ii}(:,3))) <= foci_list0_cell{ii}(:,1);
        else
            I_real4 = true(size(I_real1));
        end
        I_real5 = ismember(Imin_all0_cell{ii},find(I_real_RNA));
        I_real = I_real1 & I_real2 & I_real3 & I_real4 & I_real5;
        [~,ind00] = sort(find(I_real));
        ind_en_new = zeros(size(I_real));
        ind_en_new(I_real) = ind00;

        foci_list0_cell{ii} = foci_list0_cell{ii}(I_real,:);
        Nu_en_cell{ii} = Nu_en_cell{ii}(I_real,:);
        xyz_all0_cell{ii} = xyz_all0_cell{ii}(I_real,:);
        dxy0_cell{ii} = dxy0_cell{ii}(I_real,:);
        Imin_all0_cell{ii} = Imin_all0_cell{ii}(I_real,:);
        I_RNA0_cell{ii} = I_RNA0_cell{ii}(I_real,:);
    %     N_en0_cell{ii} = N_en0_cell{ii}(I_real,:);
        I_en0_cell{ii} = I_en0_cell{ii}(I_real,:)*r_rescale;
        EL_RNA0_cell{ii} = EL_RNA0_cell{ii}(I_real,:);
        C_RNA0_cell{ii} = C_RNA0_cell{ii}(I_real,:); 
        C_RNA0_cell{ii}(:,1:2) = C_RNA0_cell{ii}(:,1:2)*r_rescale;
        C_RNA_min = mean(C_RNA1_cell{ii}(EL_RNA1_cell{ii} > 0.8 & EL_RNA1_cell{ii} <= 1,2));
        C_RNA0_cell{ii}(:,2) = C_RNA0_cell{ii}(:,2)-C_RNA_min;
        C_RNA_minb = mean(C_RNA1_cell{ii}(EL_RNA1_cell{ii} > 0.8 & EL_RNA1_cell{ii} <= 1,3));
        C_RNA0_cell{ii}(:,3) = C_RNA0_cell{ii}(:,3)-C_RNA_minb;
        dz_en_RNA_cell{ii} = foci_list0_cell{ii}(:,8)-foci_list1_cell{ii}(Imin_all0_cell{ii},8);
        C_RNA0_cell{ii}(:,4) = C_RNA0_cell{ii}(:,4)*r_rescale;

        Imin_all0_cell{ii} = ind_RNA_new(Imin_all0_cell{ii});

        foci_list1_cell{ii} = foci_list1_cell{ii}(I_real_RNA,:);
        Nu_RNA_cell{ii} =  Nu_RNA_cell{ii}(I_real_RNA,:);
        xyz_all1_cell{ii} = xyz_all1_cell{ii}(I_real_RNA,:);
        xyz_all2_cell{ii} = xyz2(I_real_RNA,:);
        dxy1_cell{ii} =  dxy1_cell{ii}(I_real_RNA,:);
        Imin_all1_cell{ii} =  Imin_all1_cell{ii}(I_real_RNA,:);
        I_RNA1_cell{ii} =  I_RNA1_cell{ii}(I_real_RNA,:);
        I_en1_cell{ii} =  I_en1_cell{ii}(I_real_RNA,:)*r_rescale;
    %     I_en2_cell{ii} =  I_en2_cell{ii}(I_real_RNA,:)*r_rescale;
        N_en1_cell{ii} =  N_en1_cell{ii}(I_real_RNA,:);
        EL_RNA1_cell{ii} =  EL_RNA1_cell{ii}(I_real_RNA,:);
        C_RNA1_cell{ii} = C_RNA1_cell{ii}(I_real_RNA,:); 
        C_RNA1_cell{ii}(:,1:2) = C_RNA1_cell{ii}(:,1:2)*r_rescale;
        C_RNA1_cell{ii}(:,2) = C_RNA1_cell{ii}(:,2)-C_RNA_min;
        C_RNA1_cell{ii}(:,3) = C_RNA1_cell{ii}(:,3)-C_RNA_minb;
        foci_circle_area_cell{ii} =  foci_circle_area_cell{ii}(I_real_RNA,:,:);
        C_RNA1_cell{ii}(:,4) = C_RNA1_cell{ii}(:,4)*r_rescale;

        Imin_all1_cell{ii} =  ind_en_new(Imin_all1_cell{ii});


        dxy = pdist2(xyz_all1_cell{ii}(:,1:2),xyz_all0_cell{ii}(:,1:2));
        [dxymin1,I1min1] = min(dxy,[],2);
        dxy1_cell{ii} = dxymin1;
        Imin_all1_cell{ii} = I1min1;
        I_en1_cell{ii} = I_en0_cell{ii}(I1min1); 

        dxyz0 = xyz_all0_cell{ii}-xyz_all1_cell{ii}(Imin_all0_cell{ii}',:);
        dxyz_all0_cell{ii} = dxyz0;

        N_en00 = hist(Imin_all0_cell{ii}'.*(dxy0_cell{ii}' <= r_th),0:size(xyz_all1_cell{ii},1));
        N_en00 = N_en00(2:end);
        N_en0_cell{ii} = N_en00(Imin_all0_cell{ii}')';
        N_en1_cell{ii} = N_en00';
        N_en11 = hist(Imin_all0_cell{ii}'.*(dxy0_cell{ii}' > rout_min & dxy0_cell{ii}' <= rout_max),0:size(xyz_all1_cell{ii},1));
        N_en11 = N_en11(2:end);
        N_en2_cell{ii} = N_en11';

        N_outmin = hist(Imin_all0_cell{ii}'.*(dxy0_cell{ii}' <= rout_min),0:size(xyz_all1_cell{ii},1));
        N_outmin = N_outmin(2:end);
        N_outmin_cell{ii} = N_outmin';

        N_outmax = hist(Imin_all0_cell{ii}'.*(dxy0_cell{ii}' <= rout_max),0:size(xyz_all1_cell{ii},1));
        N_outmax = N_outmax(2:end);
        N_outmax_cell{ii} = N_outmax';

        en_circle_area_cell{ii} = foci_circle_area_cell{ii,1}(Imin_all0_cell{ii},:,:);

        enrich_RNA_cell{ii} = nan(size(xyz_all1_cell{ii},1),2);
    % %     load([path0,em_name{ii}(1:end-1),mat_tail],'foci_data2','fake_data2','max_image','h','r_size','resolution','resolutionz','quanti_p');
        foci_data0 = foci_data2;
        fake_data0 = fake_data2;
        I_foci_data = r_size == rsize0;
        h_area = sum(h{I_foci_data}(:));
        size0 = [size(max_image,1),size(max_image,2),z_range0];
        enrich_xyz = zeros(size(foci_data0{I_foci_data},1),3);
        [enrich_xyz(:,1),enrich_xyz(:,2),enrich_xyz(:,3)] = ind2sub(size0,foci_data0{I_foci_data}(:,11));
        dxyz = pdist2(enrich_xyz.*repmat([resolution,resolution,resolutionz],size(enrich_xyz,1),1),xyz_all1_cell{ii});
        [d12min,I1min] = min(dxyz);
        I1_couple = find(d12min <= epsilonxyz2);
        enrich1 = foci_data0{I_foci_data}(:,4);
    % %     enrich2 = fake_data0{I_foci_data}(:,4);
        enrich2 = foci_data0{I_foci_data}(:,3);
    %     enrich2 = h_area*foci_data0{I_foci_data}(:,9).*foci_data0{I_foci_data}(:,13);
        enrich0 = (enrich1-enrich2)./foci_data0{I_foci_data}(:,15)./foci_data0{I_foci_data}(:,16)/r0_cali*r_rescale;
        enrich_RNA_cell{ii}(I1_couple,1) = enrich0(I1min(I1_couple));
        enrich_all_cell{ii} = [enrich0,foci_data0{I_foci_data}(:,2)*r_rescale,foci_data0{I_foci_data}(:,7)];

        dxyz2 = pdist2(xyz_all1_cell{ii},xyz_all1_cell{ii});
        I2_couple = double(dxyz2 <= epsilonxyz);
        I2_couple(I2_couple ~= 1) = nan;
        enrich_RNA_cell{ii}(:,2) = max(repmat(N_en1_cell{ii},1,length(N_en1_cell{ii})).*I2_couple)';
    end

    I00 = ismember(cycle_em,I_cycle) & ismember(ind_em,I_plot);

    xyz_all0 = cat(1,xyz_all0_cell{I00});
    foci_list0 = cat(1,foci_list0_cell{I00});
    dxy0 = cat(1,dxy0_cell{I00});
    Imin_all0 = cat(1,Imin_all0_cell{I00});
    dxyz_all0 = cat(1,dxyz_all0_cell{I00});
    I_RNA0 = cat(1,I_RNA0_cell{I00});
    I_en0 = cat(1,I_en0_cell{I00})*rescale_en;
    N_en0 = cat(1,N_en0_cell{I00});
    EL_RNA0 = cat(1,EL_RNA0_cell{I00});
    C_RNA0 = cat(1,C_RNA0_cell{I00});

    xyz_all1 = cat(1,xyz_all1_cell{I00});
    xyz_all2 = cat(1,xyz_all2_cell{I00});
    dxy1 = cat(1,dxy1_cell{I00});
    Imin_all1 = cat(1,Imin_all1_cell{I00});
    I_RNA1 = cat(1,I_RNA1_cell{I00});
    I_en1 = cat(1,I_en1_cell{I00})*rescale_en;
    % I_en2 = cat(1,I_en2_cell{I00});
    N_en1{pp} = cat(1,N_en1_cell{I00});
    N_en2{pp} = cat(1,N_en2_cell{I00});
    EL_RNA1 = cat(1,EL_RNA1_cell{I00});
    C_RNA1 = cat(1,C_RNA1_cell{I00});
    enrich_RNA = cat(1,enrich_RNA_cell{I00});
    enrich_all = cat(1,enrich_all_cell{I00});

    N_outmin = cat(1,N_outmin_cell{I00});
    N_outmax = cat(1,N_outmax_cell{I00});

    dz_en_RNA = cat(1,dz_en_RNA_cell{I00});

    foci_circle_area = cat(1,foci_circle_area_cell{I00,1});
    en_circle_area = cat(1,en_circle_area_cell{I00});
    en_circle_diff = diff(en_circle_area,1,2);

    xl0 = foci_circle_area_cell{1,2};
    xl = xl0(1:end-1);
    xr = xl0(2:end);
    dx = xl0(2)-xl0(1);
    x0 = xl+dx/2;

    tl = repmat(dxy0,1,length(xl))-repmat(xl,length(dxy0),1) >= 0;
    tr = repmat(dxy0,1,length(xr))-repmat(xr,length(dxy0),1) < 0;
    [tmax,xind] = max(tl & tr,[],2);
    xind(tmax == 0) = length(xl);
    [~,Nz] = ismember(dz_en_RNA,dz);
    id1 = sub2ind(size(en_circle_diff),[1:length(dxy0)]',xind,Nz);
    en_ring_area = en_circle_diff(id1);
    en_ring_area_max = max(max(en_circle_diff,[],1),[],3)';
    en_ring_ratio = en_ring_area./en_ring_area_max(xind);

    Nz0 = ismember(dz,dz0);
    ind_th = find(abs(xl0-r_th) < (xl0(2)-xl0(1))*epsilon);
    ind_outmin = find(abs(xl0-rout_min) < (xl0(2)-xl0(1))*epsilon);
    ind_outmax = find(abs(xl0-rout_max) < (xl0(2)-xl0(1))*epsilon);
    Sa_foci = mean(foci_circle_area(:,ind_th,Nz0),3);
    Sp_foci = mean(foci_circle_area(:,ind_outmax,Nz0)-foci_circle_area(:,ind_outmin,Nz0),3);
    S_foci_max = max(max(foci_circle_area,[],1),[],3)';
    Sa_foci_ratio = Sa_foci/S_foci_max(ind_th);
    Sp_foci_ratio = Sp_foci/(S_foci_max(ind_outmax)-S_foci_max(ind_outmin));

    id1 = sub2ind(size(en_circle_area),[1:length(dxy0)]',ind_th*ones(length(dxy0),1),Nz);
    id2 = sub2ind(size(en_circle_area),[1:length(dxy0)]',ind_outmin*ones(length(dxy0),1),Nz);
    id3 = sub2ind(size(en_circle_area),[1:length(dxy0)]',ind_outmax*ones(length(dxy0),1),Nz);
    Sa_en = en_circle_area(id1);
    Sp_en = en_circle_area(id3)-en_circle_area(id2);
    S_en_max = max(max(en_circle_area,[],1),[],3)';
    Sa_en_ratio = Sa_en/S_en_max(ind_th);
    Sp_en_ratio = Sp_en/(S_en_max(ind_outmax)-S_en_max(ind_outmin));

    % % r00  = 0.1625*ones(size(ind_em));
    % % r11 = zeros(0);
    % % for ii = find(I00)'
    % %     r11 = cat(1,r11,r00(ii)*ones(size(I_en0_cell{ii})));%/mean(r00);
    % % end
    % % I_en0 = I_en0./r11;
    % % % % r11 = b2_em(index_all0)./b_em(index_all0);
    
xyz_all0_cell2{pp} = xyz_all0_cell;
Imin_all0_cell2{pp} = Imin_all0_cell;
dxy0_cell2{pp} = dxy0_cell;
% foci_circle_area_cell2{pp} = foci_circle_area_cell;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


%% Plot spot number distribution: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ccode = [1,0,0;0,0,1;0.8,0.8,0.8];
angle_bin0 = -pi:0.1:pi; angle_bin = (angle_bin0(1:end-1)+angle_bin0(2:end))/2;
xlim1 = [0,0.8];
% ylim1 = [0,0.6];
ylim1 = [0,20];
ylim1b = [0,150];
xlim1c = [0,0.6];
ylim1c = [0,20];
xlim2 = [-0.5,3.5];
ylim2 = [0,0.8];
ylim21 = [0,1];
xlim3 = [-0.5,3.5];
ylim3 = [-0.5,0.5];

%%
%%% P_N vs. EL    
N_en1_all = [N_en1{:}];
N_en2_all = [N_en2{:}];

EL0 = [0.2,0.4;0.4,0.6;0.75,0.95];
bin0 = 0:1; 
bin1 = 0:3; bin1_label = {'00','01','10','11'};
P_th0_all = zeros(size(EL0,1),length(bin0)^2);
P_out0_all = zeros(size(EL0,1),length(bin0)^2);
P_th0_fit_all = zeros(size(EL0,1),length(bin0)^2);
P_out0_fit_all = zeros(size(EL0,1),length(bin0)^2);
PN0_all = zeros(size(EL0,1),length(bin0)^2);
options = optimoptions(@simulannealbnd,'MaxFunctionEvaluations',15000000,'FunctionTolerance',1e-10);
Ntest = 1;%10;

for ii = 1:size(EL0,1)
    I111 = EL_RNA1 >= EL0(ii,1) & EL_RNA1 <= EL0(ii,2) & Sp_foci_ratio > smin2 & Sa_foci_ratio > smin2;
    P_th0 = hist3(N_en1_all(I111,:),'Ctrs',{bin0,bin0}); P_th1 = P_th0/sum(P_th0(:));
    P_out0 = hist3(N_en2_all(I111,:),'Ctrs',{bin0,bin0}); P_out1 = P_out0/sum(P_out0(:));
    
    [PN00,~] = deconvblind(P_th1,P_out1);
    test0 = [P_out1(:)',PN00(:)']; test0(test0 < 0) = 0;
    
    test_all = zeros(Ntest,length(test0));
    fval_all = zeros(Ntest,1);
    for sss = 1:Ntest
        [test_all(sss,:),fval_all(sss),~] = simulannealbnd(@(x) -sum(sum(log(reshape(x(1:size(x,2)/2),[sqrt(size(x,2)/2),sqrt(size(x,2)/2)])/sum(x(1:size(x,2)/2))+1e-12).*P_out0))-sum(sum(log(conv20(reshape(x(1:size(x,2)/2),[sqrt(size(x,2)/2),sqrt(size(x,2)/2)])/sum(x(1:size(x,2)/2)),reshape(x(size(x,2)/2+1:end),[sqrt(size(x,2)/2),sqrt(size(x,2)/2)])/sum(x(size(x,2)/2+1:end)))+1e-12).*P_th0)),rand(size(test0)),zeros(size(test0)),ones(size(test0)),options);
    end
    [~,Imin] = min(fval_all);
    
    P_out0_fit = reshape(test_all(Imin,1:size(test_all,2)/2)/sum(test_all(Imin,1:size(test_all,2)/2)),length(bin0),length(bin0));
    PN0_fit = reshape(test_all(Imin,size(test_all,2)/2+1:end)/sum(test_all(Imin,size(test_all,2)/2+1:end)),length(bin0),length(bin0));
    P_th_fit = conv20(P_out0_fit,PN0_fit);

    P_out0_all(ii,:) = P_out1(:)';
    P_th0_all(ii,:) = P_th1(:)';
    P_out0_fit_all(ii,:) = P_out0_fit(:)';
    PN0_all(ii,:) = PN0_fit(:)';
    P_th0_fit_all(ii,:) = P_th_fit(:)';
end

%%
asize = [1,0.7];
dsize = [1,1];
figsize = dsize+(dsize+dsize).*[size(EL0,1),2];

figure(1)
set(gcf,'Units','inches','Position',[1,1,figsize])

for ii = 1:size(EL0,1)
    axes('Units','inches','Position',[dsize+(asize+dsize).*[ii-1,0],asize])
        ha1 = bar(bin1,[P_th0_all(ii,:)',P_out0_all(ii,:)']);
        hold on
        plot(bin1+get(ha1(1),'XOffset'),P_th0_fit_all(ii,:),'o-','Color',ha1(1).FaceColor)
        plot(bin1+get(ha1(2),'XOffset'),P_out0_fit_all(ii,:),'o-','Color',ha1(2).FaceColor)
    % %     xlim(xlim2)
        ylim(ylim21)
    xlabel('# Protein spot/locus','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
    ylabel('Probability','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
    title(['EL: ',num2str(EL0(ii,1)),' - ',num2str(EL0(ii,2))],'Interpreter','none','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
    set(gca,'FontName',nfont,'FontSize',fsize2,'FontWeight','normal','Units','normal','XTick',bin1,'XTickLabel',bin1_label)
    % legend('show')
    box off

    axes('Units','inches','Position',[dsize+(asize+dsize).*[ii-1,1],dsize])
        h = heatmap(bin0,bin0,reshape(PN0_all(ii,:,:),[length(bin0),length(bin0)]));
        h.YDisplayData = flipud(h.YDisplayData);
        h.ColorLimits = [0,0.8];
    xlabel('Cad')
    ylabel('Bcd')
    P2 = PN0_all(ii,:); n1 = sum(P2([2,4])); n2 = sum(P2([3,4])); n12 = P2(4); rho = (n12-n1*n2)/sqrt((n1-n1^2))/sqrt((n2-n2^2));
    title(['EL: ',num2str(EL0(ii,1)),' - ',num2str(EL0(ii,2)),', \rho = ',num2str(rho)])
%     set(gca,'FontName',nfont,'FontSize',fsize2,'FontWeight','normal','Units','normal')
    % legend('show')
%     box off
end


%% Plot en-en spot distance distribution:
fmin = @(x,y) x;%min(x,y);
fmax = @(x,y) max(x,y);
% % EL0 = [0:0.1:0.8;0.2:0.1:1]';
EL0 = [0.05,0.95;0.7,0.95];
ELmean = mean(EL0,2)';
ccode0 = mycolors(size(EL0,1)+1);
dxy_en_all = cell(size(ELmean));
Nid_en_all = zeros(size(ELmean));
ind_dist = find(abs(xl0-rdist_max) < (xl0(2)-xl0(1))*epsilon);
r_th2 = 0.28;

for ii = find(I00)'
    dxy_en0 = pdist2(xyz_all0_cell2{1}{ii}(:,1:2),xyz_all0_cell2{2}{ii}(:,1:2));
    dz_en0 = pdist2(xyz_all0_cell2{1}{ii}(:,3),xyz_all0_cell2{2}{ii}(:,3));
    ind_en0 = pdist2(Imin_all0_cell2{1}{ii},Imin_all0_cell2{2}{ii},fmin);
    same_en0 = pdist2(Imin_all0_cell2{1}{ii},Imin_all0_cell2{2}{ii});
    ind_en1 = pdist2([1:size(xyz_all0_cell2{1}{ii},1)]',[1:size(xyz_all0_cell2{2}{ii},1)]',fmin);
    small_en0 = pdist2(dxy0_cell2{1}{ii},dxy0_cell2{2}{ii},fmin);
    large_en0 = pdist2(dxy0_cell2{1}{ii},dxy0_cell2{2}{ii},fmax);
    
    Sd_foci = mean(foci_circle_area_cell{ii,1}(:,ind_dist,Nz0),3);
    Sd_foci_ratio = Sd_foci/S_foci_max(ind_dist);
    
    I_true00 = same_en0 == 0 & small_en0 <= r_th2;% & dz_en0 < resolutionz*1.1 & dz_en0 > resolutionz*0.9;
    dxy_en = dxy_en0(I_true00);
    ind_en = ind_en0(I_true00);
    ind_en11 = ind_en1(I_true00);
% %     large_en1 = large_en0(I_true00);
% %     
% %     I_true11 = large_en1 <= r_th2;
% %     dxy_en = [dxy_en,dxy_en(I_true11)];
% %     ind_en = [ind_en,ind_en(I_true11)];
% %     ind_en11 = [ind_en11,ind_en11(I_true11)];
    
    for jj = 1:length(ELmean)
        I111 = find(EL_RNA1_cell{ii} >= EL0(jj,1) & EL_RNA1_cell{ii} <= EL0(jj,2) & Sd_foci_ratio > smin3);
        ind1 = ismember(ind_en,I111);
        
        if ii == 1
            dxy_en_all{jj} = dxy_en(ind1);
        else
            dxy_en_all{jj} = cat(1,dxy_en_all{jj},dxy_en(ind1));
        end
        Nid_en_all(ii) = numel(unique(ind_en11(ind1)));
    end
end

% % axes('Units','inches','Position',[4*dsize(1)+3*asize(1),2*dsize(2)+1*asize(2),asize])
nar = zeros(length(ELmean),length(x0));
for jj = 1:length(ELmean)
    nar0 = hist0(dxy_en_all{jj},xl,xr);
    nar(jj,:) = nar0./(xr.^2-xl.^2)/pi/sum(Nid_en_all)/resolutionz/length(dz0);
    plot(x0,nar(jj,:),'Color',ccode0(jj,:),'DisplayName',['EL: ',num2str(EL0(jj,1)),' - ',num2str(EL0(jj,2))])
    hold on
end
plot(x0,-diff(nar),'Color',ccode0(length(ELmean)+1,:),'DisplayName','Diff')
xlim(xlim1c)
% % ylim(ylim1c)
xlabel('r (\mum)','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
ylabel('#/\mum^3','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
title('Enriched spot distance distribution (EL)','Interpreter','none','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
legend('show')
set(gca,'FontName',nfont,'FontSize',fsize2,'FontWeight','normal','Units','normal')
box off

%% Output: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist(result_folder,'dir')
    mkdir(result_folder)
end
saveas(1,[result_folder,out_name,enrich_add,distri_add,fig_tail])
saveas(1,[result_folder,out_name,enrich_add,distri_add,fig_tail2],'epsc')
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








