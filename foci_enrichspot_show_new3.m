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
hist_folder2 = 'Histogram_en/';
fit_folder2 = 'Histogram_protein_A/';
% hist_folder_couple = 'Histogram_alignment/';
% hist_folder2_couple = 'Histogram_alignment_RNA2/';
fit_add = '_spot_fit';
fit_add2 = '_peak_fit';
hist_tail = '_raw.xls';
hist_link_tail = '_link.xls';
N_thresh = 0; %0;
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

r_th = 0.28;
rout_min = 0.28;
rout_max = 0.40;
rfit_min = 0.40;
rfit_max = 0.8;
r_ad = 1;
r_area = 1;
rxymax = [2.05,5.25];
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
% % % % % embryo_name = '22445_2_gal4-647_bcd_488_cad_555_par3_new_par3';
embryo_name = '22445_2_hb_cds_TMR_gal4-647_bcd_combine_488_par3_new2';
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
I_plot = [6,8:16];%[1,2,5,8,9,10];%[3:5,8:11,22:25,27:34];%[3:5,8:11,14:19,22:25,27:34];
I_cycle = 10:12;
epsilon = 1e-4;
smin = 0.1;
smin2 = 0.99;
rsize0 = 8;
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

result_folder = ['enrichspot_result/',embryo_name,'/analysis results_test/'];
out_name = [embryo_name,'_C',num2str(I_cycle)];
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

%% Data loading and adjusting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([ana_folder,embryo_name,'/',embryo_name,mat_tail]);
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
% %     load([ana_folder,embryo_name,'/',out_folder0,out_folder0(1:end-2),num2str(ii),mat_tail])
    load([path0,em_name{ii}(1:end-1),mat_tail],'foci_data2','fake_data2','max_image','h','r_size','resolution','resolutionz','quanti_p','nucleus_protein_profile');
    id0 = strfind(em_name{ii},out_folder(1:end-1));
    load([path0,em_name{ii}(1:id0-1),fit_folder2,em_name{ii}(id0+8:end-1),fit_add,mat_tail],'b');
    
    Inten_temp = xlsread([path0,em_name{ii}(1:id0-1),hist_folder_RNA2,em_name{ii}(id0+8:end-1),hist_tail]);
    xyz2 = Inten_temp(:,6:8).*repmat([resolution,resolution,resolutionz],size(Inten_temp,1),1);
    
    
% %     if ~exist([ana_folder,embryo_name,'/',out_folder0])
% %         mkdir([ana_folder,embryo_name,'/',out_folder0])
% %     end
% %     save([ana_folder,embryo_name,'/',out_folder0,out_folder0(1:end-2),num2str(ii),mat_tail],'foci_data2','fake_data2','max_image','h','r_size','resolution','resolutionz','quanti_p','nucleus_protein_profile','b','Inten_temp')
        
    
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
        load([path0,em_name{ii}(1:id0-1),hist_folder2,em_name{ii}(id0+8:end-1),fit_add,mat_tail],'p0');
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
N_en1 = cat(1,N_en1_cell{I00});
N_en2 = cat(1,N_en2_cell{I00});
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
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Plot spot number distribution: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
asize = [1,0.7];
dsize = [1,1];
figsize = dsize+(dsize+asize).*[4,2];
ccode = [1,0,0;0,0,1;0.8,0.8,0.8];
angle_bin0 = -pi:0.1:pi; angle_bin = (angle_bin0(1:end-1)+angle_bin0(2:end))/2;
xlim1 = [0,0.8];
% ylim1 = [0,0.6];
ylim1 = [0,20];
ylim1b = [0,80];
xlim2 = [-0.5,3.5];
ylim2 = [0,0.8];
ylim21 = [0,1];
xlim3 = [-0.5,3.5];
ylim3 = [-0.5,0.5];

figure(1)
set(gcf,'Units','inches','Position',[1,1,figsize])
axes('Units','inches','Position',[dsize(1),dsize(2),asize])
% %     id0 = C_RNA0(:,IC) > 30e-9 & en_ring_ratio > smin;
    id0 = EL_RNA0 >= 0.10 & EL_RNA0 <= 0.35 & en_ring_ratio > smin;
    id1 = EL_RNA1 >= 0.10 & EL_RNA1 <= 0.35;
    xdata = dxy0(id0);
    ydata = 1./en_ring_area(id0)/nnz(id1)/resolutionz/length(dz0);%./length(xdata);
    [~,nar,~] = histw(xdata,ydata,xl,xr);
    pa = polyfit(x0(x0 >= rfit_min & x0 <= rfit_max),nar(x0 >= rfit_min & x0 <= rfit_max),1);
    
% %     id0 = C_RNA0(:,IC) < 5e-9 & en_ring_ratio > smin;
    id0 = EL_RNA0 >= 0.65 & EL_RNA0 <= 0.9 & en_ring_ratio > smin;
    id1 = EL_RNA1 >= 0.65 & EL_RNA1 <= 0.9;
    xdata = dxy0(id0);
    ydata = 1./en_ring_area(id0)./nnz(id1)/resolutionz/length(dz0);%./length(xdata);
    [~,npr,~] = histw(xdata,ydata,xl,xr);
    pp = polyfit(x0(x0 >= rfit_min & x0 <= rfit_max),npr(x0 >= rfit_min & x0 <= rfit_max),1);
    
    plot(x0,nar,'r','DisplayName','Anterior')
    hold on
    plot(x0,polyval(pa,x0),'r-.','DisplayName',['Anterior fit (k = ',num2str(pa(1)),')'])
    plot(x0,npr,'b','DisplayName','Posterior');
    plot(x0,polyval(pp,x0),'b-.','DisplayName',['Posterior fit (k = ',num2str(pp(1)),')'])
    plot(r_th*ones(size(ylim)),ylim,'k--')
    plot(rout_max*ones(size(ylim)),ylim,'k--')
xlim(xlim1)
ylim(ylim1)
xlabel('r (\mum)','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
ylabel('#/\mum^3','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
% title('Distance distribution (EL)','Interpreter','none','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
legend('show')
set(gca,'FontName',nfont,'FontSize',fsize2,'FontWeight','normal','Units','normal')
box off

axes('Units','inches','Position',[2*dsize(1)+1*asize(1),dsize(2),asize])
    theta0 = angle(dxyz_all0(:,1)+1i*dxyz_all0(:,2));
    
    id0 = EL_RNA0 >= 0.10 & EL_RNA0 <= 0.35 & en_ring_ratio > smin;
    id1 = EL_RNA1 >= 0.10 & EL_RNA1 <= 0.35;
    
    id2 = dxy0 >= 0 & dxy0 < 1;
    xdata = theta0(id0 & id2);
    ydata = 1./Sa_en(id0 & id2)/nnz(id1)/resolutionz/length(dz0)/((angle_bin(2)-angle_bin(1))/2/pi);
    [~,nar,~] = histw(xdata,ydata,angle_bin0(1:end-1),angle_bin0(2:end));
    polarplot([angle_bin,angle_bin(1)],[nar,nar(1)],'-r')
    hold on
% % %     
% % %     id2 = dxy0 >= rout_min & dxy0 < rout_max;
% % %     xdata = theta0(id0 & id2);
% % %     ydata = 1./Sp_en(id0 & id2)/nnz(id1)/resolutionz/length(dz0)/((angle_bin(2)-angle_bin(1))/2/pi);
% % %     [~,nar,~] = histw(xdata,ydata,angle_bin0(1:end-1),angle_bin0(2:end));
% % %     polarplot([angle_bin,angle_bin(1)],[nar,nar(1)],'--r')
% % %     hold on
    
    id0 = EL_RNA0 >= 0.65 & EL_RNA0 <= 0.90 & en_ring_ratio > smin;
    id1 = EL_RNA1 >= 0.65 & EL_RNA1 <= 0.90;
    
    id2 = dxy0 >= 0 & dxy0 < 1;
    xdata = theta0(id0 & id2);
    ydata = 1./Sa_en(id0 & id2)/nnz(id1)/resolutionz/length(dz0)/((angle_bin(2)-angle_bin(1))/2/pi);
    [~,npr,~] = histw(xdata,ydata,angle_bin0(1:end-1),angle_bin0(2:end));
    polarplot([angle_bin,angle_bin(1)],[npr,npr(1)],'-b')
    hold on
    
% % %     id2 = dxy0 >= rout_min & dxy0 < rout_max;
% % %     xdata = theta0(id0 & id2);
% % %     ydata = 1./Sp_en(id0 & id2)/nnz(id1)/resolutionz/length(dz0)/((angle_bin(2)-angle_bin(1))/2/pi);
% % %     [~,npr,~] = histw(xdata,ydata,angle_bin0(1:end-1),angle_bin0(2:end));
% % %     polarplot([angle_bin,angle_bin(1)],[npr,npr(1)],'--b')
% % %     hold on

    rlim(ylim1b)
% % xlim(xlim1)
% % ylim(ylim1)
% % xlabel('r (\mum)','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
% % ylabel('#/\mum^3','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
% % % title('Distance distribution (EL)','Interpreter','none','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
% % legend('show')
set(gca,'FontName',nfont,'FontSize',fsize2,'FontWeight','normal','Units','normal')
% % box off
        
axes('Units','inches','Position',[3*dsize(1)+2*asize(1),dsize(2),asize])
    dN = 1;
    N0 = 0:dN:8;
%     id0 = C_RNA1(:,IC) > 10e-9;
    id0 = EL_RNA1 >= 0.2 & EL_RNA1 <= 0.4;
    [na,r0] = hist(N_en1(id0),N0);
%     id0 = C_RNA1(:,IC) < 5e-9;
% %     [np,r0] = hist(N_en2(id0),N0);
% % %     id0 = EL_RNA1 >= 0.10 & EL_RNA1 <= 0.35;
% % %     [np,r0] = hist(N_en2(id0),N0);
    id0 = EL_RNA1 >= 0.6 & EL_RNA1 <= 0.8;
    [np,r0] = hist(N_en1(id0),N0);
    
    b0 = bar(r0,[na'/sum(na),np'/sum(np)]);
    set(b0(1),'EdgeColor','k','FaceColor',ccode(1,:),'DisplayName','Inner')
    set(b0(2),'EdgeColor','k','FaceColor',ccode(2,:),'DisplayName','Outer')
% %     set(b0(1),'EdgeColor','k','FaceColor',ccode(1,:),'DisplayName','Anterior')
% %     set(b0(2),'EdgeColor','k','FaceColor',ccode(2,:),'DisplayName','Posterior')
%     bar(r0,na/sum(na)-np/sum(np),'FaceColor',[0.8,0.8,0.8],'LineColor','k')
%     hold on
%     plot(r0,np/sum(np),'b','DisplayName','Posterior');
    xlim(xlim2)
    ylim(ylim2)
xlabel('# Protein spot/locus','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
ylabel('Probability','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
% title('# En spots distribution (EL)','Interpreter','none','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
set(gca,'FontName',nfont,'FontSize',fsize2,'FontWeight','normal','Units','normal')
% legend('show')
box off

axes('Units','inches','Position',[4*dsize(1)+asize(1)*3,dsize(2),asize])
    b0 = bar(r0,na/sum(na)-np/sum(np),'EdgeColor','k','FaceColor',ccode(3,:));
    xlim(xlim3)
    ylim(ylim3)
xlabel('# Protein spot/locus','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
ylabel('\DeltaP','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
% title('# En spots distribution (EL)','Interpreter','none','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
set(gca,'FontName',nfont,'FontSize',fsize2,'FontWeight','normal','Units','normal')
% legend('show')
box off

% %     [n0,~] = deconv(na'/sum(na),np'/sum(np));
% %     b0 = bar(r0,n0,'EdgeColor','k','FaceColor',ccode(3,:));
% %     xlim(xlim3)
% %     ylim(ylim3)
% % xlabel('# Protein spot/locus','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
% % ylabel('Probability','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
% % % title('# En spots distribution (EL)','Interpreter','none','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
% % set(gca,'FontName',nfont,'FontSize',fsize2,'FontWeight','normal','Units','normal')
% % % legend('show')
% % box off

%%% P_N vs. EL    
EL0 = [0.2,0.4;0.6,0.8];
bin0 = 0:4;
P_th0_all = zeros(size(EL0,1),length(bin0));
P_out0_all = zeros(size(EL0,1),length(bin0));
P_th0_fit_all = zeros(size(EL0,1),length(bin0));
P_out0_fit_all = zeros(size(EL0,1),length(bin0));
PN0_all = zeros(size(EL0,1),length(bin0));
options = optimoptions(@simulannealbnd,'MaxFunctionEvaluations',15000000,'FunctionTolerance',1e-10);
Ntest = 10;
%%
for ii = 1:size(EL0,1)
    I111 = EL_RNA1 >= EL0(ii,1) & EL_RNA1 <= EL0(ii,2) & Sp_foci_ratio > smin2 & Sa_foci_ratio > smin2;
    P_th0 = hist(N_en1(I111),bin0); P_th1 = P_th0/sum(P_th0);
    P_out0 = hist(N_en2(I111),bin0); P_out1 = P_out0/sum(P_out0);
    
    [PN00,~] = deconv([P_th1,zeros(1,size(P_th1,2)-1)],P_out1);
    test0 = [P_out1,PN00(1:length(bin0))]; test0(test0 < 0) = 0;
    
    test_all = zeros(Ntest,length(test0));
    fval_all = zeros(Ntest,1);
    for sss = 1:Ntest
        [test_all(sss,:),fval_all(sss),~] = simulannealbnd(@(x) -sum(sum(log(x(1:size(x,2)/2)/sum(x(1:size(x,2)/2))+1e-12).*P_out0))-sum(sum(log(conv(x(1:size(x,2)/2)/sum(x(1:size(x,2)/2)),x(size(x,2)/2+1:end)/sum(x(size(x,2)/2+1:end)))+1e-12).*[P_th0,zeros(1,size(P_out0,2)-1)])),rand(size(test0)),zeros(size(test0)),ones(size(test0)),options);
    end
    [~,Imin] = min(fval_all);
    
    P_out0_fit = test_all(Imin,1:size(bin0,2))/sum(test_all(Imin,1:size(bin0,2)));
    PN0_fit = test_all(Imin,size(bin0,2)+1:end)/sum(test_all(Imin,size(bin0,2)+1:end));
    P_th_fit = conv(P_out0_fit,PN0_fit);

    P_out0_all(ii,:) = P_out1;
    P_th0_all(ii,:) = P_th1;
    P_out0_fit_all(ii,:) = P_out0_fit;
    PN0_all(ii,:) = PN0_fit;
    P_th0_fit_all(ii,:) = P_th_fit(1:size(bin0,2));
end

axes('Units','inches','Position',[3*dsize(1)+2*asize(1),2*dsize(2)+1*asize(2),asize])
    bar(bin0,PN0_all');
    xlim(xlim2)
    ylim(ylim21)
xlabel('# Enriched protein spot/locus','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
ylabel('Probability','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
% title('# En spots distribution (EL)','Interpreter','none','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
set(gca,'FontName',nfont,'FontSize',fsize2,'FontWeight','normal','Units','normal')
% legend('show')
box off


ii = 1;
axes('Units','inches','Position',[1*dsize(1)+0*asize(1),2*dsize(2)+1*asize(2),asize])
    ha1 = bar(bin0,[P_th0_all(ii,:)',P_out0_all(ii,:)']);
    hold on
    plot(bin0+get(ha1(1),'XOffset'),P_th0_fit_all(ii,:),'o-','Color',ha1(1).FaceColor)
    plot(bin0+get(ha1(2),'XOffset'),P_out0_fit_all(ii,:),'o-','Color',ha1(2).FaceColor)
    xlim(xlim2)
    ylim(ylim21)
xlabel('# Protein spot/locus','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
ylabel('Probability','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
title(['EL: ',num2str(EL0(ii,1)),' - ',num2str(EL0(ii,2))],'Interpreter','none','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
set(gca,'FontName',nfont,'FontSize',fsize2,'FontWeight','normal','Units','normal')
% legend('show')
box off


ii = 2;
axes('Units','inches','Position',[2*dsize(1)+1*asize(1),2*dsize(2)+1*asize(2),asize])
    ha1 = bar(bin0,[P_th0_all(ii,:)',P_out0_all(ii,:)']);
    hold on
    plot(bin0+get(ha1(1),'XOffset'),P_th0_fit_all(ii,:),'o-','Color',ha1(1).FaceColor)
    plot(bin0+get(ha1(2),'XOffset'),P_out0_fit_all(ii,:),'o-','Color',ha1(2).FaceColor)
    xlim(xlim2)
    ylim(ylim21)
xlabel('# Protein spot/locus','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
ylabel('Probability','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
title(['EL: ',num2str(EL0(ii,1)),' - ',num2str(EL0(ii,2))],'Interpreter','none','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
set(gca,'FontName',nfont,'FontSize',fsize2,'FontWeight','normal','Units','normal')
% legend('show')
box off


%% Output: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist(result_folder,'dir')
    mkdir(result_folder)
end
saveas(1,[result_folder,out_name,enrich_add,distri_add,fig_tail])
saveas(1,[result_folder,out_name,enrich_add,distri_add,fig_tail2],'epsc')
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Plot enrichment for different spot numbers: %%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nbin = 30;
rover = 0.6;
EL_range = [0.20,0.70];
% % EL_A = [0.15,0.45];
% % EL_P = [0.55,0.85];
EL0 = [0.20,0.80];
N_range = 0;
bin_min = 0:2:140;   %%% bin_min for [Bcd]
bin_max = 20:2:160;   %%% bin_max for [Bcd]

EL_min = 0:0.01:0.9;   %%% bin_min for EL
EL_max = 0.1:0.01:1;   %%% bin_max for EL

C_range_mean0 = [0,9];
C_range0 = [0,inf];

asize = [1,0.7];
asize2 = [0.4,0.4];
dsize = [1,1];
figsize = dsize+(dsize+asize).*[2,1]+(dsize+asize2).*[1,1];
% ccode = [1,0,0;0,1,0;0,0,1];
ccode = [1,0,0;0,0,0];
r = 0.2;
ccode2 = ccode*r+ones(size(ccode))*(1-r);
ccode3 = [1,0.5,0.5;0.5,0.5,1;0.8,0.8,0.8];
% % legend_name3 = {'All','N_e_n = 0'};
% % tick_label3 = {'Anterior','Posterior'};
tick_label4 = {'N_e_n = 0','All'};
xlim0 = [0,150];
xlim1 = [0,1];
xtick0 = 0:25:150;
xtick1 = 0:0.2:1;
ylim1 = [-5,15];
xlim2 = [0.5,2.5];
ylim2 = [-2,4];

figure(2)
set(gcf,'Units','inches','Position',[1,1,figsize])
axes
set(gca,'Units','inches','Position',[dsize,asize])

    for ii = 1:length(N_range)
        I_all = EL_RNA1 >= EL_range(1) & EL_RNA1 <= EL_range(2) & ~isnan(enrich_RNA(:,1)) & N_en1 == N_range(ii) & enrich_RNA(:,2) <= N_range(ii);
        [xbin,ybin,xerr,yerr] = equal_bin(C_RNA1(I_all,IC)/1e-9,enrich_RNA(I_all,1),Nbin,rover,true);
% %         [xbin,ybin,xerr,yerr] = equal_dist(C_RNA1(I_all,IC)/1e-9,enrich_RNA(I_all,1),bin_min,bin_max);
    % %     [x_N2,y_N2,xerr_N2,yerr_N2] = equal_dist(xx2,yy2,bin_min2,bin_max2);
        patch([xbin,xbin(end:-1:1)],[ybin-yerr,ybin(end:-1:1)+yerr(end:-1:1)],'k','FaceColor',ccode2(ii,:),'EdgeColor','none');
        hold on
        plot(xbin,ybin,'Color',ccode(ii,:),'LineWidth',1,'DisplayName',['N = ',num2str(N_range(ii))])
    end
    % % y = ybin;

    I_all = enrich_all(:,3) >= EL_range(1) & enrich_all(:,3) <= EL_range(2) & ~isnan(enrich_all(:,1));
    % % % [xbin,ybin,xerr,yerr] = equal_bin(enrich_all(:,2)/1e-9,enrich_all(:,1),Nbin,rover,true);
    [xbin,ybin,xerr,yerr] = equal_dist(enrich_all(:,2)/1e-9,enrich_all(:,1),bin_min,bin_max);
    % % I_all = EL_RNA1 >= EL_range(1) & EL_RNA1 <= EL_range(2) & ~isnan(enrich_RNA(:,1));
    % % [xbin,ybin,xerr,yerr] = equal_bin(C_RNA1(I_all,IC)/1e-9,enrich_RNA(I_all,1)/r0_cali,Nbin,rover,true);
    % % [xbin,ybin,xerr,yerr] = equal_dist(C_RNA1(I_all,IC)/1e-9,enrich_RNA(I_all,1),bin_min,bin_max);
    % % % [x_N2,y_N2,xerr_N2,yerr_N2] = equal_dist(xx2,yy2,bin_min2,bin_max2);
    patch([xbin,xbin(end:-1:1)],[ybin-yerr,ybin(end:-1:1)+yerr(end:-1:1)],'k','FaceColor',ccode2(length(N_range)+1,:),'EdgeColor','none');
    hold on
    plot(xbin,ybin,'Color',ccode(length(N_range)+1,:),'LineWidth',1,'DisplayName','All')

xlabel('Bcd concentration (nM)','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
ylabel('Enrichment','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
legend('show')
set(gca,'FontName',nfont,'FontSize',fsize2,'FontWeight','normal','XTick',xtick0)
xlim(xlim0)
ylim(ylim1)
box off
set(gca,'Units','normalized')


axes
set(gca,'Units','inches','Position',[dsize+(dsize+asize).*[1,0],asize])

    for ii = 1:length(N_range)
        I_all = ~isnan(enrich_RNA(:,1)) & N_en1 == N_range(ii) & enrich_RNA(:,2) <= N_range(ii);
        [xbin,ybin,xerr,yerr] = equal_dist(EL_RNA1(I_all),enrich_RNA(I_all,1),EL_min,EL_max);
        patch([xbin,xbin(end:-1:1)],[ybin-yerr,ybin(end:-1:1)+yerr(end:-1:1)],'k','FaceColor',ccode2(ii,:),'EdgeColor','none');
        hold on
        plot(xbin,ybin,'Color',ccode(ii,:),'LineWidth',1,'DisplayName',['N = ',num2str(N_range(ii))])
    end

    I_all = ~isnan(enrich_all(:,1));
    [xbin,ybin,xerr,yerr] = equal_dist(enrich_all(:,3),enrich_all(:,1),EL_min,EL_max);
    patch([xbin,xbin(end:-1:1)],[ybin-yerr,ybin(end:-1:1)+yerr(end:-1:1)],'k','FaceColor',ccode2(length(N_range)+1,:),'EdgeColor','none');
    hold on
    plot(xbin,ybin,'Color',ccode(length(N_range)+1,:),'LineWidth',1,'DisplayName','All')

xlabel('EL','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
ylabel('Enrichment','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
legend('show')
set(gca,'FontName',nfont,'FontSize',fsize2,'FontWeight','normal','XTick',xtick1)
xlim(xlim1)
ylim(ylim1)
box off
set(gca,'Units','normalized')



% % %%% Calculate the spot detection probability
% % y1 = ybin;
% % C_pick = xbin;
% % y0 = mean(y(C_pick >= C_range_mean0(1) & C_pick < C_range_mean0(2)));
% % y(y < y0) = 0;
% % y1(y1 < y0) = 0;
% % P_pick = (y-y1)./(y0-y1);
% % P_pick(isnan(P_pick)) = 1;
% % C_pick = C_pick(2:end);
% % P_pick = P_pick(2:end);
% % 
% % % % P_pick(:) = 1;

axes
set(gca,'Units','inches','Position',[dsize+(dsize+asize).*[2,0],asize2])
% %     I_all = EL_RNA1 >= EL_A(1) & EL_RNA1 <= EL_A(2) & ~isnan(enrich_RNA(:,1)) & N_en1 == N_range(ii) & enrich_RNA(:,2) <= N_range(ii);
% %     mean_En_0(1) = mean(enrich_RNA(I_all,1));
% %     std0_En_0(1) = std0(enrich_RNA(I_all,1));
% %     I_all = enrich_all(:,3) >= EL_A(1) & enrich_all(:,3) <= EL_A(2) & ~isnan(enrich_all(:,1));
% %     mean_En_all(1) = mean(enrich_all(I_all,1));
% %     std0_En_all(1) = std0(enrich_all(I_all,1));
% %     
% %     I_all = EL_RNA1 >= EL_P(1) & EL_RNA1 <= EL_P(2) & ~isnan(enrich_RNA(:,1)) & N_en1 == N_range(ii) & enrich_RNA(:,2) <= N_range(ii);
% %     mean_En_0(2) = mean(enrich_RNA(I_all,1));
% %     std0_En_0(2) = std0(enrich_RNA(I_all,1));
% %     I_all = enrich_all(:,3) >= EL_P(1) & enrich_all(:,3) <= EL_P(2) & ~isnan(enrich_all(:,1));
% %     mean_En_all(2) = mean(enrich_all(I_all,1));
% %     std0_En_all(2) = std0(enrich_all(I_all,1));

% %     b0 = bar([1:length(tick_label3)],[mean_En_all',mean_En_0']);
% %     set(b0(1),'EdgeColor','k','FaceColor',ccode3(1,:),'DisplayName',legend_name3{1})
% %     set(b0(2),'EdgeColor','k','FaceColor',ccode3(2,:),'DisplayName',legend_name3{2})
% %     hold on
% %     x1 = [1:length(tick_label3)]+get(b0(1),'XOffset');
% %     x2 = [1:length(tick_label3)]+get(b0(2),'XOffset');
% %     errorbar(x1,mean_En_all,std0_En_all,'LineStyle','none','Color','k')
% %     errorbar(x2,mean_En_0,std0_En_0,'LineStyle','none','Color','k')
% % set(gca,'FontName',nfont,'FontSize',fsize2,'FontWeight','normal','XTick',[1:length(tick_label3)],'XTickLabel',tick_label3)
% % ylabel('Enrichment','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
% % legend('show')
% % xlim(xlim2)
% % ylim(ylim2)
% % box off
% % set(gca,'Units','normalized')

    I_all = EL_RNA1 >= EL0(1) & EL_RNA1 <= EL0(2) & C_RNA1(:,IC)/1e-9 >= C_range0(1) & C_RNA1(:,IC)/1e-9 <= C_range0(2) & ~isnan(enrich_RNA(:,1)) & ismember(N_en1,N_range) & enrich_RNA(:,2) <= max(N_range);
    mean_En(1) = mean(enrich_RNA(I_all,1));
    std0_En(1) = std0(enrich_RNA(I_all,1));
    I_all = enrich_all(:,3) >= EL0(1) & enrich_all(:,3) <= EL0(2) & enrich_all(:,2)/1e-9 >= C_range0(1) & enrich_all(:,2)/1e-9 <= C_range0(2) & ~isnan(enrich_all(:,1));
    mean_En(2) = mean(enrich_all(I_all,1));
    std0_En(2) = std0(enrich_all(I_all,1));
    
    bar([1:length(tick_label4)],mean_En,'FaceColor',[0.8,0.8,0.8]);
    hold on
    errorbar([1:length(tick_label4)],mean_En,std0_En,'LineStyle','none','Color','k')
    
set(gca,'FontName',nfont,'FontSize',fsize2,'FontWeight','normal','XTick',[1:length(tick_label4)],'XTickLabel',tick_label4)
ylabel('Enrichment','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
% legend('show')
xlim(xlim2)
ylim(ylim2)
box off
set(gca,'Units','normalized')
    
set(2,'Units','normalized')

%% Output: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist(result_folder,'dir')
    mkdir(result_folder)
end
saveas(2,[result_folder,out_name,enrich_add,mean_add,fig_tail])
saveas(2,[result_folder,out_name,enrich_add,mean_add,fig_tail2],'epsc')
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Plot Intensity distribution and Pon vs. EL: %%%%%%%%%%%%%%%%%%%%%%%%%%
sub_pos2 = [5,10];
Nbin = prod(sub_pos2); Nbin2 = 5; dN0 = 2;
rover = 0.4;

% % I_all = EL_RNA1 >= EL_range(1) & EL_RNA1 <= EL_range(2);
% % [C_mean,~,~,~,~,C_min,C_max] = equal_bin(C_RNA1(I_all,IC),C_RNA1(I_all,IC),Nbin,rover,true);
EL_min = 0:0.02:0.90;
EL_max = 0.10:0.02:1;
EL_line = 0:0.01:EL_max(end);
EL_mean = (EL_min+EL_max)/2; 
C_mean = zeros(size(EL_min));
C_std = zeros(size(EL_min));
N_mean = zeros(size(EL_min));
C_line = 0:0.1:120;
C_all = cell(size(EL_min));

bin0 = 0:5;

dI = 0.1;
dw = 0.6;
Il = -dw/2:dI:30;
Ir = Il+dw;
Ic = (Il+Ir)/2;
r0 = Ic;
r1 = r0(1):0.01:r0(end);


% % Npeak = [1,2.2,3.4,4.4,5.6,7,8]';
% % Npeak0 = [0,1,2.2,3.4,4.4,5.6,7]';
% % Npeak1 = [2.2,3.4,4.4,5.6,7,8,9]';
% % % % Npeak = [1,3,6,9,12,15,18,21]';
% % % % Npeak0 = [0,1,3,6,9,12,15,18]';
% % % % Npeak1 = [3,6,9,12,15,18,21,24]';
Npeak = [1:25]';
Npeak0 = [0:24]';
Npeak1 = [2:26]';
% % Npeak = [1:8]';
% % Npeak0 = Npeak-1;
% % Npeak1 = Npeak+1;

Il0 = Npeak-0.55; Il0(1) = 0;
Ir0 = Npeak+0.45;
% % Il0 = [0,1.2942,2.3874,3.4206,4.4177,5.3752,6.4480,7.5344,8.5483,9.4781,10.4353,11.4637,12.5017,13.5360,14.5674,15.5385,16.5417,17.6551,18.4577,19.1951,20.2659,21.2768,22.3401,23.4907,24.6987];
% % Ir0 = [1.2942,2.3874,3.4206,4.4177,5.3752,6.4480,7.5344,8.5483,9.4781,10.4353,11.4637,12.5017,13.5360,14.5674,15.5385,16.5417,17.6551,18.4577,19.1951,20.2659,21.2768,22.3401,23.4907,24.6987,30];
dw0 = 1;

% % b0 = 0.4;
% gau2 = @(a1,a2,a3,a4,a5,b,c,x) a1/sqrt(2*pi)/c*exp(-(x*b).^2/2/c^2)+a2/sqrt(2*pi*2)/c*exp(-(x*2*b).^2/2/2/c^2)+a3/sqrt(2*pi*3)/c*exp(-(x*3*b).^2/2/3/c^2)+a4/sqrt(2*pi*4)/c*exp(-(x*4*b).^2/2/4/c^2)+a5/sqrt(2*pi*5)/c*exp(-(x*5*b).^2/2/5/c^2);   %%% multi-gaussian function to fit distributions
% cmin = [0,0,0,0,0,0.4,0.1];
% cmax = [1,1,1,1,1,0.9,0.5];
% % gau2 = @(b,a1,a2,a3,a4,a5,c1,c2,c3,c4,c5,x) a1/sqrt(2*pi)/c1*exp(-(x*b).^2/2/c1^2)+a2/sqrt(2*pi)/c2*exp(-(x*2*b).^2/2/c2^2)+a3/sqrt(2*pi)/c3*exp(-(x*3*b).^2/2/c3^2)+a4/sqrt(2*pi)/c4*exp(-(x*4*b).^2/2/c4^2)+a5/sqrt(2*pi)/c5*exp(-(x*5*b).^2/2/c5^2);   %%% multi-gaussian function to fit distributions
% % c0 =[10.^(-Npeak),b0/2*ones(size(Npeak))]'; c0 = [b0,c0(:)'];
% % cmin = [zeros(size(Npeak)),b0/6*ones(size(Npeak))]'; cmin = [b0/4,cmin(:)'];
% % cmax = [ones(size(Npeak)),b0*5*ones(size(Npeak))]'; cmax = [b0*2,cmax(:)'];
b0 = 1;
c0 =[10.^(-Npeak),b0*Npeak,b0/2*ones(size(Npeak))]'; c0 = c0(:)';
cmin = [zeros(size(Npeak)),b0*Npeak-b0*(Npeak-Npeak0)*0.45,b0/4*ones(size(Npeak))]'; cmin = cmin(:)';
cmax = [1*ones(size(Npeak)),b0*Npeak+b0*(Npeak1-Npeak)*0.45,b0*3*ones(size(Npeak))]'; cmax = cmax(:)';
% % cmin = [zeros(size(Npeak)),zeros(size(Npeak)),zeros(size(Npeak))]'; cmin = cmin(:)';
% % cmax = [1*ones(size(Npeak)),inf(size(Npeak)),inf(size(Npeak))]'; cmax = cmax(:)';
Nmax = 20;

% gau1_gen = @(n) ['x(',num2str(n),',1)*exp(-(x(',num2str(n),',2)*(xdata(:,3)-x(',num2str(n),',3)).^2))'];   %%% 1D single-gaussian model function text generator
gau1_gen = @(n) ['a',num2str(n),'*exp(-((x-b',num2str(n),')/c',num2str(n),').^2)'];   %%% 1D single-gaussian model function text generator
para_gen = @(n) ['a',num2str(n),',b',num2str(n),',c',num2str(n)];
textfun1 = [];
textpara1 = [];
for n = 1:length(Npeak)
    textfun1 = [textfun1,gau1_gen(n),'+'];
    textpara1 = [textpara1,para_gen(n),','];
end
textfun1 = ['gau1 = @(',textpara1(1:(end-1)),',x) ',textfun1(1:(end-1)),';'];
eval(textfun1);


Pon = zeros(size(EL_min));
Pon_err = zeros(size(EL_min));
Non = zeros(size(EL_min));
Preal = zeros(length(EL_min),length(Il));
Pob = zeros(length(EL_min),length(Il));
Pfake = zeros(length(EL_min),length(Il));
curve_coef = zeros(length(EL_min),length(Npeak)*3);
P_peak_all = zeros(length(EL_min),length(Npeak));
% curve_coef = zeros(length(EL_min),length(Npeak)*2+1);
P_th0_all = zeros(length(EL_min),length(bin0));
P_out0_all = zeros(length(EL_min),length(bin0));
PN0_all = zeros(length(EL_min),length(bin0));
N0_all = zeros(length(EL_min),1);
options = optimoptions(@simulannealbnd,'MaxFunctionEvaluations',15000000,'FunctionTolerance',1e-10);
Ntest = 10;

asize = [1.2,0.7];
asize2 = [0.6,0.9];
dsize = [1,1];
figsize = dsize+(dsize+asize).*[4,1];
figsize2 = dsize+(dsize+asize2).*[9,1];
ccode = [1,0,0;0,0,0];
r = 0.2;
ccode2 = ccode*r+ones(size(ccode))*(1-r);
ccode3 = mycolors(Nbin/Nbin2);
xlim1 = [0,0.8];
ylim1 = [0,20];
xlim3 = [0,1];
ylim3 = [0,0.5];
xlim4 = [0,30];
ylim4 = [-0.025,0.1];
ylim4b = [0,0.1];
xlim41 = [-1,5];
xlim5 = [0,1];
xlim50 = [0,100];
ylim5 = [0,0.6];
xtick2 = 0:0.2:1;
xtick2b = 0:20:100;
xlim5b = [0,1];
ylim5b = [0,20];

figure(30);maximize(30)
set(gcf,'Name',[embryo_name,': Distance distribution fitting'])
figure(40);maximize(40)
set(gcf,'Name',[embryo_name,': difference of spot intensity distribution'])
figure(41);maximize(41)
set(gcf,'Name',[embryo_name,': Spot number distribution'])
figure(50)
set(gcf,'Units','inches','Position',[1,1,figsize])
ax0 = axes;
set(ax0,'Units','inches','Position',[dsize,asize])
figure(60)
set(gcf,'Units','inches','Position',[1,1,figsize2])
for ii = 1:length(EL_min)
%%% Spot distance distribution vs. [TF]
    I000 = EL_RNA0 >= EL_min(ii) & EL_RNA0 <= EL_max(ii) & en_ring_ratio > smin;
    C_mean(ii) = mean(C_RNA0(I000,IC));
    C_std(ii) = std(C_RNA0(I000,IC));
    C_all{ii} = C_RNA0(I000,IC)/1e-9;
    xdata = dxy0(I000);
    ydata = 1./en_ring_area(I000)./length(xdata);
    [~,nar,~] = histw(xdata,ydata,xl,xr);
    figure(30)
    subplot(sub_pos2(1),sub_pos2(2),ii)
        plot(x0,nar,'k','DisplayName','Data')
        hold on
        plot(r_th*ones(size(ylim1)),ylim1,'k--')
        plot(rout_max*ones(size(ylim1)),ylim1,'k--')
    xlim(xlim3)
    ylim(ylim3)
    xlabel('r (\mum)','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
    ylabel('P/\mum^2','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
    title(['EL:',num2str(EL_min(ii)),' - ',num2str(EL_max(ii))],'Interpreter','none','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
    set(gca,'FontName',nfont,'FontSize',fsize2,'FontWeight','normal')
    
%%% Spot intensity distribution vs. [TF]
    I000 = EL_RNA0 >= EL_min(ii) & EL_RNA0 <= EL_max(ii) & dxy0 > rout_min & dxy0 <= rout_max & Sp_en_ratio > smin & Sa_en_ratio > smin;
    I111 = EL_RNA1 >= EL_min(ii) & EL_RNA1 <= EL_max(ii) & Sp_foci_ratio > smin & Sa_foci_ratio > smin;
    I222 = EL_RNA0 >= EL_min(ii) & EL_RNA0 <= EL_max(ii) & dxy0 <= r_th & Sp_en_ratio > smin & Sa_en_ratio > smin;

    xdata = I_en0(I000);
    ydata = 1./Sp_en(I000);
    [~,nf,~] = histw(xdata,ydata,Il,Ir);
    Pfake_all = nf/nnz(I111)*mean(Sa_foci(I111))/dw*r_area;
    %%%%%
    nr = hist0(I_en0(I222),Il,Ir);
    Pob_all = nr/nnz(I111)/dw;
    Preal_all = Pob_all-Pfake_all; 
    xx0 = r0;
    yy0 = Preal_all; yy0(yy0 < 0) = 0;
% %     curve0 = fit(xx0(xx0 <= r0(end))',yy0(xx0 <= r0(end))', ['gauss',num2str(length(Npeak))], 'StartPoint',c0, 'Lower',cmin, 'Upper',cmax);
    curve0 = fit(xx0(xx0 <= r0(end))',yy0(xx0 <= r0(end))', gau1, 'StartPoint',c0, 'Lower',cmin, 'Upper',cmax);
% % %     curve_coef0 = coeffvalues(curve0);
% % % %     rescale0 = curve_coef0(2)/b0;
% % %     dscale0 = curve_coef0(2)-b0;
% % % %     rescale0 = 1;
% % %     
% % % %     xdata = I_en0(I000)/rescale0;
% % %     xdata = I_en0(I000)-dscale0;
% % %     ydata = 1./Sp_en(I000);
% % %     [~,nf,~] = histw(xdata,ydata,Il,Ir);
% % %     Pfake_all = nf/nnz(I111)*mean(Sa_foci(I111))/dw*r_area;
% % %     %%%%%
% % % %     nr = hist0(I_en0(I222)/rescale0,Il,Ir);
% % %     nr = hist0(I_en0(I222)-dscale0,Il,Ir);
% % %     Pob_all = nr/nnz(I111)/dw;
% % %     Preal_all = Pob_all-Pfake_all; 
% % % %     Preal_all = (Pob_all-Pfake_all)/interp1(C_pick,P_pick,EL_mean(ii)/1e-9,'linear',1); 
% % %     Pfake(ii,:) = Pfake_all;
% % %     Pob(ii,:) = Pob_all;
% % %     Preal(ii,:) = Preal_all;
% % %     eval(['xx',num2str(ii),'= r0;'])
% % %     eval(['yy',num2str(ii),'= Preal_all;'])
% % %     xx0 = r0;
% % %     yy0 = Preal_all; yy0(yy0 < 0) = 0;
% % % % %     [P0max,I0max] = max(Preal_all);
% % % % %     b0 = r0(I0max);
% % % % % %     [~,Imax] = min(pdist2(r0',b0*Npeak'));
% % % % % %     Pmax = Preal_all(Imax);
% % % % % %     c0 = [Pmax/2,b0,b0/2];
% % % % % % %     c0 = [0.1,0,0,0,0,0.6,0.3];
% % % % % % % %     curve0 = fit(xx0(xx0 <= 6)',yy0(xx0 <= 6)', gau2, 'StartPoint',c0, 'Lower',cmin, 'Upper',cmax);
% % %     curve0 = fit(xx0(xx0 <= r0(end))',yy0(xx0 <= r0(end))', gau1, 'StartPoint',c0, 'Lower',cmin, 'Upper',cmax);
% % % % %     curve0 = fit(xx0(xx0 <= r0(end))',yy0(xx0 <= r0(end))', ['gauss',num2str(length(Npeak))], 'StartPoint',c0, 'Lower',cmin, 'Upper',cmax);
% % % % %     curve0 = fit(xx0(xx0 <= 6)',yy0(xx0 <= 6)', ['gauss',num2str(length(Npeak))]);
% % % % % % % % %     curve_coef(ii,:) = coeffvalues(curve0);
    tempcoeff = coeffvalues(curve0);
    curve_coef(ii,:) = tempcoeff(:)';
    N_mean(ii) = nnz(I111);
    
% %     I_fit0 = coeffvalues(curve0);
% %     I_fit1 = I_fit0(2:3:end);
% %     I_min = (I_fit1(1:end-1)+I_fit1(2:end))/2;
% %     Il0 = [0,I_min];
% %     Ir0 = [I_min,30];
    
    [~,nf0,~] = histw(xdata,ydata,Il0,Ir0);
    Pfake_all0 = nf0/nnz(I111)*mean(Sa_foci(I111))*r_area;
    nr0 = hist0(I_en0(I222),Il0,Ir0);
    Pob_all0 = nr0/nnz(I111);
    P_peak_all(ii,:) = Pob_all0-Pfake_all0;
    
    %%%%%
    figure(40)
    subplot(sub_pos2(1),sub_pos2(2),ii)
        plot(r0,Preal_all,'ko','MarkerSize',4,'DisplayName','Real');
        hold on
        plot(r1,feval(curve0,r1),'r','DisplayName','Fit')
    xlim(xlim4)
    ylim(ylim4)
    xlabel('En(#)','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
    ylabel('Frequency','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
    title(['EL: ',num2str(EL_min(ii)),' - ',num2str(EL_max(ii)),', C = ',num2str(C_mean(ii)/1e-9),' nM'],'Interpreter','none','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
    set(gca,'FontName',nfont,'FontSize',fsize2,'FontWeight','normal')
    
    if mod(ii,Nbin2) == dN0
        figure(50)
        axes(ax0);
            plot(r0,Preal_all,'Color',ccode3(ceil(ii/Nbin2),:),'DisplayName',['EL: ',num2str(EL_mean(ii))]);
            hold on
            
        figure(60)
        axes('Units','inches','Position',[[(ii-dN0)/Nbin2+1,1].*dsize+[(ii-dN0)/Nbin2,0].*asize2,asize2]);
%             plot(r0,Preal_all,'k.');
            bar(r0,Preal_all,'FaceColor',[0.8,0.8,0.8],'EdgeColor','none','DisplayName','Data');
            hold on
            plot(r1,feval(curve0,r1),'r','DisplayName','Fit')
        xlim(xlim4)
        ylim(ylim4b)
        xlabel('IEn(#)','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
        ylabel('Frequency','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
        title(['EL: ',num2str(EL_mean(ii)),', C = ',num2str(C_mean(ii)/1e-9),' nM, Cstd = ',num2str(C_std(ii)/1e-9),' nM'])
        set(gca,'FontName',nfont,'FontSize',fsize2,'FontWeight','normal','Units','normalized')
        box off
        view(90,-90)
    end
    
%%% Preal vs. EL    
    I111 = EL_RNA1 >= EL_min(ii) & EL_RNA1 <= EL_max(ii) & Sp_foci_ratio > smin & Sa_foci_ratio > smin;
    f_bk = mean((N_outmax(I111)-N_outmin(I111))./Sp_foci(I111))*r_area;
    f_bk_err = std0((N_outmax(I111)-N_outmin(I111))./Sp_foci(I111));
    Pon(ii) = mean(N_en1(I111)-f_bk*Sa_foci(I111));
% %     Pon(ii) = mean(N_en1(I111)-f_bk*Sa_foci(I111))/interp1(C_pick,P_pick,C_mean(ii)/1e-9,'linear',1);
    Pon_err(ii) = std0(N_en1(I111))+f_bk_err*mean(Sa_foci(I111))+f_bk*std0(Sa_foci(I111));
    
%%% P_N vs. EL    
    I111 = EL_RNA1 >= EL_min(ii) & EL_RNA1 <= EL_max(ii) & Sp_foci_ratio > smin2 & Sa_foci_ratio > smin2;
    P_th0 = hist(N_en1(I111),bin0); P_th0 = P_th0/sum(P_th0);
    P_out0 = hist(N_en2(I111),bin0); P_out0 = P_out0/sum(P_out0);
    
    [PN00,~] = deconv([P_th0,zeros(1,size(P_th0,2)-1)],P_out0);
    PN0 = PN00;
% %     test0 = PN00(1:length(bin0)); test0(test0 < 0) = 0;
% %     
% %     test_all = zeros(Ntest,length(test0));
% %     fval_all = zeros(Ntest,1);
% %     for sss = 1:Ntest
% %         [test_all(sss,:),fval_all(sss),~] = simulannealbnd(@(x) -sum(sum(log(conv(P_out0,x/sum(x))+1e-12).*[P_th0,zeros(1,size(x,2)-1)])),test0,zeros(size(test0)),ones(size(test0)),options);
% %     end
% %     [~,Imin] = min(fval_all);
% %     PN0 = test_all(Imin,:)/sum(test_all(Imin,:));

    P_th0_all(ii,:) = P_th0;
    P_out0_all(ii,:) = P_out0;
    PN0_all(ii,:) = PN0;
    N0_all(ii) = nnz(I111);
    
    figure(41)
    subplot(sub_pos2(1),sub_pos2(2),ii)
        bar(bin0,[P_th0',P_out0',PN0']);
    xlim(xlim41)
%     ylim(ylim41)
    xlabel('N(#)','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
    ylabel('Probability','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
    title(['EL: ',num2str(EL_min(ii)),' - ',num2str(EL_max(ii))],'Interpreter','none','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
%     legend('N_t_h','N_o_u_t','N_b_i_n_d')
    set(gca,'FontName',nfont,'FontSize',fsize2,'FontWeight','normal')
    
%%% <Nreal> vs. [TF]    
    I333 = r0 <= Nmax;
    x_temp = r0(I333);
    y_temp = Preal_all(I333); y_temp(y_temp < 0) = 0;
    Non(ii) = y_temp*x_temp'/sum(y_temp);
end

figure(50)
% set(gcf,'Units','inches','Position',[1,1,figsize])
axes(ax0);
xlim(xlim4)
ylim(ylim4)
xlabel('En(#)','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
ylabel('Frequency','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
title('Spot intensity distribution','Interpreter','none','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
set(gca,'FontName',nfont,'FontSize',fsize2,'FontWeight','normal','Units','normalized')
box off

axes
set(gca,'Units','inches','Position',[2*dsize(1)+asize(1),dsize(2),asize])
    plot(EL_mean,Pon,'k.','DisplayName','Data')
hold on
    errorbar(EL_mean,Pon,Pon_err,'k','LineStyle','none','DisplayName','Error')
    g0 = fittype( @(h1,c1,a1,h2,c2,a2,d, x) a1*x.^h1./(c1^h1+x.^h1)+a2*x.^h2./(c2^h2+x.^h2)+d );
    beta0 = [-6,0.5,0.1,-6,0.3,0.1,0];
    beta_min = zeros(size(beta0)); beta_min([1,4]) = 3*beta0([1,4]);
    beta_max = 3*beta0; beta_max([1,4]) = 0;
    curve = fit(EL_mean',Pon', g0, 'StartPoint',beta0,'Lower',beta_min,'Upper',beta_max);
    plot(EL_line,feval(curve,EL_line),'r-','DisplayName','Fit');
xlabel('EL','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
ylabel('Pon','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
title('Spot frequency vs. EL','Interpreter','none','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
legend('show')
set(gca,'FontName',nfont,'FontSize',fsize2,'FontWeight','normal','XTick',xtick2,'Units','normalized')
xlim(xlim5)
ylim(ylim5)
box off

axes
set(gca,'Units','inches','Position',[3*dsize(1)+2*asize(1),dsize(2),asize])
I_true0 = EL_mean >= 0.1 & EL_mean <= 0.60;
    plot(C_mean(I_true0)/1e-9,Pon(I_true0),'k.','DisplayName','Data')
hold on
    errorbar(C_mean(I_true0)/1e-9,Pon(I_true0),Pon_err(I_true0),'k','LineStyle','none','DisplayName','Error')
    g0 = fittype( @(h1,c1,a1,h2,c2,a2,d, x) a1*x.^h1./(c1^h1+x.^h1)+a2*x.^h2./(c2^h2+x.^h2)+d );
    beta0 = [6,10,0.1,6,20,0.1,0];
    beta_min = zeros(size(beta0));
    beta_max = 4*beta0;% beta_max(end) = 1;
    curve = fit(C_mean(I_true0)'/1e-9,Pon(I_true0)', g0, 'StartPoint',beta0,'Lower',beta_min,'Upper',beta_max);
    plot(C_line,feval(curve,C_line),'r-','DisplayName','Fit');
xlabel('[Bcd]','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
ylabel('Pon','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
title('Spot frequency vs. [TF]','Interpreter','none','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
legend('show')
set(gca,'FontName',nfont,'FontSize',fsize2,'FontWeight','normal','XTick',xtick2b,'Units','normalized')
xlim(xlim50)
ylim(ylim5)
box off

axes
set(gca,'Units','inches','Position',[4*dsize(1)+3*asize(1),dsize(2),asize])
    plot(EL_mean,Non)
xlabel('EL','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
ylabel('<Non>','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
title('Spot Intensity vs. [TF]','Interpreter','none','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
% legend('show')
set(gca,'FontName',nfont,'FontSize',fsize2,'FontWeight','normal','XTick',xtick2,'Units','normalized')
xlim(xlim5b)
ylim(ylim5b)
box off

set(50,'Units','normalized')
set(60,'Units','normalized')
%% Output: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist(result_folder,'dir')
    mkdir(result_folder)
end
saveas(30,[result_folder,out_name,enrich_add,distri_add,TF_add,EL_add,fig_tail])
saveas(30,[result_folder,out_name,enrich_add,distri_add,TF_add,EL_add,fig_tail2],'epsc')
saveas(40,[result_folder,out_name,Inten_add,distri_add,TF_add,EL_add,fig_tail])
saveas(40,[result_folder,out_name,Inten_add,distri_add,TF_add,EL_add,fig_tail2],'epsc')
saveas(50,[result_folder,out_name,Pon_add,TF_add,EL_add,fig_tail])
saveas(50,[result_folder,out_name,Pon_add,TF_add,EL_add,fig_tail2],'epsc')
saveas(60,[result_folder,out_name,Inten_add,distri_add,TF_add,EL_add,sel_add,fig_tail])
saveas(60,[result_folder,out_name,Inten_add,distri_add,TF_add,EL_add,sel_add,fig_tail2],'epsc')
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Plot Intensity distribution and Pon vs. [TF]: %%%%%%%%%%%%%%%%%%%%%%%%%%
xlim0 = [0,30];
ylim0 = [0,1];
xlim1 = [0,200];
xlim1b = [0,1];
ylim1 = [0,0.5];
xlim1c = [0,30];
ylim1c = [0,30];
asize = [2,0.9];
asize2 = [0.6,0.9];
asize3 = [0.5,0.5];
dsize = [1,1];
% % figsize = dsize+(dsize+asize).*[1,1]+(dsize+asize2).*[1,0]+(dsize+asize).*[1,0];
figsize = dsize+(dsize+asize).*[1,1]+(dsize+asize2).*[1,0]+(dsize+asize).*[1,0]+(dsize+asize3).*[0,1];

dI = 0.1;
dw = 0.3;
Imin = 0:dI:27;
Imax = Imin+dw;
Ic = Imin(1):0.01:Imin(end);
Nmax = 30;

% % Npeak = [1:8]';
% % a1 = [0.2,1.5,0.5,0.5,0.25,0.1,0.1,0.1]';
% % b1 = 0.1*[3,5,8,10,14,21,26,30]';
% % c1 = 0.25*ones(size(b1));
% % c0 =[a1,b1,c1]'; c0 = c0(:)';
% % cmin = [zeros(size(Npeak)),b1-diff([0;b1])/4,diff([0;b1])/10]'; cmin = cmin(:)';
% % cmax = [3*ones(size(Npeak)),b1+[diff(b1)/2;0.5],diff([0;b1])*5]'; cmax = cmax(:)';

% % Npeak = [1,2.2,3.4,4.4,5.6,7]';
% % Npeak0 = [0,1,2.2,3.4,4.4,5.6]';
% % Npeak1 = [2.2,3.4,4.4,5.6,7,8]';
Npeak = [1:25]';
Npeak0 = [0:24]';
Npeak1 = [2:26]';
% % % Npeak = [1,[1:8]*3]';
% % % Npeak0 = [0,1,[1:7]*3]';
% % % Npeak1 = [3,[2:9]*3]';
b0 = 1;
% % c0 =[10.^(-Npeak),b0*Npeak,b0/2*ones(size(Npeak))]'; c0 = c0(:)';
c0 =[0.1*ones(size(Npeak)),b0*Npeak,b0/2*ones(size(Npeak))]'; c0 = c0(:)';
cmin = [zeros(size(Npeak)),b0*Npeak-b0*(Npeak-Npeak0)*0.5,b0/10*ones(size(Npeak))]'; cmin = cmin(:)';
cmax = [2*ones(size(Npeak)),b0*Npeak+b0*(Npeak1-Npeak)*0.5,b0*3*ones(size(Npeak))]'; cmax = cmax(:)';

% gau1_gen = @(n) ['x(',num2str(n),',1)*exp(-(x(',num2str(n),',2)*(xdata(:,3)-x(',num2str(n),',3)).^2))'];   %%% 1D single-gaussian model function text generator
gau1_gen = @(n) ['a',num2str(n),'*exp(-((x-b',num2str(n),')/c',num2str(n),').^2)'];   %%% 1D single-gaussian model function text generator
para_gen = @(n) ['a',num2str(n),',b',num2str(n),',c',num2str(n)];
textfun1 = [];
textpara1 = [];
for n = 1:length(Npeak)
    textfun1 = [textfun1,gau1_gen(n),'+'];
    textpara1 = [textpara1,para_gen(n),','];
end
textfun1 = ['gau1 = @(',textpara1(1:(end-1)),',x) ',textfun1(1:(end-1)),';'];
eval(textfun1);

% % Icut_min = [4.5,10.5,16.5];
% % Icut_max = [10.5,16.5,30];
% % % % % % % Icut_min = [1:25]-0.5; Icut_min(1) = 0;
% % % % % % % Icut_max = [1:25]+0.5; Icut_min(end) = 30;
% % % Icut_min = [0,2,4.5,7.5,10.5,13.5,16.5,20];
% % % Icut_max = [2,4.5,7.5,10.5,13.5,16.5,20,30];

I_peak0 = curve_coef(:,2:3:end);
P_peak0 = curve_coef(:,1:3:end).*curve_coef(:,3:3:end)*sqrt(pi);
% [xbin,ybin,~] = histw(I_peak,ones(size(I_peak)),0:0.05:5,0.1:0.05:5.1);
[xbin,ybin,~] = histw(I_peak0(:),P_peak0(:),Imin,Imax);
ybin = ybin/max(ybin);
curve1 = fit(xbin(xbin <= Nmax)',ybin(xbin <= Nmax)', gau1,'StartPoint',c0,'Lower',cmin,'Upper',cmax);
coeffvalues(curve1);

% % I_min = Ic(islocalmin(feval(curve1,Ic)));

I_fit0 = coeffvalues(curve1);
I_fit1 = I_fit0(2:3:end);
I_min = (I_fit1(1:end-1)+I_fit1(2:end))/2;

Icut_min = [0,I_min];
Icut_max = [I_min,30];
ccode = mycolors(length(Icut_max)+1);
ccode2 = mycolors(size(curve_coef,2)/3);


P_peak = zeros(size(P_peak0,1),length(Icut_min));
for ii = 1:size(P_peak0,1)
    for jj = 1:size(P_peak0,2)
        kk = find(I_peak0(ii,jj) >= Icut_min & I_peak0(ii,jj) < Icut_max);
        if ~isempty(kk)
            P_peak(ii,kk) = P_peak(ii,kk)+P_peak0(ii,jj);
        end
    end
end
P_peak = [sum(P_peak0,2),P_peak];
P_peak0 = [sum(P_peak0,2),P_peak0];

P_peak_all(P_peak_all <= 0) = epsilon;
P_peak = [sum(P_peak_all,2),P_peak_all];


figure(70)
set(gcf,'Units','inches','Position',[1,1,figsize])
axes('Units','inches','Position',[dsize,asize2]);
    bar(xbin,ybin,'FaceColor',[0.8,0.8,0.8],'EdgeColor','none','DisplayName','Data')
    hold on
    plot(Ic,feval(curve1,Ic),'r','DisplayName','Fit')
xlim(xlim0)
ylim(ylim0)
xlabel('Peak intensity (A.U.)','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
ylabel('Normalized weight','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
legend('show')
set(gca,'FontName',nfont,'FontSize',fsize2,'FontWeight','normal')
box off
view(90,-90)
set(gca,'Units','normalized')

axes('Units','inches','Position',[dsize+(dsize+asize2).*[0,1],asize3]);
    n11 = 1:length(Npeak);
    I110 = coeffvalues(curve1);
    I11 = I110((n11-1)*3+2);
    pI = polyfit(n11,I11,1);
    plot(n11,I11,'k.','DisplayName','Data')
    hold on
    plot(xlim1c,polyval(pI,xlim1c),'r','DisplayName','Fit')
xlim(xlim1c)
ylim(ylim1c)
xlabel('N','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
ylabel('Peak intensity (A.U.)','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
title(['k = ',num2str(pI)],'FontName',nfont,'FontSize',fsize,'FontWeight','normal')
legend('show')
set(gca,'FontName',nfont,'FontSize',fsize2,'FontWeight','normal')
box off
set(gca,'Units','normalized')

axes('Units','inches','Position',[dsize.*[2,1]+asize2.*[1,0],asize]);
I_true0 = EL_mean >= 0.1 & EL_mean <= 0.60;
    for ii = 1:size(P_peak,2)
        plot(C_mean(I_true0)/1e-9,P_peak(I_true0,ii)','Color',ccode(ii,:),'DisplayName',['State ',num2str(ii-1)])
%         plot(EL_mean,P_peak(:,ii)','Color',ccode(ii,:),'DisplayName',['State ',num2str(ii-1)])
        hold on
    end
xlim(xlim1)
% xlim(xlim1b)
ylim(ylim1)
xlabel('Bcd concentration (nM)','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
% xlabel('EL','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
ylabel('Probability','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
legend('show')
set(gca,'FontName',nfont,'FontSize',fsize2,'FontWeight','normal')
box off
set(gca,'Units','normalized')

axes('Units','inches','Position',[dsize.*[3,1]+asize2.*[1,0]+asize.*[1,0],asize]);
    for ii = 1:size(I_peak0,2)
        plot(C_mean/1e-9,I_peak0(:,ii)','Color',ccode2(ii,:),'DisplayName',['State ',num2str(ii-1)])
%         plot(EL_mean,I_peak0(:,ii)','Color',ccode2(ii,:),'DisplayName',['State ',num2str(ii-1)])
        hold on
    end
xlim(xlim1)
% xlim(xlim1b)
ylim(xlim0)
xlabel('Bcd concentration (nM)','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
% xlabel('EL','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
ylabel('Peak intensity (A.U.)','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
legend('show')
set(gca,'FontName',nfont,'FontSize',fsize2,'FontWeight','normal')
box off
set(gca,'Units','normalized')

set(70,'Units','normalized')


xdata = C_mean(I_true0)'/1e-9;
ydata2 = P_peak(I_true0,1:end-1);
ydata2(:,1) = 1-ydata2(:,1);
ydata = ydata2.*repmat(N_mean(I_true0)',1,size(ydata2,2));

N1 = 1:24;
N2 = 2:25;

% % clim = [0,28;0,100;15,110;18,70;45,100;8,40;30,90];
% % clim2 = [8,16;0,0;20,25;0,0;0,0;12,18;0,0];
% % clim = repmat([0,500],length(N1),1);
clim = [0,10^1.3; 0,10^1.35; 10^0.8,10^1.6; 10^1.35,10^2.2; 10^1.35,10^2.2; 10^1.15,10^1.8;...
        10^1.15,10^2.3; 10^1.4,10^2.3; 10^1.2,10^1.8; 10^1.35,10^2.3; 10^1.05,10^2.3; 10^1.65,10^2.3;...
        repmat([0,500],length(N1)/2,1)];
clim2 = [8,16;0,0;20,25;0,0;0,0;12,18;0,0];
% clim = [0,200;0,200;0,200;0,200];
% % xlim1c = [0.5,1.5;0.5,2;1,2.5;1,2;1.6,2.1;0.5,2;1.5,2];
% % ylim1c = [-2.5,-1;-2,4;-1,1;-1.5,0.5;-1,0;-1,0.5;-1.5,0.5];
xlim1c = [0.5,1.5; 0.5,1.5; 0.5,2.5; 1,2.5; 1,2.5; 1,2;...
        0.5,2.5; 1,2.5; 1,2; 1,2.5; 10^1.05,10^2.3; 10^1.65,10^2.3;...
        repmat([-3,3],length(N1)*3/4,1)];
    
ylim1c = [-2.5,-0.5; -1,1; -1,2; -0.5,1.5; -1.5,1; -1,1;...
        -2,1; -1,0; -0.5,1;-1,1; 10^1.05,10^2.3; 10^1.65,10^2.3;...
        repmat([-3,3],length(N1)*3/4,1)];
cplot = 0:0.01:3;
axes_ij = [4,6];
asize3 = [0.45,0.45];
% % figsize2 = dsize+(dsize+asize3).*[length(N1),1];
figsize2 = dsize+(dsize+asize3).*fliplr(axes_ij);

figure(80)
set(gcf,'Units','inches','Position',[1,1,figsize2])
for ii = 1:length(N1)
    axes_ii = mod(ii-1,axes_ij(2))+1;
    axes_jj = floor((ii-1)/axes_ij(2))+1;
    I_true2 = xdata >= clim(ii,1) & xdata <= clim(ii,2);% & ~(xdata >= clim2(ii,1) & xdata <= clim2(ii,2));
    x00 = log10(xdata);
    y00 = log10(ydata2(:,N2(ii))./ydata2(:,N1(ii)));
    p00 = polyfit(x00(I_true2),y00(I_true2),1);
    axes('Units','inches','Position',[dsize+(dsize+asize3).*[axes_ii-1,axes_jj-1],asize3]);
        plot(x00,y00,'.','Color',[0.8,0.8,0.8],'DisplayName','Raw data')
        hold on
        plot(x00(I_true2),y00(I_true2),'k.','DisplayName','Data')
        plot(cplot,polyval(p00,cplot),'r','DisplayName',['Fit: k = ',num2str(p00(1))])
    xlim(xlim1c(ii,:))
    ylim(ylim1c(ii,:))
    xlabel('log[Bcd]','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
    ylabel(['log(P',num2str(N2(ii)-1),'/P',num2str(N1(ii)-1),')'],'FontName',nfont,'FontSize',fsize,'FontWeight','normal')
    title(['k = ',num2str(p00(1))],'FontName',nfont,'FontSize',fsize,'FontWeight','normal')
%     legend('show')
    set(gca,'FontName',nfont,'FontSize',fsize2,'FontWeight','normal')
    box off
    set(gca,'Units','normalized')
end

set(80,'Units','normalized')


N1 = 1:24;
N2 = 2:25;
figure(81); maximize(81)
% % set(gcf,'Units','inches','Position',[1,1,figsize2])
for ii = 1:length(N1)
% %     I_true2 = xdata >= clim(ii,1) & xdata <= clim(ii,2) & ~(xdata >= clim2(ii,1) & xdata <= clim2(ii,2));
    x00 = log10(xdata);
    y00 = log10(ydata2(:,N2(ii))./ydata2(:,N1(ii)));
% %     p00 = polyfit(x00(I_true2),y00(I_true2),1);
    p00 = polyfit(x00,y00,1);
% %     axes('Units','inches','Position',[dsize+(dsize+asize3).*[ii-1,0],asize3]);
    subplot(axes_ij(1),axes_ij(2),ii)
        plot(x00,y00,'.','Color',[0.8,0.8,0.8],'DisplayName','Raw data')
        hold on
% %         plot(x00(I_true2),y00(I_true2),'k.','DisplayName','Data')
        plot(cplot,polyval(p00,cplot),'r','DisplayName',['Fit: k = ',num2str(p00(1))])
% %     xlim(xlim1c(ii,:))
% %     ylim(ylim1c(ii,:))
    xlabel('log[Bcd]','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
    ylabel(['log(P',num2str(N2(ii)),'/P',num2str(N1(ii)),')'],'FontName',nfont,'FontSize',fsize,'FontWeight','normal')
    title(['k = ',num2str(p00(1))],'FontName',nfont,'FontSize',fsize,'FontWeight','normal')
% %     legend('show')
% %     set(gca,'FontName',nfont,'FontSize',fsize2,'FontWeight','normal')
    box off
% %     set(gca,'Units','normalized')
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Fitting the binding curve: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I_true2 = (xdata <= 120 & xdata >= 8) | (xdata >= 0 & xdata <= 4); %I_true2([9,15:17]) = false;
xdata00 = xdata(I_true2);
ydata00 = ydata(I_true2,:);
ydata02 = ydata2(I_true2,:);
cdata00 = C_all(I_true2);
for iii = 1:length(cdata00)
    [~,N_sort] = sort(rand(size(cdata00{iii})));
    cdata00{iii} = cdata00{iii}(N_sort(1:100));
end
% % % % kk0 = [20,1,1,20,2,1,25,3,1,30,3,1,50,3,1,60,3,1,70,3,1,3,2];
% % kk0 = [20,1,3,20,2,3,25,3,3,30,3,3,50,3,3,60,3,3,70,3,3,2];

% % c0 = [20,20,25,30,50,60,70];
% % h0 = [2,2,3,3,3,3,3];
% % ki = ones(size(h0));
% % ki0 = 3*ones(size(h0));
% % k01 = 2;
% % kk00 = [c0;h0;ki;ki0];
% % kk0 = [kk00(:)',k01];
c0 = [20,20,25,30,50,60,70];
h0 = [3,3,3,3,3,3,4];
ki = ones(size(h0));
k10 = 1.5;
% k10 = 0;
ki0 = 1.5;
kk00 = [c0;h0;ki];
kk0 = [kk00(:)',k10,ki0];

Nfit = 500;

% [kk,kk_all,fval_all,exitflag_all] = NSX_binding_fit1(xdata00,ydata00,kk0,Nfit);
[kk,kk_all,fval_all,exitflag_all] = NSX_binding_fit2(xdata00,ydata00,kk0,Nfit);
% [kk,kk_all,fval_all,exitflag_all] = NSX_binding_fit2b(cdata00,ydata00,kk0,Nfit);

% % c0 = kk(1:4:end-1);
% % h0 = kk(2:4:end-1);
% % ki = kk(3:4:end-1);
% % ki0 = kk(4:4:end-1);
% % k01 = kk(end);
c0f = kk(1:3:end-2);
h0f = kk(2:3:end-2);
kif = kk(3:3:end-2);
ki0f = zeros(size(kif)); ki0f(1) = kk(end-1); ki0f(2:end) = kk(end);
k01f = 1;

plotN = {1,2,3,4,5,6,7,8};
cx = [0:0.01:150]';
% % cy = NSX_binding1(kk,cx);
cy = NSX_binding2(kk,cx);
xlim2 = [0,100];
ylim2 = [0.5,1;0,0.2;0,0.2;0,0.2;0,0.2;0,0.06;0,0.06;0,0.06];

asize4 = [2,0.45];
dsize2 = [0.5,0.5];
figsize3 = dsize2+(dsize2+asize4).*[1,length(plotN)];

figure(90)
set(gcf,'Units','inches','Position',[1,1,figsize3])
for ii = 1:length(plotN)
    axes('Units','inches','Position',[dsize2+(dsize2+asize4).*[0,ii-1],asize4]);
        plot(xdata,ydata2(:,plotN{ii}),'.')
        hold on
        plot(cx,cy(:,plotN{ii}))
    xlim(xlim2)
    ylim(ylim2(ii,:))
    if any(ismember(plotN{ii},1))
%         title(['N = ',num2str(plotN{ii}),', kn0 = ',num2str(kn0),'k01 = ',num2str(k01)])
        title(['N = ',num2str(plotN{ii}),', k01 = ',num2str(k01f)])
    else
%         title(['N = ',num2str(plotN{ii}),', C0 = ',num2str(c0(plotN{ii}-1)),'h0 = ',num2str(h0(plotN{ii}-1))])
        title(['N = ',num2str(plotN{ii}),', C0 = ',num2str(c0f(plotN{ii}-1)),', h = ',num2str(h0f(plotN{ii}-1)),', ki = ',num2str(kif(plotN{ii}-1)),', ki0 = ',num2str(ki0f(plotN{ii}-1))])
    end

    xlabel('[Bcd] (nM)','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
    ylabel('Probability','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
    set(gca,'FontName',nfont,'FontSize',fsize2,'FontWeight','normal')
    box off
    set(gca,'Units','normalized')
end

set(90,'Units','normalized','Renderer','Painter')

%% Output: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist(result_folder,'dir')
    mkdir(result_folder)
end
saveas(70,[result_folder,out_name,Inten_add,peak_add,fig_tail])
saveas(70,[result_folder,out_name,Inten_add,peak_add,fig_tail2],'epsc')
%%
saveas(80,[result_folder,out_name,Inten_add,fit_add2,fig_tail])
saveas(80,[result_folder,out_name,Inten_add,fit_add2,fig_tail2],'epsc')
%%
saveas(90,[result_folder,out_name,binding_add,fit_add2,fig_tail])
saveas(90,[result_folder,out_name,binding_add,fit_add2,fig_tail2],'epsc')
save([result_folder,fit_name,mat_tail],'kk_all','fval_all','exitflag_all')
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Plot Intensity distribution and Pon vs. [TF]: %%%%%%%%%%%%%%%%%%%%%%%%%%
sub_pos2 = [5,8];
Nbin = prod(sub_pos2); Nbin2 = 6; dN0 = 2;
rover = 0.7;
EL_range = [0.1,0.70];

I_all = EL_RNA1 >= EL_range(1) & EL_RNA1 <= EL_range(2);
[C_mean,~,~,~,~,C_min,C_max] = equal_bin(C_RNA1(I_all,IC),C_RNA1(I_all,IC),Nbin,rover,true);
% % C_min = [0:3:6,6:3:114]*1e-9; 
% % C_max = [8:3:14,22:3:130]*1e-9; 
% % C_mean = (C_min+C_max)/2; 
C_line = 0:0.1:C_max(end)/1e-9;
% C_min(1) = -inf; 
% C_max(end) = inf;
N_mean = zeros(size(C_min));

dI = 0.5;
dw = 1.6;
Il = -dw/2:dI:30;
Ir = Il+dw;
Ic = (Il+Ir)/2;
r0 = Ic;
r1 = r0(1):0.01:r0(end);

% % Npeak = [1,2.2,3.4,4.4,5.6,7,8]';
% % Npeak0 = [0,1,2.2,3.4,4.4,5.6,7]';
% % Npeak1 = [2.2,3.4,4.4,5.6,7,8,9]';
Npeak = [1,3,6,9,12,15,18,21]';
Npeak0 = [0,1,3,6,9,12,15,18]';
Npeak1 = [3,6,9,12,15,18,21,24]';
% % Npeak = [1:8]';
% % Npeak0 = Npeak-1;
% % Npeak1 = Npeak+1;

% % b0 = 0.4;
% gau2 = @(a1,a2,a3,a4,a5,b,c,x) a1/sqrt(2*pi)/c*exp(-(x*b).^2/2/c^2)+a2/sqrt(2*pi*2)/c*exp(-(x*2*b).^2/2/2/c^2)+a3/sqrt(2*pi*3)/c*exp(-(x*3*b).^2/2/3/c^2)+a4/sqrt(2*pi*4)/c*exp(-(x*4*b).^2/2/4/c^2)+a5/sqrt(2*pi*5)/c*exp(-(x*5*b).^2/2/5/c^2);   %%% multi-gaussian function to fit distributions
% cmin = [0,0,0,0,0,0.4,0.1];
% cmax = [1,1,1,1,1,0.9,0.5];
% % gau2 = @(b,a1,a2,a3,a4,a5,c1,c2,c3,c4,c5,x) a1/sqrt(2*pi)/c1*exp(-(x*b).^2/2/c1^2)+a2/sqrt(2*pi)/c2*exp(-(x*2*b).^2/2/c2^2)+a3/sqrt(2*pi)/c3*exp(-(x*3*b).^2/2/c3^2)+a4/sqrt(2*pi)/c4*exp(-(x*4*b).^2/2/c4^2)+a5/sqrt(2*pi)/c5*exp(-(x*5*b).^2/2/c5^2);   %%% multi-gaussian function to fit distributions
% % c0 =[10.^(-Npeak),b0/2*ones(size(Npeak))]'; c0 = [b0,c0(:)'];
% % cmin = [zeros(size(Npeak)),b0/6*ones(size(Npeak))]'; cmin = [b0/4,cmin(:)'];
% % cmax = [ones(size(Npeak)),b0*5*ones(size(Npeak))]'; cmax = [b0*2,cmax(:)'];
b0 = 1;
c0 =[10.^(-Npeak),b0*Npeak,b0/2*ones(size(Npeak))]'; c0 = c0(:)';
cmin = [zeros(size(Npeak)),b0*Npeak-b0*(Npeak-Npeak0)*0.5,b0/5*ones(size(Npeak))]'; cmin = cmin(:)';
cmax = [1*ones(size(Npeak)),b0*Npeak+b0*(Npeak1-Npeak)*0.5,b0*2*ones(size(Npeak))]'; cmax = cmax(:)';
% cmax = [1*ones(size(Npeak)),b0*Npeak+b0*(Npeak1-Npeak)*0.4,b0*Npeak/2]'; cmax = cmax(:)';
% % cmin = [zeros(size(Npeak)),zeros(size(Npeak)),zeros(size(Npeak))]'; cmin = cmin(:)';
% % cmax = [1*ones(size(Npeak)),inf(size(Npeak)),inf(size(Npeak))]'; cmax = cmax(:)';
Nmax = 20;

Pon = zeros(size(C_min));
Non = zeros(size(C_min));
Preal = zeros(length(C_min),length(Il));
Pob = zeros(length(C_min),length(Il));
Pfake = zeros(length(C_min),length(Il));
curve_coef = zeros(length(C_min),length(Npeak)*3);
% curve_coef = zeros(length(C_min),length(Npeak)*2+1);

asize = [1.2,0.7];
asize2 = [0.6,0.9];
dsize = [1,1];
figsize = dsize+(dsize+asize).*[3,1];
figsize2 = dsize+(dsize+asize2).*[8,1];
ccode = [1,0,0;0,0,0];
r = 0.2;
ccode2 = ccode*r+ones(size(ccode))*(1-r);
ccode3 = mycolors(ceil(Nbin/Nbin2));
xlim3 = [0,1];
ylim3 = [0,0.5];
xlim4 = [0,30];
ylim4 = [-0.025,0.1];
ylim4b = [0,0.1];
xlim5 = [0,100];
ylim5 = [0,0.6];
xtick2 = 0:25:200;
xlim5b = [0,200];
ylim5b = [0,20];

figure(3);maximize(3)
set(gcf,'Name',[embryo_name,': Distance distribution fitting'])
figure(4);maximize(4)
set(gcf,'Name',[embryo_name,': difference of spot intensity distribution'])
figure(5)
set(gcf,'Units','inches','Position',[1,1,figsize])
ax0 = axes;
set(ax0,'Units','inches','Position',[dsize,asize])
figure(6)
set(gcf,'Units','inches','Position',[1,1,figsize2])
for ii = 1:length(C_min)
%%% Spot distance distribution vs. [TF]
    I000 = EL_RNA0 >= EL_range(1) & EL_RNA0 <= EL_range(2) & C_RNA0(:,IC) >= C_min(ii) & C_RNA0(:,IC) <= C_max(ii) & en_ring_ratio > smin;
    xdata = dxy0(I000);
    ydata = 1./en_ring_area(I000)./length(xdata);
    [~,nar,~] = histw(xdata,ydata,xl,xr);
    figure(3)
    subplot(sub_pos2(1),sub_pos2(2),ii)
        plot(x0,nar,'k','DisplayName','Data')
        hold on
        plot(r_th*ones(size(ylim1)),ylim1,'k--')
        plot(rout_max*ones(size(ylim1)),ylim1,'k--')
    xlim(xlim3)
    ylim(ylim3)
    xlabel('r (\mum)','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
    ylabel('P/\mum^2','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
    title([num2str(C_min(ii)/1e-9),' nM - ',num2str(C_max(ii)/1e-9),' nM'],'Interpreter','none','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
    set(gca,'FontName',nfont,'FontSize',fsize2,'FontWeight','normal')
    
%%% Spot intensity distribution vs. [TF]
    I000 = EL_RNA0 >= EL_range(1) & EL_RNA0 <= EL_range(2) & C_RNA0(:,IC) >= C_min(ii) & C_RNA0(:,IC) <= C_max(ii) & dxy0 > rout_min & dxy0 <= rout_max & Sp_en_ratio > smin & Sa_en_ratio > smin;
    I111 = EL_RNA1 >= EL_range(1) & EL_RNA1 <= EL_range(2) & C_RNA1(:,IC) >= C_min(ii) & C_RNA1(:,IC) <= C_max(ii) & Sp_foci_ratio > smin & Sa_foci_ratio > smin;
    I222 = EL_RNA0 >= EL_range(1) & EL_RNA0 <= EL_range(2) & C_RNA0(:,IC) >= C_min(ii) & C_RNA0(:,IC) <= C_max(ii) & dxy0 <= r_th & Sp_en_ratio > smin & Sa_en_ratio > smin;

    xdata = I_en0(I000);
    ydata = 1./Sp_en(I000);
    [~,nf,~] = histw(xdata,ydata,Il,Ir);
    Pfake_all = nf/nnz(I111)*mean(Sa_foci(I111))/dw*r_area;
    %%%%%
    nr = hist0(I_en0(I222),Il,Ir);
    Pob_all = nr/nnz(I111)/dw;
    Preal_all = Pob_all-Pfake_all; 
    xx0 = r0;
    yy0 = Preal_all; yy0(yy0 < 0) = 0;
    curve0 = fit(xx0(xx0 <= r0(end))',yy0(xx0 <= r0(end))', ['gauss',num2str(length(Npeak))], 'StartPoint',c0, 'Lower',cmin, 'Upper',cmax);
    curve_coef0 = coeffvalues(curve0);
%     rescale0 = curve_coef0(2)/b0;
    dscale0 = curve_coef0(2)-b0;
%     rescale0 = 1;
    
%     xdata = I_en0(I000)/rescale0;
    xdata = I_en0(I000)-dscale0;
    ydata = 1./Sp_en(I000);
    [~,nf,~] = histw(xdata,ydata,Il,Ir);
    Pfake_all = nf/nnz(I111)*mean(Sa_foci(I111))/dw*r_area;
    %%%%%
%     nr = hist0(I_en0(I222)/rescale0,Il,Ir);
    nr = hist0(I_en0(I222)-dscale0,Il,Ir);
    Pob_all = nr/nnz(I111)/dw;
    Preal_all = Pob_all-Pfake_all; 
%     Preal_all = (Pob_all-Pfake_all)/interp1(C_pick,P_pick,C_mean(ii)/1e-9,'linear',1); 
    Pfake(ii,:) = Pfake_all;
    Pob(ii,:) = Pob_all;
    Preal(ii,:) = Preal_all;
    eval(['xx',num2str(ii),'= r0;'])
    eval(['yy',num2str(ii),'= Preal_all;'])
    xx0 = r0;
    yy0 = Preal_all; yy0(yy0 < 0) = 0;
% %     [P0max,I0max] = max(Preal_all);
% %     b0 = r0(I0max);
% % %     [~,Imax] = min(pdist2(r0',b0*Npeak'));
% % %     Pmax = Preal_all(Imax);
% % %     c0 = [Pmax/2,b0,b0/2];
% % % %     c0 = [0.1,0,0,0,0,0.6,0.3];
% % % % %     curve0 = fit(xx0(xx0 <= 6)',yy0(xx0 <= 6)', gau2, 'StartPoint',c0, 'Lower',cmin, 'Upper',cmax);
    curve0 = fit(xx0(xx0 <= r0(end))',yy0(xx0 <= r0(end))', ['gauss',num2str(length(Npeak))], 'StartPoint',c0, 'Lower',cmin, 'Upper',cmax);
% %     curve0 = fit(xx0(xx0 <= 6)',yy0(xx0 <= 6)', ['gauss',num2str(length(Npeak))]);
    curve_coef(ii,:) = coeffvalues(curve0);
    N_mean(ii) = nnz(I111);
    
    %%%%%
    figure(4)
    subplot(sub_pos2(1),sub_pos2(2),ii)
        plot(r0,Preal_all,'ko','MarkerSize',4,'DisplayName','Real');
        hold on
        plot(r1,feval(curve0,r1),'r','DisplayName','Fit')
    xlim(xlim4)
    ylim(ylim4)
    xlabel('En(#)','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
    ylabel('Frequency','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
    title([num2str(C_min(ii)/1e-9),' nM - ',num2str(C_max(ii)/1e-9),' nM'],'Interpreter','none','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
    set(gca,'FontName',nfont,'FontSize',fsize2,'FontWeight','normal')
    
    if mod(ii,Nbin2) == dN0
        figure(5)
        axes(ax0);
            plot(r0,Preal_all,'Color',ccode3(ceil(ii/Nbin2),:),'DisplayName',[num2str(C_mean(ii)/1e-9),' nM']);
            hold on
            
        figure(6)
        axes('Units','inches','Position',[[(ii-dN0)/Nbin2+1,1].*dsize+[(ii-dN0)/Nbin2,0].*asize2,asize2]);
%             plot(r0,Preal_all,'k.');
            bar(r0,Preal_all,'FaceColor',[0.8,0.8,0.8],'EdgeColor','none','DisplayName','Data');
            hold on
            plot(r1,feval(curve0,r1),'r','DisplayName','Fit')
        xlim(xlim4)
        ylim(ylim4b)
        xlabel('IEn(#)','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
        ylabel('Frequency','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
        title([num2str(C_mean(ii)/1e-9),' nM'])
        set(gca,'FontName',nfont,'FontSize',fsize2,'FontWeight','normal')
        box off
        view(90,-90)
    end
    
%%% Preal vs. [TF]    
    I111 = EL_RNA1 >= EL_range(1) & EL_RNA1 <= EL_range(2) & C_RNA1(:,IC) >= C_min(ii) & C_RNA1(:,IC) <= C_max(ii) & Sp_foci_ratio > smin & Sa_foci_ratio > smin;
    f_bk = mean((N_outmax(I111)-N_outmin(I111))./Sp_foci(I111))*r_area;
    Pon(ii) = mean(N_en1(I111)-f_bk*Sa_foci(I111));
% %     Pon(ii) = mean(N_en1(I111)-f_bk*Sa_foci(I111))/interp1(C_pick,P_pick,C_mean(ii)/1e-9,'linear',1);
    
%%% <Nreal> vs. [TF]    
    I333 = r0 <= Nmax;
    x_temp = r0(I333);
    y_temp = Preal_all(I333); y_temp(y_temp < 0) = 0;
    Non(ii) = y_temp*x_temp'/sum(y_temp);
end

figure(5)
% set(gcf,'Units','inches','Position',[1,1,figsize])
axes(ax0);
xlim(xlim4)
ylim(ylim4)
xlabel('En(#)','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
ylabel('Frequency','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
title('Spot intensity distribution','Interpreter','none','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
set(gca,'FontName',nfont,'FontSize',fsize2,'FontWeight','normal','Units','normalized')
box off

axes
set(gca,'Units','inches','Position',[2*dsize(1)+asize(1),dsize(2),asize])
    plot(C_mean/1e-9,Pon,'k.','DisplayName','Data')
hold on
    g0 = fittype( @(h1,c1,a1,h2,c2,a2,d, x) a1*x.^h1./(c1^h1+x.^h1)+a2*x.^h2./(c2^h2+x.^h2)+d );
    beta0 = [3,10,0.1,6,40,0.3,0];
    beta_min = zeros(size(beta0));
    beta_max = 3*beta0;% beta_max(end) = 1;
    curve = fit(C_mean'/1e-9,Pon', g0, 'StartPoint',beta0,'Lower',beta_min,'Upper',beta_max);
    plot(C_line,feval(curve,C_line),'r-','DisplayName','Fit');

    g0 = fittype( @(h1,c1,a1,d, x) a1*x.^h1./(c1^h1+x.^h1)+d );
    beta0 = [6,30,0.4,0];
    beta_min = zeros(size(beta0));
    beta_max = 3*beta0;% beta_max(end) = 1;
    curve2 = fit(C_mean'/1e-9,Pon', g0, 'StartPoint',beta0,'Lower',beta_min,'Upper',beta_max);
    plot(C_line,feval(curve2,C_line),'b-','DisplayName','Fit2');

xlabel('Bcd concentration (nM)','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
ylabel('Pon','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
title('Spot frequency vs. [TF]','Interpreter','none','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
legend('show')
set(gca,'FontName',nfont,'FontSize',fsize2,'FontWeight','normal','XTick',xtick2,'Units','normalized')
xlim(xlim5)
ylim(ylim5)
box off

axes
set(gca,'Units','inches','Position',[3*dsize(1)+2*asize(1),dsize(2),asize])
    plot(C_mean/1e-9,Non)
xlabel('Bcd concentration (nM)','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
ylabel('<Non>','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
title('Spot Intensity vs. [TF]','Interpreter','none','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
% legend('show')
set(gca,'FontName',nfont,'FontSize',fsize2,'FontWeight','normal','XTick',xtick2,'Units','normalized')
xlim(xlim5b)
ylim(ylim5b)
box off

set(5,'Units','normalized')
set(6,'Units','normalized')
%% Output: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist(result_folder,'dir')
    mkdir(result_folder)
end
saveas(3,[result_folder,out_name,enrich_add,distri_add,TF_add,fig_tail])
saveas(3,[result_folder,out_name,enrich_add,distri_add,TF_add,fig_tail2],'epsc')
saveas(4,[result_folder,out_name,Inten_add,distri_add,TF_add,fig_tail])
saveas(4,[result_folder,out_name,Inten_add,distri_add,TF_add,fig_tail2],'epsc')
saveas(5,[result_folder,out_name,Pon_add,TF_add,fig_tail])
saveas(5,[result_folder,out_name,Pon_add,TF_add,fig_tail2],'epsc')
saveas(6,[result_folder,out_name,Inten_add,distri_add,TF_add,sel_add,fig_tail])
saveas(6,[result_folder,out_name,Inten_add,distri_add,TF_add,sel_add,fig_tail2],'epsc')
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Plot Intensity distribution and Pon vs. [TF]: %%%%%%%%%%%%%%%%%%%%%%%%%%
xlim0 = [0,30];
ylim0 = [0,1];
xlim1 = [0,200];
xlim1b = [0,1];
ylim1 = [0,0.5];
asize = [2,0.9];
asize2 = [0.6,0.9];
dsize = [1,1];
figsize = dsize+(dsize+asize).*[1,1]+(dsize+asize2).*[1,0]+(dsize+asize).*[1,0];

dI = 0.1;
dw = 0.3;
Imin = 0:dI:24;
Imax = Imin+dw;
Ic = Imin(1):0.01:Imin(end);
Nmax = 30;

% % Npeak = [1:8]';
% % a1 = [0.2,1.5,0.5,0.5,0.25,0.1,0.1,0.1]';
% % b1 = 0.1*[3,5,8,10,14,21,26,30]';
% % c1 = 0.25*ones(size(b1));
% % c0 =[a1,b1,c1]'; c0 = c0(:)';
% % cmin = [zeros(size(Npeak)),b1-diff([0;b1])/4,diff([0;b1])/10]'; cmin = cmin(:)';
% % cmax = [3*ones(size(Npeak)),b1+[diff(b1)/2;0.5],diff([0;b1])*5]'; cmax = cmax(:)';

% % % % % Npeak = [1,2.2,3.4,4.4,5.6,7]';
% % % % % Npeak0 = [0,1,2.2,3.4,4.4,5.6]';
% % % % % Npeak1 = [2.2,3.4,4.4,5.6,7,8]';
% % % Npeak = [1:8]'*3;
% % % Npeak0 = [0:7]'*3;
% % % Npeak1 = [2:9]'*3;
% % % b0 = 1;
% % % c0 =[10.^(-Npeak),b0*Npeak,b0/2*ones(size(Npeak))]'; c0 = c0(:)';
% % % cmin = [zeros(size(Npeak)),b0*Npeak-b0*(Npeak-Npeak0)*0.4,b0/10*ones(size(Npeak))]'; cmin = cmin(:)';
% % % cmax = [2*ones(size(Npeak)),b0*Npeak+b0*(Npeak1-Npeak)*0.4,b0*1*ones(size(Npeak))]'; cmax = cmax(:)';

Npeak = [1,[1:8]*3]';
Npeak0 = [0,1,[1:7]*3]';
Npeak1 = [3,[2:9]*3]';
b0 = 1;
c0 =[10.^(-Npeak),b0*Npeak,b0/2*ones(size(Npeak))]'; c0 = c0(:)';
cmin = [zeros(size(Npeak)),b0*Npeak-b0*(Npeak-Npeak0)*0.4,b0/10*ones(size(Npeak))]'; cmin = cmin(:)';
cmax = [2*ones(size(Npeak)),b0*Npeak+b0*(Npeak1-Npeak)*0.4,b0*2*ones(size(Npeak))]'; cmax = cmax(:)';

% gau1_gen = @(n) ['x(',num2str(n),',1)*exp(-(x(',num2str(n),',2)*(xdata(:,3)-x(',num2str(n),',3)).^2))'];   %%% 1D single-gaussian model function text generator
gau1_gen = @(n) ['a',num2str(n),'*exp(-((x-b',num2str(n),')/c',num2str(n),').^2)'];   %%% 1D single-gaussian model function text generator
para_gen = @(n) ['a',num2str(n),',b',num2str(n),',c',num2str(n)];
textfun1 = [];
textpara1 = [];
for n = 1:length(Npeak)
    textfun1 = [textfun1,gau1_gen(n),'+'];
    textpara1 = [textpara1,para_gen(n),','];
end
textfun1 = ['gau1 = @(',textpara1(1:(end-1)),',x) ',textfun1(1:(end-1)),';'];
eval(textfun1);

% % % % Icut_min = [4.5,10.5,16.5];
% % % % Icut_max = [10.5,16.5,30];
% % Icut_min = [0,2,4.5,7.5,10.5,16.5,20];
% % Icut_max = [2,4.5,7.5,10.5,16.5,20,30];
Icut_min = [0,2,4.5,7.5,10.5,13.5,16.5,20];
Icut_max = [2,4.5,7.5,10.5,13.5,16.5,20,30];
ccode = mycolors(length(Icut_max)+1);
ccode2 = mycolors(size(curve_coef,2)/3);

I_peak0 = curve_coef(:,2:3:end);
P_peak0 = curve_coef(:,1:3:end).*curve_coef(:,3:3:end)*sqrt(pi);
% [xbin,ybin,~] = histw(I_peak,ones(size(I_peak)),0:0.05:5,0.1:0.05:5.1);
[xbin,ybin,~] = histw(I_peak0(:),P_peak0(:),Imin,Imax);
ybin = ybin/max(ybin);
curve1 = fit(xbin(xbin <= Nmax)',ybin(xbin <= Nmax)', gau1,'StartPoint',c0,'Lower',cmin,'Upper',cmax);
% % curve1 = fit(xbin(xbin <= Nmax)',ybin(xbin <= Nmax)', ['gauss',num2str(length(Npeak))],'StartPoint',c0,'Lower',cmin,'Upper',cmax);
coeffvalues(curve1);

P_peak = zeros(size(P_peak0,1),length(Icut_min));
for ii = 1:size(P_peak0,1)
    for jj = 1:size(P_peak0,2)
        kk = find(I_peak0(ii,jj) >= Icut_min & I_peak0(ii,jj) < Icut_max);
        if ~isempty(kk)
            P_peak(ii,kk) = P_peak(ii,kk)+P_peak0(ii,jj);
        end
    end
end
P_peak = [sum(P_peak0,2),P_peak];
P_peak0 = [sum(P_peak0,2),P_peak0];

figure(7)
set(gcf,'Units','inches','Position',[1,1,figsize])
axes('Units','inches','Position',[dsize,asize2]);
    bar(xbin,ybin,'FaceColor',[0.8,0.8,0.8],'EdgeColor','none','DisplayName','Data')
    hold on
    plot(Ic,feval(curve1,Ic),'r','DisplayName','Fit')
xlim(xlim0)
ylim(ylim0)
xlabel('Peak intensity (A.U.)','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
ylabel('Normalized weight','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
legend('show')
set(gca,'FontName',nfont,'FontSize',fsize2,'FontWeight','normal')
box off
view(90,-90)
set(gca,'Units','normalized')

axes('Units','inches','Position',[dsize.*[2,1]+asize2.*[1,0],asize]);
    for ii = 1:size(P_peak,2)
        plot(C_mean/1e-9,P_peak(:,ii)','Color',ccode(ii,:),'DisplayName',['State ',num2str(ii-1)])
%         plot(EL_mean,P_peak(:,ii)','Color',ccode(ii,:),'DisplayName',['State ',num2str(ii-1)])
        hold on
    end
xlim(xlim1)
% xlim(xlim1b)
ylim(ylim1)
xlabel('Bcd concentration (nM)','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
% xlabel('EL','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
ylabel('Probability','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
legend('show')
set(gca,'FontName',nfont,'FontSize',fsize2,'FontWeight','normal')
box off
set(gca,'Units','normalized')

axes('Units','inches','Position',[dsize.*[3,1]+asize2.*[1,0]+asize.*[1,0],asize]);
    for ii = 1:size(I_peak0,2)
        plot(C_mean/1e-9,I_peak0(:,ii)','Color',ccode2(ii,:),'DisplayName',['State ',num2str(ii-1)])
%         plot(EL_mean,I_peak(:,ii)','Color',ccode2(ii,:),'DisplayName',['State ',num2str(ii-1)])
        hold on
    end
xlim(xlim1)
% xlim(xlim1b)
ylim(xlim0)
xlabel('Bcd concentration (nM)','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
% xlabel('EL','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
ylabel('Peak intensity (A.U.)','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
legend('show')
set(gca,'FontName',nfont,'FontSize',fsize2,'FontWeight','normal')
box off
set(gca,'Units','normalized')

set(7,'Units','normalized')

xdata = C_mean'/1e-9;
ydata2 = P_peak(:,1:end-1);
ydata2(:,1) = 1-ydata2(:,1);
ydata = ydata2.*repmat(N_mean',1,size(ydata2,2));

%% Output: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist(result_folder,'dir')
    mkdir(result_folder)
end
saveas(7,[result_folder,out_name,Inten_add,peak_add,fig_tail])
saveas(7,[result_folder,out_name,Inten_add,peak_add,fig_tail2],'epsc')
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Plot the fitting of P_N vs. [TF]: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C_range = [7,11];
C_rangeb = [12,19];
C_range2 = [6,8];
C_range2b = [8,9];
C_range2bb = [9,11];
C_range2bbb = [12,18];
C_range2c = [21,30];
C_range3 = [5,11];
C_range3b = [12,25];
C_range3bb = [29,30];
C_range3c = [12,36];
C0 = 0:0.1:40;
g0 = fittype( @(h,c,a,d, x) a*x.^h./(c^h+x.^h)+d );
beta1 = [5,8,0.1,0]; beta1min = [0,0,0,0]; beta1max = [10,10,0.2,0];
beta2 = [10,20,0.1,0]; beta2min = [0,0,0,0]; beta2max = [15,30,0.2,0];
beta3 = [20,30,0.05,0]; beta3min = [0,0,0,-inf]; beta3max = [30,40,0.2,inf];

xlim0 = [0,0.15];
xlim1 = [0,40];
ylim1 = [-0.02,0.2];
asize = [2,0.9];
asize2 = [0.9,0.9];
dsize = [1,1];
figsize = dsize+(dsize+asize2).*[1,1]+(dsize+asize2).*[1,0]+(dsize+asize).*[1,0];

figure(8)
set(gcf,'Units','inches','Position',[1,1,figsize])
%%% Plot P1 vs. P2
Isel = (C_mean/1e-9 >= C_range(1) & C_mean/1e-9 <= C_range(2)) | (C_mean/1e-9 >= C_rangeb(1) & C_mean/1e-9 <= C_rangeb(2));
p12 = polyfit(P_peak0(Isel,2),P_peak0(Isel,3),1);
axes('Units','inches','Position',[dsize,asize2]);
    plot(P_peak0(:,2),P_peak0(:,3),'.','Color',[0.8,0.8,0.8],'DisplayName','All data')
    hold on
    plot(P_peak0(Isel,2),P_peak0(Isel,3),'k.','DisplayName','Low data')
    plot(ylim1,polyval(p12,ylim1),'r','DisplayName','Fit')
xlim(xlim0)
ylim(xlim0)
xlabel('P_1','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
ylabel('P_2','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
legend('show')
set(gca,'FontName',nfont,'FontSize',fsize2,'FontWeight','normal')
box off
set(gca,'Units','normalized')

%%% Plot P2 vs. P3
Isel = (C_mean/1e-9 >= C_range2(1) & C_mean/1e-9 <= C_range2(2)) | (C_mean/1e-9 >= C_range2b(1) & C_mean/1e-9 <= C_range2b(2)) | (C_mean/1e-9 >= C_range2bb(1) & C_mean/1e-9 <= C_range2bb(2)) | (C_mean/1e-9 >= C_range2bbb(1) & C_mean/1e-9 <= C_range2bbb(2));
p23 = polyfit(P_peak0(Isel,3),P_peak0(Isel,4),1);
Iselb = C_mean/1e-9 >= C_range2c(1) & C_mean/1e-9 <= C_range2c(2);
p23b = polyfit(P_peak0(Iselb,3),P_peak0(Iselb,4),1);
axes('Units','inches','Position',[dsize+(dsize+asize2).*[1,0],asize2]);
    plot(P_peak0(:,3),P_peak0(:,4),'.','Color',[0.8,0.8,0.8],'DisplayName','All data')
    hold on
    plot(P_peak0(Isel,3),P_peak0(Isel,4),'k.','DisplayName','Low data')
    plot(ylim1,polyval(p23,ylim1),'r','DisplayName','Fit')
    plot(P_peak0(Iselb,3),P_peak0(Iselb,4),'.','Color',[0,0,1],'DisplayName','Low data')
    plot(ylim1,polyval(p23b,ylim1),'b','DisplayName','Fit')
xlim(xlim0)
ylim(xlim0)
xlabel('P_2','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
ylabel('P_3','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
legend('show')
set(gca,'FontName',nfont,'FontSize',fsize2,'FontWeight','normal')
box off
set(gca,'Units','normalized')

%%% Plot purified P's and fittings
% p23(1) = 0.25;
P10 = P_peak0(:,2)';
P20 = P_peak0(:,3)';
P30 = P_peak0(:,4)';
P1 = P10/(1-p12(1));
P2 = (P20-P10*p12(1))/(1-p23(1));
P3 = P30-P10*p12(1)*p23(1)-polyval(p23b,P2)*(1-p23(1));

Isel = (C_mean/1e-9 >= C_range3(1) & C_mean/1e-9 <= C_range3(2));
curve1 = fit(C_mean(Isel)'/1e-9,P1(Isel)', g0, 'StartPoint',beta1,'Lower',beta1min,'Upper',beta1max);
Isel = (C_mean/1e-9 >= C_range3b(1) & C_mean/1e-9 <= C_range3b(2)) | (C_mean/1e-9 >= C_range3bb(1) & C_mean/1e-9 <= C_range3bb(2));
curve2 = fit(C_mean(Isel)'/1e-9,P2(Isel)', g0, 'StartPoint',beta2,'Lower',beta2min,'Upper',beta2max);
Isel = (C_mean/1e-9 >= C_range3c(1) & C_mean/1e-9 <= C_range3c(2));
curve3 = fit(C_mean(Isel)'/1e-9,P3(Isel)', g0, 'StartPoint',beta3,'Lower',beta3min,'Upper',beta3max);

axes('Units','inches','Position',[dsize+(dsize+asize2).*[2,0],asize]);
    plot(C_mean/1e-9,P1,'Color',ccode2(2,:),'Marker','.','LineStyle','none','DisplayName','P1 data')
    hold on
    plot(C0,feval(curve1,C0),'Color',ccode2(2,:),'DisplayName','P1 fit')
    plot(C_mean/1e-9,P2,'Color',ccode2(3,:),'Marker','.','LineStyle','none','DisplayName','P2 data')
    plot(C0,feval(curve2,C0),'Color',ccode2(3,:),'DisplayName','P2 fit')
    plot(C_mean/1e-9,P3,'Color',ccode2(4,:),'Marker','.','LineStyle','none','DisplayName','P3 data')
    plot(C0,feval(curve3,C0),'Color',ccode2(4,:),'DisplayName','P3 fit')
xlim(xlim1)
ylim(ylim1)
xlabel('Bcd concentration (nM)','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
ylabel('Probability','FontName',nfont,'FontSize',fsize,'FontWeight','normal')
legend('show')
set(gca,'FontName',nfont,'FontSize',fsize2,'FontWeight','normal')
box off
set(gca,'Units','normalized')

set(8,'Units','normalized')

%% Output: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist(result_folder,'dir')
    mkdir(result_folder)
end
saveas(8,[result_folder,out_name,peak_add,fit_add2,fig_tail])
saveas(8,[result_folder,out_name,peak_add,fit_add2,fig_tail2],'epsc')
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%









