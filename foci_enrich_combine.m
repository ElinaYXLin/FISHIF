
clear all
close all

%% Parameter setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mat_tail = '.mat';
ana_folder = 'enrichspot_result/';
embryo_folder = '22445_2_hb_cds_TMR_gal4-647_bcd_combine_488_par3_new2/';
embryo_name = {'22445_2_hb_cds_TMR_gal4-647_bcd_combine_488_par3_new',...
               '22445_2_hb_cds_TMR_gal4-647_bcd_combine_488_par3_new_add'};
embryo_ind = {1:12,1:5};

out_name = '22445_2_hb_cds_TMR_gal4-647_bcd_combine_488_par3_new2';
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Data initiation: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xyz_all0_cell = cell(0);
dxy0_cell = cell(0);
Imin_all0_cell = cell(0);
I_RNA0_cell = cell(0);
I_en0_cell = cell(0);
N_en0_cell = cell(0);
EL_RNA0_cell = cell(0);
C_RNA0_cell = cell(0);
xyz_all1_cell = cell(0);
dxy1_cell = cell(0);
Imin_all1_cell = cell(0);
I_RNA1_cell = cell(0);
I_en1_cell = cell(0);
I_en2_cell = cell(0);
N_en1_cell = cell(0);
EL_RNA1_cell = cell(0);
C_RNA1_cell = cell(0);
ind_em = zeros(0);
b0_em = zeros(0);
b2_em = zeros(0);
cycle_em = zeros(0);
em_name = cell(0);
foci_list0_cell = cell(0);
foci_list1_cell = cell(0);
foci_circle_area_cell = cell(0);
Nu_en_cell = cell(0);
Nu_RNA_cell = cell(0);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Data loading: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ind0 = 0;
for I_file = 1:length(embryo_name)
    S = load([ana_folder,embryo_folder,embryo_name{I_file},mat_tail]);

    xyz_all0_cell = cat(1,xyz_all0_cell,S.xyz_all0_cell(embryo_ind{I_file},:));
    dxy0_cell = cat(1,dxy0_cell,S.dxy0_cell(embryo_ind{I_file},:));
    Imin_all0_cell = cat(1,Imin_all0_cell,S.Imin_all0_cell(embryo_ind{I_file},:));
    I_RNA0_cell = cat(1,I_RNA0_cell,S.I_RNA0_cell(embryo_ind{I_file},:));
    I_en0_cell = cat(1,I_en0_cell,S.I_en0_cell(embryo_ind{I_file},:));
    N_en0_cell = cat(1,N_en0_cell,S.N_en0_cell(embryo_ind{I_file},:));
    EL_RNA0_cell = cat(1,EL_RNA0_cell,S.EL_RNA0_cell(embryo_ind{I_file},:));
    C_RNA0_cell = cat(1,C_RNA0_cell,S.C_RNA0_cell(embryo_ind{I_file},:));
    xyz_all1_cell = cat(1,xyz_all1_cell,S.xyz_all1_cell(embryo_ind{I_file},:));
    dxy1_cell = cat(1,dxy1_cell,S.dxy1_cell(embryo_ind{I_file},:));
    Imin_all1_cell = cat(1,Imin_all1_cell,S.Imin_all1_cell(embryo_ind{I_file},:));
    I_RNA1_cell = cat(1,I_RNA1_cell,S.I_RNA1_cell(embryo_ind{I_file},:));
    I_en1_cell = cat(1,I_en1_cell,S.I_en1_cell(embryo_ind{I_file},:));
    I_en2_cell = cat(1,I_en2_cell,S.I_en2_cell(embryo_ind{I_file},:));
    N_en1_cell = cat(1,N_en1_cell,S.N_en1_cell(embryo_ind{I_file},:));
    EL_RNA1_cell = cat(1,EL_RNA1_cell,S.EL_RNA1_cell(embryo_ind{I_file},:));
    C_RNA1_cell = cat(1,C_RNA1_cell,S.C_RNA1_cell(embryo_ind{I_file},:));
    ind_em = cat(1,ind_em,[(ind0+1):(ind0+length(embryo_ind{I_file}))]');
    b0_em = cat(1,b0_em,S.b0_em(embryo_ind{I_file},:));
    b2_em = cat(1,b2_em,S.b2_em(embryo_ind{I_file},:));
    cycle_em = cat(1,cycle_em,S.cycle_em(embryo_ind{I_file},:));
    em_name = cat(1,em_name,S.em_name(embryo_ind{I_file},:));
    foci_list0_cell = cat(1,foci_list0_cell,S.foci_list0_cell(embryo_ind{I_file},:));
    foci_list1_cell = cat(1,foci_list1_cell,S.foci_list1_cell(embryo_ind{I_file},:));
    foci_circle_area_cell = cat(1,foci_circle_area_cell,S.foci_circle_area_cell(embryo_ind{I_file},:));
    Nu_en_cell = cat(1,Nu_en_cell,S.Nu_en_cell(embryo_ind{I_file},:));
    Nu_RNA_cell = cat(1,Nu_RNA_cell,S.Nu_RNA_cell(embryo_ind{I_file},:));
    
    ind0 = ind0+length(embryo_ind{I_file});
end
dz = S.dz;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Output: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save([ana_folder,embryo_folder,out_name,mat_tail],'xyz_all0_cell','dxy0_cell','Imin_all0_cell','I_RNA0_cell','I_en0_cell','N_en0_cell','EL_RNA0_cell','C_RNA0_cell','xyz_all1_cell','dxy1_cell','Imin_all1_cell','I_RNA1_cell','I_en1_cell','I_en2_cell','N_en1_cell','EL_RNA1_cell','C_RNA1_cell','ind_em','b0_em','b2_em','cycle_em','em_name','foci_list0_cell','foci_list1_cell','foci_circle_area_cell','Nu_en_cell','Nu_RNA_cell','dz');
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





