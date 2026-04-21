clear all
close all

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A program to summarize the result of dual-FISH labeling experiment %%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Parameter setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
list_name = 'Duallist.xls';
in_folder = 'stacks/';
input_name = 'matchlist.xls';
% out_folder = 'Results_protein2/';
out_folder = 'Results/';
mask_folder = 'masks/';
mask_name = 'mask.mat';
old_add = '_old';

sfit_tail = '.sfit';
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
nu_add = '_nucleus';
cyto_add = '_cytoplasm';
fate_add = '_fate';
RNA_add = '_RNA';
signal2_add = '_RNA2';
% signal2_add = '_gal4';
% signal2_add = '_nullo';
% signal2_add = '_hb-y';
protein_add = '_protein';
fluc_add = '_fluc';
D3_add = '_3D';
DAPI_add = '_DAPI';
noise_add = '_noise';
mismatch_add = '_mismatch';
z_add = '_z';
xy_add= '_xy';
fit_add2 = '_fit';
FISH_add = '_FISH';
single_add = '_single';
marker_add = '_marker';
all_add = '_all';
reg_add = '_regulation';
enrich_add = '_enrich';
binN_add = '_binN';
binC_add = '_binC';
para_add = '_para';
prefit_add = '_prefit';

% % embryo_type = {'hby2'};
% % embryo_name = {'hb-y2'};
% % sub_ij = [5,5];
% % cycle_range = 9:14;
% embryo_type = {'hgal4'};
% embryo_name = {'Hb-hb-gal4-16249-1-5M'};
% sub_ij = [5,5];
% cycle_range = 9:12;
% % embryo_type = {'bgal4C'};
% % embryo_name = {'Bcd-hb-gal4-16249-2-2M'};
% % sub_ij = [3,3];
% % cycle_range = 9:12;
% % % embryo_type = {'bgal4'};
% % % embryo_name = {'Bcd-hb-gal4-16249-1-5M'};
% % % sub_ij = [5,5];
% % % cycle_range = 9:12;
% % embryo_type = {'b'};
% % embryo_name = {'Bcd-hb-nullo-OreR'};
% % sub_ij = [4,5];
% % cycle_range = 10:13;
% % cycle_range2 = 11:13;
% % rsize0 = 2;
% % nM_con = 1;
% % Nr = 20;
% % % embryo_type = {'hby2'};
% % % embryo_name = {'hb-y'};
% % % % % embryo_type = {'hby','b9'};
% % % % % embryo_name = {'hb-y','BAC9'};
% % % embryo_type = {'cbgal4'};
% embryo_type = {'T'};
embryo_type = {'b3n'};
% embryo_name = {'Cad-Bcd-hb-gal4-16249-1-5M'};
embryo_name = {'10652-1-3M_FISHIF_lacZ_gal4_Bcd'};
sub_ij = [4,5];
cycle_range = 10:13;
cycle_range2 = 11:13;
rsize0 = 2;
% rsize0 = 3;
nM_con = 1;
Nr = 20;
% % % % embryo_type = {'h3','h4k5ac','h3k27ac'};
% % % % embryo_name = {'16249-1-5M_hb_gal4_H3','16249-1-5M_hb_gal4_H4K5ac','16249-1-5M_hb_gal4_H3K27ac'};
% % % % sub_ij = [3,3];
% % % % cycle_range = 9:14;
% % % % cycle_range2 = 9:14;
% % % % rsize0 = 1;
% % % % nM_con = 1;
% % % % Nr = 20;

EL_range_fluc0 = [0,0.5];
EL_range_fluc1 = [0,0.5];
% % EL_range_fluc0 = [0,0.5];
% % EL_range_fluc1 = [0,0.3];
% % EL_range_fluc0 = [0.8,1];
% % EL_range_fluc1 = [0.2,0.6];

tc_rescale = false;
r0_cali = 0.2697/3;
t_cali = true;

bin_min = 0:0.05:0.9;   %%% bin_min for EL
bin_max = 0.1:0.05:1;   %%% bin_max for EL
% % bin_min2 = 0:2:100;   %%% bin_min for [Bcd]
% % bin_max2 = 3:2:105;   %%% bin_max for [Bcd]
bin_min2 = 0:2:100;   %%% bin_min for [Bcd]
bin_max2 = 4:2:104;   %%% bin_max for [Bcd]
% % % % bin_min2 = 0:0.2:10;   %%% bin_min for [Bcd]
% % % % bin_max2 = 0.4:0.2:10.4;   %%% bin_max for [Bcd]
EL_range = [0.25,0.75];
% % EL_range2 = [0.15,1];
EL_range2 = [0,1];
z_range = [2,22];
Nbin = 30;
rover = 0.6;
% NbinE = 20;
roverE = 0.25;
NbinE = 20;
% roverE = 0.6;
% NbinE = 30;   %%% Equal population binning for enrichment
NmoveE = 250;
RNAmax = 150;
xlim0 = [0,1]; xlim_bin = 0.01;
xlim1 = [0,50]; xlim_bin1 = 1;
% xlim1 = [0,5]; xlim_bin1 = 1;
xlim2 = [0,50];
% % % xlim2 = [0,5];
ylim1 = [0,3];
ylim2 = [0,150];
ylim2b = [0,200]; 
% % ylim2b = [0,60000]; 
% % % ylim2b = [0,200];
ylim20 = [0,100];
ylim3 = [0,200];
% ylim3 = [0,40];
ylim4 = [0,20];
% ylim4 = [0,5];
% % ylim4 = [0,1.5];
% % % ylim4 = [0,0.5];
% % fitlim = [2,50];
fitlim = [1,22];
% % fitlim = [2,30];
fitlimN = [1,20];
fitlimC = [1,20];
% fitlim2 = [1,40];
% % fitlim2 = [1,50];
% % fitlim2 = [2,50];
fitlim2 = [2,30];
rclim = [0,20];
% % % fitlim = [1,5];
% % % % fitlim = [2,30];
% % % fitlimN = [1,5];
% % % fitlimC = [1,5];
% % % fitlim2 = [1,5];
% % % rclim = [0,5];

hc0 = [-6,0.5];
hp0 = [6,0.3];
% % hc1 = 6;
% % cc1 = 10;
hc1 = 6;
cc1 = 5;

ccode = [1,0,0;0,0,1;0,1,0;];
rc = 0.2;

% ana_folder = 'Dualprocess6_result/01102016_16249-1-5M_OreR/';
% ana_folder = 'Dualprocess6_result/03202016_hb-y2/';
% % % ana_folder = 'Dualprocess6_result/02072016_16249-1-5M_Cad_Bcd/';
% % ana_folder = 'Dualprocess6_result/02072016_16249-1-5M_Cad_Bcd_PBcd_original/';
% % % % ana_folder = 'Dualprocess6_result/02072016_16249-1-5M_Cad_Bcd_PCad_original/';
% ana_folder = 'Dualprocess6_result/02072016_16249-1-5M_Cad_Bcd_PCad_original_deback_01292018/';
ana_folder = 'Dualprocess6_result/11012014_10652-1-3M_FISHIF_lacZ_gal4_Bcd/';
% % % ana_folder = 'Dualprocess6_result/02072016_16249-1-5M_Cad_Bcd_PCad_original_deback_r3/';
% % % % ana_folder = 'Dualprocess6_result/02072016_16249-1-5M_Cad_Bcd_PCad_original_fake/';
% % % % ana_folder = 'Dualprocess6_result/12072016b_3_16249-1-5M_FISHIF_hb48TMR_gal4A647_HistoneA488_60X/';
xls_title = {'Name','Cycle','h','EL0','ymax','ymin'};
xls_title2 = {'Name','Cycle','h','C0','ymax','ymin'};
xls_title3 = {'Name','Cycle','Cmax','EL0','Cmin','Cmax0'};
% xls_title3 = {'Name','Cycle','h','EL0','Cmax','Cmin'};
xls_title0 = {'Name','Cycle','#nuclei','Nu/Cyto','Cmax (nM)'};
xls_title4 = {'Cycle','#nuclei_Ch1','#nuclei_Ch2','Enrich mean (Ch1)','Enrich err (Ch1)','Enrich mean (Ch2)','Enrich err (Ch2)','Fake mean (Ch1)','Fake err (Ch1)','Fake mean (Ch2)','Fake err (Ch2)'};

expf = @(beta,x) beta(1)*exp(-x/beta(2))+beta(3);
Hill1 = @(beta,x) beta(3)*x.^beta(1)./(x.^beta(1)+beta(2)^beta(1));
% Hill1 = @(beta,x) beta(3)*beta(2)^beta(1)./(x.^beta(1)+beta(2)^beta(1))+beta(4);
g0 = fittype( @(h,c,a,d, x) a*x.^h./(c^h+x.^h)+d );
% g0 = fittype( @(h,c,a,d, x) a*c^h./(c^h+x.^h)+d );

ratio0 = 0.0036;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for  iem = 1:length(embryo_type)
    
%% Data loading: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [num_list, folder_list] = xlsread(list_name);
    folder_list = folder_list(strcmpi(embryo_type{iem},folder_list(:,6)),:);
    [N1,N2] = size(folder_list);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Figure initialization: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    close all
    figure(1); maximize(1); set(1,'name',['RNA fate plot (RNA): ',embryo_name{iem}],'PaperPositionMode', 'auto')
    figure(2); maximize(2); set(2,'name',['RNA2 fate plot (RNA2): ',embryo_name{iem}],'PaperPositionMode', 'auto')
    figure(3); maximize(3); set(3,'name',['#foci/nu distribution (RNA): ',embryo_name{iem}],'PaperPositionMode', 'auto')
    figure(4); maximize(4); set(4,'name',['#foci/nu distribution (RNA2): ',embryo_name{iem}],'PaperPositionMode', 'auto')
    figure(5); maximize(5); set(5,'name',['#foci/nu vs. EL (RNA): ',embryo_name{iem}],'PaperPositionMode', 'auto')
    figure(6); maximize(6); set(6,'name',['#foci/nu vs. EL (RNA2): ',embryo_name{iem}],'PaperPositionMode', 'auto')
    figure(7); maximize(7); set(7,'name',['#RNA/nu vs. EL (RNA): ',embryo_name{iem}],'PaperPositionMode', 'auto')
    figure(8); maximize(8); set(8,'name',['#RNA/nu vs. EL (RNA2): ',embryo_name{iem}],'PaperPositionMode', 'auto')
    figure(9); maximize(9); set(9,'name',['[Bcd] fate plot: ',embryo_name{iem}],'PaperPositionMode', 'auto')
    figure(10); maximize(10); set(10,'name',['[Bcd] vs. EL: ',embryo_name{iem}],'PaperPositionMode', 'auto')
    figure(101); maximize(101); set(101,'name',['Bcd mean vs. variance: ',embryo_name{iem}],'PaperPositionMode', 'auto')
    figure(11); maximize(11); set(11,'name',['#foci/nu vs. [Bcd] (RNA): ',embryo_name{iem}],'PaperPositionMode', 'auto')
    figure(12); maximize(12); set(12,'name',['#foci/nu vs. [Bcd] (RNA2): ',embryo_name{iem}],'PaperPositionMode', 'auto')
    figure(13); maximize(13); set(13,'name',['#RNA/nu vs. [Bcd] (RNA): ',embryo_name{iem}],'PaperPositionMode', 'auto')
    figure(14); maximize(14); set(14,'name',['#RNA/nu vs. [Bcd] (RNA2): ',embryo_name{iem}],'PaperPositionMode', 'auto')
    
    figure(15); maximize(15); set(15,'name',['Binding vs. [Bcd] (binN): ',embryo_name{iem}],'PaperPositionMode', 'auto')
    figure(16); maximize(16); set(16,'name',['Binding vs. [Bcd] (binC): ',embryo_name{iem}],'PaperPositionMode', 'auto')
    figure(17); maximize(17); set(17,'name',['#foci/nu vs. [Bcd] pooled: ',embryo_name{iem}],'PaperPositionMode', 'auto')
    figure(18); maximize(18); set(18,'name',['#RNA/nu vs. [Bcd] pooled: ',embryo_name{iem}],'PaperPositionMode', 'auto')
    figure(19); maximize(19); set(19,'name',['#foci/nu vs. EL parameters: ',embryo_name{iem}],'PaperPositionMode', 'auto')
    figure(20); maximize(20); set(20,'name',['#RNA/nu vs. EL parameters: ',embryo_name{iem}],'PaperPositionMode', 'auto')
    figure(21); maximize(21); set(21,'name',['#foci/nu vs. [Bcd] parameters: ',embryo_name{iem}],'PaperPositionMode', 'auto')
    figure(22); maximize(22); set(22,'name',['#RNA/nu vs. [Bcd] parameters: ',embryo_name{iem}],'PaperPositionMode', 'auto')
    figure(23); maximize(23); set(23,'name',['[Bcd] vs. EL parameters: ',embryo_name{iem}],'PaperPositionMode', 'auto')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Image analysis: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    i00 = 0;
    beta_all1 = zeros(0);
    beta_all2 = zeros(0);
    beta_all3 = zeros(0);
    beta_all4 = zeros(0);
    beta_all5 = zeros(0);
    beta_all1b = zeros(0);
    beta_all2b = zeros(0);
    beta_all4b = zeros(0);
    beta_all5b = zeros(0);
    name_all = cell(0);
    cycle_all = zeros(0);
    Nnu_all = zeros(0);
    N_ratio = zeros(0);
    Cm_all = zeros(0);
    
    nucleus_protein_profile_all0 = cell(0);
    nucleus_RNA_profile_all0 = cell(0);
    nucleus_signal2_profile_all0 = cell(0);
    enrich_foci0 = cell(0);
    enrich_foci20 = cell(0);
    RNA_foci0 = cell(0);
    good_embryo = zeros(0);
    
    %%% Input existing fitting parameters (if any): %%%%%%%%%%%%%%%%%%%%%%%
    if exist([ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,foci_add,D3_add,prefit_add,output_tail])
        foci_prefit = xlsread([ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,foci_add,D3_add,prefit_add,output_tail]);
    elseif exist([ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,foci_add,D3_add,sfit_tail])
        v0 = load([ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,foci_add,D3_add,sfit_tail], '-mat');
        foci_prefit = zeros(length(v0.savedSession.AllFitdevsAndConfigs),5);
        for ii = 1:length(v0.savedSession.AllFitdevsAndConfigs)
            temp0 = v0.savedSession.AllFitdevsAndConfigs{ii}.Fitdev.Fit;
            temp00 = v0.savedSession.AllFitdevsAndConfigs{ii}.Fitdev.FitName;
            I0 = find(temp00 == ' ');
            tempI = str2num(temp00(I0(end)+1:end));
            foci_prefit(tempI,:) = [tempI,temp0.h,temp0.c,temp0.a,temp0.d];
        end
    else
        foci_prefit = [];
    end
    
    if exist([ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,foci_add,signal2_add,D3_add,prefit_add,output_tail])
        foci2_prefit = xlsread([ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,foci_add,signal2_add,D3_add,prefit_add,output_tail]);
    elseif exist([ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,foci_add,signal2_add,D3_add,sfit_tail])
        v0 = load([ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,foci_add,signal2_add,D3_add,sfit_tail], '-mat');
        foci2_prefit = zeros(length(v0.savedSession.AllFitdevsAndConfigs),5);
        for ii = 1:length(v0.savedSession.AllFitdevsAndConfigs)
            temp0 = v0.savedSession.AllFitdevsAndConfigs{ii}.Fitdev.Fit;
            temp00 = v0.savedSession.AllFitdevsAndConfigs{ii}.Fitdev.FitName;
            I0 = find(temp00 == ' ');
            tempI = str2num(temp00(I0(end)+1:end));
            foci2_prefit(tempI,:) = [tempI,temp0.h,temp0.c,temp0.a,temp0.d];
        end
    else
        foci2_prefit = [];
    end
    
    if exist([ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,fish_add,int_add,D3_add,prefit_add,output_tail])
        fish_prefit = xlsread([ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,fish_add,int_add,D3_add,prefit_add,output_tail]);
    elseif exist([ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,fish_add,int_add,D3_add,sfit_tail])
        v0 = load([ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,fish_add,int_add,D3_add,sfit_tail], '-mat');
        fish_prefit = zeros(length(v0.savedSession.AllFitdevsAndConfigs),5);
        for ii = 1:length(v0.savedSession.AllFitdevsAndConfigs)
            temp0 = v0.savedSession.AllFitdevsAndConfigs{ii}.Fitdev.Fit;
            temp00 = v0.savedSession.AllFitdevsAndConfigs{ii}.Fitdev.FitName;
            I0 = find(temp00 == ' ');
            tempI = str2num(temp00(I0(end)+1:end));
            fish_prefit(tempI,:) = [tempI,temp0.h,temp0.c,temp0.a,temp0.d];
        end
    else
        fish_prefit = [];
    end
    
    if exist([ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,fish_add,signal2_add,int_add,D3_add,prefit_add,output_tail])
        fish2_prefit = xlsread([ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,fish_add,signal2_add,int_add,D3_add,prefit_add,output_tail]);
    elseif exist([ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,fish_add,signal2_add,int_add,D3_add,sfit_tail])
        v0 = load([ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,fish_add,signal2_add,int_add,D3_add,sfit_tail], '-mat');
        fish2_prefit = zeros(length(v0.savedSession.AllFitdevsAndConfigs),5);
        for ii = 1:length(v0.savedSession.AllFitdevsAndConfigs)
            temp0 = v0.savedSession.AllFitdevsAndConfigs{ii}.Fitdev.Fit;
            temp00 = v0.savedSession.AllFitdevsAndConfigs{ii}.Fitdev.FitName;
            I0 = find(temp00 == ' ');
            tempI = str2num(temp00(I0(end)+1:end));
            fish2_prefit(tempI,:) = [tempI,temp0.h,temp0.c,temp0.a,temp0.d];
        end
    else
        fish2_prefit = [];
    end
    
    if exist([ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,foci_add,reg_add,D3_add,prefit_add,output_tail])
        rfoci_prefit = xlsread([ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,foci_add,reg_add,D3_add,prefit_add,output_tail]);
    elseif exist([ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,foci_add,reg_add,D3_add,sfit_tail])
        v0 = load([ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,foci_add,reg_add,D3_add,sfit_tail], '-mat');
        rfoci_prefit = zeros(length(v0.savedSession.AllFitdevsAndConfigs),5);
        for ii = 1:length(v0.savedSession.AllFitdevsAndConfigs)
            temp0 = v0.savedSession.AllFitdevsAndConfigs{ii}.Fitdev.Fit;
            temp00 = v0.savedSession.AllFitdevsAndConfigs{ii}.Fitdev.FitName;
            I0 = find(temp00 == ' ');
            tempI = str2num(temp00(I0(end)+1:end));
            rfoci_prefit(tempI,:) = [tempI,temp0.h,temp0.c,temp0.a,temp0.d];
        end
    else
        rfoci_prefit = [];
    end
    
    if exist([ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,foci_add,reg_add,signal2_add,D3_add,prefit_add,output_tail])
        rfoci2_prefit = xlsread([ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,foci_add,reg_add,signal2_add,D3_add,prefit_add,output_tail]);
    elseif exist([ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,foci_add,signal2_add,reg_add,D3_add,sfit_tail])
        v0 = load([ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,foci_add,signal2_add,reg_add,D3_add,sfit_tail], '-mat');
        rfoci2_prefit = zeros(length(v0.savedSession.AllFitdevsAndConfigs),5);
        for ii = 1:length(v0.savedSession.AllFitdevsAndConfigs)
            temp0 = v0.savedSession.AllFitdevsAndConfigs{ii}.Fitdev.Fit;
            temp00 = v0.savedSession.AllFitdevsAndConfigs{ii}.Fitdev.FitName;
            I0 = find(temp00 == ' ');
            tempI = str2num(temp00(I0(end)+1:end));
            rfoci2_prefit(tempI,:) = [tempI,temp0.h,temp0.c,temp0.a,temp0.d];
        end
    else
        rfoci2_prefit = [];
    end
    
    if exist([ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,fish_add,reg_add,D3_add,prefit_add,output_tail])
        rfish_prefit = xlsread([ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,fish_add,reg_add,D3_add,prefit_add,output_tail]);
    elseif exist([ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,fish_add,reg_add,D3_add,sfit_tail])
        v0 = load([ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,fish_add,reg_add,D3_add,sfit_tail], '-mat');
        rfish_prefit = zeros(length(v0.savedSession.AllFitdevsAndConfigs),5);
        for ii = 1:length(v0.savedSession.AllFitdevsAndConfigs)
            temp0 = v0.savedSession.AllFitdevsAndConfigs{ii}.Fitdev.Fit;
            temp00 = v0.savedSession.AllFitdevsAndConfigs{ii}.Fitdev.FitName;
            I0 = find(temp00 == ' ');
            tempI = str2num(temp00(I0(end)+1:end));
            rfish_prefit(tempI,:) = [tempI,temp0.h,temp0.c,temp0.a,temp0.d];
        end
    else
        rfish_prefit = [];
    end
    
    if exist([ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,fish_add,signal2_add,reg_add,D3_add,prefit_add,output_tail])
        rfish2_prefit = xlsread([ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,fish_add,signal2_add,reg_add,D3_add,prefit_add,output_tail]);
    elseif exist([ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,fish_add,signal2_add,reg_add,D3_add,sfit_tail])
        v0 = load([ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,fish_add,signal2_add,reg_add,D3_add,sfit_tail], '-mat');
        rfish2_prefit = zeros(length(v0.savedSession.AllFitdevsAndConfigs),5);
        for ii = 1:length(v0.savedSession.AllFitdevsAndConfigs)
            temp0 = v0.savedSession.AllFitdevsAndConfigs{ii}.Fitdev.Fit;
            temp00 = v0.savedSession.AllFitdevsAndConfigs{ii}.Fitdev.FitName;
            I0 = find(temp00 == ' ');
            tempI = str2num(temp00(I0(end)+1:end));
            rfish2_prefit(tempI,:) = [tempI,temp0.h,temp0.c,temp0.a,temp0.d];
        end
    else
        rfish2_prefit = [];
    end

    if exist([ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,protein_add,int_add,D3_add,prefit_add,output_tail])
        protein_prefit = xlsread([ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,protein_add,int_add,D3_add,prefit_add,output_tail]);
    elseif exist([ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,protein_add,int_add,D3_add,sfit_tail])
        v0 = load([ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,protein_add,int_add,D3_add,sfit_tail], '-mat');
        protein_prefit = zeros(length(v0.savedSession.AllFitdevsAndConfigs),4);
        for ii = 1:length(v0.savedSession.AllFitdevsAndConfigs)
            temp0 = v0.savedSession.AllFitdevsAndConfigs{ii}.Fitdev.Fit;
            temp00 = v0.savedSession.AllFitdevsAndConfigs{ii}.Fitdev.FitName;
            I0 = find(temp00 == ' ');
            tempI = str2num(temp00(I0(end)+1:end));
            protein_prefit(tempI,:) = [tempI,temp0.a,temp0.b,temp0.c];
        end
    else
        protein_prefit = [];
    end

    if exist([ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,protein_add,fluc_add,D3_add,prefit_add,output_tail])
        fluc_prefit = xlsread([ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,protein_add,fluc_add,D3_add,prefit_add,output_tail]);
    elseif exist([ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,protein_add,fluc_add,D3_add,sfit_tail])
        v0 = load([ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,protein_add,fluc_add,D3_add,sfit_tail], '-mat');
        fluc_prefit = zeros(length(v0.savedSession.AllFitdevsAndConfigs),3);
        for ii = 1:length(v0.savedSession.AllFitdevsAndConfigs)
            temp0 = v0.savedSession.AllFitdevsAndConfigs{ii}.Fitdev.Fit;
            temp00 = v0.savedSession.AllFitdevsAndConfigs{ii}.Fitdev.FitName;
            I0 = find(temp00 == ' ');
            tempI = str2num(temp00(I0(end)+1:end));
            fluc_prefit(tempI,:) = [tempI,temp0.p1,temp0.p2];
        end
    else
        fluc_prefit = [];
    end
    
    if exist([ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,foci_add,reg_add,all_add,D3_add,prefit_add,output_tail])
        rfoci_pool_prefit = xlsread([ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,foci_add,reg_add,all_add,D3_add,prefit_add,output_tail]);
    elseif exist([ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,foci_add,reg_add,all_add,D3_add,sfit_tail])
        v0 = load([ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,foci_add,reg_add,all_add,D3_add,sfit_tail], '-mat');
        rfoci_pool_prefit = zeros(length(v0.savedSession.AllFitdevsAndConfigs),5);
        for ii = 1:length(v0.savedSession.AllFitdevsAndConfigs)
            temp0 = v0.savedSession.AllFitdevsAndConfigs{ii}.Fitdev.Fit;
            temp00 = v0.savedSession.AllFitdevsAndConfigs{ii}.Fitdev.FitName;
            I0 = find(temp00 == ' ');
            tempI = str2num(temp00(I0(end)+1:end));
            rfoci_pool_prefit(tempI,:) = [tempI,temp0.h,temp0.c,temp0.a,temp0.d];
        end
    else
        rfoci_pool_prefit = [];
    end
    
    if exist([ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,fish_add,reg_add,all_add,D3_add,prefit_add,output_tail])
        rfish_pool_prefit = xlsread([ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,fish_add,reg_add,all_add,D3_add,prefit_add,output_tail]);
    elseif exist([ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,fish_add,reg_add,all_add,D3_add,sfit_tail])
        v0 = load([ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,fish_add,reg_add,all_add,D3_add,sfit_tail], '-mat');
        rfish_pool_prefit = zeros(length(v0.savedSession.AllFitdevsAndConfigs),5);
        for ii = 1:length(v0.savedSession.AllFitdevsAndConfigs)
            temp0 = v0.savedSession.AllFitdevsAndConfigs{ii}.Fitdev.Fit;
            temp00 = v0.savedSession.AllFitdevsAndConfigs{ii}.Fitdev.FitName;
            I0 = find(temp00 == ' ');
            tempI = str2num(temp00(I0(end)+1:end));
            rfish_pool_prefit(tempI,:) = [tempI,temp0.h,temp0.c,temp0.a,temp0.d];
        end
    else
        rfish_pool_prefit = [];
    end
    
    
    %%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for list_I = 1:N1
        [sub_num, sub_list, sub_all] = xlsread([folder_list{list_I,1},in_folder,input_name]);
        channel_name = eval(folder_list{list_I,5});
        if size(sub_all,2) >= 17
            I_se = ~strcmp('F',sub_all(:,end));
    %         I_se = strcmp('T',sub_all(:,end));
            sub_num = sub_num(I_se,:);
            sub_list = sub_list(I_se,:);
        end
        [M1,M2] = size(sub_list);

        for list_J = 1:M1
            i00 = i00+1;
            image_folder = [folder_list{list_I,1},in_folder,sub_list{list_J,3}];
            result_folder = [folder_list{list_I,1},out_folder];
            load([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mat_tail],'nucleus_protein_profile','nucleus_protein_profile_ab','nucleus_RNA_profile','foci_RNA_profile','nucleus_signal2_profile','foci_signal2_profile','foci_data','fake_data','foci_data2','fake_data2','h','h2','r_size','r_size2','em_mask','quanti_p');
            image_folder = [folder_list{list_I,1},in_folder,sub_list{list_J,3}(1:end-1),old_add,'/'];
            
            load([folder_list{list_I,1},mask_folder,sub_list{list_J,3},mask_name]);   %%% load 3D mask
            
            Iz = nucleus_protein_profile(:,4) >= z_range(1) & nucleus_protein_profile(:,4) <= z_range(2);
            Itrue0 = (nucleus_protein_profile(:,1) >= EL_range_fluc0(1)) & (nucleus_protein_profile(:,1) <= EL_range_fluc0(2)) & Iz;
            Itrue = (nucleus_protein_profile(:,1) >= EL_range_fluc1(1)) & (nucleus_protein_profile(:,1) <= EL_range_fluc1(2)) & Iz;
            p0_protein = polyfit(nucleus_protein_profile(Itrue0,2),nucleus_protein_profile(Itrue0,3),1);
            if isempty(fluc_prefit) || size(fluc_prefit,1) < i00 || fluc_prefit(i00,1) <= 0
                p_protein = polyfit(nucleus_protein_profile(Itrue,2),nucleus_protein_profile(Itrue,3),1);
            else
                p_protein = fluc_prefit(i00,2:end);
            end
            
            if t_cali
                r_cali = r0_cali/(p0_protein(1)/p_protein(1));
            else
                r_cali = r0_cali;
            end
            
            if tc_rescale
                rc_cali = 1/(p0_protein(1)/p_protein(1));
            else
                rc_cali = 1;
            end
            
            r_cali = r0_cali*rc_cali;
            nucleus_protein_profile_ab(:,2) = nucleus_protein_profile_ab(:,2)/rc_cali;
            
            nucleus_protein_profile_ab(:,2) = nucleus_protein_profile_ab(:,2)/1e-9;
            foci_RNA_profile = foci_RNA_profile(foci_RNA_profile(:,2) > 0,:);
            N_cycle = sub_num(list_J,13);   %%% Calculate the nuclear cycle number
            name_all = cat(1,name_all,{image_folder});
            cycle_all = [cycle_all;N_cycle];
            Nnu_all = [Nnu_all;size(nucleus_protein_profile,1)];
            N_ratio = [N_ratio;protein_profile3D_short2(mask_stack)];
            Cm_all = [Cm_all;max(nucleus_protein_profile_ab(Iz,2))];

            
%%% Foci number analysis: %%%==============================================
            %%% Dostatni plot for RNA channel
            figure(1)
%             ha1 = subplot(sub_ij(1),sub_ij(2),i00);
            ha1 = subaxis(sub_ij(1),sub_ij(2),i00, 'Spacing', 0, 'PaddingRight',0.03, 'PaddingLeft',0.05, 'PaddingTop',0.03, 'PaddingBottom',0.03, 'MarginLeft', 0.01, 'MarginRight', 0.001, 'MarginTop', 0.015, 'MarginBottom', 0.025);
                h0 = openfig([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),nu_add,fate_add,D3_add,figure_tail]);
                ha0 = get(h0,'Children');
                copyaxes(ha0,ha1,true);
                close(h0)
                axes(ha1);
%                 title([sub_list{list_J,3},': Cycle ',num2str(N_cycle)],'interpreter','none')
                title(['Sample ',num2str(i00),': Cycle ',num2str(N_cycle)])

            %%% Dostatni plot for RNA2 channel
            figure(2)
%             ha1 = subplot(sub_ij(1),sub_ij(2),i00);
            ha1 = subaxis(sub_ij(1),sub_ij(2),i00, 'Spacing', 0, 'PaddingRight',0.03, 'PaddingLeft',0.05, 'PaddingTop',0.03, 'PaddingBottom',0.03, 'MarginLeft', 0.01, 'MarginRight', 0.001, 'MarginTop', 0.015, 'MarginBottom', 0.025);
                h0 = openfig([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),nu_add,signal2_add,fate_add,D3_add,figure_tail]);
                ha0 = get(h0,'Children');
                copyaxes(ha0,ha1,true);
                close(h0)
                axes(ha1);
%                 title([sub_list{list_J,3},': Cycle ',num2str(N_cycle)],'interpreter','none')
                title(['Sample ',num2str(i00),': Cycle ',num2str(N_cycle)])

            %%% #foci/nu distribution for RNA channel
            figure(3)
%             ha1 = subplot(sub_ij(1),sub_ij(2),i00);
            ha1 = subaxis(sub_ij(1),sub_ij(2),i00, 'Spacing', 0, 'PaddingRight',0.03, 'PaddingLeft',0.05, 'PaddingTop',0.03, 'PaddingBottom',0.03, 'MarginLeft', 0.01, 'MarginRight', 0.001, 'MarginTop', 0.015, 'MarginBottom', 0.025);
                h0 = openfig([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),fish_add,num_add,D3_add,figure_tail]);
                ha0 = get(h0,'Children');
                copyaxes(ha0(2),ha1,true);
                close(h0)
                axes(ha1);
%                 title([sub_list{list_J,3},': Cycle ',num2str(N_cycle)],'interpreter','none')
                title(['Sample ',num2str(i00),': Cycle ',num2str(N_cycle)])

            %%% #foci/nu distribution for RNA2 channel
            figure(4)
%             ha1 = subplot(sub_ij(1),sub_ij(2),i00);
            ha1 = subaxis(sub_ij(1),sub_ij(2),i00, 'Spacing', 0, 'PaddingRight',0.03, 'PaddingLeft',0.05, 'PaddingTop',0.03, 'PaddingBottom',0.03, 'MarginLeft', 0.01, 'MarginRight', 0.001, 'MarginTop', 0.015, 'MarginBottom', 0.025);
                h0 = openfig([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),fish_add,signal2_add,num_add,D3_add,figure_tail]);
                ha0 = get(h0,'Children');
                copyaxes(ha0(2),ha1,true);
                close(h0)
                axes(ha1);
%                 title([sub_list{list_J,3},': Cycle ',num2str(N_cycle)],'interpreter','none')
                title(['Sample ',num2str(i00),': Cycle ',num2str(N_cycle)])

            %%% #foci/nu vs. EL for RNA channel
            figure(5)
%             ha1 = subplot(sub_ij(1),sub_ij(2),i00);
            ha1 = subaxis(sub_ij(1),sub_ij(2),i00, 'Spacing', 0, 'PaddingRight',0.03, 'PaddingLeft',0.05, 'PaddingTop',0.03, 'PaddingBottom',0.03, 'MarginLeft', 0.01, 'MarginRight', 0.001, 'MarginTop', 0.015, 'MarginBottom', 0.025);
                [xbin,ybin,xerr,yerr,N_in] = equal_dist(nucleus_RNA_profile(Iz,1),nucleus_RNA_profile(Iz,3),bin_min,bin_max);
                Itrue = (xbin >= EL_range(1)) & (xbin <= EL_range(2)) & ~isnan(ybin);
                xbin0 = xbin(Itrue);
                ybin0 = ybin(Itrue);
                if isempty(foci_prefit) || size(foci_prefit,1) < i00 || foci_prefit(i00,1) <= 0
                    [rmax,I_max] = max(ybin0);
                    beta0 = [hc0,rmax,0];
                    try
                        [beta1,r,~,~,~] = nlinfit(xbin0,ybin0,@Hill,beta0);
                    catch err
                        beta1 = nan(size(beta0));
                    end
                else
                    beta1 = foci_prefit(i00,2:end);
                end
                beta_all1 = [beta_all1;beta1];
                errorbarxy(xbin,ybin,xerr,yerr,xerr,yerr,'k','k','.','DisplayName','Binned data');
                hold on
                plot([xlim0(1):xlim_bin:xlim0(2)],Hill(beta1,[xlim0(1):xlim_bin:xlim0(2)]),'r-','DisplayName','Hill fit');
                text(mean(xlim0),mean(ylim1),['h = ',num2str(beta1(1),'%5.2f'),char(10),'EL0 = ',num2str(beta1(2),'%5.2f'),char(10),'ymax = ',num2str(beta1(3),'%5.2f')])
                xlim(xlim0)
                ylim(ylim1)
                xlabel('EL')
                ylabel('#foci/nu')
%                 title([sub_list{list_J,3},': Cycle ',num2str(N_cycle)],'interpreter','none')
                title(['Sample ',num2str(i00),': Cycle ',num2str(N_cycle)])
                hold off

            %%% #foci/nu vs. EL for RNA2 channel
            figure(6)
%             ha1 = subplot(sub_ij(1),sub_ij(2),i00);
            ha1 = subaxis(sub_ij(1),sub_ij(2),i00, 'Spacing', 0, 'PaddingRight',0.03, 'PaddingLeft',0.05, 'PaddingTop',0.03, 'PaddingBottom',0.03, 'MarginLeft', 0.01, 'MarginRight', 0.001, 'MarginTop', 0.015, 'MarginBottom', 0.025);
                [xbin,ybin,xerr,yerr,N_in] = equal_dist(nucleus_signal2_profile(Iz,1),nucleus_signal2_profile(Iz,3),bin_min,bin_max);
                Itrue = (xbin >= EL_range(1)) & (xbin <= EL_range(2)) & ~isnan(ybin);
                xbin0 = xbin(Itrue);
                ybin0 = ybin(Itrue);
                if isempty(foci2_prefit) || size(foci2_prefit,1) < i00 || foci2_prefit(i00,1) <= 0
                    [rmax,I_max] = max(ybin0);
                    beta0 = [hc0,rmax,0];
                    try
                        [beta1,r,~,~,~] = nlinfit(xbin0,ybin0,@Hill,beta0);
                    catch err
                        beta1 = nan(size(beta0));
                    end
                else
                    beta1 = foci2_prefit(i00,2:end);
                end
                beta_all1b = [beta_all1b;beta1];
                errorbarxy(xbin,ybin,xerr,yerr,xerr,yerr,'k','k','.','DisplayName','Binned data');
                hold on
                plot([xlim0(1):xlim_bin:xlim0(2)],Hill(beta1,[xlim0(1):xlim_bin:xlim0(2)]),'r-','DisplayName','Hill fit');
                text(mean(xlim0),mean(ylim1),['h = ',num2str(beta1(1),'%5.2f'),char(10),'EL0 = ',num2str(beta1(2),'%5.2f'),char(10),'ymax = ',num2str(beta1(3),'%5.2f')])
                xlim(xlim0)
                ylim(ylim1)
                xlabel('EL')
                ylabel('#foci/nu')
%                 title([sub_list{list_J,3},': Cycle ',num2str(N_cycle)],'interpreter','none')
                title(['Sample ',num2str(i00),': Cycle ',num2str(N_cycle)])

            %%% #RNA/nu vs. EL for RNA channel
            figure(7)
%             ha1 = subplot(sub_ij(1),sub_ij(2),i00);
            ha1 = subaxis(sub_ij(1),sub_ij(2),i00, 'Spacing', 0, 'PaddingRight',0.03, 'PaddingLeft',0.05, 'PaddingTop',0.03, 'PaddingBottom',0.03, 'MarginLeft', 0.01, 'MarginRight', 0.001, 'MarginTop', 0.015, 'MarginBottom', 0.025);
                [xbin,ybin,xerr,yerr,N_in] = equal_dist(nucleus_RNA_profile(Iz,1),nucleus_RNA_profile(Iz,4),bin_min,bin_max);
                Itrue = (xbin >= EL_range(1)) & (xbin <= EL_range(2)) & ~isnan(ybin);
                xbin0 = xbin(Itrue);
                ybin0 = ybin(Itrue);
                if isempty(fish_prefit) || size(fish_prefit,1) < i00 || fish_prefit(i00,1) <= 0
                    [rmax,I_max] = max(ybin0);
                    beta0 = [hc0,rmax,0];
                    try
                        [beta2,r,~,~,~] = nlinfit(xbin0,ybin0,@Hill,beta0);
                    catch err
                        beta2 = nan(size(beta0));
                    end
                else
                    beta2 = fish_prefit(i00,2:end);
                end
                beta_all2 = [beta_all2;beta2];
                plot(nucleus_RNA_profile(Iz,1),nucleus_RNA_profile(Iz,4),'.','Color',0.7*ones(1,3),'DisplayName','Data')
                hold on
                errorbarxy(xbin,ybin,xerr,yerr,xerr,yerr,'k','k','.','DisplayName','Binned data');
                plot([xlim0(1):xlim_bin:xlim0(2)],Hill(beta2,[xlim0(1):xlim_bin:xlim0(2)]),'r-','DisplayName','Hill fit');
                text(mean(xlim0),mean(ylim2),['h = ',num2str(beta2(1),'%5.2f'),char(10),'EL0 = ',num2str(beta2(2),'%5.2f'),char(10),'ymax = ',num2str(beta2(3),'%5.2f')])
                xlim(xlim0)
                ylim(ylim2)
                xlabel('EL')
                ylabel('#RNA/nu')
%                 title([sub_list{list_J,3},': Cycle ',num2str(N_cycle)],'interpreter','none')
                title(['Sample ',num2str(i00),': Cycle ',num2str(N_cycle)])
                hold off

            %%% #RNA/nu vs. EL for RNA2 channel
            figure(8)
%             ha1 = subplot(sub_ij(1),sub_ij(2),i00);
            ha1 = subaxis(sub_ij(1),sub_ij(2),i00, 'Spacing', 0, 'PaddingRight',0.03, 'PaddingLeft',0.05, 'PaddingTop',0.03, 'PaddingBottom',0.03, 'MarginLeft', 0.01, 'MarginRight', 0.001, 'MarginTop', 0.015, 'MarginBottom', 0.025);
                [xbin,ybin,xerr,yerr,N_in] = equal_dist(nucleus_signal2_profile(Iz,1),nucleus_signal2_profile(Iz,4),bin_min,bin_max);
                Itrue = (xbin >= EL_range(1)) & (xbin <= EL_range(2)) & ~isnan(ybin);
                xbin0 = xbin(Itrue);
                ybin0 = ybin(Itrue);
                if isempty(fish2_prefit) || size(fish2_prefit,1) < i00 || fish2_prefit(i00,1) <= 0
                    [rmax,I_max] = max(ybin0);
                    beta0 = [hc0,rmax,0];
                    try
                        [beta2,r,~,~,~] = nlinfit(xbin0,ybin0,@Hill,beta0);
                    catch err
                        beta2 = nan(size(beta0));
                    end
                else
                    beta2 = fish2_prefit(i00,2:end);
                end
                beta_all2b = [beta_all2b;beta2];
                plot(nucleus_signal2_profile(Iz,1),nucleus_signal2_profile(Iz,4),'.','Color',0.7*ones(1,3),'DisplayName','Data')
                hold on
                errorbarxy(xbin,ybin,xerr,yerr,xerr,yerr,'k','k','.','DisplayName','Binned data');
                plot([xlim0(1):xlim_bin:xlim0(2)],Hill(beta2,[xlim0(1):xlim_bin:xlim0(2)]),'r-','DisplayName','Hill fit');
                text(mean(xlim0),mean(ylim2b),['h = ',num2str(beta2(1),'%5.2f'),char(10),'EL0 = ',num2str(beta2(2),'%5.2f'),char(10),'ymax = ',num2str(beta2(3),'%5.2f')])
                xlim(xlim0)
                ylim(ylim2b)
                xlabel('EL')
                ylabel('#RNA/nu')
%                 title([sub_list{list_J,3},': Cycle ',num2str(N_cycle)],'interpreter','none')
                title(['Sample ',num2str(i00),': Cycle ',num2str(N_cycle)])
                hold off

                
%%% Bcd concentration analysis: %%%========================================
            %%% Bcd fate plot
            figure(9)
%             ha1 = subplot(sub_ij(1),sub_ij(2),i00);
            ha1 = subaxis(sub_ij(1),sub_ij(2),i00, 'Spacing', 0, 'PaddingRight',0.03, 'PaddingLeft',0.05, 'PaddingTop',0.03, 'PaddingBottom',0.03, 'MarginLeft', 0.01, 'MarginRight', 0.001, 'MarginTop', 0.015, 'MarginBottom', 0.025);
                h0 = openfig([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,fate_add,D3_add,figure_tail]);
                ha0 = get(h0,'Children');
                copyaxes(ha0(end),ha1,true);
                ha1xy = get(ha1,'Position');
                ha2 = axes('Position',[ha1xy(1)+ha1xy(3)/Nr,ha1xy(2),ha1xy(3)/Nr,ha1xy(4)]);
                copyaxes(ha0(1),ha2,true);
                close(h0)
                axes(ha1);
%                 title([sub_list{list_J,3},': Cycle ',num2str(N_cycle)],'interpreter','none')
                title(['Sample ',num2str(i00),': Cycle ',num2str(N_cycle)])

            %%% [Bcd] vs. EL
            figure(10)
%             ha1 = subplot(sub_ij(1),sub_ij(2),i00);
            ha1 = subaxis(sub_ij(1),sub_ij(2),i00, 'Spacing', 0, 'PaddingRight',0.03, 'PaddingLeft',0.05, 'PaddingTop',0.03, 'PaddingBottom',0.03, 'MarginLeft', 0.01, 'MarginRight', 0.001, 'MarginTop', 0.015, 'MarginBottom', 0.025);
                [xbin,ybin,xerr,yerr,N_in] = equal_dist(nucleus_protein_profile_ab(Iz,1),nucleus_protein_profile_ab(Iz,2),bin_min,bin_max);
                Itrue = (xbin >= EL_range2(1)) & (xbin <= EL_range2(2)) & ~isnan(ybin);
                xbin0 = xbin(Itrue);
                ybin0 = ybin(Itrue);
                if isempty(protein_prefit) || size(protein_prefit,1) < i00 || protein_prefit(i00,1) <= 0
                    beta0 = [max(ybin0),0.2,0];
                    try
                        [beta1,r,~,~,~] = nlinfit(xbin0,ybin0,expf,beta0);
                        H_I = beta1;
                    catch err
                        H_I = nan(size(beta0));
                    end
% % %                     [rmax,I_max] = max(ybin0);
% % %                     beta0 = [hp0,rmax,0];
% % %                     try
% % %                         [beta1,r,~,~,~] = nlinfit(xbin0,ybin0,@Hill,beta0);
% % %                     catch err
% % %                         beta1 = nan(size(beta0));
% % %                     end
                else
                    beta1 = protein_prefit(i00,2:end);
                end
                beta_all3 = [beta_all3;[beta1,max(ybin0)]];
%                 beta_all3 = [beta_all3;beta1];
                plot(nucleus_protein_profile_ab(Iz,1),nucleus_protein_profile_ab(Iz,2),'.','Color',0.7*ones(1,3),'DisplayName','Data')
                hold on
                errorbarxy(xbin,ybin,xerr,yerr,xerr,yerr,'k','k','.','DisplayName','Binned data');
                plot([xlim0(1):xlim_bin:xlim0(2)],expf(beta1,[xlim0(1):xlim_bin:xlim0(2)]),'r-','DisplayName','Exp fit');
% %                 text(mean(xlim0),mean(ylim3),['EL0 = ',num2str(beta1(2),'%5.2f'),char(10),'Cmax = ',num2str(beta1(1),'%6.2f')])
% %                 plot([xlim0(1):xlim_bin:xlim0(2)],Hill(beta1,[xlim0(1):xlim_bin:xlim0(2)]),'r-','DisplayName','Exp fit');
                text(mean(xlim0),mean(ylim2b),['h = ',num2str(beta1(1),'%5.2f'),char(10),'EL0 = ',num2str(beta1(2),'%5.2f'),char(10),'ymax = ',num2str(beta1(3),'%5.2f')])
                xlim(xlim0)
                ylim(ylim3)
                xlabel('EL')
                ylabel('[Protein] (nM)')
%                 title([sub_list{list_J,3},': Cycle ',num2str(N_cycle)],'interpreter','none')
                title(['Sample ',num2str(i00),': Cycle ',num2str(N_cycle)])
                hold off

            %%% Bcd intensity variance vs. mean
            figure(101)
%             ha1 = subplot(sub_ij(1),sub_ij(2),i00);
            ha1 = subaxis(sub_ij(1),sub_ij(2),i00, 'Spacing', 0, 'PaddingRight',0.03, 'PaddingLeft',0.05, 'PaddingTop',0.03, 'PaddingBottom',0.03, 'MarginLeft', 0.01, 'MarginRight', 0.001, 'MarginTop', 0.015, 'MarginBottom', 0.025);
                Itrue0 = (nucleus_protein_profile(:,1) >= EL_range_fluc0(1)) & (nucleus_protein_profile(:,1) <= EL_range_fluc0(2)) & Iz;
                Itrue = (nucleus_protein_profile(:,1) >= EL_range_fluc1(1)) & (nucleus_protein_profile(:,1) <= EL_range_fluc1(2)) & Iz;
% %                 p0 = polyfit(nucleus_protein_profile(Itrue0,2),nucleus_protein_profile(Itrue0,3),1);
% %                 p = polyfit(nucleus_protein_profile(Itrue,2),nucleus_protein_profile(Itrue,3),1);
                
                plot(nucleus_protein_profile(:,2),nucleus_protein_profile(:,3),'.','Color',0.7*ones(1,3),'DisplayName','Data')
                hold on
                plot(nucleus_protein_profile(Itrue,2),nucleus_protein_profile(Itrue,3),'.','Color',zeros(1,3),'DisplayName','Used')
                xlimt = xlim;
                ylimt = ylim;
                plot(xlimt,polyval(p_protein,xlimt),'r-','DisplayName','Fit');
                
% %                 text(mean(xlimt),mean(ylimt),['C_m = ',num2str(max(nucleus_protein_profile_ab(Iz,2))/p_protein(1)*p0_protein(1),'%5.2f'),char(10),'C_m_0 = ',num2str(max(nucleus_protein_profile_ab(Iz,2)),'%5.2f')])
                text(mean(xlimt),mean(ylimt),['C_m = ',num2str(max(nucleus_protein_profile_ab(Iz,2)),'%5.2f')])
                xlim(xlimt)
                ylim(ylimt)
                xlabel('EL')
                ylabel('[Protein] (nM)')
%                 title([sub_list{list_J,3},': Cycle ',num2str(N_cycle)],'interpreter','none')
                title(['Sample ',num2str(i00),': Cycle ',num2str(N_cycle)])
                hold off

                
%%% TX regulation analysis: %%%============================================
            %%% #foci/nu vs. [Bcd] for RNA channel
            figure(11)
%             ha1 = subplot(sub_ij(1),sub_ij(2),i00);
            ha1 = subaxis(sub_ij(1),sub_ij(2),i00, 'Spacing', 0, 'PaddingRight',0.03, 'PaddingLeft',0.05, 'PaddingTop',0.03, 'PaddingBottom',0.03, 'MarginLeft', 0.01, 'MarginRight', 0.001, 'MarginTop', 0.015, 'MarginBottom', 0.025);
                [xbin,ybin,xerr,yerr,N_in] = equal_bin(nucleus_protein_profile_ab(Iz,2),nucleus_RNA_profile(Iz,3),Nbin,rover,1);
                I00 = ~isnan(xbin) & ~isnan(ybin) & N_in >= 5;
                xbin = xbin(I00);
                ybin = ybin(I00);
                xerr = xerr(I00);
                yerr = yerr(I00);
                N_in = N_in(I00);
                
                I_fit = xbin >= fitlim(1) & xbin <= fitlim(2);
                [rmax,I_max] = max(ybin(~isnan(ybin)));
                [~,I_middle] = min(abs(ybin(1:I_max)-rmax/2));
%                     [~,I_middle] = min(abs(ybin(I_max:end)-rmax/2));
                beta0 = [hc1,xbin(I_middle),rmax,0];
                if isempty(rfoci_prefit) || size(rfoci_prefit,1) < i00 || rfoci_prefit(i00,1) <= 0
                    try
% %                         [beta1,r,~,~,~] = nlinfit(xbin(I_fit),ybin(I_fit),@Hill,beta0);
%                         [beta1,r,~,~,~] = nlinfit(xbin(I_fit),ybin(I_fit),Hill1,beta0);
                        curve = fit(xbin(I_fit)',ybin(I_fit)', g0, 'StartPoint',beta0,'Lower',[0,0,0,0] );
                        beta1 = [curve.h,curve.c,curve.a,curve.d];
                    catch err
                        beta1 = nan(size(beta0));
                        curve.h = nan; curve.c = nan; curve.a = nan; curve.d = nan;
                    end
                else
                    beta1 = rfoci_prefit(i00,2:end);
                    curve = fit(xbin(I_fit)',ybin(I_fit)', g0, 'StartPoint',beta0,'Lower',[0,0,0,0] );
                    curve.h = beta1(1); curve.c = beta1(2); curve.a = beta1(3); curve.d = beta1(4);
                end
                beta_all4 = [beta_all4;beta1];
                errorbarxy(xbin,ybin,xerr,yerr,xerr,yerr,'k','k','.','DisplayName','Binned data');
                hold on
% %                 plot([xlim1(1):xlim_bin1:xlim1(2)],Hill(beta1,[xlim1(1):xlim_bin1:xlim1(2)]),'r-','DisplayName','Hill fit');
%                 plot([xlim1(1):xlim_bin1:xlim1(2)],Hill1(beta1,[xlim1(1):xlim_bin1:xlim1(2)]),'r-','DisplayName','Hill fit');
                plot([xlim1(1):xlim_bin1:xlim1(2)],feval(curve,[xlim1(1):xlim_bin1:xlim1(2)]),'r-','DisplayName','Hill fit');
                text(mean(xlim1),mean(ylim1),['h = ',num2str(beta1(1),'%5.2f'),char(10),'C0 = ',num2str(beta1(2),'%5.2f'),' nM',char(10),'ymax = ',num2str(beta1(3),'%5.2f')])
                xlim(xlim1)
                ylim(ylim1)
                xlabel('[Protein] (nM)')
                ylabel('#foci/nu')
%                 title([sub_list{list_J,3},': Cycle ',num2str(N_cycle)],'interpreter','none')
                title(['Sample ',num2str(i00),': Cycle ',num2str(N_cycle)])
                hold off

            %%% #foci/nu vs. [Bcd] for RNA2 channel
            figure(12)
%             ha1 = subplot(sub_ij(1),sub_ij(2),i00);
            ha1 = subaxis(sub_ij(1),sub_ij(2),i00, 'Spacing', 0, 'PaddingRight',0.03, 'PaddingLeft',0.05, 'PaddingTop',0.03, 'PaddingBottom',0.03, 'MarginLeft', 0.01, 'MarginRight', 0.001, 'MarginTop', 0.015, 'MarginBottom', 0.025);
                [xbin,ybin,xerr,yerr,N_in] = equal_bin(nucleus_protein_profile_ab(Iz,2),nucleus_signal2_profile(Iz,3),Nbin,rover,1);
                I00 = ~isnan(xbin) & ~isnan(ybin) & N_in >= 5;
                xbin = xbin(I00);
                ybin = ybin(I00);
                xerr = xerr(I00);
                yerr = yerr(I00);
                N_in = N_in(I00);
                
                I_fit = xbin >= fitlim(1) & xbin <= fitlim(2);
                [rmax,I_max] = max(ybin);
                [~,I_middle] = min(abs(ybin(1:I_max)-rmax/2));
%                     [~,I_middle] = min(abs(ybin(I_max:end)-rmax/2));
                beta0 = [hc1,xbin(I_middle),rmax,0];
                if isempty(rfoci2_prefit) || size(rfoci2_prefit,1) < i00 || rfoci2_prefit(i00,1) <= 0
                    try
% %                         [beta1,r,~,~,~] = nlinfit(xbin(I_fit),ybin(I_fit),@Hill,beta0);
%                         [beta1,r,~,~,~] = nlinfit(xbin(I_fit),ybin(I_fit),Hill1,beta0);
                        curve = fit(xbin(I_fit)',ybin(I_fit)', g0, 'StartPoint',beta0,'Lower',[0,0,0,0] );
                        beta1 = [curve.h,curve.c,curve.a,curve.d];
                    catch err
                        beta1 = nan(size(beta0));
                        curve.h = nan; curve.c = nan; curve.a = nan; curve.d = nan;
                    end
                else
                    beta1 = rfoci2_prefit(i00,2:end);
                    curve = fit(xbin(I_fit)',ybin(I_fit)', g0, 'StartPoint',beta0,'Lower',[0,0,0,0] );
                    curve.h = beta1(1); curve.c = beta1(2); curve.a = beta1(3); curve.d = beta1(4);
                end
                beta_all4b = [beta_all4b;beta1];
                errorbarxy(xbin,ybin,xerr,yerr,xerr,yerr,'k','k','.','DisplayName','Binned data');
                hold on
% %                 plot([xlim1(1):xlim_bin1:xlim1(2)],Hill(beta1,[xlim1(1):xlim_bin1:xlim1(2)]),'r-','DisplayName','Hill fit');
%                 plot([xlim1(1):xlim_bin1:xlim1(2)],Hill1(beta1,[xlim1(1):xlim_bin1:xlim1(2)]),'r-','DisplayName','Hill fit');
                plot([xlim1(1):xlim_bin1:xlim1(2)],feval(curve,[xlim1(1):xlim_bin1:xlim1(2)]),'r-','DisplayName','Hill fit');
               text(mean(xlim1),mean(ylim1),['h = ',num2str(beta1(1),'%5.2f'),char(10),'C0 = ',num2str(beta1(2),'%5.2f'),' nM',char(10),'ymax = ',num2str(beta1(3),'%5.2f')])
                xlim(xlim1)
                ylim(ylim1)
                xlabel('[Protein] (nM)')
                ylabel('#foci/nu')
%                 title([sub_list{list_J,3},': Cycle ',num2str(N_cycle)],'interpreter','none')
                title(['Sample ',num2str(i00),': Cycle ',num2str(N_cycle)])
                hold off

            %%% #RNA/nu vs. [Bcd] for RNA channel
            figure(13)
%             ha1 = subplot(sub_ij(1),sub_ij(2),i00);
            ha1 = subaxis(sub_ij(1),sub_ij(2),i00, 'Spacing', 0, 'PaddingRight',0.03, 'PaddingLeft',0.05, 'PaddingTop',0.03, 'PaddingBottom',0.03, 'MarginLeft', 0.01, 'MarginRight', 0.001, 'MarginTop', 0.015, 'MarginBottom', 0.025);
                [xbin,ybin,xerr,yerr,N_in] = equal_bin(nucleus_protein_profile_ab(Iz,2),nucleus_RNA_profile(Iz,4),Nbin,rover,1);
                I00 = ~isnan(xbin) & ~isnan(ybin) & N_in >= 5;
                xbin = xbin(I00);
                ybin = ybin(I00);
                xerr = xerr(I00);
                yerr = yerr(I00);
                N_in = N_in(I00);
                I_fit = xbin >= fitlim(1) & xbin <= fitlim(2);
                if isempty(rfish_prefit) || size(rfish_prefit,1) < i00 || rfish_prefit(i00,1) <= 0
                    [rmax,I_max] = max(ybin);
                    [~,I_middle] = min(abs(ybin(1:I_max)-rmax/2));
%                     [~,I_middle] = min(abs(ybin(I_max:end)-rmax/2));
                    beta0 = [hc1,xbin(I_middle),rmax,0];
                    try
% %                         [beta2,r,~,~,~] = nlinfit(xbin(I_fit),ybin(I_fit),@Hill,beta0);
%                         [beta2,r,~,~,~] = nlinfit(xbin(I_fit),ybin(I_fit),Hill1,beta0);
                        curve = fit(xbin(I_fit)',ybin(I_fit)', g0, 'StartPoint',beta0,'Lower',[0,0,0,0] );
                        beta2 = [curve.h,curve.c,curve.a,curve.d];
                    catch err
                        beta2 = nan(size(beta0));
                        curve.h = nan; curve.c = nan; curve.a = nan; curve.d = nan;
                    end
                else
                    beta2 = rfish_prefit(i00,2:end);
                    try
                        curve = fit(xbin(I_fit)',ybin(I_fit)', g0, 'StartPoint',beta0,'Lower',[0,0,0,0] );
                    catch err
                    end
                    curve.h = beta2(1); curve.c = beta2(2); curve.a = beta2(3); curve.d = beta2(4); 
                end
                beta_all5 = [beta_all5;beta2];
                plot(nucleus_protein_profile_ab(Iz,2),nucleus_RNA_profile(Iz,4),'.','Color',0.7*ones(1,3),'DisplayName','Data')
                hold on
                errorbarxy(xbin,ybin,xerr,yerr,xerr,yerr,'k','k','.','DisplayName','Binned data');
% %                 plot([xlim1(1):xlim_bin1:xlim1(2)],Hill(beta2,[xlim1(1):xlim_bin1:xlim1(2)]),'r-','DisplayName','Hill fit');
%                 plot([xlim1(1):xlim_bin1:xlim1(2)],Hill1(beta2,[xlim1(1):xlim_bin1:xlim1(2)]),'r-','DisplayName','Hill fit');
                plot([xlim1(1):xlim_bin1:xlim1(2)],feval(curve,[xlim1(1):xlim_bin1:xlim1(2)]),'r-','DisplayName','Hill fit');
                text(mean(xlim1),mean(ylim2),['h = ',num2str(beta2(1),'%5.2f'),char(10),'C0 = ',num2str(beta2(2),'%5.2f'),' nM',char(10),'ymax = ',num2str(beta2(3),'%5.2f')])
                xlim(xlim1)
                ylim(ylim2)
                xlabel('[Protein] (nM)')
                ylabel('#RNA/nu')
%                 title([sub_list{list_J,3},': Cycle ',num2str(N_cycle)],'interpreter','none')
                title(['Sample ',num2str(i00),': Cycle ',num2str(N_cycle)])
                hold off

            %%% #RNA/nu vs. [Bcd] for RNA2 channel
            figure(14)
%             ha1 = subplot(sub_ij(1),sub_ij(2),i00);
            ha1 = subaxis(sub_ij(1),sub_ij(2),i00, 'Spacing', 0, 'PaddingRight',0.03, 'PaddingLeft',0.05, 'PaddingTop',0.03, 'PaddingBottom',0.03, 'MarginLeft', 0.01, 'MarginRight', 0.001, 'MarginTop', 0.015, 'MarginBottom', 0.025);
                [xbin,ybin,xerr,yerr,N_in] = equal_bin(nucleus_protein_profile_ab(Iz,2),nucleus_signal2_profile(Iz,4),Nbin,rover,1);
                I00 = ~isnan(xbin) & ~isnan(ybin) & N_in >= 5;
                xbin = xbin(I00);
                ybin = ybin(I00);
                xerr = xerr(I00);
                yerr = yerr(I00);
                N_in = N_in(I00);
                I_fit = xbin >= fitlim(1) & xbin <= fitlim(2);
                if isempty(rfish2_prefit) || size(rfish2_prefit,1) < i00 || rfish2_prefit(i00,1) <= 0
                    [rmax,I_max] = max(ybin);
                    [~,I_middle] = min(abs(ybin(1:I_max)-rmax/2));
%                     [~,I_middle] = min(abs(ybin(I_max:end)-rmax/2));
                    beta0 = [hc1,xbin(I_middle),rmax,0];
                    try
% %                         [beta2,r,~,~,~] = nlinfit(xbin(I_fit),ybin(I_fit),@Hill,beta0);
%                         [beta2,r,~,~,~] = nlinfit(xbin(I_fit),ybin(I_fit),Hill1,beta0);
                        curve = fit(xbin(I_fit)',ybin(I_fit)', g0, 'StartPoint',beta0,'Lower',[0,0,0,0] );
                        beta2 = [curve.h,curve.c,curve.a,curve.d];
                    catch err
                        beta2 = nan(size(beta0));
                        curve.h = nan; curve.c = nan; curve.a = nan; curve.d = nan;
                    end
                else
                    beta2 = rfish2_prefit(i00,2:end);
                    try
                        curve = fit(xbin(I_fit)',ybin(I_fit)', g0, 'StartPoint',beta0,'Lower',[0,0,0,0] );
                    catch err
                    end
                    curve.h = beta2(1); curve.c = beta2(2); curve.a = beta2(3); curve.d = beta2(4); 
                end
                beta_all5b = [beta_all5b;beta2];
                plot(nucleus_protein_profile_ab(Iz,2),nucleus_signal2_profile(Iz,4),'.','Color',0.7*ones(1,3),'DisplayName','Data')
                hold on
                errorbarxy(xbin,ybin,xerr,yerr,xerr,yerr,'k','k','.','DisplayName','Binned data');
% %                 plot([xlim1(1):xlim_bin1:xlim1(2)],Hill(beta2,[xlim1(1):xlim_bin1:xlim1(2)]),'r-','DisplayName','Hill fit');
%                 plot([xlim1(1):xlim_bin1:xlim1(2)],Hill1(beta2,[xlim1(1):xlim_bin1:xlim1(2)]),'r-','DisplayName','Hill fit');
                plot([xlim1(1):xlim_bin1:xlim1(2)],feval(curve,[xlim1(1):xlim_bin1:xlim1(2)]),'r-','DisplayName','Hill fit');
                text(mean(xlim1),mean(ylim2b),['h = ',num2str(beta2(1),'%5.2f'),char(10),'C0 = ',num2str(beta2(2),'%5.2f'),' nM',char(10),'ymax = ',num2str(beta2(3),'%5.2f')])
                xlim(xlim1)
                ylim(ylim2b)
                xlabel('[Protein] (nM)')
                ylabel('#RNA/nu')
%                 title([sub_list{list_J,3},': Cycle ',num2str(N_cycle)],'interpreter','none')
                title(['Sample ',num2str(i00),': Cycle ',num2str(N_cycle)])
                hold off
                
                
%%% Bcd enrichment data pooling: %%%=======================================
%             if any(N_cycle == cycle_range)
                foci_raw = foci_data{r_size == rsize0};
                fake_raw = fake_data{r_size == rsize0};
                foci_raw2 = foci_data2{r_size2 == rsize0};
                fake_raw2 = fake_data2{r_size2 == rsize0};

                h_area = sum(h{r_size == rsize0}(:));
                foci_raw(:,3) = nM_con*h_area*foci_raw(:,9).*foci_raw(:,13);
                fake_raw(:,3) = nM_con*h_area*fake_raw(:,9).*fake_raw(:,13);
                foci_raw(:,4) = foci_raw(:,4)-foci_raw(:,5)/quanti_p(1)*ratio0;
                fake_raw(:,4) = fake_raw(:,4)-fake_raw(:,5)/quanti_p(1)*ratio0;
                
                h_area2 = sum(h2{r_size2 == rsize0}(:));
                foci_raw2(:,3) = nM_con*h_area2*foci_raw2(:,9).*foci_raw2(:,13);
                fake_raw2(:,3) = nM_con*h_area2*fake_raw2(:,9).*fake_raw2(:,13);
                foci_raw2(:,4) = foci_raw2(:,4)-foci_raw2(:,18)/quanti_p(1)*ratio0;
                fake_raw2(:,4) = fake_raw2(:,4)-fake_raw2(:,18)/quanti_p(1)*ratio0;

                p_profile = nucleus_protein_profile_ab(:,2);
                reg_I = (nucleus_protein_profile_ab(:,3) >= z_range(1)) & (nucleus_protein_profile_ab(:,3) <= z_range(2));
                [nucleus_RNA_profile0,ind_foci] = foci_info(nucleus_RNA_profile,foci_RNA_profile);
                p_profile0 = p_profile(ind_foci);
                reg_I0 = reg_I(ind_foci);
                f_ob = nucleus_RNA_profile0(reg_I0,4);
                p_ob = p_profile0(reg_I0);
                EL_ob = nucleus_RNA_profile0(reg_I0,1);
                [p_ob,IX] = sort(p_ob);
                f_ob = f_ob(IX);
                EL_ob = EL_ob(IX);
                RNA_foci0 = cat(1,RNA_foci0,{[p_ob,f_ob,EL_ob]});
                
                enrich_foci0 = cat(1,enrich_foci0,{[foci_raw(:,[2,4,3]),fake_raw(:,[2,4,3]),foci_raw(:,[7:8,15])]});
                enrich_foci20 = cat(1,enrich_foci20,{[foci_raw2(:,[2,4,3]),fake_raw2(:,[2,4,3]),foci_raw2(:,[7:8,15])]});
% %                 enrich_foci0 = cat(1,enrich_foci0,{[fake_raw(:,[2,4,3]),fake_raw(:,[2,4,3]),foci_raw(:,[7:8,15])]});
% %                 enrich_foci20 = cat(1,enrich_foci20,{[fake_raw2(:,[2,4,3]),fake_raw2(:,[2,4,3]),foci_raw2(:,[7:8,15])]});
                
                nucleus_protein_profile_all0 = cat(1,nucleus_protein_profile_all0,{nucleus_protein_profile_ab(Iz,:)});
                nucleus_RNA_profile_all0 = cat(1,nucleus_RNA_profile_all0,{nucleus_RNA_profile(Iz,:)});
                nucleus_signal2_profile_all0 = cat(1,nucleus_signal2_profile_all0,{nucleus_signal2_profile(Iz,:)});

%             end
            
        end
    end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    N_cycle = length(cycle_range)+1;
    En_cycle = nan(N_cycle,4);
    En_cycle_fake = nan(N_cycle,4);
    En_nu = nan(N_cycle,2);
    En_name = cell(N_cycle,1);
    for I_cycle = 1:N_cycle
        if I_cycle < N_cycle
            cycle_range_temp = cycle_range(I_cycle);
            cycle_name = num2str(cycle_range(I_cycle));
        else
%             cycle_range_temp = cycle_range;
%             cycle_name = 'All';
            cycle_range_temp = cycle_range2;
            cycle_name = ['Cycle: ',num2str(cycle_range2(1)),' - ',num2str(cycle_range2(end))];
        end
        if any(ismember(cycle_all,cycle_range_temp))
            RNA_foci = cell2mat(RNA_foci0(ismember(cycle_all,cycle_range_temp)));
            enrich_foci = cell2mat(enrich_foci0(ismember(cycle_all,cycle_range_temp)));
            enrich_foci2 = cell2mat(enrich_foci20(ismember(cycle_all,cycle_range_temp)));
            nucleus_protein_profile_all = cell2mat(nucleus_protein_profile_all0(ismember(cycle_all,cycle_range_temp)));
            nucleus_RNA_profile_all = cell2mat(nucleus_RNA_profile_all0(ismember(cycle_all,cycle_range_temp)));
            nucleus_signal2_profile_all = cell2mat(nucleus_signal2_profile_all0(ismember(cycle_all,cycle_range_temp)));
        else
            RNA_foci = nan(1,size(RNA_foci0{1},2));
            enrich_foci = nan(1,size(enrich_foci0{1},2));
            enrich_foci2 = nan(1,size(enrich_foci20{1},2));
            nucleus_protein_profile_all = nan(1,size(nucleus_protein_profile_all0{1},2));
            nucleus_RNA_profile_all = nan(1,size(nucleus_RNA_profile_all0{1},2));
            nucleus_signal2_profile_all = nan(1,size(nucleus_signal2_profile_all0{1},2));
        end
        En_name{I_cycle} = cycle_name;
        
%% Enrichment vs. [Bcd] for RNA and RNA2 channels
        %%% Equal population bin:
        rz = foci_data{r_size == rsize0}(1,16);
        Itrue = (~isnan(enrich_foci(:,2))) & (~isnan(enrich_foci(:,3))) & enrich_foci(:,7) >= EL_range(1) & enrich_foci(:,7) <= EL_range(2) & (enrich_foci(:,8) >= z_range(1)) & (enrich_foci(:,8) <= z_range(2));
        Itrue0 = RNA_foci(:,3) >= EL_range(1) & RNA_foci(:,3) <= EL_range(2);
        Itrue000 = (~isnan(enrich_foci(:,2))) & (~isnan(enrich_foci(:,3)));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        enrich_foci000 = enrich_foci(Itrue000,:);
        enrich_foci = enrich_foci(Itrue,:);
        xx = enrich_foci(:,1)/1e-9/rc_cali; %xx(xx<0) = 0;
        yy = (enrich_foci(:,2)-enrich_foci(:,3))./enrich_foci(:,9)/rz/r_cali;
        zz = (enrich_foci(:,5)-enrich_foci(:,6))./enrich_foci(:,9)/rz/r_cali;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        rrx = RNA_foci(:,1);
        rrr = RNA_foci(:,2);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        yy000 = (enrich_foci000(:,2)-enrich_foci000(:,3))./enrich_foci000(:,9)/rz/r_cali;
        zz000 = (enrich_foci000(:,5)-enrich_foci000(:,6))./enrich_foci000(:,9)/rz/r_cali;
        En_cycle(I_cycle,1) = mean(yy000);
        En_cycle(I_cycle,2) = std0(yy000);
        En_cycle_fake(I_cycle,1) = mean(zz000);
        En_cycle_fake(I_cycle,2) = std0(zz000);
        En_nu(I_cycle,1) = length(yy000);
        
%         [x_N,y_N,xerr_N,yerr_N,~,x_start,x_end] = equal_bin(xx,yy,NbinE,NmoveE);
        [x_N,y_N,xerr_N,yerr_N,~,x_start,x_end] = equal_bin(xx,yy,NbinE,roverE,1);
        [rmax,I_max] = max(y_N);
        if ~isempty(I_max)
            [~,I_middle] = min(abs(y_N(1:I_max)-rmax/2));
            beta0 = [hc1,x_N(I_middle),rmax,0];
        else
            beta0 = [hc1,cc1,10,0];
        end
        I_fit = x_N >= fitlimN(1) & x_N <= fitlimN(2);
        try
            [beta1,r,~,~,~] = nlinfit(x_N(I_fit),y_N(I_fit),@Hill,beta0);
        catch err
            beta1 = nan(size(beta0));
        end
        [rx_N,ry_N,rxerr_N,ryerr_N] = equal_dist(rrx,rrr,x_start,x_end);
        Itruer = x_N >= rclim(1) & x_N <= rclim(2);

        rz = foci_data2{r_size2 == rsize0}(1,16);
        Itrue2 = (~isnan(enrich_foci2(:,2))) & (~isnan(enrich_foci2(:,3))) & enrich_foci2(:,7) >= EL_range(1) & enrich_foci2(:,7) <= EL_range(2) & (enrich_foci2(:,8) >= z_range(1)) & (enrich_foci2(:,8) <= z_range(2));
        Itrue2000 = (~isnan(enrich_foci2(:,2))) & (~isnan(enrich_foci2(:,3)));
        
        enrich_foci2000 = enrich_foci2(Itrue2000,:);
        enrich_foci2 = enrich_foci2(Itrue2,:);
        
        xx2 = enrich_foci2(:,1)/1e-9/rc_cali; %xx(xx<0) = 0;
        yy2 = (enrich_foci2(:,2)-enrich_foci2(:,3))./enrich_foci2(:,9)/rz/r_cali;
        zz2 = (enrich_foci2(:,5)-enrich_foci2(:,6))./enrich_foci2(:,9)/rz/r_cali;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        yy2000 = (enrich_foci2000(:,2)-enrich_foci2000(:,3))./enrich_foci2000(:,9)/rz/r_cali;
        zz2000 = (enrich_foci2000(:,5)-enrich_foci2000(:,6))./enrich_foci2000(:,9)/rz/r_cali;
        En_cycle(I_cycle,3) = mean(yy2000);
        En_cycle(I_cycle,4) = std0(yy2000);
        En_cycle_fake(I_cycle,3) = mean(zz2000);
        En_cycle_fake(I_cycle,4) = std0(zz2000);
        En_nu(I_cycle,2) = length(yy2000);

%         [x_N2,y_N2,xerr_N2,yerr_N2,~,x_start2,x_end2] = equal_bin(xx2,yy2,NbinE,NmoveE);
        [x_N2,y_N2,xerr_N2,yerr_N2,~,x_start2,x_end2] = equal_bin(xx2,yy2,NbinE,roverE,1);
        [rmax,I_max] = max(y_N2);
        if ~isempty(I_max)
            [~,I_middle] = min(abs(y_N2(1:I_max)-rmax/2));
            beta0 = [hc1,x_N2(I_middle),rmax,0];
        else
            beta0 = [hc1,cc1,10,0];
        end
        I_fit = x_N2 >= fitlimN(1) & x_N2 <= fitlimN(2);
        try
            [beta2,r,~,~,~] = nlinfit(x_N2(I_fit),y_N2(I_fit),@Hill,beta0);
        catch err
            beta2 = nan(size(beta0));
        end
        [rx_N2,ry_N2,rxerr_N2,ryerr_N2] = equal_dist(rrx,rrr,x_start2,x_end2);
        Itruer2 = x_N2 >= rclim(1) & x_N2 <= rclim(2);

        figure(15)
        ha1 = subaxis(4,N_cycle,I_cycle,1, 'Spacing', 0, 'PaddingRight',0.03, 'PaddingLeft',0.05, 'PaddingTop',0.03, 'PaddingBottom',0.03, 'MarginLeft', 0.01, 'MarginRight', 0.001, 'MarginTop', 0.015, 'MarginBottom', 0.025);
            errorbarxy(x_N,y_N,xerr_N,yerr_N,xerr_N,yerr_N,ccode(1,:),ccode(1,:),'.');
            hold on
            plot([xlim1(1):xlim_bin1:xlim1(2)],Hill(beta1,[xlim1(1):xlim_bin1:xlim1(2)]),'k-','DisplayName','Hill fit');
            text(mean(xlim1),mean(ylim4),['h = ',num2str(beta1(1),'%5.2f'),char(10),'C0 = ',num2str(beta1(2),'%5.2f'),' nM',char(10),'ymax = ',num2str(beta1(3),'%5.2f')])
            xlabel('[Protein] (nM)')
            ylabel('Protein bound')
            xlim(xlim2)
            ylim(ylim4)
            title(['Protein enrichment vs. [Protein] for Channel 1, Cycle: ',cycle_name])
        ha1 = subaxis(4,N_cycle,I_cycle,2, 'Spacing', 0, 'PaddingRight',0.03, 'PaddingLeft',0.05, 'PaddingTop',0.03, 'PaddingBottom',0.03, 'MarginLeft', 0.01, 'MarginRight', 0.001, 'MarginTop', 0.015, 'MarginBottom', 0.025);
            errorbarxy(x_N2,y_N2,xerr_N2,yerr_N2,xerr_N2,yerr_N2,ccode(2,:),ccode(2,:),'.');
            hold on
            plot([xlim1(1):xlim_bin1:xlim1(2)],Hill(beta2,[xlim1(1):xlim_bin1:xlim1(2)]),'k-','DisplayName','Hill fit');
            text(mean(xlim1),mean(ylim4),['h = ',num2str(beta2(1),'%5.2f'),char(10),'C0 = ',num2str(beta2(2),'%5.2f'),' nM',char(10),'ymax = ',num2str(beta2(3),'%5.2f')])
            xlabel('[Protein] (nM)')
            ylabel('Protein bound')
            xlim(xlim2)
            ylim(ylim4)
            title(['Protein enrichment vs. [Protein] for Channel 2, Cycle: ',cycle_name])
        ha1 = subaxis(4,N_cycle,I_cycle,3, 'Spacing', 0, 'PaddingRight',0.03, 'PaddingLeft',0.05, 'PaddingTop',0.03, 'PaddingBottom',0.03, 'MarginLeft', 0.01, 'MarginRight', 0.001, 'MarginTop', 0.015, 'MarginBottom', 0.025);
            errorbarxy(y_N,ry_N,yerr_N,ryerr_N,yerr_N,ryerr_N,ccode(1,:)*rc+(1-rc),ccode(1,:)*rc+(1-rc),'.','DisplayName','Raw Data');
            hold on
            errorbarxy(y_N(Itruer),ry_N(Itruer),yerr_N(Itruer),ryerr_N(Itruer),yerr_N(Itruer),ryerr_N(Itruer),ccode(1,:),ccode(1,:),'.','DisplayName','Gated Data');
            p = polyfit(y_N(Itruer),ry_N(Itruer),1);
            plot(ylim4,ylim4*p(1)+p(2),'k-','DisplayName','Fit');
            xlabel('Protein bound')
            ylabel('#RNA/foci')
            xlim(ylim4)
            ylim(ylim20)
            legend('show')
            title(['TX vs. Protein enrichment for Channel 1, Cycle: ',cycle_name])
        ha1 = subaxis(4,N_cycle,I_cycle,4, 'Spacing', 0, 'PaddingRight',0.03, 'PaddingLeft',0.05, 'PaddingTop',0.03, 'PaddingBottom',0.03, 'MarginLeft', 0.01, 'MarginRight', 0.001, 'MarginTop', 0.015, 'MarginBottom', 0.025);
            errorbarxy(y_N2,ry_N2,yerr_N2,ryerr_N2,yerr_N2,ryerr_N2,ccode(2,:)*rc+(1-rc),ccode(2,:)*rc+(1-rc),'.','DisplayName','Raw Data');
            hold on
            errorbarxy(y_N2(Itruer2),ry_N2(Itruer2),yerr_N2(Itruer2),ryerr_N2(Itruer2),yerr_N2(Itruer2),ryerr_N2(Itruer2),ccode(2,:),ccode(2,:),'.','DisplayName','Gated Data');
            p = polyfit(y_N2(Itruer2),ry_N2(Itruer2),1);
            plot(ylim4,ylim4*p(1)+p(2),'k-','DisplayName','Fit');
            xlabel('Protein bound')
            ylabel('#RNA/foci')
            xlim(ylim4)
            ylim(ylim20)
            legend('show')
            title(['TX vs. Protein enrichment for Channel 2, Cycle: ',cycle_name])

        %%% Equal distance bin
        [x_N,y_N,xerr_N,yerr_N] = equal_dist(xx,yy,bin_min2,bin_max2);
        [rmax,I_max] = max(y_N);
        if ~isempty(I_max)
            [~,I_middle] = min(abs(y_N(1:I_max)-rmax/2));
            beta0 = [hc1,x_N(I_middle),rmax,0];
        else
            beta0 = [hc1,cc1,10,0];
        end
        I_fit = x_N >= fitlimC(1) & x_N <= fitlimC(2);
        try
            [beta1,r,~,~,~] = nlinfit(x_N(I_fit),y_N(I_fit),@Hill,beta0);
        catch err
            beta1 = nan(size(beta0));
        end
        [rx_N,ry_N,rxerr_N,ryerr_N] = equal_dist(rrx,rrr,bin_min2,bin_max2);
        Itruer = x_N >= rclim(1) & x_N <= rclim(2);

        [x_N2,y_N2,xerr_N2,yerr_N2] = equal_dist(xx2,yy2,bin_min2,bin_max2);
        [rmax,I_max] = max(y_N2);
        if ~isempty(I_max)
            [~,I_middle] = min(abs(y_N2(1:I_max)-rmax/2));
            beta0 = [hc1,x_N2(I_middle),rmax,0];
        else
            beta0 = [hc1,cc1,10,0];
        end
        I_fit = x_N2 >= fitlimC(1) & x_N2 <= fitlimC(2);
        try
            [beta2,r,~,~,~] = nlinfit(x_N2(I_fit),y_N2(I_fit),@Hill,beta0);
        catch err
            beta2 = nan(size(beta0));
        end
        [rx_N2,ry_N2,rxerr_N2,ryerr_N2] = equal_dist(rrx,rrr,bin_min2,bin_max2);
        Itruer2 = x_N2 >= rclim(1) & x_N2 <= rclim(2);

        figure(16)
        ha1 = subaxis(4,N_cycle,I_cycle,1, 'Spacing', 0, 'PaddingRight',0.03, 'PaddingLeft',0.05, 'PaddingTop',0.03, 'PaddingBottom',0.03, 'MarginLeft', 0.01, 'MarginRight', 0.001, 'MarginTop', 0.015, 'MarginBottom', 0.025);
            errorbarxy(x_N,y_N,xerr_N,yerr_N,xerr_N,yerr_N,ccode(1,:),ccode(1,:),'.');
            hold on
            plot([xlim1(1):xlim_bin1:xlim1(2)],Hill(beta1,[xlim1(1):xlim_bin1:xlim1(2)]),'k-','DisplayName','Hill fit');
            text(mean(xlim1),mean(ylim4),['h = ',num2str(beta1(1),'%5.2f'),char(10),'C0 = ',num2str(beta1(2),'%5.2f'),' nM',char(10),'ymax = ',num2str(beta1(3),'%5.2f')])
            xlabel('[Protein] (nM)')
            ylabel('Protein bound')
            xlim(xlim2)
            ylim(ylim4)
            title(['Protein enrichment vs. [Protein] for Channel 1, Cycle: ',cycle_name])
        ha1 = subaxis(4,N_cycle,I_cycle,2, 'Spacing', 0, 'PaddingRight',0.03, 'PaddingLeft',0.05, 'PaddingTop',0.03, 'PaddingBottom',0.03, 'MarginLeft', 0.01, 'MarginRight', 0.001, 'MarginTop', 0.015, 'MarginBottom', 0.025);
            errorbarxy(x_N2,y_N2,xerr_N2,yerr_N2,xerr_N2,yerr_N2,ccode(2,:),ccode(2,:),'.');
            hold on
            plot([xlim1(1):xlim_bin1:xlim1(2)],Hill(beta2,[xlim1(1):xlim_bin1:xlim1(2)]),'k-','DisplayName','Hill fit');
            text(mean(xlim1),mean(ylim4),['h = ',num2str(beta2(1),'%5.2f'),char(10),'C0 = ',num2str(beta2(2),'%5.2f'),' nM',char(10),'ymax = ',num2str(beta2(3),'%5.2f')])
            xlabel('[Protein] (nM)')
            ylabel('Protein bound')
            xlim(xlim2)
            ylim(ylim4)
            title(['Protein enrichment vs. [Protein] for Channel 2, Cycle: ',cycle_name])
        ha1 = subaxis(4,N_cycle,I_cycle,3, 'Spacing', 0, 'PaddingRight',0.03, 'PaddingLeft',0.05, 'PaddingTop',0.03, 'PaddingBottom',0.03, 'MarginLeft', 0.01, 'MarginRight', 0.001, 'MarginTop', 0.015, 'MarginBottom', 0.025);
            errorbarxy(y_N,ry_N,yerr_N,ryerr_N,yerr_N,ryerr_N,ccode(1,:)*rc+(1-rc),ccode(1,:)*rc+(1-rc),'.','DisplayName','Raw Data');
            hold on
            errorbarxy(y_N(Itruer),ry_N(Itruer),yerr_N(Itruer),ryerr_N(Itruer),yerr_N(Itruer),ryerr_N(Itruer),ccode(1,:),ccode(1,:),'.','DisplayName','Gated Data');
            p = polyfit(y_N(Itruer),ry_N(Itruer),1);
            plot(ylim4,ylim4*p(1)+p(2),'k-','DisplayName','Fit');
            xlabel('Protein bound')
            ylabel('#RNA/foci')
            xlim(ylim4)
            ylim(ylim20)
            legend('show')
            title(['TX vs. Protein enrichment for Channel 1, Cycle: ',cycle_name])
        ha1 = subaxis(4,N_cycle,I_cycle,4, 'Spacing', 0, 'PaddingRight',0.03, 'PaddingLeft',0.05, 'PaddingTop',0.03, 'PaddingBottom',0.03, 'MarginLeft', 0.01, 'MarginRight', 0.001, 'MarginTop', 0.015, 'MarginBottom', 0.025);
            errorbarxy(y_N2,ry_N2,yerr_N2,ryerr_N2,yerr_N2,ryerr_N2,ccode(2,:)*rc+(1-rc),ccode(2,:)*rc+(1-rc),'.','DisplayName','Raw Data');
            hold on
            errorbarxy(y_N2(Itruer2),ry_N2(Itruer2),yerr_N2(Itruer2),ryerr_N2(Itruer2),yerr_N2(Itruer2),ryerr_N2(Itruer2),ccode(2,:),ccode(2,:),'.','DisplayName','Gated Data');
            p = polyfit(y_N2(Itruer2),ry_N2(Itruer2),1);
            plot(ylim4,ylim4*p(1)+p(2),'k-','DisplayName','Fit');
            xlabel('Protein bound')
            ylabel('#RNA/foci')
            xlim(ylim4)
            ylim(ylim20)
            legend('show')
            title(['TX vs. Protein enrichment for Channel 2, Cycle: ',cycle_name])
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        


%% TX regulation analysis (all embryos pooled together): %%%%%%%%%%%%%%%%%%
        figure(17)
        %%% #foci/nu vs. [Bcd] for RNA channel
        ha1 = subaxis(2,N_cycle,I_cycle,1, 'Spacing', 0, 'PaddingRight',0.03, 'PaddingLeft',0.05, 'PaddingTop',0.03, 'PaddingBottom',0.03, 'MarginLeft', 0.01, 'MarginRight', 0.001, 'MarginTop', 0.015, 'MarginBottom', 0.025);
            [xbin,ybin,xerr,yerr,N_in] = equal_bin(nucleus_protein_profile_all(:,2),nucleus_RNA_profile_all(:,3),Nbin,rover,1);
            I00 = ~isnan(xbin) & ~isnan(ybin) & N_in >= 5;
            xbin = xbin(I00);
            ybin = ybin(I00);
            xerr = xerr(I00);
            yerr = yerr(I00);
            N_in = N_in(I00);
            [rmax,I_max] = max(ybin);
            [rmin,I_min] = min(ybin);
            if ~isempty(I_max)
                [~,I_middle] = min(abs(ybin(1:I_max)-rmax/2));
%                 [~,I_middle] = min(abs(ybin(I_max:end)-rmax/2));
                beta0 = [hc1,xbin(I_middle),rmax,rmin];
            else
                beta0 = [hc1,cc1,2,rmin];
            end
            I_fit = xbin >= fitlim2(1) & xbin <= fitlim2(2);
            
            jj0 = I_cycle*2-1;
            if isempty(rfoci_pool_prefit) || size(rfoci_pool_prefit,1) < jj0 || rfoci_pool_prefit(jj0,1) <= 0
                try
    % %                 [beta1,r,~,~,~] = nlinfit(xbin(I_fit),ybin(I_fit),@Hill,beta0);
    %                 [beta1,r,~,~,~] = nlinfit(xbin(I_fit),ybin(I_fit),Hill1,beta0);
                    curve = fit(xbin(I_fit)',ybin(I_fit)', g0, 'StartPoint',beta0,'Lower',[0,0,0,0] );
                    beta1 = [curve.h,curve.c,curve.a,curve.d];
                catch err
                    beta1 = nan(size(beta0));
                    curve.h = nan; curve.c = nan; curve.a = nan; curve.d = nan;
                end
            else
                beta1 = rfoci_pool_prefit(jj0,2:end);
                curve = fit(xbin(I_fit)',ybin(I_fit)', g0, 'StartPoint',beta0,'Lower',[0,0,0,0] );
                curve.h = beta1(1); curve.c = beta1(2); curve.a = beta1(3); curve.d = beta1(4);
            end
            
            errorbarxy(xbin,ybin,xerr,yerr,xerr,yerr,'k','k','.','DisplayName','Binned data');
            hold on
% %             plot([xlim1(1):xlim_bin1:xlim1(2)],Hill(beta1,[xlim1(1):xlim_bin1:xlim1(2)]),'r-','DisplayName','Hill fit');
%             plot([xlim1(1):xlim_bin1:xlim1(2)],Hill1(beta1,[xlim1(1):xlim_bin1:xlim1(2)]),'r-','DisplayName','Hill fit');
            plot([xlim1(1):xlim_bin1:xlim1(2)],feval(curve,[xlim1(1):xlim_bin1:xlim1(2)]),'r-','DisplayName','Hill fit');
            text(mean(xlim1),mean(ylim1),['h = ',num2str(beta1(1),'%5.2f'),char(10),'C0 = ',num2str(beta1(2),'%5.2f'),' nM',char(10),'ymax = ',num2str(beta1(3),'%5.2f')])
            xlim(xlim1)
            ylim(ylim1)
            xlabel('[Protein] (nM)')
            ylabel('#foci/nu')
            title(['#foci/nu vs. [Bcd] for RNA (pooled), Cycle: ',cycle_name],'interpreter','none')
            hold off

        %%% #foci/nu vs. [Bcd] for RNA2 channel
        ha1 = subaxis(2,N_cycle,I_cycle,2, 'Spacing', 0, 'PaddingRight',0.03, 'PaddingLeft',0.05, 'PaddingTop',0.03, 'PaddingBottom',0.03, 'MarginLeft', 0.01, 'MarginRight', 0.001, 'MarginTop', 0.015, 'MarginBottom', 0.025);
            [xbin,ybin,xerr,yerr,N_in] = equal_bin(nucleus_protein_profile_all(:,2),nucleus_signal2_profile_all(:,3),Nbin,rover,1);
            I00 = ~isnan(xbin) & ~isnan(ybin) & N_in >= 5;
            xbin = xbin(I00);
            ybin = ybin(I00);
            xerr = xerr(I00);
            yerr = yerr(I00);
            N_in = N_in(I00);
            [rmax,I_max] = max(ybin);
            if ~isempty(I_max)
                [~,I_middle] = min(abs(ybin(1:I_max)-rmax/2));
%                 [~,I_middle] = min(abs(ybin(I_max:end)-rmax/2));
                beta0 = [hc1,xbin(I_middle),rmax,0];
            else
                beta0 = [hc1,cc1,10,0];
            end
            I_fit = xbin >= fitlim2(1) & xbin <= fitlim2(2);
            
            jj0 = I_cycle*2;
            if isempty(rfoci_pool_prefit) || size(rfoci_pool_prefit,1) < jj0 || rfoci_pool_prefit(jj0,1) <= 0
                try
    % %                 [beta1,r,~,~,~] = nlinfit(xbin(I_fit),ybin(I_fit),@Hill,beta0);
    %                 [beta1,r,~,~,~] = nlinfit(xbin(I_fit),ybin(I_fit),Hill1,beta0);
                    curve = fit(xbin(I_fit)',ybin(I_fit)', g0, 'StartPoint',beta0,'Lower',[0,0,0,0] );
                    beta1 = [curve.h,curve.c,curve.a,curve.d];
                catch err
                    beta1 = nan(size(beta0));
                    curve.h = nan; curve.c = nan; curve.a = nan; curve.d = nan;
                end
            else
                beta1 = rfoci_pool_prefit(jj0,2:end);
                curve = fit(xbin(I_fit)',ybin(I_fit)', g0, 'StartPoint',beta0,'Lower',[0,0,0,0] );
                curve.h = beta1(1); curve.c = beta1(2); curve.a = beta1(3); curve.d = beta1(4);
            end
                
            errorbarxy(xbin,ybin,xerr,yerr,xerr,yerr,'k','k','.','DisplayName','Binned data');
            hold on
% %             plot([xlim1(1):xlim_bin1:xlim1(2)],Hill(beta1,[xlim1(1):xlim_bin1:xlim1(2)]),'r-','DisplayName','Hill fit');
%             plot([xlim1(1):xlim_bin1:xlim1(2)],Hill1(beta1,[xlim1(1):xlim_bin1:xlim1(2)]),'r-','DisplayName','Hill fit');
            plot([xlim1(1):xlim_bin1:xlim1(2)],feval(curve,[xlim1(1):xlim_bin1:xlim1(2)]),'r-','DisplayName','Hill fit');
            text(mean(xlim1),mean(ylim1),['h = ',num2str(beta1(1),'%5.2f'),char(10),'C0 = ',num2str(beta1(2),'%5.2f'),' nM',char(10),'ymax = ',num2str(beta1(3),'%5.2f')])
            xlim(xlim1)
            ylim(ylim1)
            xlabel('[Protein] (nM)')
            ylabel('#foci/nu')
            title(['#foci/nu vs. [Bcd] for RNA2 (pooled), Cycle: ',cycle_name],'interpreter','none')
            hold off

        figure(18)
        %%% #RNA/nu vs. [Bcd] for RNA channel
        ha1 = subaxis(2,N_cycle,I_cycle,1, 'Spacing', 0, 'PaddingRight',0.03, 'PaddingLeft',0.05, 'PaddingTop',0.03, 'PaddingBottom',0.03, 'MarginLeft', 0.01, 'MarginRight', 0.001, 'MarginTop', 0.015, 'MarginBottom', 0.025);
            [xbin,ybin,xerr,yerr,N_in] = equal_bin(nucleus_protein_profile_all(:,2),nucleus_RNA_profile_all(:,4),Nbin,rover,1);
            I00 = ~isnan(xbin) & ~isnan(ybin) & N_in >= 5;
            xbin = xbin(I00);
            ybin = ybin(I00);
            xerr = xerr(I00);
            yerr = yerr(I00);
            N_in = N_in(I00);
            [rmax,I_max] = max(ybin);
            if ~isempty(I_max)
                [~,I_middle] = min(abs(ybin(1:I_max)-rmax/2));
%                 [~,I_middle] = min(abs(ybin(I_max:end)-rmax/2));
                beta0 = [hc1,xbin(I_middle),rmax,0];
            else
                beta0 = [hc1,cc1,10,0];
            end
            I_fit = xbin >= fitlim2(1) & xbin <= fitlim2(2);
            
            jj0 = I_cycle*2-1;
            if isempty(rfish_pool_prefit) || size(rfish_pool_prefit,1) < jj0 || rfish_pool_prefit(jj0,1) <= 0
                try
    % %                 [beta2,r,~,~,~] = nlinfit(xbin(I_fit),ybin(I_fit),@Hill,beta0);
    %                 [beta2,r,~,~,~] = nlinfit(xbin(I_fit),ybin(I_fit),Hill1,beta0);
                    curve = fit(xbin(I_fit)',ybin(I_fit)', g0, 'StartPoint',beta0,'Lower',[0,0,0,0] );
                    beta2 = [curve.h,curve.c,curve.a,curve.d];
                catch err
                    beta2 = nan(size(beta0));
                    curve.h = nan; curve.c = nan; curve.a = nan; curve.d = nan;
                end
            else
                beta2 = rfish_pool_prefit(jj0,2:end);
                curve = fit(xbin(I_fit)',ybin(I_fit)', g0, 'StartPoint',beta0,'Lower',[0,0,0,0] );
                curve.h = beta2(1); curve.c = beta2(2); curve.a = beta2(3); curve.d = beta2(4);
            end
            
            errorbarxy(xbin,ybin,xerr,yerr,xerr,yerr,'k','k','.','DisplayName','Binned data');
            hold on
%             plot([xlim1(1):xlim_bin1:xlim1(2)],Hill(beta2,[xlim1(1):xlim_bin1:xlim1(2)]),'r-','DisplayName','Hill fit');
            plot([xlim1(1):xlim_bin1:xlim1(2)],feval(curve,[xlim1(1):xlim_bin1:xlim1(2)]),'r-','DisplayName','Hill fit');
            text(mean(xlim1),mean(ylim2),['h = ',num2str(beta2(1),'%5.2f'),char(10),'C0 = ',num2str(beta2(2),'%5.2f'),' nM',char(10),'ymax = ',num2str(beta2(3),'%5.2f')])
            xlim(xlim1)
            ylim(ylim2)
            xlabel('[Protein] (nM)')
            ylabel('#RNA/nu')
            title(['#RNA/nu vs. [Bcd] for RNA (pooled), Cycle: ',cycle_name],'interpreter','none')
            hold off

        %%% #RNA/nu vs. [Bcd] for RNA2 channel
        ha1 = subaxis(2,N_cycle,I_cycle,2, 'Spacing', 0, 'PaddingRight',0.03, 'PaddingLeft',0.05, 'PaddingTop',0.03, 'PaddingBottom',0.03, 'MarginLeft', 0.01, 'MarginRight', 0.001, 'MarginTop', 0.015, 'MarginBottom', 0.025);
            [xbin,ybin,xerr,yerr,N_in] = equal_bin(nucleus_protein_profile_all(:,2),nucleus_signal2_profile_all(:,4),Nbin,rover,1);
            I00 = ~isnan(xbin) & ~isnan(ybin) & N_in >= 5;
            xbin = xbin(I00);
            ybin = ybin(I00);
            xerr = xerr(I00);
            yerr = yerr(I00);
            N_in = N_in(I00);
            [rmax,I_max] = max(ybin);
            if ~isempty(I_max)
                [~,I_middle] = min(abs(ybin(1:I_max)-rmax/2));
%                 [~,I_middle] = min(abs(ybin(I_max:end)-rmax/2));
                beta0 = [hc1,xbin(I_middle),rmax,0];
            else
                beta0 = [hc1,cc1,10,0];
            end
            I_fit = xbin >= fitlim2(1) & xbin <= fitlim2(2);
            
            jj0 = I_cycle*2;
            if isempty(rfish_pool_prefit) || size(rfish_pool_prefit,1) < jj0 || rfish_pool_prefit(jj0,1) <= 0
                try
    % %                 [beta2,r,~,~,~] = nlinfit(xbin(I_fit),ybin(I_fit),@Hill,beta0);
    %                 [beta2,r,~,~,~] = nlinfit(xbin(I_fit),ybin(I_fit),Hill1,beta0);
                    curve = fit(xbin(I_fit)',ybin(I_fit)', g0, 'StartPoint',beta0,'Lower',[0,0,0,0] );
                    beta2 = [curve.h,curve.c,curve.a,curve.d];
                catch err
                    beta2 = nan(size(beta0));
                    curve.h = nan; curve.c = nan; curve.a = nan; curve.d = nan;
                end
            else
                beta2 = rfish_pool_prefit(jj0,2:end);
                curve = fit(xbin(I_fit)',ybin(I_fit)', g0, 'StartPoint',beta0,'Lower',[0,0,0,0] );
                curve.h = beta2(1); curve.c = beta2(2); curve.a = beta2(3); curve.d = beta2(4);
            end
            
            errorbarxy(xbin,ybin,xerr,yerr,xerr,yerr,'k','k','.','DisplayName','Binned data');
            hold on
%             plot([xlim1(1):xlim_bin1:xlim1(2)],Hill(beta2,[xlim1(1):xlim_bin1:xlim1(2)]),'r-','DisplayName','Hill fit');
            plot([xlim1(1):xlim_bin1:xlim1(2)],feval(curve,[xlim1(1):xlim_bin1:xlim1(2)]),'r-','DisplayName','Hill fit');
            text(mean(xlim1),mean(ylim2b),['h = ',num2str(beta2(1),'%5.2f'),char(10),'C0 = ',num2str(beta2(2),'%5.2f'),' nM',char(10),'ymax = ',num2str(beta2(3),'%5.2f')])
            xlim(xlim1)
            ylim(ylim2b)
            xlabel('[Protein] (nM)')
            ylabel('#RNA/nu')
            title(['#RNA/nu vs. [Bcd] for RNA2 (pooled), Cycle: ',cycle_name],'interpreter','none')
            hold off
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    end


%% Plots of fitting parameters: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% #foci/nu vs. EL:
    figure(19)
    for ii = 1:size(beta_all1,2)
        ha1 = subaxis(2,4,ii,1, 'Spacing', 0, 'PaddingRight',0.03, 'PaddingLeft',0.05, 'PaddingTop',0.03, 'PaddingBottom',0.03, 'MarginLeft', 0.01, 'MarginRight', 0.001, 'MarginTop', 0.015, 'MarginBottom', 0.025);
        error_plot_all(cycle_all,beta_all1(:,ii),cycle_range);
        xlim([cycle_range(1)-0.5,cycle_range(end)+0.5])
        ylim00 = ylim; [~,Ia] = min(abs(ylim00)); ylim00(Ia) = 0; ylim(ylim00);
        xlabel('Cycle')
        ylabel(xls_title{2+ii})
        title('#foci/nu vs. EL (RNA channel)')
    end
    for ii = 1:size(beta_all1b,2)
        ha1 = subaxis(2,4,ii,2, 'Spacing', 0, 'PaddingRight',0.03, 'PaddingLeft',0.05, 'PaddingTop',0.03, 'PaddingBottom',0.03, 'MarginLeft', 0.01, 'MarginRight', 0.001, 'MarginTop', 0.015, 'MarginBottom', 0.025);
        error_plot_all(cycle_all,beta_all1b(:,ii),cycle_range);
        xlim([cycle_range(1)-0.5,cycle_range(end)+0.5])
        ylim00 = ylim; [~,Ia] = min(abs(ylim00)); ylim00(Ia) = 0; ylim(ylim00);
        xlabel('Cycle')
        ylabel(xls_title{2+ii})
        title('#foci/nu vs. EL (RNA2 channel)')
    end

    %%% #RNA/nu vs. EL:
    figure(20)
    for ii = 1:size(beta_all2,2)
        ha1 = subaxis(2,4,ii,1, 'Spacing', 0, 'PaddingRight',0.03, 'PaddingLeft',0.05, 'PaddingTop',0.03, 'PaddingBottom',0.03, 'MarginLeft', 0.01, 'MarginRight', 0.001, 'MarginTop', 0.015, 'MarginBottom', 0.025);
        error_plot_all(cycle_all,beta_all2(:,ii),cycle_range);
        xlim([cycle_range(1)-0.5,cycle_range(end)+0.5])
        ylim00 = ylim; [~,Ia] = min(abs(ylim00)); ylim00(Ia) = 0; ylim(ylim00);
        xlabel('Cycle')
        ylabel(xls_title{2+ii})
        title('#RNA/nu vs. EL (RNA channel)')
    end
    for ii = 1:size(beta_all2b,2)
        ha1 = subaxis(2,4,ii,2, 'Spacing', 0, 'PaddingRight',0.03, 'PaddingLeft',0.05, 'PaddingTop',0.03, 'PaddingBottom',0.03, 'MarginLeft', 0.01, 'MarginRight', 0.001, 'MarginTop', 0.015, 'MarginBottom', 0.025);
        error_plot_all(cycle_all,beta_all2b(:,ii),cycle_range);
        xlim([cycle_range(1)-0.5,cycle_range(end)+0.5])
        ylim00 = ylim; [~,Ia] = min(abs(ylim00)); ylim00(Ia) = 0; ylim(ylim00);
        xlabel('Cycle')
        ylabel(xls_title{2+ii})
        title('#RNA/nu vs. EL (RNA2 channel)')
    end

    %%% #foci/nu vs. [Bcd]:
    figure(21)
    for ii = 1:size(beta_all4,2)
        ha1 = subaxis(2,4,ii,1, 'Spacing', 0, 'PaddingRight',0.03, 'PaddingLeft',0.05, 'PaddingTop',0.03, 'PaddingBottom',0.03, 'MarginLeft', 0.01, 'MarginRight', 0.001, 'MarginTop', 0.015, 'MarginBottom', 0.025);
        error_plot_all(cycle_all,beta_all4(:,ii),cycle_range);
        xlim([cycle_range(1)-0.5,cycle_range(end)+0.5])
        ylim00 = ylim; [~,Ia] = min(abs(ylim00)); ylim00(Ia) = 0; ylim(ylim00);
        xlabel('Cycle')
        ylabel(xls_title2{2+ii})
        title('#foci/nu vs. [Bcd] (RNA channel)')
    end
    for ii = 1:size(beta_all4b,2)
        ha1 = subaxis(2,4,ii,2, 'Spacing', 0, 'PaddingRight',0.03, 'PaddingLeft',0.05, 'PaddingTop',0.03, 'PaddingBottom',0.03, 'MarginLeft', 0.01, 'MarginRight', 0.001, 'MarginTop', 0.015, 'MarginBottom', 0.025);
        error_plot_all(cycle_all,beta_all4b(:,ii),cycle_range);
        xlim([cycle_range(1)-0.5,cycle_range(end)+0.5])
        ylim00 = ylim; [~,Ia] = min(abs(ylim00)); ylim00(Ia) = 0; ylim(ylim00);
        xlabel('Cycle')
        ylabel(xls_title2{2+ii})
        title('#foci/nu vs. [Bcd] (RNA2 channel)')
    end

    %%% #RNA/nu vs. [Bcd]:
    figure(22)
    for ii = 1:size(beta_all5,2)
        ha1 = subaxis(2,4,ii,1, 'Spacing', 0, 'PaddingRight',0.03, 'PaddingLeft',0.05, 'PaddingTop',0.03, 'PaddingBottom',0.03, 'MarginLeft', 0.01, 'MarginRight', 0.001, 'MarginTop', 0.015, 'MarginBottom', 0.025);
        error_plot_all(cycle_all,beta_all5(:,ii),cycle_range);
        xlim([cycle_range(1)-0.5,cycle_range(end)+0.5])
        ylim00 = ylim; [~,Ia] = min(abs(ylim00)); ylim00(Ia) = 0; ylim(ylim00);
        xlabel('Cycle')
        ylabel(xls_title2{2+ii})
        title('#RNA/nu vs. [Bcd] (RNA channel)')
    end
    for ii = 1:size(beta_all5b,2)
        ha1 = subaxis(2,4,ii,2, 'Spacing', 0, 'PaddingRight',0.03, 'PaddingLeft',0.05, 'PaddingTop',0.03, 'PaddingBottom',0.03, 'MarginLeft', 0.01, 'MarginRight', 0.001, 'MarginTop', 0.015, 'MarginBottom', 0.025);
        error_plot_all(cycle_all,beta_all5b(:,ii),cycle_range);
        xlim([cycle_range(1)-0.5,cycle_range(end)+0.5])
        ylim00 = ylim; [~,Ia] = min(abs(ylim00)); ylim00(Ia) = 0; ylim(ylim00);
        xlabel('Cycle')
        ylabel(xls_title2{2+ii})
        title('#RNA/nu vs. [Bcd] (RNA2 channel)')
    end

    %%% [Bcd] vs. EL:
    figure(23)
    for ii = 1:size(beta_all3,2)
        ha1 = subaxis(1,4,ii,1, 'Spacing', 0, 'PaddingRight',0.03, 'PaddingLeft',0.05, 'PaddingTop',0.03, 'PaddingBottom',0.03, 'MarginLeft', 0.01, 'MarginRight', 0.001, 'MarginTop', 0.015, 'MarginBottom', 0.025);
        error_plot_all(cycle_all,beta_all3(:,ii),cycle_range);
        xlim([cycle_range(1)-0.5,cycle_range(end)+0.5])
        ylim00 = ylim; [~,Ia] = min(abs(ylim00)); ylim00(Ia) = 0; ylim(ylim00);
        xlabel('Cycle')
        ylabel(xls_title3{2+ii})
        title('[Bcd] vs. EL')
    end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    

%% Data output: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if exist([ana_folder,embryo_name{iem},'/']) ~= 7
        mkdir([ana_folder,embryo_name{iem},'/']);
    end
    
    saveas(1,[ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,fate_add,D3_add,figure_tail]);
    saveas(2,[ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,signal2_add,fate_add,D3_add,figure_tail]);
    saveas(3,[ana_folder,embryo_name{iem},'/',embryo_name{iem},fish_add,num_add,D3_add,figure_tail]);
    saveas(4,[ana_folder,embryo_name{iem},'/',embryo_name{iem},fish_add,signal2_add,num_add,D3_add,figure_tail]);
    saveas(5,[ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,foci_add,D3_add,figure_tail]);
    saveas(6,[ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,foci_add,signal2_add,D3_add,figure_tail]);
    saveas(7,[ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,fish_add,int_add,D3_add,figure_tail]);
    saveas(8,[ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,fish_add,signal2_add,int_add,D3_add,figure_tail]);
    saveas(9,[ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,protein_add,fate_add,D3_add,figure_tail]);
    saveas(10,[ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,protein_add,int_add,D3_add,figure_tail]);
    saveas(101,[ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,protein_add,fluc_add,D3_add,figure_tail]);
    saveas(11,[ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,foci_add,reg_add,D3_add,figure_tail]);
    saveas(12,[ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,foci_add,signal2_add,reg_add,D3_add,figure_tail]);
    saveas(13,[ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,fish_add,reg_add,D3_add,figure_tail]);
    saveas(14,[ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,fish_add,signal2_add,reg_add,D3_add,figure_tail]);
    saveas(15,[ana_folder,embryo_name{iem},'/',embryo_name{iem},protein_add,enrich_add,binN_add,D3_add,figure_tail]);
    saveas(16,[ana_folder,embryo_name{iem},'/',embryo_name{iem},protein_add,enrich_add,binC_add,D3_add,figure_tail]);
    saveas(17,[ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,foci_add,reg_add,all_add,D3_add,figure_tail]);
    saveas(18,[ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,fish_add,reg_add,all_add,D3_add,figure_tail]);
    saveas(19,[ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,foci_add,para_add,D3_add,figure_tail]);
    saveas(20,[ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,fish_add,para_add,D3_add,figure_tail]);
    saveas(21,[ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,foci_add,reg_add,para_add,D3_add,figure_tail]);
    saveas(22,[ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,fish_add,reg_add,para_add,D3_add,figure_tail]);
    saveas(23,[ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,protein_add,int_add,para_add,D3_add,figure_tail]);

    saveas(1,[ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,fate_add,D3_add,figure_tail2]);
    saveas(2,[ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,signal2_add,fate_add,D3_add,figure_tail2]);
    saveas(3,[ana_folder,embryo_name{iem},'/',embryo_name{iem},fish_add,num_add,D3_add,figure_tail2]);
    saveas(4,[ana_folder,embryo_name{iem},'/',embryo_name{iem},fish_add,signal2_add,num_add,D3_add,figure_tail2]);
    saveas(5,[ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,foci_add,D3_add,figure_tail2]);
    saveas(6,[ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,foci_add,signal2_add,D3_add,figure_tail2]);
    saveas(7,[ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,fish_add,int_add,D3_add,figure_tail2]);
    saveas(8,[ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,fish_add,signal2_add,int_add,D3_add,figure_tail2]);
    saveas(9,[ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,protein_add,fate_add,D3_add,figure_tail2]);
    saveas(10,[ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,protein_add,int_add,D3_add,figure_tail2]);
    saveas(101,[ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,protein_add,fluc_add,D3_add,figure_tail2]);
    saveas(11,[ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,foci_add,reg_add,D3_add,figure_tail2]);
    saveas(12,[ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,foci_add,signal2_add,reg_add,D3_add,figure_tail2]);
    saveas(13,[ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,fish_add,reg_add,D3_add,figure_tail2]);
    saveas(14,[ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,fish_add,signal2_add,reg_add,D3_add,figure_tail2]);
    saveas(15,[ana_folder,embryo_name{iem},'/',embryo_name{iem},protein_add,enrich_add,binN_add,D3_add,figure_tail2]);
    saveas(16,[ana_folder,embryo_name{iem},'/',embryo_name{iem},protein_add,enrich_add,binC_add,D3_add,figure_tail2]);
    saveas(17,[ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,foci_add,reg_add,all_add,D3_add,figure_tail2]);
    saveas(18,[ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,fish_add,reg_add,all_add,D3_add,figure_tail2]);
    saveas(19,[ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,foci_add,para_add,D3_add,figure_tail2]);
    saveas(20,[ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,fish_add,para_add,D3_add,figure_tail2]);
    saveas(21,[ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,foci_add,reg_add,para_add,D3_add,figure_tail2]);
    saveas(22,[ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,fish_add,reg_add,para_add,D3_add,figure_tail2]);
    saveas(23,[ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,protein_add,int_add,para_add,D3_add,figure_tail2]);

    %%% Save the Hill fitting results to Excel files:
    out_ex0 = cat(1,xls_title,cat(2,name_all,num2cell(cycle_all),num2cell(beta_all1)));
    xlswrite([ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,foci_add,D3_add,output_tail],out_ex0);

    out_ex0 = cat(1,xls_title,cat(2,name_all,num2cell(cycle_all),num2cell(beta_all1b)));
    xlswrite([ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,foci_add,signal2_add,D3_add,output_tail],out_ex0);

    out_ex0 = cat(1,xls_title,cat(2,name_all,num2cell(cycle_all),num2cell(beta_all2)));
    xlswrite([ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,fish_add,int_add,D3_add,output_tail],out_ex0);

    out_ex0 = cat(1,xls_title,cat(2,name_all,num2cell(cycle_all),num2cell(beta_all2b)));
    xlswrite([ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,fish_add,signal2_add,int_add,D3_add,output_tail],out_ex0);

    out_ex0 = cat(1,xls_title2,cat(2,name_all,num2cell(cycle_all),num2cell(beta_all4)));
    xlswrite([ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,foci_add,reg_add,D3_add,output_tail],out_ex0);

    out_ex0 = cat(1,xls_title2,cat(2,name_all,num2cell(cycle_all),num2cell(beta_all4b)));
    xlswrite([ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,foci_add,reg_add,signal2_add,D3_add,output_tail],out_ex0);

    out_ex0 = cat(1,xls_title2,cat(2,name_all,num2cell(cycle_all),num2cell(beta_all5)));
    xlswrite([ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,fish_add,reg_add,D3_add,output_tail],out_ex0);

    out_ex0 = cat(1,xls_title2,cat(2,name_all,num2cell(cycle_all),num2cell(beta_all5b)));
    xlswrite([ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,fish_add,signal2_add,reg_add,D3_add,output_tail],out_ex0);

    out_ex0 = cat(1,xls_title3,cat(2,name_all,num2cell(cycle_all),num2cell(beta_all3)));
    xlswrite([ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,protein_add,int_add,D3_add,output_tail],out_ex0);

    out_ex0 = cat(1,xls_title0,cat(2,name_all,num2cell(cycle_all),num2cell(Nnu_all),num2cell(N_ratio),num2cell(Cm_all)));
    xlswrite([ana_folder,embryo_name{iem},'/',embryo_name{iem},nu_add,D3_add,output_tail],out_ex0);

    out_ex0 = cat(1,xls_title4,cat(2,En_name,num2cell(En_nu),num2cell(En_cycle),num2cell(En_cycle_fake)));
    xlswrite([ana_folder,embryo_name{iem},'/',embryo_name{iem},enrich_add,D3_add,output_tail],out_ex0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        

end




