clear all
close all

tic
xy_mismatch_check = true;
epsilonxy = 0.55;
epsilonxyz = 0.55;
epsilonz = 1;

%% Parameter setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
% out_folder = 'Results/';
out_folder = 'Results_protein2/';
% out_folder0 = 'Results_3Dz1/';
% out_folder0 = 'Results_decross/';
% out_folder0 = 'Results_original/';
% % out_folder0 = 'Results0/';
% hist_folder = 'Histogram/';
hist_folder = 'Histogram_alignment/';
fit_folder = 'Histogram_A/';
hist_folder2 = 'Histogram_en2/';
fit_folder2 = 'Histogram_protein2_A/';
% hist_folder_couple = 'Histogram_alignment/';
% hist_folder2_couple = 'Histogram_alignment_RNA2/';
fit_add = '_spot_fit';
hist_tail = '_raw.xls';
hist_tail2 = '_raw.xlsx';
hist_link_tail = '_link.xls';
N_thresh = 3;
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
enrich_add = '_enrich';
spot_add = '_spot';
hist_add = '_hist';

pro_ind = 14;
sub_pos = [6,7];
ana_folder = 'enrichspot_result/';
embryo_name = '02072016_16249-1-5M_FISHIF_hb_gal4_Cad_Bcd_PBcd_par2';
embryo_type = {'cbgal4'};
% embryo_name = '11012014_10652-1-3M_FISHIF_lacZ_gal4_Bcd_par2';
% embryo_type = {'b3n'};
% % embryo_name = '09302015_16249-1-5M_FISHIF_hb_gal4_Bcd_par2';
% % embryo_type = {'bgal4'};
var_list = {'nucleus_protein_profile','nucleus_protein_profile_ab','quanti_p','max_image','em_mask'};
dzmax = 1;
dz = -dzmax:dzmax;
inf0 = 1e8;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Data loading: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[num_list, folder_list] = xlsread(list_name);
folder_list = folder_list(strcmpi(embryo_type,folder_list(:,6)),:);
[N1,N2] = size(folder_list);
% run_list = {'[1]','[2]'};
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Index loading: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
list_ind = zeros(0);
em_index_all = 0;
for list_I = 1:N1
    [sub_num, sub_list] = xlsread([folder_list{list_I,1},in_folder,input_name]);
    [M1,M2] = size(sub_list);
    
    for list_J = 1:M1%eval(run_list{list_I})
        em_index_all = em_index_all+1;
        list_ind = cat(1,list_ind,[list_I,list_J]);
    end
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Output initialization: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xyz_all0_cell = cell(em_index_all,1);
dxy0_cell = cell(em_index_all,1);
Imin_all0_cell = cell(em_index_all,1);
I_RNA0_cell = cell(em_index_all,1);
I_en0_cell = cell(em_index_all,1);
N_en0_cell = cell(em_index_all,1);
EL_RNA0_cell = cell(em_index_all,1);
C_RNA0_cell = cell(em_index_all,1);
foci_list0_cell = cell(em_index_all,1);
Nu_en_cell = cell(em_index_all,1);

xyz_all1_cell = cell(em_index_all,1);
dxy1_cell = cell(em_index_all,1);
Imin_all1_cell = cell(em_index_all,1);
I_RNA1_cell = cell(em_index_all,1);
I_en1_cell = cell(em_index_all,1);
I_en2_cell = cell(em_index_all,1);
N_en1_cell = cell(em_index_all,1);
EL_RNA1_cell = cell(em_index_all,1);
C_RNA1_cell = cell(em_index_all,1);
foci_list1_cell = cell(em_index_all,1);
foci_circle_area_cell = cell(em_index_all,2);
Nu_RNA_cell = cell(em_index_all,1);

ind_em =zeros(em_index_all,1);
b0_em = zeros(em_index_all,1);
b2_em = zeros(em_index_all,1);
cycle_em = zeros(em_index_all,1);
em_name = cell(em_index_all,1);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Data analysis: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pool_name = parpool(21);
% matlabpool(9)
parfor em_index = 1:size(list_ind,1)
    list_I = list_ind(em_index,1);
    list_J = list_ind(em_index,2);
    [sub_num, sub_list] = xlsread([folder_list{list_I,1},in_folder,input_name]);
        
    if isempty(strfind(sub_list{list_J,3},'_60X'))
        [~,~,mismatch_matrix] = xlsread(mismatch_name);
        xymismatch_name0 = xymismatch_name;
    else
        [~,~,mismatch_matrix] = xlsread(mismatch_name2);
        xymismatch_name0 = xymismatch_name2;
    end

%         image_folder = [folder_list{list_I,1},in_folder,sub_list{list_J,3}];
        image_folder = [folder_list{list_I,1},in_folder,sub_list{list_J,3}(1:end-1),old_add,'/'];
        
        if isempty(flip0)
            flip_axis = any(cellfun(@(x) ~isempty(strfind(image_folder,x)),flip_list));   %%% Check whether the protein profile has an opposite direction (to Bcd)
        else
            flip_axis = flip0;
        end
        
        result_folder = [folder_list{list_I,1},out_folder];
        resolution = sub_num(list_J,9);
        resolutionz = sub_num(list_J,11);
        RNA_channel = sub_num(list_J,10);
        signal2_channel = sub_num(list_J,pro_ind);
        Nbin = sub_num(list_J,2);
        Mdim = sub_num(list_J,3);
%         load([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mat_tail]);
        embryo_var = cell(1,length(var_list));
        for ss = 1:length(var_list)
            embryo_var{ss} = var_load([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mat_tail],var_list{ss});
        end
        
        mask_stack = var_load([folder_list{list_I,1},mask_folder,sub_list{list_J,3},mask_name],'mask_stack');   %%% load 3D mask

        z_size = size(mask_stack,3);
        b = var_load([folder_list{list_I,1},fit_folder,sub_list{list_J,3}(1:end-1),fit_add,mat_tail],'b');   %%% load spot intensity fitting result
        Inten_thresh = b*N_thresh;   %%% set foci intensity threshold
        single_Inten = b;
%         b2 = var_load([folder_list{list_I,1},fit_folder2,sub_list{list_J,3}(1:end-1),fit_add,mat_tail],'b');   %%% load spot intensity fitting result
        b2 = embryo_var{3}(1)/resolution^2/resolutionz/6.02e8;
        
        Inten_thresh2 = b2*N_thresh;   %%% set foci intensity threshold
        single_Inten2 = b2;
%         N_cycle = round(log2(max(mask_stack(:))))+2;   %%% Calculate the nuclear cycle number
        N_cycle = sub_num(list_J,13);

        all_color = eval(folder_list{list_I,5});
        RNA_color = all_color{RNA_channel};
        signal2_color = all_color{signal2_channel};
        RNA_signal2_mismatch = mismatch_matrix{strcmp(RNA_color,mismatch_matrix(1,:)),strcmp(signal2_color,mismatch_matrix(:,1))};
        if xy_mismatch_check
            RNA_signal2_xymismatch = cell(1,7);
            RNA_signal2_xymismatch{1} = var_load(xymismatch_name0,[RNA_color,'_',signal2_color]);
            RNA_signal2_xymismatch{2} = var_load(xymismatch_name0,[RNA_color,'_',signal2_color,'_con']);
            RNA_signal2_xymismatch{3} = var_load(xymismatch_name0,[RNA_color,'_',signal2_color,'_x0']);
            resolution_mismatch = var_load(xymismatch_name0,'resolution_mismatch');
            RNA_signal2_xymismatch(4:end) = {Nbin,Mdim,resolution,resolution_mismatch};
        else
            RNA_signal2_xymismatch = [];
        end
        signal2_RNA_mismatch = mismatch_matrix{strcmp(signal2_color,mismatch_matrix(:,1)),strcmp(RNA_color,mismatch_matrix(1,:))};
        if xy_mismatch_check
            signal2_RNA_xymismatch = cell(1,7);
            signal2_RNA_xymismatch{1} = var_load(xymismatch_name0,[signal2_color,'_',RNA_color]);
            signal2_RNA_xymismatch{2} = var_load(xymismatch_name0,[signal2_color,'_',RNA_color,'_con']);
            signal2_RNA_xymismatch{3} = var_load(xymismatch_name0,[signal2_color,'_',RNA_color,'_x0']);
            signal2_RNA_xymismatch(4:end) = {Nbin,Mdim,resolution,resolution_mismatch};
        else
            signal2_RNA_xymismatch = [];
        end
        
        if exist([folder_list{list_I,1},hist_folder,sub_list{list_J,3}(1:(end-1)),hist_tail],'file')
            [foci_list,~,~] = xlsread([folder_list{list_I,1},hist_folder,sub_list{list_J,3}(1:(end-1)),hist_tail]);
        else
            [foci_list,~,~] = xlsread([folder_list{list_I,1},hist_folder,sub_list{list_J,3}(1:(end-1)),hist_tail2]);
        end
        if exist([folder_list{list_I,1},hist_folder2,sub_list{list_J,3}(1:(end-1)),hist_tail],'file')
            [foci_list2,~,~] = xlsread([folder_list{list_I,1},hist_folder2,sub_list{list_J,3}(1:(end-1)),hist_tail]);
        else
            [foci_list2,~,~] = xlsread([folder_list{list_I,1},hist_folder2,sub_list{list_J,3}(1:(end-1)),hist_tail2]);
        end
        
        foci_listb = foci_list;
        foci_list2b = foci_list2;
        foci_list2(:,6:8) = position_adjust(foci_list2(:,6:8),embryo_var{4},RNA_signal2_mismatch,RNA_signal2_xymismatch);
        foci_listb(:,6:8) = position_adjust(foci_list(:,6:8),embryo_var{4},signal2_RNA_mismatch,signal2_RNA_xymismatch);
        
        %%% Determine which nuclei the foci belong to:
        mask_expand = dilate3D(mask_stack,8);
% % %         mask_expand = dilate3D(mask_stack,20);
% % %         mask_expand = mask_stack;
        Nuclei_foci = mask_expand(sub2ind(size(mask_stack),round(foci_list(:,6)),round(foci_list(:,7)),round(foci_list(:,8))));
        foci_list = foci_list(Nuclei_foci > 0,:);
        Nuclei_foci = Nuclei_foci(Nuclei_foci > 0);
        
        x00 = min(max(round(foci_list2(:,6)),1),size(mask_stack,1)); Nx00 = round(foci_list2(:,6)) >= 1 & round(foci_list2(:,6)) <= size(mask_stack,1);
        y00 = min(max(round(foci_list2(:,7)),1),size(mask_stack,2)); Ny00 = round(foci_list2(:,7)) >= 1 & round(foci_list2(:,7)) <= size(mask_stack,2);
        z00 = min(max(round(foci_list2(:,8)),1),size(mask_stack,3)); Nz00 = round(foci_list2(:,8)) >= 1 & round(foci_list2(:,8)) <= size(mask_stack,3);
        Nuclei_foci2 = mask_expand(sub2ind(size(mask_stack),x00,y00,z00));
        N00 = Nuclei_foci2 > 0 & Nx00 & Ny00 & Nz00;
        foci_list2 = foci_list2(N00,:);
        Nuclei_foci2 = Nuclei_foci2(N00);
%         C_foci = C_foci(N00);
        
        C_nuclei = embryo_var{1}(:,2)/embryo_var{3}(1);
        C_foci = foci_list2(:,5)/embryo_var{3}(1);
        
        I_foci = prod(foci_list(:,1:3),2)*2*pi/single_Inten;
        I_foci2 = prod(foci_list2(:,1:3),2)*2*pi/single_Inten2;
        
        foci_list0_cell(em_index,:) = {foci_list2};
        foci_list1_cell(em_index,:) = {foci_list};
        
        %%% Calculate the EL coordinate for each foci:
        EL_info = get_EL(embryo_var{5});   %%% get EL information (extreme points, EL length) from the embryo mask
        x0 = EL_info(1);
        y0 = EL_info(2);
        x1 = EL_info(3);
        y1 = EL_info(4);
        L2_extreme = EL_info(5);
        foci_xy = foci_list(:,[7,6]);
        foci_EL = 1-dot((foci_xy-repmat([x0,y0],size(foci_xy,1),1)),repmat([(x1-x0),(y1-y0)],size(foci_xy,1),1),2)/L2_extreme;
        if xor(mean(I_foci(foci_EL >= 0.5)) > mean(I_foci(foci_EL <= 0.5)),flip_axis)
            foci_EL = 1-foci_EL;
        end
        
        %%% Calculate the area of foci centered circles:
        dx = 0.02;
        xl0 = 0:dx:(3+dx);
        xl = xl0/resolution;
        foci_circle_area = zeros(size(foci_list,1),length(xl),length(dz));

        warning('off')
        for jj = 1:size(foci_list,1)
            for zz = 1:length(dz)
                if (round(foci_list(jj,8))+dz(zz) > 0) && (round(foci_list(jj,8))+dz(zz) <= size(mask_expand,3))
                    nuclei0 = mask_expand(:,:,round(foci_list(jj,8))+dz(zz)) == Nuclei_foci(jj);
                    if nnz(nuclei0)
                        nuclei_bound0 = bwboundaries(nuclei0);
                        nuclei_bound = nuclei_bound0{1};

                        foci_same_nuclei = (Nuclei_foci == Nuclei_foci(jj)) & (foci_list(:,8)-dzmax <= foci_list(jj,8)+dz(zz)) & (foci_list(:,8)+dzmax >= foci_list(jj,8)+dz(zz));
                        foci_same_nuclei(jj) = false;
                        id_foci = find(foci_same_nuclei);

                        for kk = id_foci'
                            xf = foci_list([jj,kk],6);
                            yf = foci_list([jj,kk],7);
                            xp = [(xf(2)+xf(1))/2-(yf(2)-yf(1))/2;(xf(2)+xf(1))/2+(yf(2)-yf(1))/2]; 
                            yp = [(yf(2)+yf(1))/2+(xf(2)-xf(1))/2;(yf(2)+yf(1))/2-(xf(2)-xf(1))/2];
                            if xf(1) < xf(2)
                                s = 'R';
                            elseif xf(1) > xf(2)
                                s = 'L';
                            elseif yf(1) < yf(2)
                                s = 'T';
                            else
                                s = 'B';
                            end
                            nuclei_bound = cutpolygon(nuclei_bound,[xp,yp],s);
                        end
                        if ~isempty(nuclei_bound)
                            Pnb = polyshape(nuclei_bound);
                            for kk = 1:length(xl)
                                if xl(kk) > 0
                                    Pcircle = nsidedpoly(128,'Center',foci_list(jj,6:7),'Radius',xl(kk));
                                    foci_circle_area(jj,kk,zz) = area(intersect(Pnb,Pcircle))*resolution^2;
                                end
                            end
                        end
                    end
                end
            end
        end
        warning('on')
        foci_circle_area_cell(em_index,:) = {foci_circle_area,xl0};

%%% foci coupling:
        xyz1 = foci_list(:,6:8).*repmat([resolution,resolution,resolutionz],size(foci_list,1),1);
        xyz2 = foci_list2(:,6:8).*repmat([resolution,resolution,resolutionz],size(foci_list2,1),1);
        xyz_all0_cell(em_index,:) = {xyz2};
        xyz_all1_cell(em_index,:) = {xyz1};
        
%         dxyz = pdist2(xyz1,xyz2);
        dxy0 = pdist2(xyz1(:,1:2),xyz2(:,1:2));
        dxy1 = (pdist2(Nuclei_foci,Nuclei_foci2) > 0)*inf0;
        dxy2 = (pdist2(xyz1(:,3),xyz2(:,3)) > dzmax)*inf0;
        dxy = dxy0+dxy1+dxy2;
        
        [dxymin,I1min] = min(dxy);
        [dxymin1,I1min1] = min(dxy,[],2);
        Imin_all0_cell(em_index,:) = {I1min'};
        Imin_all1_cell(em_index,:) = {I1min1};
        
        N_en00 = hist(I1min.*(dxymin <= 0.35),0:size(foci_list,1));
        N_en00 = N_en00(2:end);
        
%         dxyz0 = cat(1,dxyz0,dxyz(:));
        dxy0_cell(em_index,:) = {dxymin'};
        I_RNA0_cell(em_index,:) = {I_foci(I1min')};
        I_en0_cell(em_index,:) = {I_foci2};
        N_en0_cell(em_index,:) = {N_en00(I1min)'};
        EL_RNA0_cell(em_index,:) = {foci_EL(I1min')};
        C_RNA0_cell(em_index,:) = {[embryo_var{2}(Nuclei_foci(I1min'),2),C_nuclei(Nuclei_foci(I1min')),C_foci]};
        Nu_en_cell(em_index,:) = {Nuclei_foci2};
 
        dxy1_cell(em_index,:) = {dxymin1};
        I_RNA1_cell(em_index,:) = {I_foci};
        I_en1_cell(em_index,:) = {I_foci2(I1min1)};
        N_en1_cell(em_index,:) = {N_en00'};
        EL_RNA1_cell(em_index,:) = {foci_EL};
        C_RNA1_cell(em_index,:) = {[embryo_var{2}(Nuclei_foci,2),C_nuclei(Nuclei_foci),C_foci(I1min1)]};
        Nu_RNA_cell(em_index,:) = {Nuclei_foci};
        
        ind_em(em_index,1) = em_index;
        cycle_em(em_index,1) = N_cycle;
        b0_em(em_index,1) = b;
        b2_em(em_index,1) = b2;
        em_name(em_index,1) = {[result_folder,sub_list{list_J,3}]};
%         disp(['Embryo ',num2str(em_index),': t = ',num2str(toc)])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
end
delete(pool_name)
% matlabpool close
toc


%% Output: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if exist([ana_folder,embryo_name,'/']) ~= 7
    mkdir([ana_folder,embryo_name,'/']);
end
save([ana_folder,embryo_name,'/',embryo_name,mat_tail],'xyz_all0_cell','dxy0_cell','Imin_all0_cell','I_RNA0_cell','I_en0_cell','N_en0_cell','EL_RNA0_cell','C_RNA0_cell','xyz_all1_cell','dxy1_cell','Imin_all1_cell','I_RNA1_cell','I_en1_cell','I_en2_cell','N_en1_cell','EL_RNA1_cell','C_RNA1_cell','ind_em','b0_em','b2_em','cycle_em','em_name','foci_list0_cell','foci_list1_cell','foci_circle_area_cell','Nu_en_cell','Nu_RNA_cell','dz');

%%% Plots:        
figure(1);maximize(1)
figure(2);maximize(2)
for em_index = 1:size(list_ind,1)
    figure(1)
    subplot(sub_pos(1),sub_pos(2),em_index)
        plot(xyz_all0_cell{em_index}(:,1)-xyz_all1_cell{em_index}(Imin_all0_cell{em_index}',1),xyz_all0_cell{em_index}(:,2)-xyz_all1_cell{em_index}(Imin_all0_cell{em_index}',2),'k.')
        axis equal
        xlim([-1.5,1.5])
        ylim([-1.5,1.5])
        xlabel('x (um)')
        ylabel('y (um)')
        title(['Embryo ',num2str(em_index),', Cycle ',num2str(cycle_em(em_index))])

    figure(2)
    subplot(sub_pos(1),sub_pos(2),em_index)
        theta = angle(xyz_all0_cell{em_index}(:,1)-xyz_all1_cell{em_index}(Imin_all0_cell{em_index}',1)+1i*(xyz_all0_cell{em_index}(:,2)-xyz_all1_cell{em_index}(Imin_all0_cell{em_index}',2)));
        polarhistogram(theta,50)
        title(['Embryo ',num2str(em_index),', Cycle ',num2str(cycle_em(em_index))])
end
saveas(1,[ana_folder,embryo_name,'/',embryo_name,enrich_add,spot_add,D3_add,figure_tail]);
saveas(1,[ana_folder,embryo_name,'/',embryo_name,enrich_add,spot_add,D3_add,figure_tail2]);
saveas(2,[ana_folder,embryo_name,'/',embryo_name,enrich_add,spot_add,hist_add,D3_add,figure_tail]);
saveas(2,[ana_folder,embryo_name,'/',embryo_name,enrich_add,spot_add,hist_add,D3_add,figure_tail2]);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





