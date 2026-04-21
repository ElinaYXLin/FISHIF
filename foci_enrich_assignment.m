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
out_folder = 'Results/';
% out_folder0 = 'Results_3Dz1/';
% out_folder0 = 'Results_decross/';
% out_folder0 = 'Results_original/';
out_folder0 = 'Results0/';
% hist_folder = 'Histogram/';
hist_folder = 'Histogram_alignment/';
fit_folder = 'Histogram_A/';
hist_folder2 = 'Histogram_en/';
fit_folder2 = 'Histogram_protein_A/';
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

sub_pos = [5,7];
ana_folder = 'enrichspot_result/';
% % embryo_name = '11012014_10652-1-3M_FISHIF_lacZ_gal4_Bcd';
% % embryo_type = {'b3n'};
embryo_name = '09302015_16249-1-5M_FISHIF_hb_gal4_Bcd';
embryo_type = {'bgal4'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Output initialization: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xyz_all0_cell = cell(0);
dxy0_cell = cell(0);
Imin_all0_cell = cell(0);
I_RNA0_cell = cell(0);
I_en0_cell = cell(0);
N_en0_cell = cell(0);
EL_RNA0_cell = cell(0);
C_RNA0_cell = cell(0);
foci_list0_cell = cell(0);

xyz_all1_cell = cell(0);
dxy1_cell = cell(0);
Imin_all1_cell = cell(0);
I_RNA1_cell = cell(0);
I_en1_cell = cell(0);
I_en2_cell = cell(0);
N_en1_cell = cell(0);
EL_RNA1_cell = cell(0);
C_RNA1_cell = cell(0);
foci_list1_cell = cell(0);
foci_circle_area_cell = cell(0);

em_index = 0;
ind_em =[];
b0_em = [];
b2_em = [];
cycle_em = [];
em_name = cell(0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Data loading: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[num_list, folder_list] = xlsread(list_name);
folder_list = folder_list(strcmpi(embryo_type,folder_list(:,6)),:);
[N1,N2] = size(folder_list);
% run_list = {'[1]','[2]'};
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1);maximize(1)
figure(2);maximize(2)
warning('off')
%% Data analysis: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for list_I = 1:N1
    [sub_num, sub_list] = xlsread([folder_list{list_I,1},in_folder,input_name]);
    [M1,M2] = size(sub_list);
% %     if ~isempty(out_folder0)
% %         copyfile([folder_list{list_I,1},out_folder],[folder_list{list_I,1},out_folder0]);
% %     end
    
    for list_J = 1:M1%eval(run_list{list_I})
        if isempty(strfind(sub_list{list_J,3},'_60X'))
            [~,~,mismatch_matrix] = xlsread(mismatch_name);
            load(xymismatch_name);
        else
            [~,~,mismatch_matrix] = xlsread(mismatch_name2);
            load(xymismatch_name2);
        end
        
        em_index = em_index+1;
        
%         image_folder = [folder_list{list_I,1},in_folder,sub_list{list_J,3}];
        image_folder = [folder_list{list_I,1},in_folder,sub_list{list_J,3}(1:end-1),old_add,'/'];
        
        if isempty(flip0)
            flip_axis = any(cellfun(@(x) ~isempty(strfind(image_folder,x)),flip_list));   %%% Check whether the protein profile has an opposite direction (to Bcd)
        else
            flip_axis = flip0;
        end
        
        result_folder = [folder_list{list_I,1},out_folder];
        resolutionz = sub_num(list_J,11);
        signal2_channel = sub_num(list_J,8);
        Nbin = sub_num(list_J,2);
        Mdim = sub_num(list_J,3);
        load([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mat_tail]);
        image_folder = [folder_list{list_I,1},in_folder,sub_list{list_J,3}(1:end-1),old_add,'/'];

        load([folder_list{list_I,1},mask_folder,sub_list{list_J,3},mask_name]);   %%% load 3D mask
        z_size = size(mask_stack,3);
        load([folder_list{list_I,1},fit_folder,sub_list{list_J,3}(1:end-1),fit_add,mat_tail]);   %%% load spot intensity fitting result
        Inten_thresh = b*N_thresh;   %%% set foci intensity threshold
        single_Inten = b;
%         load([folder_list{list_I,1},fit_folder2,sub_list{list_J,3}(1:end-1),fit_add,mat_tail]);   %%% load spot intensity fitting result
        b2 = quanti_p(1)/resolution^2/resolutionz/6.02e8;
        
        Inten_thresh2 = b2*N_thresh;   %%% set foci intensity threshold
        single_Inten2 = b2;
%         N_cycle = round(log2(max(mask_stack(:))))+2;   %%% Calculate the nuclear cycle number
        N_cycle = sub_num(list_J,13);

        all_color = eval(folder_list{list_I,5});
        RNA_color = all_color{RNA_channel};
        signal2_color = all_color{signal2_channel};
        RNA_signal2_mismatch = mismatch_matrix{strcmp(RNA_color,mismatch_matrix(1,:)),strcmp(signal2_color,mismatch_matrix(:,1))};
        if xy_mismatch_check
            RNA_signal2_xymismatch = {eval([RNA_color,'_',signal2_color]),eval([RNA_color,'_',signal2_color,'_con']),eval([RNA_color,'_',signal2_color,'_x0']),Nbin,Mdim,resolution,resolution_mismatch};
        else
            RNA_signal2_xymismatch = [];
        end
        signal2_RNA_mismatch = mismatch_matrix{strcmp(signal2_color,mismatch_matrix(:,1)),strcmp(RNA_color,mismatch_matrix(1,:))};
        if xy_mismatch_check
            signal2_RNA_xymismatch = {eval([signal2_color,'_',RNA_color]),eval([signal2_color,'_',RNA_color,'_con']),eval([signal2_color,'_',RNA_color,'_x0']),Nbin,Mdim,resolution,resolution_mismatch};
        else
            signal2_RNA_xymismatch = [];
        end
        
        if exist([folder_list{list_I,1},hist_folder,sub_list{list_J,3}(1:(end-1)),hist_tail])
            [foci_list,~,~] = xlsread([folder_list{list_I,1},hist_folder,sub_list{list_J,3}(1:(end-1)),hist_tail]);
        else
            [foci_list,~,~] = xlsread([folder_list{list_I,1},hist_folder,sub_list{list_J,3}(1:(end-1)),hist_tail2]);
        end
        if exist([folder_list{list_I,1},hist_folder2,sub_list{list_J,3}(1:(end-1)),hist_tail])
            [foci_list2,~,~] = xlsread([folder_list{list_I,1},hist_folder2,sub_list{list_J,3}(1:(end-1)),hist_tail]);
        else
            [foci_list2,~,~] = xlsread([folder_list{list_I,1},hist_folder2,sub_list{list_J,3}(1:(end-1)),hist_tail2]);
        end
        
        C_nuclei = nucleus_protein_profile(:,2)/quanti_p(1);
        C_foci = foci_list2(:,5)/quanti_p(1);
        
        %%% Determine which nuclei the foci belong to:
% % %         mask_expand = dilate3D(mask_stack,20);
        mask_expand = mask_stack;
        Nuclei_foci = mask_expand(sub2ind(size(mask_stack),round(foci_list(:,6)),round(foci_list(:,7)),round(foci_list(:,8))));
        foci_list = foci_list(Nuclei_foci > 0,:);
        Nuclei_foci = Nuclei_foci(Nuclei_foci > 0);
        
        
        foci_listb = foci_list;
        foci_list2b = foci_list2;
        foci_list2(:,6:8) = position_adjust(foci_list2(:,6:8),max_image,RNA_signal2_mismatch,RNA_signal2_xymismatch);
        foci_listb(:,6:8) = position_adjust(foci_list(:,6:8),max_image,signal2_RNA_mismatch,signal2_RNA_xymismatch);
        I_foci = prod(foci_list(:,1:3),2)*2*pi/single_Inten;
        I_foci2 = prod(foci_list2(:,1:3),2)*2*pi/single_Inten2;
        
        foci_list0_cell = cat(1,foci_list0_cell,{foci_list2});
        foci_list1_cell = cat(1,foci_list1_cell,{foci_list});
        
        %%% Calculate the EL coordinate for each foci:
        EL_info = get_EL(em_mask);   %%% get EL information (extreme points, EL length) from the embryo mask
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
        xl = (0:dx:(3+dx))/resolution;
        foci_circle_area = zeros(size(foci_list,1),length(xl));
        
        for jj = 1:size(foci_list,1)
            nuclei0 = mask_stack(:,:,round(foci_list(jj,8))) == Nuclei_foci(jj);
            nuclei_bound0 = bwboundaries(nuclei0);
            nuclei_bound = nuclei_bound0{1};
            
            foci_same_nuclei = Nuclei_foci == Nuclei_foci(jj);
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
                        foci_circle_area(jj,kk) = area(intersect(Pnb,Pcircle))*resolution^2;
                    end
                end
            end
        end
        foci_circle_area_cell = cat(1,foci_circle_area_cell,{foci_circle_area});

%%% foci coupling:
        xyz1 = foci_list(:,6:8).*repmat([resolution,resolution,resolutionz],size(foci_list,1),1);
        xyz2 = foci_list2(:,6:8).*repmat([resolution,resolution,resolutionz],size(foci_list2,1),1);
        xyz_all0_cell = cat(1,xyz_all0_cell,{xyz2});
        xyz_all1_cell = cat(1,xyz_all1_cell,{xyz1});
        
%         dxyz = pdist2(xyz1,xyz2);
        dxy = pdist2(xyz1(:,1:2),xyz2(:,1:2));
        
        [dxymin,I1min] = min(dxy);
        [dxymin1,I1min1] = min(dxy,[],2);
        Imin_all0_cell = cat(1,Imin_all0_cell,{I1min'});
        Imin_all1_cell = cat(1,Imin_all1_cell,{I1min1});
        
        N_en00 = hist(I1min.*(dxymin <= 0.35),0:size(foci_list,1));
        N_en00 = N_en00(2:end);
        
%         dxyz0 = cat(1,dxyz0,dxyz(:));
        dxy0_cell = cat(1,dxy0_cell,{dxymin'});
        I_RNA0_cell = cat(1,I_RNA0_cell,{I_foci(I1min')});
        I_en0_cell = cat(1,I_en0_cell,{I_foci2});
        N_en0_cell = cat(1,N_en0_cell,{N_en00(I1min)'});
        EL_RNA0_cell = cat(1,EL_RNA0_cell,{foci_EL(I1min')});
        C_RNA0_cell = cat(1,C_RNA0_cell,{[nucleus_protein_profile_ab(Nuclei_foci(I1min'),2),C_nuclei(Nuclei_foci(I1min')),C_foci]});
 
        dxy1_cell = cat(1,dxy1_cell,{dxymin1});
        I_RNA1_cell = cat(1,I_RNA1_cell,{I_foci});
        I_en1_cell = cat(1,I_en1_cell,{I_foci2(I1min1)});
        N_en1_cell = cat(1,N_en1_cell,{N_en00'});
        EL_RNA1_cell = cat(1,EL_RNA1_cell,{foci_EL});
        C_RNA1_cell = cat(1,C_RNA1_cell,{[nucleus_protein_profile_ab(Nuclei_foci,2),C_nuclei(Nuclei_foci),C_foci(I1min1)]});
        
        ind_em = cat(1,ind_em,em_index);
        cycle_em = cat(1,cycle_em,N_cycle);
        b0_em = cat(1,b0_em,b);
        b2_em = cat(1,b2_em,b2);
        em_name = cat(1,em_name,{[result_folder,sub_list{list_J,3}]});
        
%%% Plots:        
        figure(1)
        subplot(sub_pos(1),sub_pos(2),em_index)
            plot(xyz2(:,1)-xyz1(I1min,1),xyz2(:,2)-xyz1(I1min,2),'k.')
            axis equal
            xlim([-1.5,1.5])
            ylim([-1.5,1.5])
            xlabel('x (um)')
            ylabel('y (um)')
            title(['Embryo ',num2str(em_index),', Cycle ',num2str(N_cycle)])

        figure(2)
        subplot(sub_pos(1),sub_pos(2),em_index)
            theta = angle(xyz2(:,1)-xyz1(I1min,1)+1i*(xyz2(:,2)-xyz1(I1min,2)));
            polarhistogram(theta,50)
            title(['Embryo ',num2str(em_index),', Cycle ',num2str(N_cycle)])
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        
    end

end
toc


%% Output: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if exist([ana_folder,embryo_name,'/']) ~= 7
    mkdir([ana_folder,embryo_name,'/']);
end
saveas(1,[ana_folder,embryo_name,'/',embryo_name,enrich_add,spot_add,D3_add,figure_tail]);
saveas(1,[ana_folder,embryo_name,'/',embryo_name,enrich_add,spot_add,D3_add,figure_tail2]);
saveas(2,[ana_folder,embryo_name,'/',embryo_name,enrich_add,spot_add,hist_add,D3_add,figure_tail]);
saveas(2,[ana_folder,embryo_name,'/',embryo_name,enrich_add,spot_add,hist_add,D3_add,figure_tail2]);
save([ana_folder,embryo_name,'/',embryo_name,mat_tail],'xyz_all0_cell','dxy0_cell','Imin_all0_cell','I_RNA0_cell','I_en0_cell','N_en0_cell','EL_RNA0_cell','C_RNA0_cell','xyz_all1_cell','dxy1_cell','Imin_all1_cell','I_RNA1_cell','I_en1_cell','I_en2_cell','N_en1_cell','EL_RNA1_cell','C_RNA1_cell','ind_em','b0_em','b2_em','cycle_em','em_name','foci_list0_cell','foci_list1_cell','foci_circle_area_cell','xl');
warning('on')
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
% % % % %%
% % % % I_plot = [2:26,29:30,32];%[1:2,4:7,10:14,16,20];
% % % % I_cycle = 12;
% % % % cycle_all0 = cycle_em(index_all0);
% % % % cycle_all1 = cycle_em(index_all1);
% % % % fit_initial = [0.1,0.5,0.5,0.1,1,1,0.1,2,2];
% % % % fit_lower = [0,0,0.2,0,0,0.2,0,0,0.2];
% % % % fit_upper = [1,1,1,1,2,2,1,3,3];
% % % % r00  = zeros(em_index,1);
% % % % % % r00 = [0.3149,0.3953,0.4212,0.4021,0.2401,0.3934,0.3607,0.3054,0.3183,0.274,0.3*ones(1,10)]'/0.3;
% % % % 
% % % % %%% Fitting for rescaling
% % % % % % % figure
% % % % % % % for ii = 1:em_index
% % % % % % %     dx = 0.1;
% % % % % % %     xl = 0:dx:6;
% % % % % % %     xr = xl+dx;
% % % % % % %     x0 = xl+dx/2;
% % % % % % %     [n0,r0] = hist(I_en0(ismember(index_all0,ii) ),x0);
% % % % % % %     fitobject = fit(r0(1:end-1)',n0(1:end-1)'/sum(n0),'gauss3','StartPoint',fit_initial,'Lower',fit_lower,'Upper',fit_upper);
% % % % % % %     r00(ii) = min([fitobject.b1,fitobject.b2,fitobject.b3]);
% % % % % % %     
% % % % % % %     subplot(5,7,ii)
% % % % % % %     r0b = xl(1):0.01:xl(end);
% % % % % % %     plot(r0,n0/sum(n0),'ko')
% % % % % % %     hold on
% % % % % % %     plot(r0b,feval(fitobject,r0b),'r')
% % % % % % %     xlim([0,xl(end)-1])
% % % % % % %     title(['Embryo ',num2str(ii)])
% % % % % % % end
% % % % 
% % % % r00(:) = 0.13;
% % % % r11 = r00(index_all0);%/mean(r00);
% % % % % r11 = b2_em(index_all0)./b_em(index_all0);
% % % % 
% % % %     figure 
% % % %     set(gcf,'Name',result_folder)
% % % %     subplot(2,2,1)
% % % %         dx = 0.02;
% % % %         xl = 0:dx:3;
% % % %         xr = xl+dx;
% % % %         x0 = xl+dx/2;
% % % %         [n0,r0] = hist(dxy0,x0);
% % % %         nr = n0./(xr.^2-xl.^2)/pi;
% % % % %         [na,r0] = hist(dxy0(ismember(index_all0,I_plot) & EL_RNA0 < 0.4),x0);
% % % %         [na,r0] = hist(dxy0(ismember(index_all0,I_plot) & C_RNA0 > 10e-9),x0);
% % % %         nar = na./(xr.^2-xl.^2)/pi;
% % % % %         [np,r0] = hist(dxy0(ismember(index_all0,I_plot) & EL_RNA0 >= 0.7),x0);
% % % %         [np,r0] = hist(dxy0(ismember(index_all0,I_plot) & C_RNA0 < 5e-9),x0);
% % % %         npr = np./(xr.^2-xl.^2)/pi;
% % % % % %         plot(r0,nr,'k','DisplayName','All');
% % % % % %         hold on
% % % %         plot(r0,nar/sum(na),'r','DisplayName','Anterior')
% % % %         hold on
% % % %         plot(r0,npr/sum(np),'b','DisplayName','Posterior');
% % % %         xlim([0,xl(end)-1])
% % % %         xlabel('r (um)')
% % % %         ylabel('%/um^2')
% % % %         title('Distance distribution (EL)','Interpreter','none')
% % % %         legend('show')
% % % %     subplot(2,2,3)
% % % %         [nH,r0] = hist(dxy0(ismember(index_all0,I_plot) & I_RNA0 > 7),x0);
% % % %         nHr = nH./(xr.^2-xl.^2)/pi;
% % % %         [nL,r0] = hist(dxy0(ismember(index_all0,I_plot) & I_RNA0 <= 2),x0);
% % % %         nLr = nL./(xr.^2-xl.^2)/pi;
% % % % % %         plot(r0,nr,'k','DisplayName','All');
% % % % % %         hold on
% % % %         plot(r0,nHr/sum(nH),'r','DisplayName','Active')
% % % %         hold on
% % % %         plot(r0,nLr/sum(nL),'b','DisplayName','Inactive');
% % % %         xlim([0,xl(end)-1])
% % % %         xlabel('r (um)')
% % % %         ylabel('%/um^2')
% % % %         title('Distance distribution (TX)','Interpreter','none')
% % % %         legend('show')
% % % %     subplot(2,2,2)
% % % %         dx = 0.5;
% % % %         dw = 0.75;
% % % %         xl = 0:dx:31;
% % % %         xr = xl+dw;
% % % %         x0 = (xl+xr)/2;
% % % %         r0 = x0;
% % % % %         I000 = ismember(index_all0,I_plot) & EL_RNA0 <= 0.45 & EL_RNA0 >= 0.15;
% % % %         I000 = N_en0 == 2 & ismember(cycle_all0,I_cycle) & ismember(index_all0,I_plot) & EL_RNA0 <= 0.7 & EL_RNA0 >= 0.25 & C_RNA0 > 15e-9 & C_RNA0 <= 20e-9;
% % % % % % %         [n0,r0] = hist(I_en0,x0);
% % % % % % %         [na,r0] = hist(I_en0(I000 & dxy0 <= 0.3 & dxy0 > 0.05)./r11(I000 & dxy0 <= 0.3 & dxy0 > 0.05),x0);
% % % % % % %         [np,r0] = hist(I_en0(I000 & dxy0 > 0.3 & dxy0 <= 0.6)./r11(I000 & dxy0 > 0.3 & dxy0 <= 0.6),x0);
% % % %         n0 = hist0(I_en0,xl,xr);
% % % %         na = hist0(I_en0(I000 & dxy0 <= 0.35 & dxy0 >= 0)./r11(I000 & dxy0 <= 0.35 & dxy0 >= 0),xl,xr);
% % % %         np = hist0(I_en0(I000 & dxy0 > 0.35 & dxy0 <= 0.6)./r11(I000 & dxy0 > 0.35 & dxy0 <= 0.6),xl,xr);
% % % % % %         plot(r0,n0,'k','DisplayName','All');
% % % % % %         hold on
% % % %         plot(r0,na,'r','DisplayName','Anterior inner')
% % % %         hold on
% % % %         plot(r0,np,'b','DisplayName','Anterior outer');
% % % %         plot(r0,na-np*(0.35^2-0^2)/(0.6^2-0.35^2),'k','DisplayName','Anterior pure');
% % % %         xlim([0,xl(end)-1])
% % % %         xlabel('En(#)')
% % % %         ylabel('#')
% % % %         title('En distribution (Anterior)','Interpreter','none')
% % % %         legend('show')
% % % %     subplot(2,2,4)
% % % % %         I000 = ismember(index_all0,I_plot) & EL_RNA0 > 0.7;
% % % %         I000 = N_en0 == 2 & ismember(cycle_all0,I_cycle) & ismember(index_all0,I_plot) & EL_RNA0 > 0.7 & C_RNA0 < 5e-9 & C_RNA0 >= 1e-9;
% % % % % % %         [na,r0] = hist(I_en0(I000 & dxy0 <= 0.3 & dxy0 > 0.05)./r11(I000 & dxy0 <= 0.3 & dxy0 > 0.05),x0);
% % % % % % %         [np,r0] = hist(I_en0(I000 & dxy0 > 0.3 & dxy0 <= 0.6)./r11(I000 & dxy0 > 0.3 & dxy0 <= 0.6),x0);
% % % %         na = hist0(I_en0(I000 & dxy0 <= 0.35 & dxy0 >= 0)./r11(I000 & dxy0 <= 0.35 & dxy0 >= 0),xl,xr);
% % % %         np = hist0(I_en0(I000 & dxy0 > 0.35 & dxy0 <= 0.6)./r11(I000 & dxy0 > 0.35 & dxy0 <= 0.6),xl,xr);
% % % % % %         plot(r0,n0,'k','DisplayName','All');
% % % % % %         hold on
% % % %         plot(r0,na,'r','DisplayName','Posterior inner')
% % % %         hold on
% % % %         plot(r0,np,'b','DisplayName','Posterior outer');
% % % %         plot(r0,na-np*(0.35^2-0^2)/(0.6^2-0.35^2),'k','DisplayName','Posterior pure');
% % % %         xlim([0,xl(end)-1])
% % % %         xlabel('En(#)')
% % % %         ylabel('#')
% % % %         title('En distribution (Posterior)','Interpreter','none')
% % % %         legend('show')
% % % %         
% % % % % %     subplot(2,2,4)
% % % % % %         [nH,r0] = hist(I_en0(I_RNA0 > 7 & dxy0 <= 0.3),x0);
% % % % % %         [nL,r0] = hist(I_en0(I_RNA0 <= 2 & dxy0 <= 0.3),x0);
% % % % % % % %         plot(r0,n0,'k','DisplayName','All');
% % % % % % % %         hold on
% % % % % %         plot(r0,nH/sum(nH),'r','DisplayName','Active')
% % % % % %         hold on
% % % % % %         plot(r0,nL/sum(nL),'b','DisplayName','Inactive');
% % % % % %         xlim([0,xl(end)-1])
% % % % % %         xlabel('En(#)')
% % % % % %         ylabel('%')
% % % % % %         title('En distribution (TX)','Interpreter','none')
% % % % % %         legend('show')
% % % %     
% % % %         
% % % %     figure 
% % % %     set(gcf,'Name',result_folder)
% % % %     subplot(2,3,1)
% % % %         dx = 0.02;
% % % %         xl = 0:dx:3;
% % % %         xr = xl+dx;
% % % %         x0 = xl+dx/2;
% % % %         [n0,r0] = hist(dxy1,x0);
% % % %         nr = n0./(xr.^2-xl.^2)/pi;
% % % %         [na,r0] = hist(dxy1(EL_RNA1 < 0.4),x0);
% % % %         nar = na./(xr.^2-xl.^2)/pi;
% % % %         [np,r0] = hist(dxy1(EL_RNA1 >= 0.7),x0);
% % % %         npr = np./(xr.^2-xl.^2)/pi;
% % % % % %         plot(r0,nr,'k','DisplayName','All');
% % % % % %         hold on
% % % %         plot(r0,nar/sum(na),'r','DisplayName','Anterior')
% % % %         hold on
% % % %         plot(r0,npr/sum(np),'b','DisplayName','Posterior');
% % % %         xlim([0,xl(end)-1])
% % % %         xlabel('r (um)')
% % % %         ylabel('%/um^2')
% % % %         title('Distance distribution (EL)','Interpreter','none')
% % % %         legend('show')
% % % %     subplot(2,3,4)
% % % %         [nH,r0] = hist(dxy1(I_RNA1 > 7),x0);
% % % %         nHr = nH./(xr.^2-xl.^2)/pi;
% % % %         [nL,r0] = hist(dxy1(I_RNA1 <= 2),x0);
% % % %         nLr = nL./(xr.^2-xl.^2)/pi;
% % % % % %         plot(r0,nr,'k','DisplayName','All');
% % % % % %         hold on
% % % %         plot(r0,nHr/sum(nH),'r','DisplayName','Active')
% % % %         hold on
% % % %         plot(r0,nLr/sum(nL),'b','DisplayName','Inactive');
% % % %         xlim([0,xl(end)-1])
% % % %         xlabel('r (um)')
% % % %         ylabel('%/um^2')
% % % %         title('Distance distribution (TX)','Interpreter','none')
% % % %         legend('show')
% % % %     subplot(2,3,2)
% % % %         dx = 0.05;
% % % %         xl = 0:dx:4;
% % % %         xr = xl+dx;
% % % %         x0 = xl+dx/2;
% % % %         [n0,r0] = hist(I_en1,x0);
% % % %         [na,r0] = hist(I_en1(EL_RNA1 < 0.3),x0);
% % % %         [np,r0] = hist(I_en1(EL_RNA1 >= 0.7),x0);
% % % % % %         plot(r0,n0,'k','DisplayName','All');
% % % % % %         hold on
% % % %         plot(r0,na/sum(na),'r','DisplayName','Anterior')
% % % %         hold on
% % % %         plot(r0,np/sum(np),'b','DisplayName','Posterior');
% % % %         xlim([0,xl(end)-1])
% % % %         xlabel('En(#)')
% % % %         ylabel('%')
% % % %         title('En distribution (EL)','Interpreter','none')
% % % %         legend('show')
% % % %     subplot(2,3,5)
% % % %         [nH,r0] = hist(I_en1(I_RNA1 > 7),x0);
% % % %         [nL,r0] = hist(I_en1(I_RNA1 <= 2),x0);
% % % % % %         plot(r0,n0,'k','DisplayName','All');
% % % % % %         hold on
% % % %         plot(r0,nH/sum(nH),'r','DisplayName','Active')
% % % %         hold on
% % % %         plot(r0,nL/sum(nL),'b','DisplayName','Inactive');
% % % %         xlim([0,xl(end)-1])
% % % %         xlabel('En(#)')
% % % %         ylabel('%')
% % % %         title('En distribution (TX)','Interpreter','none')
% % % %         legend('show')
% % % %     subplot(2,3,3)
% % % %         dx = 1;
% % % %         x0 = 0:dx:8;
% % % %         xl = x0;
% % % %         [n0,r0] = hist(N_en1,x0);
% % % %         [na,r0] = hist(N_en1(EL_RNA1 < 0.3),x0);
% % % %         [np,r0] = hist(N_en1(EL_RNA1 >= 0.7),x0);
% % % % % %         plot(r0,n0,'k','DisplayName','All');
% % % % % %         hold on
% % % %         plot(r0,na/sum(na),'r','DisplayName','Anterior')
% % % %         hold on
% % % %         plot(r0,np/sum(np),'b','DisplayName','Posterior');
% % % %         xlim([0,xl(end)-1])
% % % %         xlabel('En spot #')
% % % %         ylabel('%')
% % % %         title('# En spots distribution (EL)','Interpreter','none')
% % % %         legend('show')
% % % %     subplot(2,3,6)
% % % %         [nH,r0] = hist(N_en1(I_RNA1 > 7),x0);
% % % %         [nL,r0] = hist(N_en1(I_RNA1 <= 2),x0);
% % % % % %         plot(r0,n0,'k','DisplayName','All');
% % % % % %         hold on
% % % %         plot(r0,nH/sum(nH),'r','DisplayName','Active')
% % % %         hold on
% % % %         plot(r0,nL/sum(nL),'b','DisplayName','Inactive');
% % % %         xlim([0,xl(end)-1])
% % % %         xlabel('En spot #')
% % % %         ylabel('%')
% % % %         title('# En spots distribution (TX)','Interpreter','none')
% % % %         legend('show')
% % % % 
        