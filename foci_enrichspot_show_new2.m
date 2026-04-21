clear all
close all

tic

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
spot_add = '_spot';
hist_add = '_hist';
fit_tail = '.sfit';

r_th = 0.3;
rout_min = 0.3;
rout_max = 0.5;
rfit_min = 0.4;
rfit_max = 1;
r_ad = 1;
rxmax = 3;
rymax = 3;
relip = 1;
Imin0 = 3e3;
sub_pos = [5,7];
ana_folder = 'enrichspot_result/';
% embryo_name = 'test';
% embryo_name = '09302015_16249-1-5M_FISHIF_hb_gal4_Bcd_par2';
% embryo_name = '11012014_10652-1-3M_FISHIF_lacZ_gal4_Bcd_par2';
embryo_name = '02072016_16249-1-5M_FISHIF_hb_gal4_Cad_Bcd_PBcd_par2';
fit_name = 'enrich_hist_fit';
I_plot = 1:42;%[1:2,4:7,10:14,16,20];%[2:26,29:30,32];%[1:2,4:7,10:14,16,20];
I_cycle = 13;
epsilon = 1e-4;
smin = 1e-2;

gau2_gen = @(n) ['m(',num2str(n),',1)/sqrt(2*pi)/m(',num2str(n),',3)*exp(-(x-m(',num2str(n),',2)).^2/2/m(',num2str(n),',3)^2)'];   %%% single-gaussian model function text generator
gau3_gen = @(n) ['a',num2str(n),'/sqrt(2*pi)/c*exp(-(x-',num2str(n),'*b).^2/2/c^2)'];   %%% single-gaussian model function text generator
lb = [0,0,0.4];
ub = [1,inf,10];
ib = [1,1,1];

IC = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Data loading and adjusting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([ana_folder,embryo_name,'/',embryo_name,mat_tail]);
N_outmin_cell = cell(size(dxy0_cell));
N_outmax_cell = cell(size(dxy0_cell));
en_circle_area_cell = cell(size(dxy0_cell));
dz_en_RNA_cell = cell(size(dxy0_cell));

for ii = 1:length(dxy0_cell)
    [~,I_real1] = spfilter(foci_list0_cell{ii},rxmax,rymax,relip,Imin0);
    I_real2 = Nu_en_cell{ii} == Nu_RNA_cell{ii}(Imin_all0_cell{ii});
    I_real3 = ismember(foci_list0_cell{ii}(:,8)-foci_list1_cell{ii}(Imin_all0_cell{ii},8),dz);
    I_real = I_real1 & I_real2 & I_real3;
    foci_list0_cell{ii} = foci_list0_cell{ii}(I_real,:);
    Nu_en_cell{ii} = Nu_en_cell{ii}(I_real,:);
    xyz_all0_cell{ii} = xyz_all0_cell{ii}(I_real,:);
    dxy0_cell{ii} = dxy0_cell{ii}(I_real,:);
    Imin_all0_cell{ii} = Imin_all0_cell{ii}(I_real,:);
    I_RNA0_cell{ii} = I_RNA0_cell{ii}(I_real,:);
    I_en0_cell{ii} = I_en0_cell{ii}(I_real,:);
    EL_RNA0_cell{ii} = EL_RNA0_cell{ii}(I_real,:);
    C_RNA0_cell{ii} = C_RNA0_cell{ii}(I_real,:);
    C_RNA_min = mean(C_RNA1_cell{ii}(EL_RNA1_cell{ii} > 0.8 & EL_RNA1_cell{ii} <= 1,2));
    C_RNA0_cell{ii}(:,2) = C_RNA0_cell{ii}(:,2)-C_RNA_min;
    C_RNA1_cell{ii}(:,2) = C_RNA1_cell{ii}(:,2)-C_RNA_min;
    C_RNA_minb = mean(C_RNA1_cell{ii}(EL_RNA1_cell{ii} > 0.8 & EL_RNA1_cell{ii} <= 1,3));
    C_RNA0_cell{ii}(:,3) = C_RNA0_cell{ii}(:,3)-C_RNA_minb;
    C_RNA1_cell{ii}(:,3) = C_RNA1_cell{ii}(:,3)-C_RNA_minb;
    dz_en_RNA_cell{ii} = foci_list0_cell{ii}(:,8)-foci_list1_cell{ii}(Imin_all0_cell{ii},8);
    
    dxy = pdist2(xyz_all1_cell{ii}(:,1:2),xyz_all0_cell{ii}(:,1:2));
    [dxymin1,I1min1] = min(dxy,[],2);
    dxy1_cell{ii} = dxymin1;
    Imin_all1_cell{ii} = I1min1;
    I_en1_cell{ii} = I_en0_cell{ii}(I1min1); 
    
    N_en00 = hist(Imin_all0_cell{ii}'.*(dxy0_cell{ii}' <= r_th),0:size(xyz_all1_cell{ii},1));
    N_en00 = N_en00(2:end);
    N_en0_cell{ii} = N_en00(Imin_all0_cell{ii}')';
    N_en1_cell{ii} = N_en00';
    
    N_outmin = hist(Imin_all0_cell{ii}'.*(dxy0_cell{ii}' <= rout_min),0:size(xyz_all1_cell{ii},1));
    N_outmin = N_outmin(2:end);
    N_outmin_cell{ii} = N_outmin';
    
    N_outmax = hist(Imin_all0_cell{ii}'.*(dxy0_cell{ii}' <= rout_max),0:size(xyz_all1_cell{ii},1));
    N_outmax = N_outmax(2:end);
    N_outmax_cell{ii} = N_outmax';
    
    en_circle_area_cell{ii} = foci_circle_area_cell{ii,1}(Imin_all0_cell{ii},:,:);
end

I00 = ismember(cycle_em,I_cycle) & ismember(ind_em,I_plot);

xyz_all0 = cat(1,xyz_all0_cell{I00});
dxy0 = cat(1,dxy0_cell{I00});
Imin_all0 = cat(1,Imin_all0_cell{I00});
I_RNA0 = cat(1,I_RNA0_cell{I00});
I_en0 = cat(1,I_en0_cell{I00});
N_en0 = cat(1,N_en0_cell{I00});
EL_RNA0 = cat(1,EL_RNA0_cell{I00});
C_RNA0 = cat(1,C_RNA0_cell{I00});

xyz_all1 = cat(1,xyz_all1_cell{I00});
dxy1 = cat(1,dxy1_cell{I00});
Imin_all1 = cat(1,Imin_all1_cell{I00});
I_RNA1 = cat(1,I_RNA1_cell{I00});
I_en1 = cat(1,I_en1_cell{I00});
% I_en2 = cat(1,I_en2_cell{I00});
N_en1 = cat(1,N_en1_cell{I00});
EL_RNA1 = cat(1,EL_RNA1_cell{I00});
C_RNA1 = cat(1,C_RNA1_cell{I00});

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

ind_th = find(abs(xl0-r_th) < (xl0(2)-xl0(1))*epsilon);
ind_outmin = find(abs(xl0-rout_min) < (xl0(2)-xl0(1))*epsilon);
ind_outmax = find(abs(xl0-rout_max) < (xl0(2)-xl0(1))*epsilon);
Sa_foci = mean(foci_circle_area(:,ind_th,:),3);
Sp_foci = mean(foci_circle_area(:,ind_outmax,:)-foci_circle_area(:,ind_outmin,:),3);
id1 = sub2ind(size(en_circle_area),[1:length(dxy0)]',ind_th*ones(length(dxy0),1),Nz);
id2 = sub2ind(size(en_circle_area),[1:length(dxy0)]',ind_outmin*ones(length(dxy0),1),Nz);
id3 = sub2ind(size(en_circle_area),[1:length(dxy0)]',ind_outmax*ones(length(dxy0),1),Nz);
Sa_en = en_circle_area(id1);
Sp_en = en_circle_area(id3)-en_circle_area(id2);
S_en_max = max(max(en_circle_area,[],1),[],3)';
Sa_en_ratio = Sa_en/S_en_max(ind_th);
Sp_en_ratio = Sp_en/(S_en_max(ind_outmax)-S_en_max(ind_outmin));

r00  = 0.1625*ones(size(ind_em));

r11 = zeros(0);
for ii = find(I00)'
    r11 = cat(1,r11,r00(ii)*ones(size(I_en0_cell{ii})));%/mean(r00);
end
I_en0 = I_en0./r11;
% r11 = b2_em(index_all0)./b_em(index_all0);


%% Data analysis: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot general properties of spots:
figure(1);maximize(1)
set(gcf,'Name',embryo_name)
subplot(2,3,1)
    id0 = C_RNA0(:,IC) > 10e-9 & en_ring_ratio > smin;
    xdata = dxy0(id0);
    ydata = 1./en_ring_area(id0)./length(xdata);
    [~,nar,~] = histw(xdata,ydata,xl,xr);
    pa = polyfit(x0(x0 >= rfit_min & x0 <= rfit_max),nar(x0 >= rfit_min & x0 <= rfit_max),1);
    
    id0 = C_RNA0(:,IC) < 5e-9 & en_ring_ratio > smin;
    xdata = dxy0(id0);
    ydata = 1./en_ring_area(id0)./length(xdata);
    [~,npr,~] = histw(xdata,ydata,xl,xr);
    pp = polyfit(x0(x0 >= rfit_min & x0 <= rfit_max),npr(x0 >= rfit_min & x0 <= rfit_max),1);
% %         plot(x0,nr,'k','DisplayName','All');
% %         hold on
    plot(x0,nar,'r','DisplayName','Anterior')
    hold on
    plot(x0,polyval(pa,x0),'r-.','DisplayName',['Anterior fit (k = ',num2str(pa(1)),')'])
    plot(x0,npr,'b','DisplayName','Posterior');
    plot(x0,polyval(pp,x0),'b-.','DisplayName',['Posterior fit (k = ',num2str(pp(1)),')'])
    xlim([0,xl(end)-1])
    xlabel('r (um)')
    ylabel('%/um^2')
    title('Distance distribution (EL)','Interpreter','none')
    legend('show')
    
subplot(2,3,4)
    id0 = I_RNA0 > 7 & en_ring_ratio > smin;
    xdata = dxy0(id0);
    ydata = 1./en_ring_area(id0)./length(xdata);
    [~,nHr,~] = histw(xdata,ydata,xl,xr);
    pH = polyfit(x0(x0 >= rfit_min & x0 <= rfit_max),nHr(x0 >= rfit_min & x0 <= rfit_max),1);
    id0 = I_RNA0 <= 2 & en_ring_ratio > smin;
    xdata = dxy0(id0);
    ydata = 1./en_ring_area(id0)./length(xdata);
    [~,nLr,~] = histw(xdata,ydata,xl,xr);
    pL = polyfit(x0(x0 >= rfit_min & x0 <= rfit_max),nLr(x0 >= rfit_min & x0 <= rfit_max),1);
% %         plot(x0,nr,'k','DisplayName','All');
% %         hold on
    plot(x0,nHr,'r','DisplayName','Active')
    hold on
    plot(x0,polyval(pH,x0),'r-.','DisplayName',['Active fit (k = ',num2str(pH(1)),')'])
    plot(x0,nLr,'b','DisplayName','Inactive');
    plot(x0,polyval(pL,x0),'b-.','DisplayName',['Inactive fit (k = ',num2str(pL(1)),')'])
    xlim([0,xl(end)-1])
    xlabel('r (um)')
    ylabel('%/um^2')
    title('Distance distribution (TX)','Interpreter','none')
    legend('show')
    
    
subplot(2,3,2)
    id0 = C_RNA0(:,IC) > 10e-9 & en_ring_ratio > smin;
    xdata = dxy0(id0);
    ydata = I_en0(id0);
    wdata = 1./en_ring_area(id0)./length(xdata);
    [~,en_a,~,~] = equal_distw(xdata,ydata,wdata,xl,xr);
    id0 = C_RNA0(:,IC) < 5e-9 & en_ring_ratio > smin;
    xdata = dxy0(id0);
    ydata = I_en0(id0);
    wdata = 1./en_ring_area(id0)./length(xdata);
    [~,en_p,~,~] = equal_distw(xdata,ydata,wdata,xl,xr);
    plot(x0,en_a,'r','DisplayName','Anterior')
    hold on
    plot(x0,en_p,'b','DisplayName','Posterior');
    xlim([0,xl(end)-1])
    xlabel('r (um)')
    ylabel('<#En>')
    title('#En vs. distance (EL)','Interpreter','none')
    legend('show')
    
subplot(2,3,5)
    id0 = I_RNA0 > 7 & en_ring_ratio > smin;
    xdata = dxy0(id0);
    ydata = I_en0(id0);
    wdata = 1./en_ring_area(id0)./length(xdata);
    [~,en_a,~,~] = equal_distw(xdata,ydata,wdata,xl,xr);
    id0 = I_RNA0 <= 2 & en_ring_ratio > smin;
    xdata = dxy0(id0);
    ydata = I_en0(id0);
    wdata = 1./en_ring_area(id0)./length(xdata);
    [~,en_p,~,~] = equal_distw(xdata,ydata,wdata,xl,xr);
    plot(x0,en_a,'r','DisplayName','Active')
    hold on
    plot(x0,en_p,'b','DisplayName','Inactive');
    xlim([0,xl(end)-1])
    xlabel('r (um)')
    ylabel('<#En>')
    title('#En vs. distance (TX)','Interpreter','none')
    legend('show')
    
    
subplot(2,3,3)
    dN = 1;
    N0 = 0:dN:8;
    id0 = C_RNA1(:,IC) > 10e-9;
    [na,r0] = hist(N_en1(id0),N0);
    id0 = C_RNA1(:,IC) < 5e-9;
    [np,r0] = hist(N_en1(id0),N0);
% %         plot(r0,n0,'k','DisplayName','All');
% %         hold on
    plot(r0,na/sum(na),'r','DisplayName','Anterior')
    hold on
    plot(r0,np/sum(np),'b','DisplayName','Posterior');
    xlim([0,N0(end)-1])
    xlabel('En spot #')
    ylabel('%')
    title('# En spots distribution (EL)','Interpreter','none')
    legend('show')
    
subplot(2,3,6)
    id0 = I_RNA1 > 7;
    [nH,r0] = hist(N_en1(id0),N0);
    id0 = I_RNA1 <= 2;
    [nL,r0] = hist(N_en1(id0),N0);
% %         plot(r0,n0,'k','DisplayName','All');
% %         hold on
    plot(r0,nH/sum(nH),'r','DisplayName','Active')
    hold on
    plot(r0,nL/sum(nL),'b','DisplayName','Inactive');
    xlim([0,N0(end)-1])
    xlabel('En spot #')
    ylabel('%')
    title('# En spots distribution (TX)','Interpreter','none')
    legend('show')
    
    

%%% Plot the intensity distributions of spots (A/P, N = 1,2...):
dI = 0.5;
dw = 0.75;
Il = 0:dI:31;
Ir = Il+dw;
Ic = (Il+Ir)/2;
r0 = Ic;
N_range = 1:2;
N_range0 = 0:max(N_range);
EL_name = {'[Bcd] > 25 nM','[Bcd] > 15 nM','[Bcd] > 10 nM','[Bcd] > 5 nM','[Bcd] < 4 nM'};
EL_range = {[0.15,0.75],[0.15,0.75],[0.15,0.75],[0.15,0.75],[0.15,0.75]};
C_range = {[25e-9,50e-9],[15e-9,25e-9],[10e-9,15e-9],[5e-9,10e-9],[0e-9,4e-9]};
    
PNob = cell(length(EL_name),length(N_range0));
Pob = cell(length(EL_name),length(N_range0));
PNreal = cell(length(EL_name),length(N_range0));
Preal = cell(length(EL_name),length(N_range0));
freal = cell(length(EL_name),length(N_range0));

figure(2);maximize(2)
set(gcf,'Name',[embryo_name,': Multi-spot model'])
for ii = 1:length(EL_name)
    I000 = EL_RNA0 >= EL_range{ii}(1) & EL_RNA0 <= EL_range{ii}(2) & C_RNA0(:,IC) >= C_range{ii}(1) & C_RNA0(:,IC) <= C_range{ii}(2) & dxy0 > rout_min & dxy0 <= rout_max;
    I111 = EL_RNA1 >= EL_range{ii}(1) & EL_RNA1 <= EL_range{ii}(2) & C_RNA1(:,IC) >= C_range{ii}(1) & C_RNA1(:,IC) <= C_range{ii}(2);
    
    Sa0 = pi*r_th^2;
    Sp0 = pi*(rout_max^2-rout_min^2);
    
    xdata = I_en0(I000);
    ydata = 1./Sp_en(I000);
    [~,nf,~] = histw(xdata,ydata,Il,Ir);
    Pfake = nf/nnz(I111)*Sa0;
    ffake = nf/sum(nf);
    kf = sum(nf)/nnz(I111)*r_ad;
    
    for jj = 1:length(N_range0)
        I000 = N_en0 == N_range0(jj) & EL_RNA0 >= EL_range{ii}(1) & EL_RNA0 <= EL_range{ii}(2) & C_RNA0(:,IC) >= C_range{ii}(1) & C_RNA0(:,IC) <= C_range{ii}(2) & dxy0 <= r_th;
        N111 = N_en1 == N_range0(jj) & EL_RNA1 >= EL_range{ii}(1) & EL_RNA1 <= EL_range{ii}(2) & C_RNA1(:,IC) >= C_range{ii}(1) & C_RNA1(:,IC) <= C_range{ii}(2);
        if jj == 1
            PNob{ii,jj} = nnz(N111)/nnz(I111);
            Pob{ii,jj} = zeros(size(Il));
            PNreal{ii,jj} = sum(1./exp(-kf*Sa_foci(N111)))/nnz(I111);
            Preal{ii,jj} = zeros(size(Il));
            freal{ii,jj} = zeros(size(Il));
        else
            xdata = I_en0(I000);
            ydata = 1./Sa_en(I000);
            nob = hist0(xdata,Il,Ir,true);
            PNob{ii,jj} = nnz(N111)/nnz(I111);
            Pob{ii,jj} = nob/nnz(I111);
            
            PNreal{ii,jj} = PNob{ii,jj};
            Preal{ii,jj} = Pob{ii,jj};
            for kk = 1:jj-1
                PNreal{ii,jj} = PNreal{ii,jj}-PNreal{ii,kk}*mean(exp(-kf*Sa_foci(N111)).*(kf*Sa_foci(N111)).^(N_range0(jj)-N_range0(kk))/factorial(N_range0(jj)-N_range0(kk)));
                Preal{ii,jj} = Preal{ii,jj}-(N_range0(kk)*freal{ii,kk}+(N_range0(jj)-N_range0(kk))*ffake)*PNreal{ii,kk}*mean(exp(-kf*Sa_foci(N111)).*(kf*Sa_foci(N111)).^(N_range0(jj)-N_range0(kk))/factorial(N_range0(jj)-N_range0(kk)));
            end
            PNreal{ii,jj} = PNreal{ii,jj}/exp(-kf*Sa0);
            Preal{ii,jj} = Preal{ii,jj}/exp(-kf*Sa0);
            freal{ii,jj} = Preal{ii,jj}/PNreal{ii,jj};
        end

        jj0 = find(N_range0(jj) == N_range);
        if jj0 
            subplot(length(EL_name),length(N_range),(ii-1)*length(N_range)+jj0)
                plot(r0,Pob{ii,jj},'r','DisplayName','Observed')
                hold on
                plot(r0,Pfake,'b','DisplayName','Fake');
                plot(r0,Preal{ii,jj},'k','DisplayName','Real');
                xlim([0,Il(end)-1])
                ylim([-0.1,0.2])
                xlabel('En(#)')
                ylabel('%')
                title(['En distribution (',EL_name{ii},', N = ',num2str(N_range(jj0)),')'],'Interpreter','none')
                legend('show')
        end
    end
end



%%% Plot the intensity distributions of spots (A/P, Assuming N = 1):
dI = 1;
dw = 1;
Il = -0.5:dI:31;
Ir = Il+dw;
Ic = (Il+Ir)/2;
r0 = Ic;
r1 = r0(1):0.01:r0(end);
EL_range = [0.15,0.75];
% % C_min = [0:2.5e-9:22.5e-9,27.5e-9:5e-9:42.5e-9];
% % C_max = [5e-9:2.5e-9:27.5e-9,37.5e-9:5e-9:52.5e-9];
% % C_mean = (C_min+C_max)/2;

Nbin = 16;
rover = 0.6;
I_all = EL_RNA1 >= EL_range(1) & EL_RNA1 <= EL_range(2);
[C_mean,~,~,~,~,C_min,C_max] = equal_bin(C_RNA1(I_all,IC),C_RNA1(I_all,IC),Nbin,rover,true);

sub_pos2 = [4,4];

N_range = 1:7;
dN = 4;

textfun = 'gau2 = @(m,x) ';
for n = 1:length(N_range)
    textfun = [textfun,gau2_gen(n),'+'];
end
textfun = [textfun(1:(end-1)),';'];
eval(textfun);

lb0 = repmat(lb,length(N_range),1);
ub0 = repmat(ub,length(N_range),1);
ib0 = repmat(ib,length(N_range),1);
lb0(:,2) = dN*N_range'*(1-0.4);
ub0(:,2) = dN*N_range'*(1+0.4);
ib0(:,2) = dN*N_range';
options = optimset('Display','off');


N_min = dN*N_range'*(1-0.5);
N_max = dN*N_range'*(1+0.5);


% % % textfun = '';
% % % textfun0 = '';
% % % for n = 1:length(N_range)
% % %     textfun = [textfun,gau3_gen(n),'+'];
% % %     textfun0 = [textfun0,'a',num2str(n),','];
% % % end
% % % textfun = ['gau2 = @(',textfun0,'b,c,x) ',textfun(1:(end-1)),';'];
% % % eval(textfun);
% % % 
% % % lb0 = lb([ones(size(N_range)),2,3]);
% % % ub0 = ub([ones(size(N_range)),2,3]);
% % % ib0 = ib([ones(size(N_range)),2,3]);
% % % lb0(end-1) = dN*(1-0.5);
% % % ub0(end-1) = dN*(1+0.5);
% % % ib0(end-1) = dN;


% % % PNob = cell(length(EL_name),length(N_range0));
% % % Pob = cell(length(EL_name),length(N_range0));
% % % PNreal = cell(length(EL_name),length(N_range0));
Preal = zeros(length(C_min),length(Il));
Preal0 = zeros(length(C_min),length(N_max));
% % % freal = cell(length(EL_name),length(N_range0));
P00 = zeros(size(C_min));
mEn = zeros(size(C_min));
en_fit_cell = cell(length(C_min),3);
en_fit_ob_cell = cell(length(C_min),3);
en_fit_fake_cell = cell(length(C_min),3);
r_ad_all = zeros(size(C_min));

figure(3);maximize(3)
set(gcf,'Name',[embryo_name,': Distance distribution fitting'])
figure(4);maximize(4)
set(gcf,'Name',[embryo_name,': Single-spot model final'])

for ii = 1:length(C_min)
    I000 = EL_RNA0 >= EL_range(1) & EL_RNA0 <= EL_range(2) & C_RNA0(:,IC) >= C_min(ii) & C_RNA0(:,IC) <= C_max(ii) & en_ring_ratio > smin;
    I111 = EL_RNA1 >= EL_range(1) & EL_RNA1 <= EL_range(2) & C_RNA1(:,IC) >= C_min(ii) & C_RNA1(:,IC) <= C_max(ii);
%     N111 = N_en1 == 1 & EL_RNA1 >= EL_range{ii}(1) & EL_RNA1 <= EL_range{ii}(2) & C_RNA1(:,IC) >= C_range{ii}(1) & C_RNA1(:,IC) <= C_range{ii}(2);

    xdata = dxy0(I000);
    ydata = 1./en_ring_area(I000)./length(xdata);
    [~,nar,~] = histw(xdata,ydata,xl,xr);
    pa = polyfit(x0(x0 >= rfit_min & x0 <= rfit_max),nar(x0 >= rfit_min & x0 <= rfit_max),1);
    
    figure(3)
    subplot(sub_pos2(1),sub_pos2(2),ii)
        plot(x0,nar,'k','DisplayName','Data')
        hold on
        plot(x0,polyval(pa,x0),'r-.','DisplayName',['Fit (k = ',num2str(pa(1)),')'])
        xlim([0,xl(end)-1])
        ylim([0,0.6])
        xlabel('r (um)')
        ylabel('%/um^2')
        legend('show')
        title([num2str(C_min(ii)/1e-9),' nM - ',num2str(C_max(ii)/1e-9),' nM'],'Interpreter','none')


    I000 = EL_RNA0 >= EL_range(1) & EL_RNA0 <= EL_range(2) & C_RNA0(:,IC) >= C_min(ii) & C_RNA0(:,IC) <= C_max(ii) & dxy0 > rout_min & dxy0 <= rout_max & Sp_en_ratio > smin;
    I111 = EL_RNA1 >= EL_range(1) & EL_RNA1 <= EL_range(2) & C_RNA1(:,IC) >= C_min(ii) & C_RNA1(:,IC) <= C_max(ii);
    I222 = EL_RNA0 >= EL_range(1) & EL_RNA0 <= EL_range(2) & C_RNA0(:,IC) >= C_min(ii) & C_RNA0(:,IC) <= C_max(ii) & dxy0 <= r_th;
    xdata = I_en0(I000);
    ydata = 1./Sp_en(I000);
    [~,nf,~] = histw(xdata,ydata,Il,Ir);
    Pfake = nf/nnz(I111)*mean(Sa_foci(I111))/dw;
    
    nr = hist0(I_en0(I222),Il,Ir,true);
    Pob_all = nr/nnz(I111)/dw;
    Preal_all = Pob_all-Pfake; 
%     Preal_all = max(Preal_all,0);
    P00(ii) = 1-sum(Preal_all);
    mEn(ii) = Preal_all*r0'/sum(Preal_all);
    Preal(ii,:) = Preal_all;
    
% % %     nr0 = hist0(I_en0(I000 & dxy0 <= r_th)./r11(I000 & dxy0 <= r_th),N_min,N_max,true);
% % %     Pob_all0 = nr0/nnz(I111);
% % %     Preal_all0 = Pob_all0-Pfake0;
% % %     Preal_all0 = max(Preal_all0,0);
% % %     Preal0(ii,:) = Preal_all0;
    
    eval(['xxx',num2str(ii),' = r0;'])
    eval(['yyy',num2str(ii),' = Preal_all;'])

% %     N111 = N_en1 == 0 & EL_RNA1 >= EL_range{ii}(1) & EL_RNA1 <= EL_range{ii}(2) & C_RNA1(:,IC) >= C_range{ii}(1) & C_RNA1(:,IC) <= C_range{ii}(2);
% %     N_en1_data = N_en1(I111);
% %     custpdf = @(data,P0) P0*exp(-kf*Sa)*(kf*Sa).^data./factorial(data)+(1-P0)*exp(-kf*Sa)*(kf*Sa).^max(data-1,0)./factorial(max(data-1,0)).*(data > 0);
% %     P00 = mle(N_en1_data,'pdf',custpdf,'start',nnz(N111)/nnz(I111));
% % % % %     P00(ii) = 1-sum(Preal_all);
% %     mEn(ii) = Preal_all*r0'/sum(Preal_all);

    dNN = pdist2(r0',dN*N_range');
    [~,Imin] = min(dNN);
    ib0(:,1) = Preal_all(Imin)';
% % %     ib0(N_range) = Preal_all(Imin);
% % %     en_fit0 = fit(r0', Preal_all', gau2, 'StartPoint', ib0, 'Upper', ub0, 'Lower', lb0);
% % %     en_fit = zeros(0);
% % %     for ss = 1:length(N_range)
% % %         en_fit = cat(1,en_fit,[eval(['en_fit0.a',num2str(ss)]),en_fit0.b,en_fit0.c]);
% % %     end
% % %     en_fit_cell(ii,:) = {en_fit,[],[]};
    [en_fit,resnorm,~,exitflag,~] = lsqcurvefit(gau2,ib0,r0',Preal_all',lb0,ub0,options);
    en_fit_cell(ii,:) = {en_fit,resnorm,exitflag};

    
    figure(4)
    subplot(sub_pos2(1),sub_pos2(2),ii)
% %         plot(r0,Pob_all,'r','DisplayName','Observed')
% %         hold on
% %         plot(r0,Pfake,'b','DisplayName','Fake');
% %         plot(r0,Preal_all,'k','DisplayName','Real');
        plot(r0,Preal_all,'ko','DisplayName','Real');
        hold on
        plot(r1,gau2(en_fit,r1),'r','DisplayName','Fit');
% %         plot(r1,feval(en_fit0,r1),'r','DisplayName','Fit');
        xlim([0,Il(end)-1])
        ylim([-0.01,0.05])
        xlabel('En(#)')
        ylabel('%')
% % %         title(['[Bcd] = ',num2str(C_mean(ii)/1e-9),' nM'],'Interpreter','none')
        title([num2str(C_min(ii)/1e-9),' nM - ',num2str(C_max(ii)/1e-9),' nM'],'Interpreter','none')
        legend('show')
end
% P00

%%
% % % en_fit_cell2 = cell(size(en_fit_cell,1),1);
% % % v0 = load([ana_folder,embryo_name,'\',fit_name,fit_tail], '-mat');
% % % for ii = 1:length(v0.savedSession.AllFitdevsAndConfigs)
% % %     temp0 = v0.savedSession.AllFitdevsAndConfigs{ii}.Fitdev.Fit;
% % %     temp00 = v0.savedSession.AllFitdevsAndConfigs{ii}.Fitdev.FitName;
% % %     I0 = find(temp00 == ' ');
% % %     tempI = str2num(temp00(I0(end)+1:end));
% % %     
% % %     fit_temp = zeros(0);
% % %     for jj = 1:N_range(end)
% % %         fit_temp = cat(1,fit_temp,[eval(['temp0.a',num2str(jj)]),eval(['temp0.b',num2str(jj)]),eval(['temp0.c',num2str(jj)])]);
% % %     end
% % %     en_fit_cell2{tempI} = fit_temp;
% % % end

en_fit_all = cat(3,en_fit_cell{:,1});
P_all = reshape(en_fit_all(:,1,:),[size(en_fit_all,1),size(en_fit_all,3)]);
P_all = [1/dw-sum(P_all);P_all];
K0 = P_all(2:end,:)./P_all(1:end-1,:);

% % Preal0 = [1-sum(Preal0,2),Preal0]';
% % K0 = Preal0(2:end,:)./Preal0(1:end-1,:);

sub_pos3 = [2,4];

figure(5)
for Ik = 1:size(K0,1)
    subplot(sub_pos3(1),sub_pos3(2),Ik)
        I0 = find(K0(Ik,:) <= 5 & ~isnan(log(K0(Ik,:))) & ~isnan(log(C_mean)) & isfinite(log(C_mean)));
        plot(log10(C_mean(I0)),log10(K0(Ik,I0)),'.')
        xlabel('log10([Bcd])')
        ylabel('log10(k)')
        title(['k',num2str(Ik),num2str(Ik-1)],'Interpreter','none')
end

figure(6)
for Ik = 1:size(K0,1)
    subplot(sub_pos3(1),sub_pos3(2),Ik)
        I0 = find(K0(Ik,:) <= 5);
        plot(C_mean(I0)/1e-9,K0(Ik,I0),'.')
        xlabel('[Bcd] (nM)')
        ylabel('k')
        title(['k',num2str(Ik),num2str(Ik-1)],'Interpreter','none')
end

figure(7)
subplot(1,2,1)
    plot(C_mean,1-P00)
    xlabel('[Bcd] (nM)')
    ylabel('P_O_N')
    title('P_O_N vs. [Bcd]')
subplot(1,2,2)
    plot(C_mean,mEn,'o')
    xlabel('[Bcd] (nM)')
    ylabel('<I_s_p_o_t>')
    title('I_s_p_o_t vs. [Bcd]')



