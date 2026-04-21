function [dist_all,sisterI,rr1] = TXcorr_plot(foci_RNA_profile,EL_nu,ind_nu,dim_xyz,N_cycle,image_folder,sub_pos,varargin)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TXcorr_plot: A function to plot the correlation between paired foci. %%%
%% 1. Calculate the distances between foci from the same nuclei;
%% 2. plot the distance between paired foci, and determine the threshold for sister-foci
%% 3. Calculate/plot correlation between paired foci.
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% foci_RNA_profile (input): foci RNA profile(:,[2,3,5]) = [ind_nu,RNA#,linear_ind];
%% EL_nu (input): list of nuclear EL;
%% ind_nu (input): indices for good nuclei;
%% dim_xyz (input): xyz dimensions;
%% N_cycle (input): cycle of embryo;
%% image_folder (input): image folder;
%% varargin{1} (input): output figure numbers
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Initialization: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Itrue = ismember(foci_RNA_profile(:,2),ind_nu);
foci_RNA_profile = foci_RNA_profile(Itrue,:);

if isempty(varargin{1})
    h1 = 103;
    h2 = 104;
    h3 = 105;
    h4 = 106;
else
    h1 = varargin{1}(1);
    h2 = varargin{1}(2);
    h3 = varargin{1}(3);
    h4 = varargin{1}(4);
end

EL_min0 = 0.1;
EL_max0 = 0.5;
EL_bin = 0.05;
EL_delta = 0.025;
EL_min = EL_min0:EL_delta:(EL_max0-EL_bin);
EL_max = (EL_min0+EL_bin):EL_delta:EL_max0;

dth = 9;

xbin = [0:100];
ccode = 'crgbky';
icolor = N_cycle-9;
icolor(icolor < 1) = 1;icolor(icolor > 6) = 6;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% 1. Calculate the distances between foci from the same nuclei: %%%%%%%%%%
dist_all = zeros(0);
sisterI = zeros(0);
S_nu = intersect(find(EL_nu >= EL_min0 & EL_nu <= EL_max0),ind_nu);
for ii = 1:length(S_nu)
    ind_foci = (foci_RNA_profile(:,2) == S_nu(ii));
    N0 = nnz(ind_foci);
    if N0 >= 2
        linear_ind = foci_RNA_profile(ind_foci,5);
        [x_ind,y_ind,z_ind] = ind2sub(dim_xyz,linear_ind);
        dist0 = pdist([x_ind,y_ind]);
        dist_all = cat(2,dist_all,dist0);
        
        if any(dist0 <= dth)
            ind0 = find(ind_foci);
            ind2 = tril(squareform(dist0 <= dth),-1);
            [ind2x,ind2y] = find(ind2);
            sisterI = cat(1,sisterI,[foci_RNA_profile(ind0(ind2x),3),foci_RNA_profile(ind0(ind2y),3)]);
        end
    end
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% 2. plot the distance between paired foci, and determine the threshold for sister-foci
figure(h1)
clf
    hist(dist_all,xbin)
    set(gca,'FontName','Arial','FontSize',16)
    xlabel('Distance (pixel)','FontName','Arial','FontSize',16,'FontWeight','bold')
    ylabel('Count (#)','FontName','Arial','FontSize',16,'FontWeight','bold')
    title([image_folder,': Cycle ',num2str(N_cycle)],'FontName','Arial','FontSize',16,'FontWeight','bold','Interpreter','none')
figure(h3)
subaxis(sub_pos(1),sub_pos(2),sub_pos(3), 'Spacing', 0, 'PaddingRight',0.01, 'PaddingLeft',0.02, 'PaddingTop',0.045, 'PaddingBottom',0.04, 'Margin', 0.01);
    hist(dist_all,xbin)
    set(gca,'FontName','Arial','FontSize',8)
    xlabel('Distance (pixel)','FontName','Arial','FontSize',8)
    ylabel('Count (#)','FontName','Arial','FontSize',8)
    xlim([min(xbin),max(xbin)])
    title([image_folder,': Cycle ',num2str(N_cycle)],'FontName','Arial','FontSize',8,'Interpreter','none')
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% 3. Calculate/plot correlation between paired foci
if size(sisterI,1) >= 2
    rr = corrcoef(sisterI);
    rr1 = rr(1,2);
else
    rr1 = nan;
end
figure(h2)
clf
    plot(sisterI(:,1),sisterI(:,2),'b.','DisplayName',['R = ',num2str(rr1)])
    set(gca,'FontName','Arial','FontSize',16)
    xlabel('TX1 (#)','FontName','Arial','FontSize',16,'FontWeight','bold')
    ylabel('TX2 (#)','FontName','Arial','FontSize',16,'FontWeight','bold')
    title([image_folder,': Cycle ',num2str(N_cycle)],'FontName','Arial','FontSize',16,'FontWeight','bold','Interpreter','none')
    legend('show')
figure(h4)
subaxis(sub_pos(1),sub_pos(2),sub_pos(3), 'Spacing', 0, 'PaddingRight',0.01, 'PaddingLeft',0.02, 'PaddingTop',0.045, 'PaddingBottom',0.04, 'Margin', 0.01);
    plot(sisterI(:,1),sisterI(:,2),'b.','DisplayName',['R = ',num2str(rr1)])
    set(gca,'FontName','Arial','FontSize',8)
    xlabel('TX1 (#)','FontName','Arial','FontSize',8)
    ylabel('TX2 (#)','FontName','Arial','FontSize',8)
    title(['R = ',num2str(rr1),', N = ',num2str(size(sisterI,1))],'FontName','Arial','FontSize',8,'Interpreter','none')
%     legend('show')

