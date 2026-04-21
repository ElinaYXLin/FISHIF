function [dual_foci_out,dual_fake_out,N_foci] = enrichment2_compare(foci_data,fake_data,h,r_size,Inten_thresh,foci_data2,fake_data2,control_name,image_folder,N_cycle,sub_pos,k_DAPI,N_bin2,N_bin3,varargin)

% N_bin = 25;
Lmin = 0;%0.15;
Lmax = 1;%0.85;
cover_lim = 0;
r_optimize = 4;
N_grid = 1000;
nM_con = 1;%6.02e8;

if isempty(varargin)
    unit1 = 'A.U.';
    unit2 = 'A.U.';
else
    unit1 = '#';
    unit2 = 'M';
    if length(varargin) >= 2
        r_optimize = varargin{2};
    end
    if length(varargin) >= 3
        z_size = varargin{3};
    end
    if length(varargin) >= 4
        real_name = varargin{4};
    end
    if length(varargin) >= 5
        sigmax = varargin{5}(1);
        sigmay = varargin{5}(2);
        sigmaz = varargin{5}(3);
        sigmax2 = varargin{5}(4);
        sigmay2 = varargin{5}(5);
        sigmaz2 = varargin{5}(6);
    else
        sigmax = 1.35;
        sigmay = 1.35;
        sigmaz = 0.8;
        sigmax2 = 1.8;
        sigmay2 = 1.8;
        sigmaz2 = 1.0;
    end
end

op_ratio0 = cell(size(r_size));
ofp_ratio0 = cell(size(r_size));

op_ratio02 = cell(size(r_size));
ofp_ratio02 = cell(size(r_size));

Gm = zeros(1);
Gm2 = zeros(1);


ir = find(r_size == r_optimize);
ind0 = find(h{ir});
[xx,yy,zz] = ind2sub(size(h{ir}),ind0);
G0 = 0;
G02 = 0;
for ii = 1:length(xx)
    for jj = 1:length(xx)
        G0 = G0 + exp(-(xx(ii)-xx(jj))^2/4/sigmax^2-(yy(ii)-yy(jj))^2/4/sigmay^2-(zz(ii)-zz(jj))^2/4/sigmaz^2);
        G02 = G02 + exp(-(xx(ii)-xx(jj))^2/4/sigmax2^2-(yy(ii)-yy(jj))^2/4/sigmay2^2-(zz(ii)-zz(jj))^2/4/sigmaz2^2);
    end
end
Gm = G0/length(xx);
Gm2 = G02/length(xx);

h_area = sum(h{ir}(:));
foci_data{ir}(:,3) = nM_con*h_area*foci_data{ir}(:,9).*foci_data{ir}(:,13);
fake_data{ir}(:,3) = nM_con*h_area*fake_data{ir}(:,9).*fake_data{ir}(:,13);
foci_data2{ir}(:,3) = nM_con*h_area*foci_data2{ir}(:,9).*foci_data2{ir}(:,13);
fake_data2{ir}(:,3) = nM_con*h_area*fake_data2{ir}(:,9).*fake_data2{ir}(:,13);

idx = (foci_data{ir}(:,7)>=Lmin) & (foci_data{ir}(:,7)<=Lmax) & (~isnan(foci_data{ir}(:,4))) & (~isnan(foci_data{ir}(:,3))) & (foci_data{ir}(:,8) > 1) & (foci_data{ir}(:,8) < z_size) & (foci_data{ir}(:,9) >= cover_lim) & (foci_data{ir}(:,1) >= Inten_thresh);
idxf = (foci_data{ir}(:,7)>=Lmin) & (foci_data{ir}(:,7)<=Lmax) & (~isnan(fake_data{ir}(:,4))) & (~isnan(fake_data{ir}(:,3))) & (fake_data{ir}(:,8) > 1) & (fake_data{ir}(:,8) < z_size) & (foci_data{ir}(:,9) >= cover_lim) & (foci_data{ir}(:,1) >= Inten_thresh);
idx2 = (foci_data2{ir}(:,7)>=Lmin) & (foci_data2{ir}(:,7)<=Lmax) & (~isnan(foci_data2{ir}(:,4))) & (~isnan(foci_data2{ir}(:,3))) & (foci_data2{ir}(:,8) > 1) & (foci_data2{ir}(:,8) < z_size) & (foci_data2{ir}(:,9) >= cover_lim);
idx2f = (foci_data2{ir}(:,7)>=Lmin) & (foci_data2{ir}(:,7)<=Lmax) & (~isnan(fake_data2{ir}(:,4))) & (~isnan(fake_data2{ir}(:,3))) & (fake_data2{ir}(:,8) > 1) & (fake_data2{ir}(:,8) < z_size) & (foci_data2{ir}(:,9) >= cover_lim);
id0 = idx & idxf & idx2 & idx2f;

foci_data{ir} = foci_data{ir}(id0,:);
fake_data{ir} = fake_data{ir}(id0,:);
foci_data2{ir} = foci_data2{ir}(id0,:);
fake_data2{ir} = fake_data2{ir}(id0,:);
EL_foci = foci_data{ir}(:,7);
EL_fake = fake_data{ir}(:,7);

op_ratio0{ir} = foci_data{ir}(:,4)-foci_data{ir}(:,3);
ofp_ratio0{ir} = fake_data{ir}(:,4)-fake_data{ir}(:,3);
op_ratio02{ir} = foci_data2{ir}(:,4)-foci_data2{ir}(:,3);
ofp_ratio02{ir} = fake_data2{ir}(:,4)-fake_data2{ir}(:,3);

dual_foci_out = [foci_data{ir}(:,2),foci_data2{ir}(:,2),op_ratio0{ir},op_ratio02{ir},foci_data{ir}(:,1),foci_data{ir}(:,4),foci_data2{ir}(:,4),foci_data{ir}(:,7),EL_foci];
dual_fake_out = [fake_data{ir}(:,2),fake_data2{ir}(:,2),ofp_ratio0{ir},ofp_ratio02{ir},fake_data{ir}(:,1),fake_data{ir}(:,4),fake_data2{ir}(:,4),fake_data{ir}(:,7),EL_fake];
N_foci = size(dual_foci_out,1);

%% Figure plotting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Data sorting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plim1 = [0,max(foci_data{ir}(:,2))];
plim2 = [0,max(foci_data2{ir}(:,2))];
d1 = (plim1(2)-plim1(1))/(N_bin2+1);
d2 = (plim2(2)-plim2(1))/(N_bin2+1);
w1 = d1;
w2 = d2;
pctr1 = (plim1(1)+d1):d1:(plim1(2)-d1);
pctr2 = (plim2(1)+d2):d2:(plim2(2)-d2);

id1 = zeros(1,size(foci_data{ir},1));
id2 = zeros(1,size(foci_data2{ir},1));

for ii = 1:N_bin2
    id1((foci_data{ir}(:,2) >= pctr1(ii)-w1) & (foci_data{ir}(:,2) <= pctr1(ii)+w1)) = ii;
    id2((foci_data2{ir}(:,2) >= pctr2(ii)-w2) & (foci_data2{ir}(:,2) <= pctr2(ii)+w2)) = ii;
end
%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% Data binning:  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
En1_mean = nan(N_bin2);
En1_err = nan(N_bin2);
En1f_mean = nan(N_bin2);
En1f_err = nan(N_bin2);
En2_mean = nan(N_bin2);
En2_err = nan(N_bin2);
En2f_mean = nan(N_bin2);
En2f_err = nan(N_bin2);
Corr12 = nan(N_bin2);
Corr12f = nan(N_bin2);
TX_mean = nan(N_bin2);
TX_std = nan(N_bin2);

for ii = 1:N_bin2
    for jj = 1:N_bin2
        ind12 = (id1 == ii) & (id2 == jj);
        En1_mean(ii,jj) = mean(op_ratio0{ir}(ind12));
        En1_err(ii,jj) = std0(op_ratio0{ir}(ind12));
        En1f_mean(ii,jj) = mean(ofp_ratio0{ir}(ind12));
        En1f_err(ii,jj) = std0(ofp_ratio0{ir}(ind12));
        
        En2_mean(ii,jj) = mean(op_ratio02{ir}(ind12));
        En2_err(ii,jj) = std0(op_ratio02{ir}(ind12));
        En2f_mean(ii,jj) = mean(ofp_ratio02{ir}(ind12));
        En2f_err(ii,jj) = std0(ofp_ratio02{ir}(ind12));
        
        if nnz(ind12) >= 10
%             c0 = corrcoef(op_ratio0{ir}(ind12),op_ratio02{ir}(ind12));
%             Corr12(ii,jj) = c0(1,2);
%             c0 = corrcoef(ofp_ratio0{ir}(ind12),ofp_ratio02{ir}(ind12));
%             Corr12f(ii,jj) = c0(1,2);
            c1 = corr([foci_data{ir}(ind12,4),foci_data2{ir}(ind12,4)]);
            c0 = corr([fake_data{ir}(ind12,4),fake_data2{ir}(ind12,4)]);
            Corr12(ii,jj) = (c1(1,2)-c0(1,2));%/sqrt((c1(1,1)-c0(1,1))*(c1(2,2)-c0(2,2)));
            Corr12f(ii,jj) = c0(1,2);%/sqrt(c0(1,1)*c0(2,2));
        end
        
        TX_mean(ii,jj) = mean(foci_data{ir}(ind12,1));
        TX_std(ii,jj) = std(foci_data{ir}(ind12,1));
    end
end
%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% Data binning for 5766:  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% N_bin3 = 7;
En1m = En1_mean(:);
En2m = En2_mean(:);
TXm = TX_mean(:);
TXs = TX_std(:);
Itrue = ~isnan(En1m) & ~isnan(En2m) & ~isnan(TXm);
En1m = En1m(Itrue);
En2m = En2m(Itrue);
TXm = TXm(Itrue);
TXs = TXs(Itrue);

elim1 = [0,max(En1m(~isnan(En1m)))];
elim2 = [0,max(En2m(~isnan(En2m)))];
ed1 = (elim1(2)-elim1(1))/(N_bin3+1);
ed2 = (elim2(2)-elim2(1))/(N_bin3+1);
ew1 = ed1;
ew2 = ed2;
ectr1 = (elim1(1)+ed1):ed1:(elim1(2)-ed1);
ectr2 = (elim2(1)+ed2):ed2:(elim2(2)-ed2);

eid1 = zeros(1,size(En1m,1));
eid2 = zeros(1,size(En2m,1));

for ii = 1:N_bin3
    eid1((En1m >= ectr1(ii)-ew1) & (En1m <= ectr1(ii)+ew1)) = ii;
    eid2((En2m >= ectr2(ii)-ew2) & (En2m <= ectr2(ii)+ew2)) = ii;
end

eTX_mean = nan(N_bin3);
eTX_std = nan(N_bin3);
for ii = 1:N_bin3
    for jj = 1:N_bin3
        ind12 = (eid1 == ii) & (eid2 == jj);
        eTX_mean(ii,jj) = mean(TXm(ind12));
        eTX_std(ii,jj) = sqrt(mean(TXs(ind12).^2));
    end
end
%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% Figure 5763: Hb@hb(Hb,Bcd) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(5763)
clf
    subplot(2,2,1)
    imagesc(pctr2,pctr1,En1_mean)
    axis xy
    title([image_folder,': ',real_name,' (cycle ',num2str(N_cycle),') enrichment'],'Interpreter','none')
    ylabel([real_name,' level (',unit2,')'])
    xlabel([control_name,' level (',unit2,')'])
    hc = colorbar;
    ylabel(hc,[real_name,' level (',unit1,')'])

    subplot(2,2,2)
    imagesc(pctr2,pctr1,En1_err)
    axis xy
    title([image_folder,': ',real_name,' (cycle ',num2str(N_cycle),') enrichment error'],'Interpreter','none')
    ylabel([real_name,' level (',unit2,')'])
    xlabel([control_name,' level (',unit2,')'])
    hc = colorbar;
    ylabel(hc,[real_name,' error (',unit1,')'])

    subplot(2,2,3)
    imagesc(pctr2,pctr1,En1f_mean)
    axis xy
    title([image_folder,': ',real_name,' (cycle ',num2str(N_cycle),') fake enrichment'],'Interpreter','none')
    ylabel([real_name,' level (',unit2,')'])
    xlabel([control_name,' level (',unit2,')'])
    hc = colorbar;
    ylabel(hc,[real_name,' fake level (',unit1,')'])

    subplot(2,2,4)
    imagesc(pctr2,pctr1,En1f_err)
    axis xy
    title([image_folder,': ',real_name,' (cycle ',num2str(N_cycle),') fake enrichment error'],'Interpreter','none')
    ylabel([real_name,' level (',unit2,')'])
    xlabel([control_name,' level (',unit2,')'])
    hc = colorbar;
    ylabel(hc,[real_name,' fake error (',unit1,')'])
    
figure(50763)
    subaxis(sub_pos(1),sub_pos(2),sub_pos(3), 'Spacing', 0, 'PaddingRight',0.01, 'PaddingLeft',0.02, 'PaddingTop',0.045, 'PaddingBottom',0.04, 'Margin', 0.01);
    imagesc(pctr2,pctr1,En1f_err,'DisplayName',[real_name,' fake error'])
    hold on
    imagesc(pctr2,pctr1,En1f_mean,'DisplayName',[real_name,' fake level'])
    imagesc(pctr2,pctr1,En1_err,'DisplayName',[real_name,' error'])
    imagesc(pctr2,pctr1,En1_mean,'DisplayName',[real_name,' level'])
    axis xy
    title([image_folder,': ',real_name,' (cycle ',num2str(N_cycle),') enrichment'],'Interpreter','none')
    ylabel([real_name,' level (',unit2,')'])
    xlabel([control_name,' level (',unit2,')'])
    hc = colorbar;
    ylabel(hc,[real_name,' level (',unit1,')'])
%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% Figure 5764: Bcd@hb(Hb,Bcd) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(5764)
clf
    subplot(2,2,1)
    imagesc(pctr2,pctr1,En2_mean)
    axis xy
    title([image_folder,': ',control_name,' (cycle ',num2str(N_cycle),') enrichment'],'Interpreter','none')
    ylabel([real_name,' level (',unit2,')'])
    xlabel([control_name,' level (',unit2,')'])
    hc = colorbar;
    ylabel(hc,[control_name,' level (',unit1,')'])

    subplot(2,2,2)
    imagesc(pctr2,pctr1,En2_err)
    axis xy
    title([image_folder,': ',control_name,' (cycle ',num2str(N_cycle),') enrichment error'],'Interpreter','none')
    ylabel([real_name,' level (',unit2,')'])
    xlabel([control_name,' level (',unit2,')'])
    hc = colorbar;
    ylabel(hc,[control_name,' error (',unit1,')'])

    subplot(2,2,3)
    imagesc(pctr2,pctr1,En2f_mean)
    axis xy
    title([image_folder,': ',control_name,' (cycle ',num2str(N_cycle),') fake enrichment'],'Interpreter','none')
    ylabel([real_name,' level (',unit2,')'])
    xlabel([control_name,' level (',unit2,')'])
    hc = colorbar;
    ylabel(hc,[control_name,' fake level (',unit1,')'])

    subplot(2,2,4)
    imagesc(pctr2,pctr1,En2f_err)
    axis xy
    title([image_folder,': ',control_name,' (cycle ',num2str(N_cycle),') fake enrichment error'],'Interpreter','none')
    ylabel([real_name,' level (',unit2,')'])
    xlabel([control_name,' level (',unit2,')'])
    hc = colorbar;
    ylabel(hc,[control_name,' fake error (',unit1,')'])

figure(50764)
    subaxis(sub_pos(1),sub_pos(2),sub_pos(3), 'Spacing', 0, 'PaddingRight',0.01, 'PaddingLeft',0.02, 'PaddingTop',0.045, 'PaddingBottom',0.04, 'Margin', 0.01);
    imagesc(pctr2,pctr1,En2f_err,'DisplayName',[control_name,' fake error'])
    hold on
    imagesc(pctr2,pctr1,En2f_mean,'DisplayName',[control_name,' fake level'])
    imagesc(pctr2,pctr1,En2_err,'DisplayName',[control_name,' error'])
    imagesc(pctr2,pctr1,En2_mean,'DisplayName',[control_name,' level'])
    axis xy
    title([image_folder,': ',control_name,' (cycle ',num2str(N_cycle),') enrichment'],'Interpreter','none')
    ylabel([real_name,' level (',unit2,')'])
    xlabel([control_name,' level (',unit2,')'])
    hc = colorbar;
    ylabel(hc,[control_name,' level (',unit1,')'])
%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% Figure 5765: corr(Hb@hb,Bcd@hb) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(5765)
clf
    subplot(1,2,1)
    imagesc(pctr2,pctr1,Corr12)
    axis xy
    title([image_folder,' (cycle ',num2str(N_cycle),'): dual enrichment correlation'],'Interpreter','none')
    ylabel([real_name,' level (',unit2,')'])
    xlabel([control_name,' level (',unit2,')'])
    hc = colorbar;
    ylabel(hc,'Correlation coefficient')

    subplot(1,2,2)
    imagesc(pctr2,pctr1,Corr12f)
    axis xy
    title([image_folder,' (cycle ',num2str(N_cycle),'): fake dual enrichment correlation'],'Interpreter','none')
    ylabel([real_name,' level (',unit2,')'])
    xlabel([control_name,' level (',unit2,')'])
    hc = colorbar;
    ylabel(hc,'Fake correlation coefficient')

figure(50765)
    subaxis(sub_pos(1),sub_pos(2),sub_pos(3), 'Spacing', 0, 'PaddingRight',0.01, 'PaddingLeft',0.02, 'PaddingTop',0.045, 'PaddingBottom',0.04, 'Margin', 0.01);
    imagesc(pctr2,pctr1,Corr12f,'DisplayName','Fake enrichment correlation')
    hold on
    imagesc(pctr2,pctr1,Corr12,'DisplayName','Enrichment correlation')
    axis xy
    title([image_folder,' (cycle ',num2str(N_cycle),'): enrichment correlation'],'Interpreter','none')
    ylabel([real_name,' level (',unit2,')'])
    xlabel([control_name,' level (',unit2,')'])
    hc = colorbar;
    ylabel(hc,'Correlation coefficient')
%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% Figure 5766: mean corr(Hb@hb,Bcd@hb) bar %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(5766)
clf
    bar([1,2],[mean(Corr12(isfinite(Corr12))),mean(Corr12f(isfinite(Corr12f)))])
    set(gca,'XTick',[1,2],'XTickLabel',{'Foci','Fake'})
    title([image_folder,' (cycle ',num2str(N_cycle),'): dual enrichment correlation mean'],'Interpreter','none')
    ylabel('Correlation coefficient')

figure(50766)
    subaxis(sub_pos(1),sub_pos(2),sub_pos(3), 'Spacing', 0, 'PaddingRight',0.01, 'PaddingLeft',0.02, 'PaddingTop',0.045, 'PaddingBottom',0.04, 'Margin', 0.01);
    bar([1,2],[mean(Corr12(~isnan(Corr12))),mean(Corr12f(~isnan(Corr12f)))])
    set(gca,'XTick',[1,2],'XTickLabel',{'Foci','Fake'})
    title([image_folder,' (cycle ',num2str(N_cycle),'): dual enrichment correlation mean'],'Interpreter','none')
    ylabel('Correlation coefficient')
%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% Figure 5767: TX(Hb@hb,Bcd@hb) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TXmmax = 50;
TXsmax = 40;

figure(5767)
clf
    subplot(2,2,1)
    imagesc(ectr2,ectr1,eTX_mean)
    axis xy
    colormap([[0,0,0];colormap]);
    set(gca,'CLim',[0,TXmmax]);
    title([image_folder,' (cycle ',num2str(N_cycle),'): dual TF transcription activation relation'],'Interpreter','none')
    ylabel([real_name,' mean enrichment (',unit1,')'])
    xlabel([control_name,' mean enrichment (',unit1,')'])
    hc = colorbar;
    ylabel(hc,['Mean RNA level (',unit1,')'])

    subplot(2,2,2)
    imagesc(ectr2,ectr1,eTX_std)
    axis xy
    colormap([[0,0,0];colormap]);
    set(gca,'CLim',[0,TXsmax]);
    title([image_folder,' (cycle ',num2str(N_cycle),'): dual TF transcription activation std'],'Interpreter','none')
    ylabel([real_name,' mean enrichment (',unit1,')'])
    xlabel([control_name,' mean enrichment (',unit1,')'])
    hc = colorbar;
    ylabel(hc,['RNA level std (',unit1,')'])

cmap = colormap;
    subplot(2,2,3)
    dTXm = max(TXm)/(size(cmap,1)-1);
    cTXm = round(TXm/dTXm+1);
    cTXm(cTXm < 1) = 1;
    for Ip = 1:length(En1m)
        plot(En2m(Ip),En1m(Ip),'Marker','o','Color',cmap(cTXm(Ip),:))
        hold on
    end
    title([image_folder,' (cycle ',num2str(N_cycle),'): dual TF transcription activation relation (scatter)'],'Interpreter','none')
    ylabel([real_name,' mean enrichment (',unit1,')'])
    xlabel([control_name,' mean enrichment (',unit1,')'])
    set(gca,'CLim',[0,max(TXm)])
    hc = colorbar;
    ylabel(hc,['Mean RNA level (',unit1,')'])

    subplot(2,2,4)
    dTXs = max(TXs)/(size(cmap,1)-1);
    cTXs = round(TXs/dTXs+1);
    cTXs(cTXs < 1) = 1;
    for Ip = 1:length(En1m)
        plot(En2m(Ip),En1m(Ip),'Marker','o','Color',cmap(cTXs(Ip),:))
        hold on
    end
    title([image_folder,' (cycle ',num2str(N_cycle),'): dual TF transcription activation std (scatter)'],'Interpreter','none')
    ylabel([real_name,' mean enrichment (',unit1,')'])
    xlabel([control_name,' mean enrichment (',unit1,')'])
    set(gca,'CLim',[0,max(TXs)])
    hc = colorbar;
    ylabel(hc,['RNA level std (',unit1,')'])

figure(50767)
    subaxis(sub_pos(1),sub_pos(2),sub_pos(3), 'Spacing', 0, 'PaddingRight',0.01, 'PaddingLeft',0.02, 'PaddingTop',0.045, 'PaddingBottom',0.04, 'Margin', 0.01);
    imagesc(ectr2,ectr1,eTX_std,'DisplayName','TX activation std')
    colormap([[0,0,0];colormap]);
    set(gca,'CLim',[0,TXmmax]);
    hold on
    imagesc(ectr2,ectr1,eTX_mean,'DisplayName','TX activation')
    axis xy
    title([image_folder,' (cycle ',num2str(N_cycle),'): dual TF transcription activation relation'],'Interpreter','none')
    ylabel([real_name,' mean enrichment (',unit1,')'])
    xlabel([control_name,' mean enrichment (',unit1,')'])
    hc = colorbar;
    ylabel(hc,['Mean RNA level (',unit1,')'])
%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%










