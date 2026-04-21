function [H_I,H_N,H_A,para_EL,para_Bcd] = dual_all_show(total_name,out_name,info_data,HH_I,pro_out,FISHI_out,FISHN_out,fociN_out,RNAI_out,RNAN_out,null_out,pro_center0,pro_center,p_name,h0,varargin)

%% Initialization: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(varargin) && ~isempty(varargin{1})
    fn = varargin{1};
else
    fn = [122,105,107,112,115,117,131,191,101];
end
if length(varargin) > 1
    color_code = varargin{2};
else
	color_code = mycolors(6);
end

good_data = logical(info_data(:,1));
cycle_data = info_data(:,2);
Nmin = 10;
Nbin = 20;
fth = @(alpha1,x) 2/3+8/3/alpha1(1).*x-(8/alpha1(1)^2/alpha1(2)+16/3/alpha1(1)/alpha1(2)).*x.^2+16/alpha1(1)^2/alpha1(2)^2.*x.^3+64/alpha1(1)^4/alpha1(2)^3.*x.^4-128/alpha1(1)^4/alpha1(2)^4.*x.^5+(-32/alpha1(1)^3/alpha1(2)^2.*x.^3+64*(1/alpha1(1)^3/alpha1(2)^3-1/alpha1(1)^4/alpha1(2)^3).*x.^4+128/alpha1(1)^4/alpha1(2)^4.*x.^5).*exp(-1/2*alpha1(1)*alpha1(2)./x);
alpha0 = [0.2,65];
EL_min = 0.20;%0.25;
EL_max = 0.7;%0.75;
pro_th = [1.5e-9,1e-9];
cycle_diff = 14;
N01 = 4;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Mean plot of protein profile: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_data = find(good_data);
    figure(fn(1))
for ii = 1:length(N_data)
    i_data = N_data(ii);
    plot(pro_out{i_data}(:,1),pro_out{i_data}(:,2),'Color',color_code(mod(ii-1,size(color_code,1))+1,:),'DisplayName',[out_name{i_data},', cycle ',num2str(cycle_data(i_data))])
    hold on
end

xy_out0 = cell2mat(pro_out(good_data));
nu_L0 = xy_out0(:,1);
nu_mean0 = xy_out0(:,2);

nu_L = sort(unique(nu_L0));
nu_mean = zeros(size(nu_L));
nu_err = zeros(size(nu_L));
for I_nu = 1:length(nu_L)
    mean_raw = nu_mean0(nu_L0 == nu_L(I_nu));
    I_good = ~isnan(mean_raw);
    nu_mean(I_nu) = mean(mean_raw(I_good));
    nu_err(I_nu) = std0(mean_raw(I_good));
end

figure(fn(1))
    errorbar(nu_L,nu_mean,nu_err,'r','DisplayName','Mean profile')
    title([total_name,': ',p_name,' profile (from ',num2str(nnz(good_data)),' embryos)'],'Interpreter','none')
    ylabel([p_name,' concentration (M)'])
    xlabel('AP axis (normalized)')
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Mean plot of FISH intensity vs EL: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_data = find(good_data);
for ii = 1:length(N_data)
    i_data = N_data(ii);
    figure(fn(2))
        plot(FISHI_out{i_data}(:,1),FISHI_out{i_data}(:,2),'Color',color_code(mod(ii-1,size(color_code,1))+1,:),'DisplayName',[out_name{i_data},', cycle ',num2str(cycle_data(i_data))])
        hold on
    figure(fn(7))
        Itrue = FISHI_out{i_data}(:,1) >= EL_min & FISHI_out{i_data}(:,1) <= EL_max & FISHI_out{i_data}(:,4) >= Nmin;
    subplot(1,2,1)
        plot(FISHI_out{i_data}(Itrue,2),FISHI_out{i_data}(Itrue,3).^2./FISHI_out{i_data}(Itrue,2),'.','Color',color_code(mod(ii-1,size(color_code,1))+1,:),'DisplayName',[out_name{i_data},', cycle ',num2str(cycle_data(i_data))])
        hold on
        plot(FISHI_out{i_data}(Itrue,2),FISHI_out{i_data}(Itrue,7)./FISHI_out{i_data}(Itrue,2),'.','Color',color_code(mod(ii-1,size(color_code,1))+1,:)/2,'DisplayName',[out_name{i_data},', cycle ',num2str(cycle_data(i_data)),' EX'])
    subplot(1,2,2)
        plot(FISHI_out{i_data}(Itrue,2),FISHI_out{i_data}(Itrue,5),'.','Color',color_code(mod(ii-1,size(color_code,1))+1,:),'DisplayName',[out_name{i_data},', cycle ',num2str(cycle_data(i_data))])
        hold on
        plot(FISHI_out{i_data}(Itrue,2),FISHI_out{i_data}(Itrue,6),'.','Color',color_code(mod(ii-1,size(color_code,1))+1,:)/2,'DisplayName',[out_name{i_data},', cycle ',num2str(cycle_data(i_data)),' control'])
        plot(FISHI_out{i_data}(Itrue,2),FISHI_out{i_data}(Itrue,7)./FISHI_out{i_data}(Itrue,3).^2,'.','Color',color_code(mod(ii-1,size(color_code,1))+1,:)/4,'DisplayName',[out_name{i_data},', cycle ',num2str(cycle_data(i_data)),' EX'])
        plot(FISHI_out{i_data}(Itrue,2),FISHI_out{i_data}(Itrue,8),'.','Color',color_code(mod(ii-1,size(color_code,1))+1,:),'DisplayName',[out_name{i_data},', cycle ',num2str(cycle_data(i_data)),' digital'])
        plot(FISHI_out{i_data}(Itrue,2),FISHI_out{i_data}(Itrue,9),'.','Color',color_code(mod(ii-1,size(color_code,1))+1,:)/2,'DisplayName',[out_name{i_data},', cycle ',num2str(cycle_data(i_data)),' digital control'])
        plot(FISHI_out{i_data}(Itrue,2),FISHI_out{i_data}(Itrue,10),'-','Color',color_code(mod(ii-1,size(color_code,1))+1,:)/4,'DisplayName',[out_name{i_data},', cycle ',num2str(cycle_data(i_data)),' p(0)'])
end

xy_out0 = cell2mat(FISHI_out(good_data));
nu_L0 = xy_out0(:,1);
nu_mean0 = xy_out0(:,2);
nu_var0 = xy_out0(:,3).^2./xy_out0(:,2);
nu_n0 = xy_out0(:,4);
nu_corr0 = xy_out0(:,5);
nu_corrn0 = xy_out0(:,6);
nu_corr0(isnan(nu_corr0)) = 0;
nu_corrn0(isnan(nu_corrn0)) = 0;
nu_vex0 = xy_out0(:,7)./xy_out0(:,2);
nu_cex0 = xy_out0(:,7)./xy_out0(:,3).^2;
nu_corrb0 = xy_out0(:,8);
nu_corrnb0 = xy_out0(:,9);
nu_corrb0(isnan(nu_corrb0)) = 0;
nu_corrnb0(isnan(nu_corrnb0)) = 0;
nu_p0 = xy_out0(:,10);

[nu_L,EL_min0,EL_max0] = seg_eqN(nu_L0,Nbin,1./Nbin);
% nu_L = sort(unique(nu_L0));
EL_min0 = EL_min0(nu_L >= EL_min & nu_L <= EL_max);
EL_max0 = EL_max0(nu_L >= EL_min & nu_L <= EL_max);
nu_L = nu_L(nu_L >= EL_min & nu_L <= EL_max);
nu_mean = zeros(size(nu_L));
nu_err = zeros(size(nu_L));
% % var_mean = zeros(size(nu_L));
% % var_err = zeros(size(nu_L));
% % corr_mean = zeros(size(nu_L));
% % corr_err = zeros(size(nu_L));
for I_nu = 1:length(nu_L)
    mean_raw = nu_mean0(nu_L0 >= EL_min0(I_nu) & nu_L0 <= EL_max0(I_nu));
% %     var_raw = nu_var0(nu_L0 >= EL_min0(I_nu) & nu_L0 <= EL_max0(I_nu));
% %     corr_raw = nu_corr0(nu_L0 >= EL_min0(I_nu) & nu_L0 <= EL_max0(I_nu));
    I_good = nu_n0(nu_L0 >= EL_min0(I_nu) & nu_L0 <= EL_max0(I_nu)) >= Nmin;
    nu_mean(I_nu) = mean(mean_raw(I_good));
    nu_err(I_nu) = std0(mean_raw(I_good));
% %     var_mean(I_nu) = mean(var_raw(I_good));
% %     var_err(I_nu) = std0(var_raw(I_good));
% %     corr_mean(I_nu) = mean(corr_raw(I_good));
% %     corr_err(I_nu) = std0(corr_raw(I_good));
end

I_good = nu_n0 >= Nmin & nu_L0 >= EL_min & nu_L0 <= EL_max;
[nu_mean1,var_mean1,nu_err1,var_err1] = equal_bin(nu_mean0(I_good),nu_var0(I_good),Nbin);
[nu_mean1,corr_mean1,nu_err1,corr_err1] = equal_bin(nu_mean0(I_good),nu_corr0(I_good),Nbin);
[nu_mean1,corrn_mean1,nu_err1,corrn_err1] = equal_bin(nu_mean0(I_good),nu_corrn0(I_good),Nbin);
[nu_mean1,vex_mean1,nu_err1,vex_err1] = equal_bin(nu_mean0(I_good),nu_vex0(I_good),Nbin);
[nu_mean1,cex_mean1,nu_err1,cex_err1] = equal_bin(nu_mean0(I_good),nu_cex0(I_good),Nbin);
[nu_mean1,corrb_mean1,nu_err1,corrb_err1] = equal_bin(nu_mean0(I_good),nu_corrb0(I_good),Nbin);
[nu_mean1,corrnb_mean1,nu_err1,corrnb_err1] = equal_bin(nu_mean0(I_good),nu_corrnb0(I_good),Nbin);
[nu_mean1,p0_mean1,nu_err1,p0_err1] = equal_bin(nu_mean0(I_good),nu_p0(I_good),Nbin);


try
    [alpha1,r,J,cov,~] = nlinfit(nu_mean1,var_mean1,fth,alpha0);
    koff0 = alpha1(1)*(alpha1(2)./2./nu_mean1-1);
    fp0th = P0_2S(alpha1(1),koff0,alpha1(2));
catch
    alpha1 = nan(size(alpha0));
    fp0th = nan(size(nu_mean1));
end


figure(fn(2))
    errorbar(nu_L,nu_mean,nu_err,'r','DisplayName','Mean profile')
    title([total_name,': TX level profile (from ',num2str(nnz(good_data)),' embryos)'],'Interpreter','none')
    ylabel('# of mRNA')
    xlabel('AP axis (normalized)')
    
figure(fn(7))
subplot(1,2,1)
    errorbarxy(nu_mean1,var_mean1,nu_err1,var_err1,nu_err1,var_err1,'r','r','.','DisplayName','Mean curve')
    hold on
    if all(~isnan(alpha1))
        plot(nu_mean1,fth(alpha1,nu_mean1),'m-','DisplayName',['Fit: k_o_n*T = ',num2str(alpha1(1),'%4.3f'),', k_T_X*T = ',num2str(alpha1(2),'%4.3f')]);
    else
        plot([],[],'m-','DisplayName',['Fit: k_o_n*T = ',num2str(alpha1(1),'%4.3f'),', k_T_X*T = ',num2str(alpha1(2),'%4.3f')]);
    end
    errorbarxy(nu_mean1,vex_mean1,nu_err1,vex_err1,nu_err1,vex_err1,'b','b','.','DisplayName','Mean ex curve')
    title([total_name,': TX mean vs variance (from ',num2str(nnz(good_data)),' embryos)',char(10),'Fit: kon*T = ',num2str(alpha1(1),'%4.3f'),', kTX*T = ',num2str(alpha1(2),'%4.3f')],'Interpreter','none')
    xlabel('Mean # mRNA (EL binned)')
    ylabel('Var/mean (EL binned)')
    
subplot(1,2,2)
    errorbarxy(nu_mean1,corr_mean1,nu_err1,corr_err1,nu_err1,corr_err1,'r','r','.','DisplayName','Mean curve')
    hold on
    errorbarxy(nu_mean1,corrn_mean1,nu_err1,corrn_err1,nu_err1,corrn_err1,'k','k','.','DisplayName','Mean control curve')
    errorbarxy(nu_mean1,cex_mean1,nu_err1,cex_err1,nu_err1,cex_err1,'b','b','.','DisplayName','Mean ex curve')
    errorbarxy(nu_mean1,corrb_mean1,nu_err1,corrb_err1,nu_err1,corrb_err1,[0.5,0,0],[0.5,0,0],'.','DisplayName','Mean digital curve')
    errorbarxy(nu_mean1,corrnb_mean1,nu_err1,corrnb_err1,nu_err1,corrnb_err1,[0.5,0.5,0.5],[0.5,0.5,0.5],'.','DisplayName','Mean digital control curve')
    errorbarxy(nu_mean1,p0_mean1,nu_err1,p0_err1,nu_err1,p0_err1,[0.5,0,0.5],[0.5,0,0.5],'-','DisplayName','Mean p(0) curve')
    errorbarxy(nu_mean1,fp0th,nu_err1,zeros(size(nu_err1)),nu_err1,zeros(size(nu_err1)),'g','g','-','DisplayName','Theoretical mean p(0) curve')
    title([total_name,': TX mean vs correlation (from ',num2str(nnz(good_data)),' embryos)'],'Interpreter','none')
    xlabel('Mean # mRNA (EL binned)')
    ylabel('Correlation coefficient (EL binned)')
    
if all(~isnan(alpha1))
    ud_alpha1 = nlparci(alpha1,r,'covar',cov);
    err_alpha1 = alpha1-ud_alpha1(:,1)';
else
    err_alpha1 = nan(size(alpha0));
end
para_EL = [nan(1,2),alpha1,nan(1,2);nan(1,2),err_alpha1,nan(1,2)];
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Mean plot of FISH foci # vs EL: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(fn(3))
N_data = find(good_data);
for ii = 1:length(N_data)
    i_data = N_data(ii);
    plot(FISHN_out{i_data}(:,1),FISHN_out{i_data}(:,2),'Color',color_code(mod(ii-1,size(color_code,1))+1,:),'DisplayName',[out_name{i_data},', cycle ',num2str(cycle_data(i_data))])
    hold on
end

xy_out0 = cell2mat(FISHN_out(good_data));
nu_L0 = xy_out0(:,1);
nu_mean0 = xy_out0(:,2);

nu_L = sort(unique(nu_L0));
nu_mean = zeros(size(nu_L));
nu_err = zeros(size(nu_L));
for I_nu = 1:length(nu_L)
    mean_raw = nu_mean0(nu_L0 == nu_L(I_nu));
    I_good = ~isnan(mean_raw);
    nu_mean(I_nu) = mean(mean_raw(I_good));
    nu_err(I_nu) = std0(mean_raw(I_good));
end

figure(fn(3))
    errorbar(nu_L,nu_mean,nu_err,'r','DisplayName','Mean profile')
    title([total_name,': foci # profile (from ',num2str(nnz(good_data)),' embryos)'],'Interpreter','none')
    ylabel('# of foci')
    xlabel('AP axis (normalized)')
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Mean plot of active nuclei # vs EL: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(fn(9))
N_data = find(good_data);
for ii = 1:length(N_data)
    i_data = N_data(ii);
    subplot(1,2,1)
        plot(fociN_out{i_data}(:,1),fociN_out{i_data}(:,2),'Color',[0.5,0,0.5],'DisplayName',[out_name{i_data},', cycle ',num2str(cycle_data(i_data)),' P(foci active)'])
        hold on
        
    subplot(1,2,2)
        Itrue = FISHI_out{i_data}(:,1) >= EL_min & FISHI_out{i_data}(:,1) <= EL_max & FISHI_out{i_data}(:,4) >= Nmin;
        temp_all = sum(FISHI_out{i_data}(:,10:13),2);
        plot(FISHI_out{i_data}(Itrue,1),FISHI_out{i_data}(Itrue,10)./temp_all(Itrue),'Color',[0.5,0.5,0.5],'DisplayName',[out_name{i_data},', cycle ',num2str(cycle_data(i_data)),' P(n = 0)'])
        hold on
        plot(FISHI_out{i_data}(Itrue,1),FISHI_out{i_data}(Itrue,11)./temp_all(Itrue),'Color',[0,0,0.5],'DisplayName',[out_name{i_data},', cycle ',num2str(cycle_data(i_data)),' P(n = 1)'])
        plot(FISHI_out{i_data}(Itrue,1),FISHI_out{i_data}(Itrue,12)./temp_all(Itrue),'Color',[0,0.5,0],'DisplayName',[out_name{i_data},', cycle ',num2str(cycle_data(i_data)),' P(n = 2)'])
        plot(FISHI_out{i_data}(Itrue,1),FISHI_out{i_data}(Itrue,13)./temp_all(Itrue),'Color',[0.5,0,0],'DisplayName',[out_name{i_data},', cycle ',num2str(cycle_data(i_data)),' P(n > 2)'])
end

xy_out0 = cell2mat(fociN_out(good_data));
nu_L0 = xy_out0(:,1);
nu_mean0 = xy_out0(:,2);

xy_out0 = cell2mat(FISHI_out(good_data));
nn0 = xy_out0(:,10);
nn1 = xy_out0(:,11);
nn2 = xy_out0(:,12);
nn3 = xy_out0(:,13);
nn_all = nn0+nn1+nn2+nn3;

nu_L = sort(unique(nu_L0));
nu_mean = zeros(size(nu_L));
nu_err = zeros(size(nu_L));
mp0 = zeros(size(nu_L));
mp1 = zeros(size(nu_L));
mp2 = zeros(size(nu_L));
mp3 = zeros(size(nu_L));
ep0 = zeros(size(nu_L));
ep1 = zeros(size(nu_L));
ep2 = zeros(size(nu_L));
ep3 = zeros(size(nu_L));

for I_nu = 1:length(nu_L)
    mean_raw = nu_mean0(nu_L0 == nu_L(I_nu));
    I_good = ~isnan(mean_raw);
    nu_mean(I_nu) = mean(mean_raw(I_good));
    nu_err(I_nu) = std0(mean_raw(I_good));
    
    I_good = (nu_L0 == nu_L(I_nu)) & (nn_all >= Nmin);
    mp0(I_nu) = mean(nn0(I_good)./nn_all(I_good));
    mp1(I_nu) = mean(nn1(I_good)./nn_all(I_good));
    mp2(I_nu) = mean(nn2(I_good)./nn_all(I_good));
    mp3(I_nu) = mean(nn3(I_good)./nn_all(I_good));
    ep0(I_nu) = std0(nn0(I_good)./nn_all(I_good));
    ep1(I_nu) = std0(nn1(I_good)./nn_all(I_good));
    ep2(I_nu) = std0(nn2(I_good)./nn_all(I_good));
    ep3(I_nu) = std0(nn3(I_good)./nn_all(I_good));
end

figure(fn(9))
subplot(1,2,1)
    errorbar(nu_L,nu_mean,nu_err,'m','DisplayName','P(foci active) Mean profile')
    title([total_name,': Active foci profile (from ',num2str(nnz(good_data)),' embryos)'],'Interpreter','none')
    ylabel('Frequency')
    xlabel('AP axis (normalized)')
subplot(1,2,2)
    errorbar(nu_L,mp0,ep0,'k','DisplayName','P(n = 0) Mean profile')
    errorbar(nu_L,mp1,ep1,'b','DisplayName','P(n = 1) Mean profile')
    errorbar(nu_L,mp2,ep2,'g','DisplayName','P(n = 2) Mean profile')
    errorbar(nu_L,mp3,ep3,'r','DisplayName','P(n > 2) Mean profile')
    title([total_name,': Active nuclei profile (from ',num2str(nnz(good_data)),' embryos)'],'Interpreter','none')
    ylabel('Frequency')
    xlabel('AP axis (normalized)')
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Mean plot of I regulation: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_data = find(good_data);
pp_out0 = zeros(length(N_data),length(pro_center0));
yk_out0 = zeros(length(N_data),length(pro_center0));
RNAIn_out = zeros(length(N_data),length(pro_center0));
for ii = 1:length(N_data)
    i_data = N_data(ii);
    Itrue = RNAI_out(i_data,1+2*length(pro_center0):3*length(pro_center0)) >= Nmin;
    r00 = RNAI_out(i_data,1:length(pro_center0));
    figure(fn(4))
    subplot(1,2,1)
        RNAIn_out(ii,:) = RNAI_out(i_data,1:length(pro_center0));%./max(r00(Itrue));
        plot(RNAI_out(i_data,1+9*length(pro_center0):10*length(pro_center0)),RNAIn_out(ii,:),'Color',color_code(mod(ii-1,size(color_code,1))+1,:),'DisplayName',[out_name{i_data},', cycle ',num2str(cycle_data(i_data))])
        hold on
    subplot(1,2,2)
        pp = RNAI_out(i_data,1+9*length(pro_center0):10*length(pro_center0));
        yk = RNAI_out(i_data,1+10*length(pro_center0):11*length(pro_center0));
        
        if ismember(str2num(total_name(end-2:end)),cycle_diff)
            pro_th0 = pro_th(2);
        else
            pro_th0 = pro_th(1);
        end
    
        II0 = yk > 0 & isfinite(yk) & pp >= pro_th0;
%         ppc = geomean(pp(II0));
%         ykc = geomean(yk(II0));
        ppc = HH_I(i_data,2);
        ykc = 1;
        pp = pp/ppc*1e-8;
        yk = yk/ykc;
        loglog(pp,yk,'Color',color_code(mod(ii-1,size(color_code,1))+1,:),'DisplayName',[out_name{i_data},', cycle ',num2str(cycle_data(i_data))])
%         pp(~II0) = nan;
        yk(~II0) = nan;
        pp_out0(ii,:) = pp;
        yk_out0(ii,:) = yk;
        hold on
    figure(fn(8))
        mean_temp = RNAI_out(i_data,1:length(pro_center0));
        std_temp = RNAI_out(i_data,1+length(pro_center0):2*length(pro_center0));
        corr_temp = RNAI_out(i_data,1+3*length(pro_center0):4*length(pro_center0));
        corrn_temp = RNAI_out(i_data,1+4*length(pro_center0):5*length(pro_center0));
        vex_temp = RNAI_out(i_data,1+5*length(pro_center0):6*length(pro_center0));
        corrb_temp = RNAI_out(i_data,1+6*length(pro_center0):7*length(pro_center0));
        corrnb_temp = RNAI_out(i_data,1+7*length(pro_center0):8*length(pro_center0));
        p0_temp = RNAI_out(i_data,1+8*length(pro_center0):9*length(pro_center0));
    subplot(1,2,1)
        plot(mean_temp(Itrue),std_temp(Itrue).^2./mean_temp(Itrue),'.','Color',color_code(mod(ii-1,size(color_code,1))+1,:),'DisplayName',[out_name{i_data},', cycle ',num2str(cycle_data(i_data))])
        hold on
        plot(mean_temp(Itrue),vex_temp(Itrue)./mean_temp(Itrue),'.','Color',color_code(mod(ii-1,size(color_code,1))+1,:),'DisplayName',[out_name{i_data},', cycle ',num2str(cycle_data(i_data)),' EX'])
    subplot(1,2,2)
        plot(mean_temp(Itrue),corr_temp(Itrue),'.','Color',color_code(mod(ii-1,size(color_code,1))+1,:),'DisplayName',[out_name{i_data},', cycle ',num2str(cycle_data(i_data))])
        hold on
        plot(mean_temp(Itrue),corrn_temp(Itrue),'.','Color',color_code(mod(ii-1,size(color_code,1))+1,:)/2,'DisplayName',[out_name{i_data},', cycle ',num2str(cycle_data(i_data)),' control'])
        plot(mean_temp(Itrue),vex_temp(Itrue)./std_temp(Itrue).^2,'.','Color',color_code(mod(ii-1,size(color_code,1))+1,:)/4,'DisplayName',[out_name{i_data},', cycle ',num2str(cycle_data(i_data)),' EX'])
        plot(mean_temp(Itrue),corrb_temp(Itrue),'.','Color',color_code(mod(ii-1,size(color_code,1))+1,:),'DisplayName',[out_name{i_data},', cycle ',num2str(cycle_data(i_data)),' digital'])
        plot(mean_temp(Itrue),corrnb_temp(Itrue),'.','Color',color_code(mod(ii-1,size(color_code,1))+1,:)/2,'DisplayName',[out_name{i_data},', cycle ',num2str(cycle_data(i_data)),' digital control'])
        plot(mean_temp(Itrue),p0_temp(Itrue),'-','Color',color_code(mod(ii-1,size(color_code,1))+1,:)/4,'DisplayName',[out_name{i_data},', cycle ',num2str(cycle_data(i_data)),' p(0)'])
end

RNAI_true = RNAI_out(good_data,1+2*length(pro_center0):3*length(pro_center0)) >= Nmin;
pro_bin0 = RNAI_out(good_data,1+9*length(pro_center0):10*length(pro_center0));
pro_bin0 = pro_bin0(RNAI_true);
pro_bin0 = pro_bin0(:);
RNAI_out0 = RNAI_out(good_data,1:length(pro_center0));
RNAI_value = RNAIn_out(RNAI_true);
RNAI_value = RNAI_value(:);
pp_value = pp_out0(RNAI_true);
pp_value = pp_value(:);
yk_value = yk_out0(RNAI_true);
yk_value = yk_value(:);
var_out0 = RNAI_out(good_data,1+length(pro_center0):2*length(pro_center0)).^2./RNAI_out(good_data,1:length(pro_center0));
var_value = RNAI_out(good_data,1+length(pro_center0):2*length(pro_center0)).^2;
var_value = var_value(RNAI_true);
var_value = var_value(:);
corr_out0 = RNAI_out(good_data,1+3*length(pro_center0):4*length(pro_center0));
% % corr_value = corr_out0;
% % corr_value(~RNAI_true) = 0;
corrn_out0 = RNAI_out(good_data,1+4*length(pro_center0):5*length(pro_center0));
% % corrn_value = corrn_out0;
% % corrn_value(~RNAI_true) = 0;
vex_out0 = RNAI_out(good_data,1+5*length(pro_center0):6*length(pro_center0))./RNAI_out(good_data,1:length(pro_center0));
cex_out0 = RNAI_out(good_data,1+5*length(pro_center0):6*length(pro_center0))./RNAI_out(good_data,1+length(pro_center0):2*length(pro_center0)).^2;
corrb_out0 = RNAI_out(good_data,1+6*length(pro_center0):7*length(pro_center0));
corrnb_out0 = RNAI_out(good_data,1+7*length(pro_center0):8*length(pro_center0));
p0_out0 = RNAI_out(good_data,1+8*length(pro_center0):9*length(pro_center0));

[pro_mean,RNAI_mean,pro_err,RNAI_std0] = equal_bin(pro_bin0,RNAI_value,Nbin);
[pro_mean,var_mean,pro_err,var_std0] = equal_bin(pro_bin0,var_value,Nbin);
[pro_mean0,yk_mean,pro_err0,yk_std0] = equal_bin(pp_value(yk_value > 0 & isfinite(yk_value)),yk_value(yk_value > 0 & isfinite(yk_value)),Nbin);
% % % RNAI_mean = sum(RNAI_value,1)./sum(RNAI_true,1);
% % % RNAI_std0 = sqrt((sum(RNAI_value.^2,1)./sum(RNAI_true,1)-RNAI_mean.^2)./sum(RNAI_true,1));
% % % var_mean = sum(var_value,1)./sum(RNAI_true,1);
% % % var_std0 = sqrt((sum(var_value.^2,1)./sum(RNAI_true,1)-var_mean.^2)./sum(RNAI_true,1));
% % corr_mean = sum(corr_value,1)./sum(RNAI_true,1);
% % corr_std0 = sqrt((sum(corr_value.^2,1)./sum(RNAI_true,1)-corr_mean.^2)./sum(RNAI_true,1));
good_bin = (~isnan(RNAI_std0)) & (RNAI_std0 > 0) & (~isnan(var_std0)) & (var_std0 > 0);

RNAI_value1 = RNAI_out0(RNAI_true);
var_value1 = var_out0(RNAI_true);
corr_value1 = corr_out0(RNAI_true);
corrn_value1 = corrn_out0(RNAI_true);
corr_value1(isnan(corr_value1)) = 0;
corrn_value1(isnan(corrn_value1)) = 0;
vex_value1 = vex_out0(RNAI_true);
cex_value1 = cex_out0(RNAI_true);
corrb_value1 = corrb_out0(RNAI_true);
corrnb_value1 = corrnb_out0(RNAI_true);
corrb_value1(isnan(corrb_value1)) = 0;
corrnb_value1(isnan(corrnb_value1)) = 0;
p0_value1 = p0_out0(RNAI_true);
[RNAI_mean1,var_mean1,RNAI_err1,var_err1] = equal_bin(RNAI_value1(:),var_value1(:),Nbin);
[RNAI_mean1,corr_mean1,RNAI_err1,corr_err1] = equal_bin(RNAI_value1(:),corr_value1(:),Nbin);
[RNAI_mean1,corrn_mean1,RNAI_err1,corrn_err1] = equal_bin(RNAI_value1(:),corrn_value1(:),Nbin);
[RNAI_mean1,vex_mean1,RNAI_err1,vex_err1] = equal_bin(RNAI_value1(:),vex_value1(:),Nbin);
[RNAI_mean1,cex_mean1,RNAI_err1,cex_err1] = equal_bin(RNAI_value1(:),cex_value1(:),Nbin);
[RNAI_mean1,corrb_mean1,RNAI_err1,corrb_err1] = equal_bin(RNAI_value1(:),corrb_value1(:),Nbin);
[RNAI_mean1,corrnb_mean1,RNAI_err1,corrnb_err1] = equal_bin(RNAI_value1(:),corrnb_value1(:),Nbin);
[RNAI_mean1,p0_mean1,RNAI_err1,p0_err1] = equal_bin(RNAI_value1(:),p0_value1(:),Nbin);

try
    [alpha1,r,J,cov,~] = nlinfit(RNAI_mean1,var_mean1,fth,alpha0);
    koff0 = alpha1(1)*(alpha1(2)./2./RNAI_mean1-1);
    fp0th = P0_2S(alpha1(1),koff0,alpha1(2));
catch
    alpha1 = nan(size(alpha0));
    fp0th = nan(size(RNAI_mean1));
end

figure(fn(4))
subplot(1,2,1)
    errorbar(pro_mean,RNAI_mean,RNAI_std0,'.k','DisplayName','Mean profile')
    hold on

    %%% Nonlinear fitting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [v_max,I_max] = max(RNAI_mean);
    pro_I = find(pro_mean <= pro_mean(I_max));
    [~,I_middle] = min(abs(RNAI_mean(pro_I)-max(RNAI_mean)/2));
    beta0 = [h0,pro_mean(pro_I(I_middle)),v_max,0];
    try
        [beta,r,~,~,~] = nlinfit(pro_mean(good_bin),RNAI_mean(good_bin),@Hill,beta0);
        plot(pro_mean,Hill(beta,pro_mean),'r','DisplayName','Fitted curve for Mean profile')
        H_I = beta;
    catch err
        plot([],[],'r','DisplayName','Fitted curve for Mean profile')
        H_I = nan(size(beta0));
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    hold off
    title([total_name,': ',p_name,' transcription (I) Regulation Curve (from ',num2str(nnz(good_data)),' embryos), h = ',num2str(H_I(1))],'Interpreter','none')
    xlabel([p_name,' concentration (M)'])
    ylabel('RNA level (#)')
%     legend('show')
%     legend('hide')

subplot(1,2,2)
    errorbar(pro_mean0,yk_mean,yk_std0,'.k','DisplayName','Mean profile')
    hold on

    %%% Nonlinear fitting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    pro_mean00 = pro_mean0(N01:end-N01+1);
    yk_mean0 = yk_mean(N01:end-N01+1);
%     if ismember(str2num(total_name(end-2:end)),cycle_diff)
%         pro_th0 = pro_th(2);
%     else
%         pro_th0 = pro_th(1);
%     end
%     pk = polyfit(log(pro_mean00(pro_mean00 > pro_th0 & yk_mean0 > 0 & isfinite(yk_mean0))),log(yk_mean0(pro_mean00 > pro_th0 & yk_mean0 > 0 & isfinite(yk_mean0))),1);
    pk = polyfit(log(pro_mean00(yk_mean0 > 0 & isfinite(yk_mean0))),log(yk_mean0(yk_mean0 > 0 & isfinite(yk_mean0))),1);
%     pk = polyfit(log(pro_mean0(pro_mean0 > pro_th & yk_mean > 0 & isfinite(yk_mean))),log(yk_mean(pro_mean0 > pro_th & yk_mean > 0 & isfinite(yk_mean))),1);
    loglog(pro_mean0,pro_mean0.^pk(1).*exp(pk(2)),'r','DisplayName','Fit')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    hold off
    title([total_name,': ',p_name,' vs koff/kon (from ',num2str(nnz(good_data)),' embryos), h = ',num2str(-pk(1))],'Interpreter','none')
    xlabel([p_name,' concentration (M)'])
    ylabel('koff/kon')
%     legend('show')


figure(fn(8))
subplot(1,2,1)
    errorbarxy(RNAI_mean1,var_mean1,RNAI_err1,var_err1,RNAI_err1,var_err1,'r','r',[],'DisplayName','Mean curve')
    hold on
    if all(~isnan(alpha1))
        plot(RNAI_mean1,fth(alpha1,RNAI_mean1),'m-','DisplayName',['Fit: k_o_n*T = ',num2str(alpha1(1),'%4.3f'),', k_T_X*T = ',num2str(alpha1(2),'%4.3f')]);
    else
        plot([],[],'m-','DisplayName',['Fit: k_o_n*T = ',num2str(alpha1(1),'%4.3f'),', k_T_X*T = ',num2str(alpha1(2),'%4.3f')]);
    end
    errorbarxy(RNAI_mean1,vex_mean1,RNAI_err1,vex_err1,RNAI_err1,vex_err1,'b','b',[],'DisplayName','Mean ex curve')
    title([total_name,': TX mean vs variance (from ',num2str(nnz(good_data)),' embryos)',char(10),'Fit: kon*T = ',num2str(alpha1(1),'%4.3f'),', kTX*T = ',num2str(alpha1(2),'%4.3f')],'Interpreter','none')
    xlabel('Mean # mRNA (concentration binned)')
    ylabel('Var/mean (concentration binned)')
    
subplot(1,2,2)
    errorbarxy(RNAI_mean1,corr_mean1,RNAI_err1,corr_err1,RNAI_err1,corr_err1,'r','r','.','DisplayName','Mean curve')
    hold on
    errorbarxy(RNAI_mean1,corrn_mean1,RNAI_err1,corrn_err1,RNAI_err1,corrn_err1,'k','k','.','DisplayName','Mean control curve')
    errorbarxy(RNAI_mean1,cex_mean1,RNAI_err1,cex_err1,RNAI_err1,cex_err1,'b','b','.','DisplayName','Mean ex curve')
    errorbarxy(RNAI_mean1,corrb_mean1,RNAI_err1,corrb_err1,RNAI_err1,corrb_err1,[0.5,0,0],[0.5,0,0],'.','DisplayName','Mean digital curve')
    errorbarxy(RNAI_mean1,corrnb_mean1,RNAI_err1,corrnb_err1,RNAI_err1,corrnb_err1,[0.5,0.5,0.5],[0.5,0.5,0.5],'.','DisplayName','Mean digital control curve')
    errorbarxy(RNAI_mean1,p0_mean1,RNAI_err1,p0_err1,RNAI_err1,p0_err1,[0.5,0,0.5],[0.5,0,0.5],'-','DisplayName','Mean p(0) curve')
    errorbarxy(RNAI_mean1,fp0th,RNAI_err1,zeros(size(RNAI_err1)),RNAI_err1,zeros(size(RNAI_err1)),'g','g','-','DisplayName','Theoretical mean p(0) curve')
    title([total_name,': TX mean vs correlation (from ',num2str(nnz(good_data)),' embryos)'],'Interpreter','none')
    xlabel('Mean # mRNA (concentration binned)')
    ylabel('Correlation coefficient (concentration binned)')
    
if all(~isnan(alpha1))
    ud_alpha1 = nlparci(alpha1,r,'covar',cov);
    err_alpha1 = alpha1-ud_alpha1(:,1)';
else
    err_alpha1 = nan(size(alpha0));
end
para_Bcd = [nan(1,2),alpha1,nan(1,2),-pk(1),(exp(pk(2))).^(-1./pk(1));nan(1,2),err_alpha1,nan(1,2),nan(1,2)];
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Mean plot of N regulation: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(fn(5))
N_data = find(good_data);
for ii = 1:length(N_data)
    i_data = N_data(ii);
    plot(pro_center0,RNAN_out(i_data,:),'Color',color_code(mod(ii-1,size(color_code,1))+1,:),'DisplayName',[out_name{i_data},', cycle ',num2str(cycle_data(i_data))])
    hold on
end

RNAN_out0 = RNAN_out(good_data,:);
RNAN_value = RNAN_out0;
RNAN_value(isnan(RNAN_out0)) = 0;
RNAN_mean = sum(RNAN_value,1)./sum(~isnan(RNAN_out0),1);
RNAN_std0 = sqrt((sum(RNAN_value.^2,1)./sum(~isnan(RNAN_out0),1)-RNAN_mean.^2)./sum(~isnan(RNAN_out0),1));
good_bin = (~isnan(RNAN_std0)) & (RNAN_std0 > 0);

figure(fn(5))
    errorbar(pro_center0,RNAN_mean,RNAN_std0,'.k','DisplayName','Mean profile')
    hold on

    %%% Nonlinear fitting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [v_max,I_max] = max(RNAN_mean);
    pro_I = find(pro_center0 <= pro_center0(I_max));
    [~,I_middle] = min(abs(RNAN_mean(pro_I)-max(RNAN_mean)/2));
    beta0 = [h0,pro_center0(pro_I(I_middle)),v_max,0];
    try
        [beta,r,~,~,~] = nlinfit(pro_center0(good_bin),RNAN_mean(good_bin),@Hill,beta0);
        plot(pro_center0,Hill(beta,pro_center0),'r','DisplayName','Fitted curve for Mean profile')
        H_N = beta;
    catch err
        plot([],[],'r','DisplayName','Fitted curve for Mean profile')
        H_N = nan(size(beta0));
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    hold off
    title([total_name,': ',p_name,' Transcription (N) Regulation Curve (from ',num2str(nnz(good_data)),' embryos), h = ',num2str(H_N(1))],'Interpreter','none')
    xlabel([p_name,' concentration (M)'])
    ylabel('# of foci per nucleus')
%     legend('show')
%     legend('hide')
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Mean plot of null rate: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
null_out0 = null_out(good_data,:);
null_value = null_out0;
null_value(isnan(null_out0)) = 0;
null_mean = sum(null_value,1)./sum(~isnan(null_out0),1);
null_std = sqrt(sum(null_value.^2,1)./sum(~isnan(null_out0),1)-null_mean.^2);
dnull_mean = [diff(null_mean),nan]/(pro_center(2)-pro_center(1));
dnull_std = [sqrt(null_std(1:end-1).^2+null_std(2:end).^2)/2,nan];
pshow_max = 5e-8;

figure(fn(6))
h_axes1 = subplot(1,2,1);
h_axes3 = subplot(1,2,2);
axes(h_axes1)
%     h_axes1 = axes;
    h_err1 = errorbar(pro_center,null_mean,null_std,'b');
    box off
    hold on

    %%% Nonlinear fitting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    v_max = max(null_mean);
    pro_I = find(pro_center <= pshow_max);
    [~,I_middle] = min(abs(null_mean(pro_I)-max(null_mean)/2));
    beta0 = [h0,pro_center(pro_I(I_middle)),v_max,0];
    try
        [beta,r,~,~,~] = nlinfit(pro_center(pro_I),null_mean(pro_I),@Hill,beta0);
        h_err11 = plot(pro_center(pro_I),Hill(beta,pro_center(pro_I)),'r');
        H_A = beta;
    catch err
        h_err11 = plot([],[],'r');
        H_A = nan(size(beta0));
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    xlabel([p_name,' concentration (M)']);
    ylabel('Active rate')
    hold off

    h_axes2 = axes;
    h_err2 = errorbar(pro_center,dnull_mean,dnull_std,'g');
    set(gca,'XAxisLocation','top','YAxisLocation','right','Color','none','Position',get(h_axes1,'Position'),'XTick',[])
    ylabel('Diff active rate (1/M)')
    box off

    legend([h_err1,h_err11,h_err2],{'TX site active rate',['Fitted curve: h = ',num2str(H_A(1))],'diff(TX site silence rate)'})
    title([total_name,': ',p_name,' Transcription site active rate'],'Interpreter','none')

    
% subplot(1,2,2)
axes(h_axes3)
    N_data = find(good_data);
    for ii = 1:length(N_data)
        i_data = N_data(ii);
        
        Itrue = RNAI_out(i_data,1+2*length(pro_center0):3*length(pro_center0)) >= Nmin;
        nn0_temp = RNAI_out(i_data,1+11*length(pro_center0):12*length(pro_center0));
        nn1_temp = RNAI_out(i_data,1+12*length(pro_center0):13*length(pro_center0));
        nn2_temp = RNAI_out(i_data,1+13*length(pro_center0):14*length(pro_center0));
        nn3_temp = RNAI_out(i_data,1+14*length(pro_center0):15*length(pro_center0));
        temp_all = nn0_temp+nn1_temp+nn2_temp+nn3_temp;

        plot(pro_center0(Itrue),nn0_temp(Itrue)./temp_all(Itrue),'Color',[0.5,0.5,0.5],'DisplayName',[out_name{i_data},', cycle ',num2str(cycle_data(i_data)),' P(n = 0)'])
        hold on
        plot(pro_center0(Itrue),nn1_temp(Itrue)./temp_all(Itrue),'Color',[0,0,0.5],'DisplayName',[out_name{i_data},', cycle ',num2str(cycle_data(i_data)),' P(n = 1)'])
        plot(pro_center0(Itrue),nn2_temp(Itrue)./temp_all(Itrue),'Color',[0,0.5,0],'DisplayName',[out_name{i_data},', cycle ',num2str(cycle_data(i_data)),' P(n = 2)'])
        plot(pro_center0(Itrue),nn3_temp(Itrue)./temp_all(Itrue),'Color',[0.5,0,0],'DisplayName',[out_name{i_data},', cycle ',num2str(cycle_data(i_data)),' P(n > 2)'])
    end

    RNAI_true = RNAI_out(good_data,1+2*length(pro_center0):3*length(pro_center0)) >= Nmin;
    RNAI_all = RNAI_out(good_data,1+11*length(pro_center0):12*length(pro_center0))+RNAI_out(good_data,1+12*length(pro_center0):13*length(pro_center0))+RNAI_out(good_data,1+13*length(pro_center0):14*length(pro_center0))+RNAI_out(good_data,1+14*length(pro_center0):15*length(pro_center0));

    RNAI_p0 = RNAI_out(good_data,1+11*length(pro_center0):12*length(pro_center0))./RNAI_all;
    RNAI_p0(~RNAI_true) = 0;
    RNAI_p1 = RNAI_out(good_data,1+12*length(pro_center0):13*length(pro_center0))./RNAI_all;
    RNAI_p1(~RNAI_true) = 0;
    RNAI_p2 = RNAI_out(good_data,1+13*length(pro_center0):14*length(pro_center0))./RNAI_all;
    RNAI_p2(~RNAI_true) = 0;
    RNAI_p3 = RNAI_out(good_data,1+14*length(pro_center0):15*length(pro_center0))./RNAI_all;
    RNAI_p3(~RNAI_true) = 0;

    mp0 = sum(RNAI_p0,1)./sum(RNAI_true,1);
    ep0 = sqrt((sum(RNAI_p0.^2,1)./sum(RNAI_true,1)-mp0.^2)./sum(RNAI_true,1));
    mp1 = sum(RNAI_p1,1)./sum(RNAI_true,1);
    ep1 = sqrt((sum(RNAI_p1.^2,1)./sum(RNAI_true,1)-mp1.^2)./sum(RNAI_true,1));
    mp2 = sum(RNAI_p2,1)./sum(RNAI_true,1);
    ep2 = sqrt((sum(RNAI_p2.^2,1)./sum(RNAI_true,1)-mp2.^2)./sum(RNAI_true,1));
    mp3 = sum(RNAI_p3,1)./sum(RNAI_true,1);
    ep3 = sqrt((sum(RNAI_p3.^2,1)./sum(RNAI_true,1)-mp3.^2)./sum(RNAI_true,1));

    errorbar(pro_center,mp0,ep0,'k','DisplayName','P(n = 0) Mean profile')
    errorbar(pro_center,mp1,ep1,'b','DisplayName','P(n = 1) Mean profile')
    errorbar(pro_center,mp2,ep2,'g','DisplayName','P(n = 2) Mean profile')
    errorbar(pro_center,mp3,ep3,'r','DisplayName','P(n > 2) Mean profile')
    title([total_name,': Active nuclei profile (from ',num2str(nnz(good_data)),' embryos)'],'Interpreter','none')
    ylabel('Frequency')
    xlabel([p_name,' concentration (M)'])
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

function y = Hill(beta,x)
 y = beta(3)*x.^beta(1)./(x.^beta(1)+beta(2)^beta(1))+beta(4);
end