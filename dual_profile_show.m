function [H_I,H_N,pro_center0,RNAI_out,RNAN_out,H_A,pro_center,null_rate,fit_out] = dual_profile_show(nucleus_protein_profile_ab,nucleus_RNA_profile,foci_RNA_profile,image_folder,N_cycle,sub_pos,color_code,varargin)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A function to re-plot transcription regulation curve and fit it to a hill function %%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Parameter setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_pro = 50;
w_pro = 0.05;
reg_min = 0.20;%0.25;
reg_max = 0.7;%0.75;
single_add = '_new';
pmin = 0;
pmax = 5e-8;%max(nucleus_protein_profile_ab(:,2));
pbin = (pmax-pmin)/N_pro;
pwindow = (pmax-pmin)*w_pro;
pro_center0 = (pmin+pbin/2):pbin:(pmax-pbin/2);
pro_th = [1.5e-9,1e-9];
cycle_diff = 14;
N01 = 4;

p_profile = nucleus_protein_profile_ab(:,2);
reg_I = (nucleus_protein_profile_ab(:,1) >= reg_min) & (nucleus_protein_profile_ab(:,1) <= reg_max) & nucleus_protein_profile_ab(:,3);

if ~isempty(varargin)
    p_name = varargin{1};
else
    p_name = 'Protein';
end

if length(varargin) > 1
    fn = varargin{2};
else
    fn = [12,120,15,150,17,170,18,91,910,34,340,51,510];
end

if length(varargin) > 2
    pmin2 = varargin{3}(1);
    pmax2 = varargin{3}(2);
else
    pmin2 = 1e-8;
    pmax2 = 1.5e-8;
end
Cbin = 0:5:500;
Nmin = 10;
Nmean = 4;
d_fit = 3;
fth = @(alpha1,x) 2/3+8/3/alpha1(1).*x-(8/alpha1(1)^2/alpha1(2)+16/3/alpha1(1)/alpha1(2)).*x.^2+16/alpha1(1)^2/alpha1(2)^2.*x.^3+64/alpha1(1)^4/alpha1(2)^3.*x.^4-128/alpha1(1)^4/alpha1(2)^4.*x.^5+(-32/alpha1(1)^3/alpha1(2)^2.*x.^3+64*(1/alpha1(1)^3/alpha1(2)^3-1/alpha1(1)^4/alpha1(2)^3).*x.^4+128/alpha1(1)^4/alpha1(2)^4.*x.^5).*exp(-1/2*alpha1(1)*alpha1(2)./x);
alpha0 = [0.2,65];
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Re-plot the transcription regulaton curve using total nuclear foci intensity %%
RNA_n = nan(size(pro_center0));   %%% number of nuclei in the bin
RNA_Bcd = nan(size(pro_center0));   %%% mean [Bcd] in the bin
RNA_mean = nan(size(pro_center0));   %%% mean foci intensity in the bin
RNA_std = nan(size(pro_center0));   %%% err foci intensity in the bin
RNA_std2 = nan(size(pro_center0));   %%% std foci intensity in the bin
RNA_corr = nan(size(pro_center0));   %%% correlation coefficient of foci intensity in the bin
RNA_corrn = nan(size(pro_center0));   %%% correlation coefficient control of foci intensity in the bin
RNA_corrb = nan(size(pro_center0));   %%% digital correlation coefficient of foci intensity in the bin
RNA_corrnb = nan(size(pro_center0));   %%% digital correlation coefficient control of foci intensity in the bin
RNA_p0 = nan(size(pro_center0));   %%% probability of silent foci in the bin
[nucleus_RNA_profile0,ind_foci] = foci_info(nucleus_RNA_profile,foci_RNA_profile);
p_profile0 = p_profile(ind_foci);
reg_I0 = (nucleus_protein_profile_ab(ind_foci,1) >= reg_min) & (nucleus_protein_profile_ab(ind_foci,1) <= reg_max) & nucleus_protein_profile_ab(ind_foci,3) & (nucleus_protein_profile_ab(ind_foci,2) >= pmin) & (nucleus_protein_profile_ab(ind_foci,2) <= pmax);

[pro_bin0,bin_min0,bin_max0] = seg_eqN(p_profile0,length(pro_center0),w_pro/2);

for Ip = 1:length(pro_center0)
    Itrue = reg_I0&(p_profile0 >= bin_min0(Ip))&(p_profile0 <= bin_max0(Ip));
    RNA_Bcd(Ip) = mean(p_profile0(Itrue));
    RNA_n(Ip) = nnz(unique(nucleus_RNA_profile0(Itrue,2)));
    if RNA_n(Ip) >= Nmin
        RNA_mean(Ip) = mean(nucleus_RNA_profile0(Itrue,4));
        RNA_std(Ip) = std0(nucleus_RNA_profile0(Itrue,4));
        RNA_std2(Ip) = std(nucleus_RNA_profile0(Itrue,4));
        [RNA_corr(Ip),RNA_corrn(Ip),RNA_corrb(Ip),RNA_corrnb(Ip)] = inu_corr(ind_foci(Itrue),nucleus_RNA_profile0(Itrue,:));
        RNA_p0(Ip) = mean(nucleus_RNA_profile0(Itrue,4) == 0);
%         RNA_p0(Ip) = mean(nucleus_RNA_profile(unique(ind_foci(Itrue)),4) == 0);
    end
end
good_bin = (~isnan(RNA_std)) & (RNA_std > 0);

% RNAI_out = nan(size(pro_center0));
% RNAI_out(good_bin) = RNA_mean(good_bin);
% RNAI_out2 = nan(size(pro_center0));
% RNAI_out2(good_bin) = RNA_std2(good_bin);
% RNAI_out3 = nan(size(pro_center0));
% RNAI_out3(good_bin) = RNA_corr(good_bin);

RNA_vex = nan(size(pro_center0));   %%% Deviation of mean foci intensity in the bin
for I_bin = 1:length(pro_center0)
    Itrue = reg_I0&(p_profile0 >= bin_min0(I_bin))&(p_profile0 <= bin_max0(I_bin));
    if RNA_n(I_bin) >= Nmin
        IImin = max(1,I_bin-d_fit);
        IImax = min(I_bin+d_fit,length(pro_center0));
        p0 = polyfit(RNA_Bcd(IImin:IImax),RNA_mean(IImin:IImax),2);
        dev0 = p0(1)*p_profile0(Itrue).^2+p0(2)*p_profile0(Itrue)-p0(1)*RNA_Bcd(I_bin)^2-p0(2)*RNA_Bcd(I_bin);
        RNA_vex(I_bin) = sum(dev0.^2)/(nnz(Itrue)-1);
    end
end


if nnz(good_bin) > 5% && sum(nucleus_RNA_profile0((reg_I),3)) > 20 && pro_center0(end) > 10000 && isempty(strfind(image_folder((end-4):end), single_add))
    %%% Nonlinear fitting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [v_max,I_max] = max(RNA_mean);
    pro_I = find(pro_bin0 <= pro_bin0(I_max));
    [~,I_middle] = min(abs(RNA_mean(pro_I)-max(RNA_mean)/2));
    beta0 = [4,pro_bin0(pro_I(I_middle)),v_max,0];
    try
        [beta,r,~,~,~] = nlinfit(pro_bin0(good_bin),RNA_mean(good_bin),@Hill,beta0);
        H_I = beta;
        yk = (beta(3))./(RNA_mean-beta(4))-1;
    catch err
        H_I = nan(size(beta0));
        RNA_sort = sort(RNA_mean(good_bin));
        yk = (mean(RNA_sort(end-Nmean+1:end))-mean(RNA_sort(1:Nmean)))./(RNA_mean-mean(RNA_sort(1:Nmean)))-1;
    end
    
    pro_bin00 = pro_bin0(N01:end-N01+1);
    yk0 = yk(N01:end-N01+1);
    good_bin0 = good_bin(N01:end-N01+1);
    if ismember(N_cycle,cycle_diff)
        pro_th0 = pro_th(2);
    else
        pro_th0 = pro_th(1);
    end
%     yk = (max(RNA_mean(good_bin))-min(RNA_mean(good_bin)))./(RNA_mean-min(RNA_mean(good_bin)))-1;
    pk = polyfit(log(pro_bin00(good_bin0 & pro_bin00 > pro_th0 & yk0 > 0 & isfinite(yk0))),log(yk0(good_bin0 & pro_bin00 > pro_th0 & yk0 > 0 & isfinite(yk0))),1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    plot([],[],'r')
    H_I = nan(1,4);
    
    yk = nan(size(RNA_mean));
    pk = nan(1,2);
end
RNAI_out = [RNA_mean,RNA_std2,RNA_n,RNA_corr,RNA_corrn,RNA_vex,RNA_corrb,RNA_corrnb,RNA_p0,pro_bin0,yk];


if fn(1)
    figure(fn(1))
    subaxis(sub_pos(1),sub_pos(2),sub_pos(3), 'Spacing', 0, 'PaddingRight',0.01, 'PaddingLeft',0.02, 'PaddingTop',0.045, 'PaddingBottom',0.04, 'Margin', 0.01);
        plot(p_profile0(reg_I0),nucleus_RNA_profile0((reg_I0),4),'b.')
        hold on
        errorbar(pro_bin0(good_bin),RNA_mean(good_bin),RNA_std(good_bin),'.k')
        if all(~isnan(H_I))
            plot(pro_bin0(good_bin),Hill(H_I,pro_bin0(good_bin)),'r')
        else
            plot([],[],'r')
        end
        hold off
        title([image_folder,', cycle = ',num2str(N_cycle),char(10),', h = ',num2str(H_I(1))],'Interpreter','none')
        xlabel([p_name,' concentration (M)'])
        ylabel('RNA level (#)')
        legend('Individual nuclei','Averaged profile','Fitting curve for averaged data')
        legend('hide')
end


if fn(12)
    figure(fn(12))
    subaxis(sub_pos(1),sub_pos(2),sub_pos(3), 'Spacing', 0, 'PaddingRight',0.01, 'PaddingLeft',0.02, 'PaddingTop',0.045, 'PaddingBottom',0.04, 'Margin', 0.01);
        loglog(pro_bin0(good_bin),yk(good_bin),'b.')
        hold on
        if all(~isnan(pk))
            loglog(pro_bin0(good_bin),pro_bin0(good_bin).^pk(1).*exp(pk(2)),'r')
        else
            plot([],[],'r')
        end
        hold off
        title([image_folder,', cycle = ',num2str(N_cycle),char(10),', h = ',num2str(-pk(1))],'Interpreter','none')
        xlabel([p_name,' concentration (M)'])
        ylabel('koff/kon')
        legend('Averaged profile','Fitting curve for averaged data')
        legend('hide')
end

fmean = RNA_mean(RNA_n >= Nmin)';
fvar = RNA_std2(RNA_n >= Nmin)';
fcorr = RNA_corr(RNA_n >= Nmin)';
fcorrn = RNA_corrn(RNA_n >= Nmin)';
fcorr(isnan(fcorr)) = 0;
fcorrn(isnan(fcorrn)) = 0;
fvar_ex = RNA_vex(RNA_n >= Nmin)';
fcorrb = RNA_corrb(RNA_n >= Nmin)';
fcorrnb = RNA_corrnb(RNA_n >= Nmin)';
fcorrb(isnan(fcorrb)) = 0;
fcorrnb(isnan(fcorrnb)) = 0;
fp0 = RNA_p0(RNA_n >= Nmin)';
try
    alpha1 = nlinfit(fmean,fvar.^2./fmean,fth,alpha0);
    koff0 = alpha1(1)*(alpha1(2)./2./fmean-1);
    fp0th = P0_2S(alpha1(1),koff0,alpha1(2));
catch
    alpha1 = nan(size(alpha0));
    fp0th = nan(size(fmean));
end

if fn(8)
    figure(fn(8))
    subaxis(sub_pos(1),sub_pos(2),sub_pos(3), 'Spacing', 0, 'PaddingRight',0.01, 'PaddingLeft',0.02, 'PaddingTop',0.045, 'PaddingBottom',0.04, 'Margin', 0.01);
        if all(~isnan(alpha1))
            [AX,H1,H2] = plotyy(fmean,[fvar.^2./fmean,fth(alpha1,fmean),fvar_ex./fmean],fmean,[fcorr,fcorrn,fvar_ex./fvar.^2,fcorrb,fcorrnb,fp0,fp0th]);
        else
            [AX,H1,H2] = plotyy(fmean,[fvar.^2./fmean,nan(size(fmean)),fvar_ex./fmean],fmean,[fcorr,fcorrn,fvar_ex./fvar.^2,fcorrb,fcorrnb,fp0,fp0th]);
        end
        set(H1(1),'LineStyle','none','Marker','.','Color','b','DisplayName','\sigma^2/n vs n data')
        set(H2(1),'LineStyle','none','Marker','o','Color','g','DisplayName','correlation data')
        set(H1(2),'LineStyle','-','Color','r','DisplayName',['Fit: k_o_n*T = ',num2str(alpha1(1),'%4.3f'),', k_T_X*T = ',num2str(alpha1(2),'%4.3f')])
        set(H2(2),'LineStyle','none','Marker','*','Color','c','DisplayName','correlation data control')
        set(H1(3),'LineStyle','none','Marker','.','Color','k','DisplayName','ex\_var/n vs n data')
        set(H2(3),'LineStyle','none','Marker','o','Color',[0.5,0.5,0.5],'DisplayName','ex\_corr data')
        set(H2(4),'LineStyle','none','Marker','*','Color',[0,0.5,0],'DisplayName','digital correlation data')
        set(H2(5),'LineStyle','none','Marker','*','Color',[0,0.5,0.5],'DisplayName','digital correlation data control')
        set(H2(6),'LineStyle','-','Color',[0.5,0,0.5],'DisplayName','probability of silent foci')
        set(H2(7),'LineStyle','-','Color','r','DisplayName','theoretical probability of silent foci')
        legend('show')
        xlabel(AX(1),'Mean # mRNA (concentration binned)')
        ylabel(AX(1),'Var/mean (concentration binned)')
        ylabel(AX(2),'Correlation coefficient (concentration binned)')
        title([image_folder,', cycle = ',num2str(N_cycle)],'Interpreter','none')
end
    
temp_data = nucleus_RNA_profile0(reg_I0&(p_profile0 >= pmin2)&(p_profile0 <= pmax2),4);
nRNA = hist(temp_data,Cbin);
phat = gamfit(temp_data);
    
if fn(10)
    figure(fn(10))
    subaxis(sub_pos(1),sub_pos(2),sub_pos(3), 'Spacing', 0, 'PaddingRight',0.01, 'PaddingLeft',0.02, 'PaddingTop',0.045, 'PaddingBottom',0.04, 'Margin', 0.01);
        bar(Cbin,nRNA/sum(nRNA)/(Cbin(2)-Cbin(1)));
        hold on
        plot(Cbin,gampdf(Cbin,phat(1),phat(2)),'r-')

        title([image_folder,', cycle = ',num2str(N_cycle),char(10),'([Bcd] = ',num2str(pmin2/1e-9,'%4.2f'),'nM - ',num2str(pmax2/1e-9,'%4.2f'),'nM)'],'Interpreter','none')
        xlabel('# mRNA')
        ylabel('Frequency')
        legend('Raw data',['Fit: a=',num2str(phat(1),'%4.2f'),', b=',num2str(phat(2),'%4.2f')])
end

    %%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
if fn(2)
    figure(fn(2))
        clf
        plot(p_profile0(reg_I0),nucleus_RNA_profile0((reg_I0),4),'b.')
        hold on
        errorbar(pro_bin0(good_bin),RNA_mean(good_bin),RNA_std(good_bin),'.k')
        if ~any(isnan(H_I))
            plot(pro_bin0(good_bin),Hill(H_I,pro_bin0(good_bin)),'r')
        else
            plot([],[],'r')
        end
        hold off
        title([image_folder,', cycle = ',num2str(N_cycle),char(10),', h = ',num2str(H_I(1))],'Interpreter','none')
        xlabel([p_name,' concentration (M)'])
        ylabel('RNA level (#)')
        legend('Individual nuclei','Averaged profile','Fitting curve for averaged data')
        legend('hide')
end

if fn(13)
    figure(fn(13))
        clf
        loglog(pro_bin0(good_bin),yk(good_bin),'b.')
        hold on
        if all(~isnan(pk))
            loglog(pro_bin0(good_bin),pro_bin0(good_bin).^pk(1).*exp(pk(2)),'r')
        else
            plot([],[],'r')
        end
        hold off
        title([image_folder,', cycle = ',num2str(N_cycle),char(10),', h = ',num2str(-pk(1))],'Interpreter','none')
        xlabel([p_name,' concentration (M)'])
        ylabel('koff/kon')
        legend('Averaged profile','Fitting curve for averaged data')
        legend('hide')
end


if fn(9)
    figure(fn(9))
        clf
        if all(~isnan(alpha1))
            [AX,H1,H2] = plotyy(fmean,[fvar.^2./fmean,fth(alpha1,fmean),fvar_ex./fmean],fmean,[fcorr,fcorrn,fvar_ex./fvar.^2,fcorrb,fcorrnb,fp0,fp0th]);
        else
            [AX,H1,H2] = plotyy(fmean,[fvar.^2./fmean,nan(size(fmean)),fvar_ex./fmean],fmean,[fcorr,fcorrn,fvar_ex./fvar.^2,fcorrb,fcorrnb,fp0,fp0th]);
        end
        set(H1(1),'LineStyle','none','Marker','.','Color','b','DisplayName','\sigma^2/n vs n data')
        set(H2(1),'LineStyle','none','Marker','o','Color','g','DisplayName','correlation data')
        set(H1(2),'LineStyle','-','Color','r','DisplayName',['Fit: k_o_n*T = ',num2str(alpha1(1),'%4.3f'),', k_T_X*T = ',num2str(alpha1(2),'%4.3f')])
        set(H2(2),'LineStyle','none','Marker','*','Color','c','DisplayName','correlation data control')
        set(H1(3),'LineStyle','none','Marker','.','Color','k','DisplayName','ex\_var/n vs n data')
        set(H2(3),'LineStyle','none','Marker','o','Color',[0.5,0.5,0.5],'DisplayName','ex\_corr data')
        set(H2(4),'LineStyle','none','Marker','*','Color',[0,0.5,0],'DisplayName','digital correlation data')
        set(H2(5),'LineStyle','none','Marker','*','Color',[0,0.5,0.5],'DisplayName','digital correlation data control')
        set(H2(6),'LineStyle','-','Color',[0.5,0,0.5],'DisplayName','probability of silent foci')
        set(H2(7),'LineStyle','-','Color','r','DisplayName','theoretical probability of silent foci')
        legend('show')
        xlabel(AX(1),'Mean # mRNA (concentration binned)')
        ylabel(AX(1),'Var/mean (concentration binned)')
        ylabel(AX(2),'Correlation coefficient (concentration binned)')
        title([image_folder,', cycle = ',num2str(N_cycle)],'Interpreter','none')
end
    
if fn(11)
    figure(fn(11))
        clf
        bar(Cbin,nRNA/sum(nRNA)/(Cbin(2)-Cbin(1)));
        hold on
        plot(Cbin,gampdf(Cbin,phat(1),phat(2)),'r-')

        title([image_folder,', cycle = ',num2str(N_cycle),char(10),'([Bcd] = ',num2str(pmin2/1e-9,'%4.2f'),'nM - ',num2str(pmax2/1e-9,'%4.2f'),'nM)'],'Interpreter','none')
        xlabel('# mRNA')
        ylabel('Frequency')
        legend('Raw data',['Fit: a=',num2str(phat(1),'%4.2f'),', b=',num2str(phat(2),'%4.2f')])
end

fit_out = [phat,alpha1,mean(fcorr),mean(fcorrn),-pk(1),(exp(pk(2))).^(-1./pk(1))];
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Re-plot the transcription regulaton curve using nuclear foci # %%%%%%%%%%%%%%%%%
RNA_mean = zeros(size(pro_center0));
RNA_std = zeros(size(pro_center0));

for Ip = 1:length(pro_center0)
    RNA_mean(Ip) = mean(nucleus_RNA_profile(reg_I&(p_profile >= (pro_center0(Ip)-pwindow/2))&(p_profile <= (pro_center0(Ip)+pwindow/2)),3));
    RNA_std(Ip) = std0(nucleus_RNA_profile(reg_I&(p_profile >= (pro_center0(Ip)-pwindow/2))&(p_profile <= (pro_center0(Ip)+pwindow/2)),3));
end
good_bin = (~isnan(RNA_std)) & (RNA_std > 0);
RNAN_out = nan(size(pro_center0));
RNAN_out(good_bin) = RNA_mean(good_bin);

if nnz(good_bin) > 5% && sum(nucleus_RNA_profile((reg_I),3)) > 20 && pro_center0(end) > 10000 && isempty(strfind(image_folder((end-4):end), single_add))
    %%% Nonlinear fitting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [v_max,I_max] = max(RNA_mean);
    pro_I = find(pro_center0 <= pro_center0(I_max));
    [~,I_middle] = min(abs(RNA_mean(pro_I)-max(RNA_mean)/2));
    beta0 = [4,pro_center0(pro_I(I_middle)),v_max,0];
    try
        [beta,r,~,~,~] = nlinfit(pro_center0(good_bin),RNA_mean(good_bin),@Hill,beta0);
        H_N = beta;
    catch err
        H_N = nan(size(beta0));
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    plot([],[],'r')
    H_N = nan(1,4);
end

if fn(3)
    figure(fn(3))
    subaxis(sub_pos(1),sub_pos(2),sub_pos(3), 'Spacing', 0, 'PaddingRight',0.01, 'PaddingLeft',0.02, 'PaddingTop',0.045, 'PaddingBottom',0.04, 'Margin', 0.01);
        plot(p_profile(reg_I),nucleus_RNA_profile(reg_I,3),'b.')
        hold on
        errorbar(pro_center0(good_bin),RNA_mean(good_bin),RNA_std(good_bin),'.k')
        if all(~isnan(H_N))
            plot(pro_center0(good_bin),Hill(H_N,pro_center0(good_bin)),'r')
        else
            plot([],[],'r')
        end
        hold off
        title([image_folder,', cycle = ',num2str(N_cycle),char(10),', h = ',num2str(H_N(1))],'Interpreter','none')
        xlabel([p_name,' concentration (M)'])
        ylabel('# of foci per nucleus')
        legend('Individual nuclei','Averaged profile','Fitting curve')
        legend('hide')
end
    
if fn(4)
    figure(fn(4))
        clf
        plot(p_profile(reg_I),nucleus_RNA_profile(reg_I,3),'b.')
        hold on
        errorbar(pro_center0(good_bin),RNA_mean(good_bin),RNA_std(good_bin),'.k')
        if ~any(isnan(H_N))
            plot(pro_center0(good_bin),Hill(H_N,pro_center0(good_bin)),'r')
        else
            plot([],[],'r')
        end
        hold off
        title([image_folder,', cycle = ',num2str(N_cycle),char(10),', h = ',num2str(H_N(1))],'Interpreter','none')
        xlabel([p_name,' concentration (M)'])
        ylabel('# of foci per nucleus')
        legend('Individual nuclei','Averaged profile','Fitting curve')
        legend('hide')
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Plot the silent foci % vs protein concentration %%%%%%%%%%%%%%%%%%%%%%%%
% N_pro = 100;
% w_pro = 0.02;
% pmin = 0;
% pmax = 5e-8;
pshow_max = 5e-8;
% pbin = (pmax-pmin)/N_pro;
% pwindow = (pmax-pmin)*w_pro;
% pro_center = (pmin+pbin/2):pbin:(pmax-pbin/2);
pro_center = pro_center0;

null_rate = nan(size(pro_center));
fn0 = nan(size(pro_center));
fn1 = nan(size(pro_center));
fn2 = nan(size(pro_center));
fn3 = nan(size(pro_center));

for Ip = 1:length(pro_center)
    N_temp = nucleus_RNA_profile(reg_I&(p_profile >= (pro_center(Ip)-pwindow/2))&(p_profile <= (pro_center(Ip)+pwindow/2)),3);
    fn0(Ip) = nnz(N_temp == 0);
    fn1(Ip) = nnz(N_temp == 1);
    fn2(Ip) = nnz(N_temp == 2);
    fn3(Ip) = nnz(N_temp > 2);
    
    if length(N_temp) >= Nmin
        null_rate(Ip) = (nnz(N_temp == 0)*2+nnz(N_temp == 1))/length(N_temp)/2;
    end
end
null_rate = 1-null_rate;
RNAI_out = [RNAI_out,fn0,fn1,fn2,fn3];
fn_all = fn0+fn1+fn2+fn3;
fn0 = fn0./fn_all;
fn1 = fn1./fn_all;
fn2 = fn2./fn_all;
fn3 = fn3./fn_all;

%%% Nonlinear fitting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
v_max = max(null_rate);
pro_I = find(pro_center <= pshow_max);
[~,I_middle] = min(abs(null_rate(pro_I)-max(null_rate)/2));
beta0 = [4,pro_center(pro_I(I_middle)),v_max,0];
try
    [beta,r,~,~,~] = nlinfit(pro_center(pro_I),null_rate(pro_I),@Hill,beta0);
    H_A = beta;
catch err
    H_A = nan(size(beta0));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if fn(5)
    figure(fn(5))
    subaxis(sub_pos(1),sub_pos(2),sub_pos(3), 'Spacing', 0, 'PaddingRight',0.01, 'PaddingLeft',0.02, 'PaddingTop',0.045, 'PaddingBottom',0.04, 'Margin', 0.01);
        plot(pro_center,null_rate,'Color',[0.5,0.5,0.5])
        hold on
        if all(~isnan(H_A))
            plot(pro_center(pro_I),Hill(H_A,pro_center(pro_I)),'Color',[0.5,0,0])
        else
            plot([],[],'r')
        end
        plot(pro_center,fn0,'k',pro_center,fn1,'b',pro_center,fn2,'g',pro_center,fn3,'r')
        hold off
        title([image_folder,', cycle = ',num2str(N_cycle),char(10),', h = ',num2str(H_A(1))],'Interpreter','none')
        xlabel([p_name,' concentration (M)'])
        ylabel('Frequency')
        legend('P(foci on)','Fitted curve','P(n = 0)','P(n = 1)','P(n = 2)','P(n > 2)')
        legend('hide')
end
    
if fn(6)
    figure(fn(6))
        clf
        plot(pro_center,null_rate,'Color',[0.5,0.5,0.5])
        hold on
        if all(~isnan(H_A))
            plot(pro_center(pro_I),Hill(H_A,pro_center(pro_I)),'Color',[0.5,0,0])
        else
            plot([],[],'r')
        end
        plot(pro_center,fn0,'k',pro_center,fn1,'b',pro_center,fn2,'g',pro_center,fn3,'r')
        hold off
        title([image_folder,', cycle = ',num2str(N_cycle),char(10),', h = ',num2str(H_A(1))],'Interpreter','none')
        xlabel([p_name,' concentration (M)'])
        ylabel('Frequency')
        legend('P(foci on)','Fitted curve','P(n = 0)','P(n = 1)','P(n = 2)','P(n > 2)')
        legend('hide')
end

if fn(7)
    figure(fn(7))
        plot(pro_center,null_rate,'Color',color_code,'DisplayName',[image_folder,', cycle = ',num2str(N_cycle)])
        hold on
        title('All embryos overlapped','Interpreter','none')
        xlabel([p_name,' concentration (M)'])
        ylabel('Active rate')
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
    

function y = Hill(beta,x)
 y = beta(3)*x.^beta(1)./(x.^beta(1)+beta(2)^beta(1))+beta(4);

