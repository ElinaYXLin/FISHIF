function [FISHI_out,FISHN_out,fociN_out,FISHI_out2,FISHN_out2,fociN_out2,fit_out,S2_out,S3_out] = RNA_profile3D_show(nucleus_RNA_profile,foci_RNA_profile,cytoplasmic_RNA_profile,mask_stack,image_folder,N_cycle,sub_pos,flip_EL,nucleus_protein_profile_ab,varargin)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A function to re-plot RNA profile %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Parameter setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nucleus_bin = 0:0.02:1;%0.025:0.05:0.975;
average_radius = 0.05;
EL_min = 0.3;%0.25;
EL_max = 0.7;%0.75;
EL_A = [0.2,0.5];
EL_P = [0.7,0.9];
Nmin = 10;
d_fit = 3;
fth = @(alpha1,x) 2/3+8/3/alpha1(1).*x-(8/alpha1(1)^2/alpha1(2)+16/3/alpha1(1)/alpha1(2)).*x.^2+16/alpha1(1)^2/alpha1(2)^2.*x.^3+64/alpha1(1)^4/alpha1(2)^3.*x.^4-128/alpha1(1)^4/alpha1(2)^4.*x.^5+(-32/alpha1(1)^3/alpha1(2)^2.*x.^3+64*(1/alpha1(1)^3/alpha1(2)^3-1/alpha1(1)^4/alpha1(2)^3).*x.^4+128/alpha1(1)^4/alpha1(2)^4.*x.^5).*exp(-1/2*alpha1(1)*alpha1(2)./x);
alpha0 = [0.2,65];
% EL_min2 = 0.3;
% EL_max2 = 0.4;
Cbin = 0:5:500;
dii = 2;
intensity_Nbin = 50;
bin_max = min(nucleus_bin+average_radius,1);
bin_min = max(nucleus_bin-average_radius,0);
if flip_EL
    nucleus_RNA_profile(:,1) = 1-nucleus_RNA_profile(:,1);
    cytoplasmic_RNA_profile(:,1) = 1-cytoplasmic_RNA_profile(:,1);
end
if ~isempty(varargin)
    fn = varargin{1};
else
    fn = [5,31,33,50,310,330,7,70,6,8,11,110,10,100,13,130];
end
if length(varargin) > 1
    EL_min2 = varargin{2}(1);
    EL_max2 = varargin{2}(2);
else
    EL_min2 = 0.3;
    EL_max2 = 0.4;
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Calculate the background FISH intensity: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S2_out(1) = mean(nucleus_RNA_profile(nucleus_RNA_profile(:,1)>=EL_A(1) & nucleus_RNA_profile(:,1)<=EL_A(2),5));
S2_out(2) = mean(nucleus_RNA_profile(nucleus_RNA_profile(:,1)>=EL_P(1) & nucleus_RNA_profile(:,1)<=EL_P(2),5));
S3_out = nnz(bwconvhull(max(logical(mask_stack),[],3)));
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Re-plot FISH intensity vs EL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[nucleus_RNA_profile0,ind_foci] = foci_info(nucleus_RNA_profile,foci_RNA_profile);
nucleus_distance = nucleus_RNA_profile0(:,1);
% Ifoci_nucleus = nucleus_RNA_profile0(:,4);
[nucleus_bin0,bin_min0,bin_max0] = seg_eqN(nucleus_distance,length(nucleus_bin),average_radius);

Itrue = (nucleus_distance >= EL_min) & (nucleus_distance <= EL_max) & nucleus_protein_profile_ab(ind_foci,3);
[x_out,y_out] = corr_plot(ind_foci(Itrue),nucleus_RNA_profile0(Itrue,:));
if fn(15)
    figure(fn(15))
    subaxis(sub_pos(1),sub_pos(2),sub_pos(3), 'Spacing', 0, 'PaddingRight',0.01, 'PaddingLeft',0.02, 'PaddingTop',0.045, 'PaddingBottom',0.04, 'Margin', 0.01);
        plot(x_out,y_out,'.')
        xlabel('Foci 1 intensity')
        ylabel('Foci 2 intensity')
        title([image_folder,', cycle = ',num2str(N_cycle)],'Interpreter','none')
end

if fn(16)
    figure(fn(16))
    clf
        plot(x_out,y_out,'.')
        xlabel('Foci 1 intensity')
        ylabel('Foci 2 intensity')
        title([image_folder,', cycle = ',num2str(N_cycle)],'Interpreter','none')
end


fin = nan(size(nucleus_bin));   %%% number of nuclei in the bin
EL0 = nan(size(nucleus_bin));   %%% mean EL in the bin
fi0 = nan(size(nucleus_bin));   %%% mean foci intensity in the bin
fi1 = nan(size(nucleus_bin));   %%% std foci intensity in the bin
fi2 = nan(size(nucleus_bin));   %%% err foci intensity in the bin
fi3 = nan(size(nucleus_bin));   %%% correlation coefficient of foci intensity in the bin
fi4 = nan(size(nucleus_bin));   %%% correlation coefficient control of foci intensity in the bin
fi30 = nan(size(nucleus_bin));   %%% digital correlation coefficient of foci intensity in the bin
fi40 = nan(size(nucleus_bin));   %%% digital correlation coefficient control of foci intensity in the bin
fi00 = nan(size(nucleus_bin));   %%% probability of silent foci in the bin
for I_bin = 1:length(nucleus_bin)
    Itrue = (nucleus_distance >= bin_min0(I_bin)) & (nucleus_distance <= bin_max0(I_bin)) & nucleus_protein_profile_ab(ind_foci,3);
    EL0(I_bin) = mean(nucleus_distance(Itrue));
    fin(I_bin) = nnz(unique(nucleus_RNA_profile0(Itrue,2)));
    if fin(I_bin) >= Nmin
        fi0(I_bin) = mean(nucleus_RNA_profile0(Itrue,4));
        fi1(I_bin) = std(nucleus_RNA_profile0(Itrue,4));
        fi2(I_bin) = fi1(I_bin)/sqrt(nnz(Itrue));
        [fi3(I_bin),fi4(I_bin),fi30(I_bin),fi40(I_bin)] = inu_corr(ind_foci(Itrue),nucleus_RNA_profile0(Itrue,:));
        fi00(I_bin) = mean(nucleus_RNA_profile0(Itrue,4) == 0);
%         fi00(I_bin) = mean(nucleus_RNA_profile(unique(ind_foci(Itrue)),4) == 0);
    end
end
fi3(isnan(fi3)) = 0;
fi4(isnan(fi4)) = 0;
fi30(isnan(fi30)) = 0;
fi40(isnan(fi40)) = 0;

fi5 = nan(size(nucleus_bin));   %%% Deviation of mean foci intensity in the bin
for I_bin = 1:length(nucleus_bin)
    Itrue = (nucleus_distance >= bin_min0(I_bin))&(nucleus_distance <= bin_max0(I_bin)) & nucleus_protein_profile_ab(ind_foci,3);
    if fin(I_bin) >= Nmin
        IImin = max(1,I_bin-d_fit);
        IImax = min(I_bin+d_fit,length(nucleus_bin));
        p0 = polyfit(EL0(IImin:IImax),fi0(IImin:IImax),2);
        dev0 = p0(1)*nucleus_distance(Itrue).^2+p0(2)*nucleus_distance(Itrue)-p0(1)*EL0(I_bin)^2-p0(2)*EL0(I_bin);
        fi5(I_bin) = sum(dev0.^2)/(nnz(Itrue)-1);
    end
end


FISHI_out = [nucleus_bin0',fi0',fi1',fin',fi3',fi4',fi5',fi30',fi40',fi00'];
good_bin = ~isnan(fi1);

xpatch = [nucleus_bin0(good_bin),nucleus_bin0(good_bin(end:-1:1))];
ypatch = [fi0(good_bin)-fi1(good_bin),fi0(good_bin(end:-1:1))+fi1(good_bin(end:-1:1))];
    
if fn(1)
    figure(fn(1))
    subaxis(sub_pos(1),sub_pos(2),sub_pos(3), 'Spacing', 0, 'PaddingRight',0.01, 'PaddingLeft',0.02, 'PaddingTop',0.045, 'PaddingBottom',0.04, 'Margin', 0.01);
        plot(nucleus_RNA_profile0(:,1),nucleus_RNA_profile0(:,2),'bo'); hold on
        plot(nucleus_RNA_profile0(:,1),nucleus_RNA_profile0(:,4),'Color',[0.7,0.7,0.7],'Marker','.','LineStyle','none')
        plot(nucleus_RNA_profile0(:,1),nucleus_RNA_profile0(:,5),'r^',cytoplasmic_RNA_profile(:,1),cytoplasmic_RNA_profile(:,2),'c*')

        errorbar(nucleus_bin0(good_bin),fi0(good_bin),fi2(good_bin),'k')
        patch(xpatch,ypatch,'r','FaceAlpha',0.5,'EdgeColor','none');

        title([image_folder,', cycle = ',num2str(N_cycle)],'Interpreter','none')
        xlabel('EL')
        ylabel('# of mRNA')
        legend('Nucleus mean intensity','Nucleus foci intensity','Nucleus background mean intensity','Cytoplasmic mean intensity','Mean nucleus foci intensity','Std of nucleus foci intensity')
        legend('hide')
end

fmean = fi0(nucleus_bin0 >= EL_min & nucleus_bin0 <= EL_max & fin >= Nmin)';
fvar = fi1(nucleus_bin0 >= EL_min & nucleus_bin0 <= EL_max & fin >= Nmin)';
fcorr = fi3(nucleus_bin0 >= EL_min & nucleus_bin0 <= EL_max & fin >= Nmin)';
fcorrn = fi4(nucleus_bin0 >= EL_min & nucleus_bin0 <= EL_max & fin >= Nmin)';
fcorrb = fi30(nucleus_bin0 >= EL_min & nucleus_bin0 <= EL_max & fin >= Nmin)';
fcorrnb = fi40(nucleus_bin0 >= EL_min & nucleus_bin0 <= EL_max & fin >= Nmin)';
fvar_ex = fi5(nucleus_bin0 >= EL_min & nucleus_bin0 <= EL_max & fin >= Nmin)';
fp0 = fi00(nucleus_bin0 >= EL_min & nucleus_bin0 <= EL_max & fin >= Nmin)';
try
    alpha1 = nlinfit(fmean,fvar.^2./fmean,fth,alpha0);
    koff0 = alpha1(1)*(alpha1(2)./2./fmean-1);
    fp0th = P0_2S(alpha1(1),koff0,alpha1(2));
catch
    alpha1 = nan(size(alpha0));
    fp0th = nan(size(fmean));
end

if fn(2)
    figure(fn(2))
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
        xlabel(AX(1),'Mean # mRNA (EL binned)')
        ylabel(AX(1),'Var/mean (EL binned)')
        ylabel(AX(2),'Correlation coefficient (EL binned)')
        title([image_folder,', cycle = ',num2str(N_cycle)],'Interpreter','none')
end


temp_data = nucleus_RNA_profile0((nucleus_distance >= EL_min2)&(nucleus_distance <= EL_max2),4);
nRNA = hist(temp_data,Cbin);
phat = gamfit(temp_data);
    
if fn(3)
    figure(fn(3))
    subaxis(sub_pos(1),sub_pos(2),sub_pos(3), 'Spacing', 0, 'PaddingRight',0.01, 'PaddingLeft',0.02, 'PaddingTop',0.045, 'PaddingBottom',0.04, 'Margin', 0.01);
        bar(Cbin,nRNA/sum(nRNA)/(Cbin(2)-Cbin(1)));
        hold on
        plot(Cbin,gampdf(Cbin,phat(1),phat(2)),'r-')

        title([image_folder,', cycle = ',num2str(N_cycle),char(10),'(EL = ',num2str(EL_min2),' - ',num2str(EL_max2),')'],'Interpreter','none')
        xlabel('# mRNA')
        ylabel('Frequency')
        legend('Raw data',['Fit: a=',num2str(phat(1),'%4.2f'),', b=',num2str(phat(2),'%4.2f')])
end

    %%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ELg = nucleus_bin0(good_bin & nucleus_bin0 >= EL_min & nucleus_bin0 <= EL_max);
    fi0g = fi0(good_bin & nucleus_bin0 >= EL_min & nucleus_bin0 <= EL_max);
    fi1g = fi1(good_bin & nucleus_bin0 >= EL_min & nucleus_bin0 <= EL_max);
    fi2g = fi2(good_bin & nucleus_bin0 >= EL_min & nucleus_bin0 <= EL_max);
    
    [max_FISHI,max_I] = max(fi0g);
    [min_FISHI,min_I] = min(fi0g(max_I:end));
    min_I = max_I+min_I-1;
    max_ELI = ELg(max_I);
    max_stdI = fi1g(max_I);
    max_errI = fi2g(max_I);
%     mid_FISHI = (max_FISHI+min_FISHI)/2;
%     [~,I0] = min(abs(fi0g(max_I:min_I)-mid_FISHI));
%     p0 = polyfit(ELg((I0+max_I-2):(I0+max_I)),fi0g((I0+max_I-2):(I0+max_I)),1);
%     mid_ELI = (mid_FISHI-p0(2))/p0(1);
    
    dfi = nan(size(ELg));
    for ii = 1:length(ELg)
        iimin = max(1,ii-dii);
        iimax = min(length(ELg),ii+dii);
        pp = polyfit(ELg(iimin:iimax),fi0g(iimin:iimax),1);
        dfi(ii) = -pp(1);
    end
    [~,I0] = max(dfi(max_I:min_I));
    mid_ELI = ELg(I0+max_I-1);
    
    FISHI_out2 = [max_FISHI,max_stdI,max_errI,max_ELI,mid_ELI];
    %%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
if fn(4)
    figure(fn(4))
    clf
        plot(nucleus_RNA_profile0(:,1),nucleus_RNA_profile0(:,2),'bo'); hold on
        plot(nucleus_RNA_profile0(:,1),nucleus_RNA_profile0(:,4),'Color',[0.7,0.7,0.7],'Marker','.','LineStyle','none')
        plot(nucleus_RNA_profile0(:,1),nucleus_RNA_profile0(:,5),'r^',cytoplasmic_RNA_profile(:,1),cytoplasmic_RNA_profile(:,2),'c*')
        errorbar(nucleus_bin0(good_bin),fi0(good_bin),fi2(good_bin),'k')
        patch(xpatch,ypatch,'r','FaceAlpha',0.5,'EdgeColor','none');

        title([image_folder,', cycle = ',num2str(N_cycle)],'Interpreter','none')
        xlabel('EL')
        ylabel('# of mRNA')
        legend('Nucleus mean intensity','Nucleus foci intensity','Nucleus background mean intensity','Cytoplasmic mean intensity','Mean nucleus foci intensity','Std of nucleus foci intensity')
end

if fn(5)
    figure(fn(5))
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
        xlabel(AX(1),'Mean # mRNA (EL binned)')
        ylabel(AX(1),'Var/mean (EL binned)')
        ylabel(AX(2),'Correlation coefficient (EL binned)')
        title([image_folder,', cycle = ',num2str(N_cycle)],'Interpreter','none')
end

if fn(6)
    figure(fn(6))
    clf
        bar(Cbin,nRNA/sum(nRNA)/(Cbin(2)-Cbin(1)));
        hold on
        plot(Cbin,gampdf(Cbin,phat(1),phat(2)),'r-')

        title([image_folder,', cycle = ',num2str(N_cycle),char(10),'(EL = ',num2str(EL_min2),' - ',num2str(EL_max2),')'],'Interpreter','none')
        xlabel('# mRNA')
        ylabel('Frequency')
        legend('Raw data',['Fit: a=',num2str(phat(1),'%4.2f'),', b=',num2str(phat(2),'%4.2f')])
end

fit_out = [phat,alpha1,mean(fcorr),mean(fcorrn)];
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Re-plot FISH intensity vs EL new %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nucleus_distance = nucleus_RNA_profile(:,1);
fi0 = zeros(size(nucleus_bin));
fi1 = zeros(size(nucleus_bin));
fi2 = zeros(size(nucleus_bin));
for I_bin = 1:length(nucleus_bin)
    fi0(I_bin) = mean(nucleus_RNA_profile((nucleus_distance >= bin_min(I_bin))&(nucleus_distance <= bin_max(I_bin)) & nucleus_protein_profile_ab(:,3),3));
    fi1(I_bin) = std(nucleus_RNA_profile((nucleus_distance >= bin_min(I_bin))&(nucleus_distance <= bin_max(I_bin)) & nucleus_protein_profile_ab(:,3),3));
    fi2(I_bin) = fi1(I_bin)/sqrt(nnz((nucleus_distance >= bin_min(I_bin))&(nucleus_distance <= bin_max(I_bin)) & nucleus_protein_profile_ab(:,3)));
end
FISHN_out = [nucleus_bin',fi0',fi1'];
good_bin = find(~isnan(fi1));

xpatch = [nucleus_bin(good_bin),nucleus_bin(good_bin(end:-1:1))];
ypatch = [fi0(good_bin)-fi1(good_bin),fi0(good_bin(end:-1:1))+fi1(good_bin(end:-1:1))];
    
if fn(7)
    figure(fn(7))
    subaxis(sub_pos(1),sub_pos(2),sub_pos(3), 'Spacing', 0, 'PaddingRight',0.01, 'PaddingLeft',0.02, 'PaddingTop',0.045, 'PaddingBottom',0.04, 'Margin', 0.01);
        plot(nucleus_RNA_profile(:,1),nucleus_RNA_profile(:,3),'Color',[0.7,0.7,0.7],'Marker','.','LineStyle','none'); hold on
        errorbar(nucleus_bin(good_bin),fi0(good_bin),fi2(good_bin),'k')
        patch(xpatch,ypatch,'r','FaceAlpha',0.5,'EdgeColor','none');

        title([image_folder,', cycle = ',num2str(N_cycle)],'Interpreter','none')
        xlabel('EL')
        ylabel('# of foci')
        legend('Nucleus foci number','Mean nucleus foci number','Std of nucleus foci number')
        legend('hide')
end

    %%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ELg = nucleus_bin(~isnan(fi1) & nucleus_bin >= EL_min & nucleus_bin <= EL_max);
    fi0g = fi0(~isnan(fi1) & nucleus_bin >= EL_min & nucleus_bin <= EL_max);
    fi1g = fi1(~isnan(fi1) & nucleus_bin >= EL_min & nucleus_bin <= EL_max);
    fi2g = fi2(~isnan(fi1) & nucleus_bin >= EL_min & nucleus_bin <= EL_max);
    
    [max_FISHN,max_I] = max(fi0g);
    [min_FISHN,min_I] = min(fi0g(max_I:end));
    min_I = max_I+min_I-1;
    max_ELN = ELg(max_I);
    max_stdN = fi1g(max_I);
    max_errN = fi2g(max_I);
%     mid_FISHN = (max_FISHN+min_FISHN)/2;
%     [~,I0] = min(abs(fi0g(max_I:min_I)-mid_FISHN));
%     p0 = polyfit(ELg((I0+max_I-2):(I0+max_I)),fi0g((I0+max_I-2):(I0+max_I)),1);
%     mid_ELN = (mid_FISHN-p0(2))/p0(1);
    
    dfi = nan(size(ELg));
    for ii = 1:length(ELg)
        iimin = max(1,ii-dii);
        iimax = min(length(ELg),ii+dii);
        pp = polyfit(ELg(iimin:iimax),fi0g(iimin:iimax),1);
        dfi(ii) = -pp(1);
    end
    [~,I0] = max(dfi(max_I:min_I));
    mid_ELN = ELg(I0+max_I-1);
    
    
    FISHN_out2 = [max_FISHN,max_stdN,max_errN,max_ELN,mid_ELN];
    %%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
if fn(8)
    figure(fn(8))
    clf
        plot(nucleus_RNA_profile(:,1),nucleus_RNA_profile(:,3),'Color',[0.7,0.7,0.7],'Marker','.','LineStyle','none'); hold on
        errorbar(nucleus_bin(good_bin),fi0(good_bin),fi2(good_bin),'k')
        patch(xpatch,ypatch,'r','FaceAlpha',0.5,'EdgeColor','none');

        title([image_folder,', cycle = ',num2str(N_cycle)],'Interpreter','none')
        xlabel('EL')
        ylabel('# of foci')
        legend('Nucleus foci number','Mean nucleus foci number','Std of nucleus foci number')
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Re-plot FISH spot #/nuclei vs EL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nfoci_nucleus = nucleus_RNA_profile(:,3);

fn0 = zeros(size(nucleus_bin));
fn1 = zeros(size(nucleus_bin));
fn2 = zeros(size(nucleus_bin));
fn3 = zeros(size(nucleus_bin));
for I_bin = 1:length(nucleus_bin)
    fn0(I_bin) = sum((nucleus_distance(Nfoci_nucleus == 0) >= bin_min(I_bin)).*(nucleus_distance(Nfoci_nucleus == 0) <= bin_max(I_bin)));
    fn1(I_bin) = sum((nucleus_distance(Nfoci_nucleus == 1) >= bin_min(I_bin)).*(nucleus_distance(Nfoci_nucleus == 1) <= bin_max(I_bin)));
    fn2(I_bin) = sum((nucleus_distance(Nfoci_nucleus == 2) >= bin_min(I_bin)).*(nucleus_distance(Nfoci_nucleus == 2) <= bin_max(I_bin)));
    fn3(I_bin) = sum((nucleus_distance(Nfoci_nucleus  > 2) >= bin_min(I_bin)).*(nucleus_distance(Nfoci_nucleus  > 2) <= bin_max(I_bin)));
end
%fn0 = hist(nucleus_distance(Nfoci_nucleus == 0), nucleus_bin);
%fn1 = hist(nucleus_distance(Nfoci_nucleus == 1), nucleus_bin);
%fn2 = hist(nucleus_distance(Nfoci_nucleus == 2), nucleus_bin);
%fn3 = hist(nucleus_distance(Nfoci_nucleus > 2), nucleus_bin);
fn_all = fn0+fn1+fn2+fn3+(fn0+fn1+fn2+fn3 == 0);
FISHI_out = [FISHI_out,fn0',fn1',fn2',fn3'];

if fn(9)
    figure(fn(9))
    subaxis(sub_pos(1),sub_pos(2),sub_pos(3), 'Spacing', 0, 'PaddingRight',0.01, 'PaddingLeft',0.02, 'PaddingTop',0.045, 'PaddingBottom',0.04, 'Margin', 0.01);
        plot(nucleus_bin,(1-fn0./fn_all)*100,'m',nucleus_bin,fn1./fn_all*100,'b',nucleus_bin,fn2./fn_all*100,'g',nucleus_bin,fn3./fn_all*100,'r')
        title([image_folder,', cycle = ',num2str(N_cycle)],'Interpreter','none')
        xlabel('EL')
        ylabel('%')
        legend('nuclei with active foci','nuclei with 1 foci','nuclei with 2 foci','nuclei with >3 foci')
        legend('hide')
end

    %%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fn00 = 1-fn0./fn_all;
    fn01 = fn00-fn00.^2;
    fn02 = fn01./sqrt(fn_all);
    ELg = nucleus_bin(~isnan(fi1) & nucleus_bin >= EL_min & nucleus_bin <= EL_max);
    fi0g = fn00(~isnan(fi1) & nucleus_bin >= EL_min & nucleus_bin <= EL_max);
    fi1g = fn01(~isnan(fi1) & nucleus_bin >= EL_min & nucleus_bin <= EL_max);
    fi2g = fn02(~isnan(fi1) & nucleus_bin >= EL_min & nucleus_bin <= EL_max);
    
    fociN_out = [nucleus_bin',fn00',fn01'];

    [max_FISHN0,max_I] = max(fi0g);
    [min_FISHN0,min_I] = min(fi0g(max_I:end));
    min_I = max_I+min_I-1;
    max_ELN0 = ELg(max_I);
    max_stdN0 = fi1g(max_I);
    max_errN0 = fi2g(max_I);
%     mid_FISHN0 = (max_FISHN0+min_FISHN0)/2;
%     [~,I0] = min(abs(fi0g(max_I:min_I)-mid_FISHN0));
%     p0 = polyfit(ELg((I0+max_I-2):(I0+max_I)),fi0g((I0+max_I-2):(I0+max_I)),1);
%     mid_ELN0 = (mid_FISHN0-p0(2))/p0(1);
    
    dfi = nan(size(ELg));
    for ii = 1:length(ELg)
        iimin = max(1,ii-dii);
        iimax = min(length(ELg),ii+dii);
        pp = polyfit(ELg(iimin:iimax),fi0g(iimin:iimax),1);
        dfi(ii) = -pp(1);
    end
    [~,I0] = max(dfi(max_I:min_I));
    mid_ELN0 = ELg(I0+max_I-1);
    
    
    fociN_out2 = [max_FISHN0,max_stdN0,max_errN0,max_ELN0,mid_ELN0];
    %%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Re-plot foci Intensity distribution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
foci_intensity = foci_RNA_profile(:,3);
[Ifoci_n,xout] = hist(foci_intensity,intensity_Nbin);
if fn(10)
    figure(fn(10))
    subaxis(sub_pos(1),sub_pos(2),sub_pos(3), 'Spacing', 0, 'PaddingRight',0.01, 'PaddingLeft',0.02, 'PaddingTop',0.045, 'PaddingBottom',0.04, 'Margin', 0.01);
        bar(xout,Ifoci_n/(sum(Ifoci_n)+(sum(Ifoci_n) == 0))*100)
        title([image_folder,', cycle = ',num2str(N_cycle)],'Interpreter','none')
        xlabel('Equavalent # of mRNA')
        ylabel('%')
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% FISH spot #/nuclei heat map %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
im0 = logical(ismember(mask_stack,find(Nfoci_nucleus == 0)));
im1 = logical(ismember(mask_stack,find(Nfoci_nucleus == 1)));
imfate(:,:,3) = max(double(im0 | im1),[],3);
clear im1
im2 = logical(ismember(mask_stack,find(Nfoci_nucleus == 2)));
imfate(:,:,2) = max(double(im0 | im2),[],3);
clear im2
im3 = logical(ismember(mask_stack,find(Nfoci_nucleus > 2)));
imfate(:,:,1) = max(double(im0 | im3),[],3);
clear im3

if fn(11)
    figure(fn(11))
    subaxis(sub_pos(1),sub_pos(2),sub_pos(3), 'Spacing', 0, 'PaddingRight',0.01, 'PaddingLeft',0.02, 'PaddingTop',0.045, 'PaddingBottom',0.04, 'Margin', 0.01);
        imshow(imfate);
        title([image_folder,' (cycle = ',num2str(N_cycle),char(10),', white: 0 foci, blue: 1 foci, green: 2 foci, red: > 2 foci)'],'Interpreter','none');
end
    
if fn(12)
    figure(fn(12))
        clf
        imshow(imfate);
        title([image_folder,' (cycle = ',num2str(N_cycle),char(10),', white: 0 foci, blue: 1 foci, green: 2 foci, red: > 2 foci)'],'Interpreter','none');
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    


%% FISH spot mRNA #/nuclei heat map %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ifoci_nucleus = nucleus_RNA_profile(:,4);
Ifoci_nucleus = [0;Ifoci_nucleus];
imfateI0 = Ifoci_nucleus(mask_stack+1);
mask2Dsum = sum(logical(mask_stack),3);
imfateI = sum(imfateI0,3)./mask2Dsum;
imfateI(~mask2Dsum) = 0;
    
if fn(13)
    figure(fn(13))
    subaxis(sub_pos(1),sub_pos(2),sub_pos(3), 'Spacing', 0, 'PaddingRight',0.01, 'PaddingLeft',0.02, 'PaddingTop',0.045, 'PaddingBottom',0.04, 'Margin', 0.01);
        patch([0.5,size(imfateI,2)+0.5,size(imfateI,2)+0.5,0.5],[0.5,0.5,size(imfateI,1)+0.5,size(imfateI,1)+0.5],'k','EdgeColor','none');
        hold on
        temp = imagesc(imfateI);
        alpha(temp,mask2Dsum)
    %     set(gca,'Color','k')
        axis off
        axis equal
        colorbar
        title([image_folder,' (cycle = ',num2str(N_cycle),')'],'Interpreter','none');
end
    
if fn(14)
    figure(fn(14))
        clf
        patch([0.5,size(imfateI,2)+0.5,size(imfateI,2)+0.5,0.5],[0.5,0.5,size(imfateI,1)+0.5,size(imfateI,1)+0.5],'k','EdgeColor','none');
        hold on
        temp = imagesc(imfateI);
        alpha(temp,mask2Dsum)
    %     set(gca,'Color','k')
        axis off
        axis equal
        colorbar
        title([image_folder,' (cycle = ',num2str(N_cycle),')'],'Interpreter','none');
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



