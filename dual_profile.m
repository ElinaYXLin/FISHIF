function dual_profile(nucleus_protein_profile,nucleus_RNA_profile,foci_RNA_profile,image_folder,N_cycle,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Protein - RNA regulation curve plot: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameter setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_pro = 100;
w_pro = 0.05;
reg_min = 0.25;
reg_max = 0.75;
subtp = 'background subtraction';
single_add = '_new';

if ~isempty(varargin) && size(nucleus_protein_profile,2) >= 4
    z_size = varargin{1};
    reg_I = (nucleus_protein_profile(:,1) >= reg_min) & (nucleus_protein_profile(:,1) <= reg_max) & (nucleus_protein_profile(:,4) > 1) & (nucleus_protein_profile(:,4) < z_size);
else
    reg_I = (nucleus_protein_profile(:,1) >= reg_min) & (nucleus_protein_profile(:,1) <= reg_max);
end
if length(varargin) > 1
    unit1 = varargin{2}{1};
    unit2 = varargin{2}{2};
else
    unit1 = 'A.U.';
    unit2 = 'A.U.';
end
if length(varargin) > 2
    figure_n = varargin{3};
else
    figure_n = [12,15,13];
end


if isempty(subtp)
    pmin = 0;
else
    pmin = min(nucleus_protein_profile(reg_I,2));
end
pmax = max(nucleus_protein_profile(reg_I,2));
pmax = pmax-pmin;
pmin = 0;
pbin = (pmax-pmin)/N_pro;
pwindow = (pmax-pmin)*w_pro;
pro_center = (pmin+pbin/2):pbin:(pmax-pbin/2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Background subtraction: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pro_center = pro_center-min(nucleus_protein_profile(:,2));
p_profile = nucleus_protein_profile(:,2)-min(nucleus_protein_profile(:,2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Regulation curve plot: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(figure_n(1))
clf
RNA_mean = zeros(size(pro_center));
RNA_std = zeros(size(pro_center));
h1 = zeros(0);
h2 = zeros(0);
    for Ip = 1:length(pro_center)
        RNA_mean(Ip) = mean(nucleus_RNA_profile(reg_I&(p_profile >= (pro_center(Ip)-pwindow/2))&(p_profile <= (pro_center(Ip)+pwindow/2)),4));
        RNA_std(Ip) = std(nucleus_RNA_profile(reg_I&(p_profile >= (pro_center(Ip)-pwindow/2))&(p_profile <= (pro_center(Ip)+pwindow/2)),4));
    end
    plot(p_profile(reg_I),nucleus_RNA_profile((reg_I),4),'bo')
    hold(gca,'on')
    errorbar(pro_center((RNA_std)>0),RNA_mean((RNA_std)>0),RNA_std((RNA_std)>0),'.k')
    if sum((RNA_std)>0) > 5% && sum(nucleus_RNA_profile((reg_I),3)) > 20 && pro_center(end) > 10000 && isempty(strfind(image_folder((end-4):end), single_add))
        %%% Nonlinear fitting 1: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [~,I_max] = max(RNA_mean);
        pro_I = find(pro_center <= pro_center(I_max));
        [~,I_middle] = min(abs(RNA_mean(pro_I)-max(RNA_mean)/2));
        beta0 = [4,pro_center(pro_I(I_middle)),max(RNA_mean)];
        try
            [beta,r,~,~,~] = nlinfit(pro_center((RNA_std)>0),RNA_mean((RNA_std)>0),@Hill,beta0);
            plot(pro_center(~isnan(RNA_std)),Hill(beta,pro_center(~isnan(RNA_std))),'r')
            h1 = beta(1);
        catch err
            plot([],[],'r')
            h1 = NaN;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%% Nonlinear fitting 2: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        try
            [beta,r,~,~,~] = nlinfit(p_profile(reg_I),nucleus_RNA_profile((reg_I),4),@Hill,beta0);
            plot(pro_center(~isnan(RNA_std)),Hill(beta,pro_center(~isnan(RNA_std))),'--g')
            h2 = beta(1);
        catch err
            plot([],[],'--g')
            h2 = NaN;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    else
        plot([],[],'r')
        plot([],[],'--g')
        h1 = NaN;
        h2 = NaN;
    end
    ax1 = gca;
    ax2=axes('position',get(ax1,'position'),'yaxislocation','right','color','none','YColor','m','XLim',get(ax1,'XLim'),'XTick',get(ax1,'XTick'),'XTickLabel','');
    hold(ax2,'on')
    plot(ax2,pro_center(~isnan(RNA_std)),RNA_std(~isnan(RNA_std))./RNA_mean(~isnan(RNA_std)),'-*m')
    
    hold(ax1,'off')
    hold(ax2,'off')
    set(ax1,'Box','off');
    title(ax1,['Nucleus regulation (RNA level) curve: ',image_folder,', cycle = ',num2str(N_cycle),', h1 = ',num2str(h1),', h2 = ',num2str(h2)],'Interpreter','none')
    xlabel(ax1,['Protein level (',unit1,')'])
    ylabel(ax1,['RNA level (',unit2,')'])
    ylabel(ax2,['Std/mean of RNA level (',unit2,')'])
    h1 = get(ax1,'Children');
    h2 = get(ax2,'Children');
    legend([h1(end:-1:1);h2(end:-1:1)],'Individual nuclei','Averaged profile','Fitting curve for averaged data','Fitting curve for raw data')

figure(figure_n(2))
clf
RNA_mean = zeros(size(pro_center));
RNA_std = zeros(size(pro_center));
    for Ip = 1:length(pro_center)
        RNA_mean(Ip) = mean(nucleus_RNA_profile(reg_I&(p_profile >= (pro_center(Ip)-pwindow/2))&(p_profile <= (pro_center(Ip)+pwindow/2)),3));
        RNA_std(Ip) = std(nucleus_RNA_profile(reg_I&(p_profile >= (pro_center(Ip)-pwindow/2))&(p_profile <= (pro_center(Ip)+pwindow/2)),3));
    end
    plot(p_profile(reg_I),nucleus_RNA_profile(reg_I,3),'bo')
    hold(gca,'on')
    errorbar(pro_center((RNA_std)>0),RNA_mean((RNA_std)>0),RNA_std((RNA_std)>0),'.k')
    if sum((RNA_std)>0) > 5% && sum(nucleus_RNA_profile((reg_I),3)) > 20 && pro_center(end) > 10000 && isempty(strfind(image_folder((end-4):end), single_add))
        %%% Nonlinear fitting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [~,I_max] = max(RNA_mean);
        pro_I = find(pro_center <= pro_center(I_max));
        [~,I_middle] = min(abs(RNA_mean(pro_I)-max(RNA_mean)/2));
        beta0 = [4,pro_center(pro_I(I_middle)),max(RNA_mean)];
        try
            [beta,r,~,~,~] = nlinfit(pro_center,RNA_mean,@Hill,beta0);
            plot(pro_center((RNA_std)>0),Hill(beta,pro_center((RNA_std)>0)),'r')
            h1 = beta(1);
        catch err
            plot([],[],'r')
            h1 = NaN;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    else
        plot([],[],'r')
        h1 = NaN;
    end
    ax1 = gca;
    ax2=axes('position',get(ax1,'position'),'yaxislocation','right','color','none','YColor','m','XLim',get(ax1,'XLim'),'XTick',get(ax1,'XTick'),'XTickLabel','');
    hold(ax2,'on')
    plot(ax2,pro_center((RNA_std)>0),RNA_std((RNA_std)>0)./RNA_mean((RNA_std)>0),'-*m')
    
    hold(ax1,'off')
    hold(ax2,'off')
    set(ax1,'Box','off');
    title(['Nucleus regulation (foci #) curve: ',image_folder,', cycle = ',num2str(N_cycle),', h = ',num2str(h1)],'Interpreter','none')
    xlabel(ax1,['Protein level (',unit1,')'])
    ylabel(ax1,'# of foci per nucleus')
    ylabel(ax2,'Std/mean of RNA #')
    h1 = get(ax1,'Children');
    h2 = get(ax2,'Children');
    legend([h1(end:-1:1);h2(end:-1:1)],'Individual nuclei','Averaged profile','Fitting curve','Std/mean of averaged profile')
    
figure(figure_n(3))
clf
RNA_mean = zeros(size(pro_center));
RNA_std = zeros(size(pro_center));
    if ~isempty(foci_RNA_profile)
        I_foci = find(foci_RNA_profile(:,2) > 0);
        I_foci = I_foci((nucleus_protein_profile(foci_RNA_profile(I_foci,2),1) >= reg_min)&(nucleus_protein_profile(foci_RNA_profile(I_foci,2),1) <= reg_max));
        for Ip = 1:length(pro_center)
            RNA_mean(Ip) = mean(foci_RNA_profile(I_foci((p_profile(foci_RNA_profile(I_foci,2)) >= (pro_center(Ip)-pwindow/2))&(p_profile(foci_RNA_profile(I_foci,2)) <= (pro_center(Ip)+pwindow/2))),3));
            RNA_std(Ip) = std(double(foci_RNA_profile(I_foci((p_profile(foci_RNA_profile(I_foci,2)) >= (pro_center(Ip)-pwindow/2))&(p_profile(foci_RNA_profile(I_foci,2)) <= (pro_center(Ip)+pwindow/2))),3)));
        end
        plot(p_profile(foci_RNA_profile(I_foci,2)),foci_RNA_profile(I_foci,3),'o')
        hold on
        errorbar(pro_center((RNA_std)>0),RNA_mean((RNA_std)>0),RNA_std((RNA_std)>0),'k')
        hold off
        legend('Individual nuclei','Averaged profile')
    else
        I_foci = [];
        plot([],[],'o')
    end
    title(['Foci regulation curve: ',image_folder,', cycle = ',num2str(N_cycle)],'Interpreter','none')
    xlabel(['Protein level (',unit1,')'])
    ylabel(['RNA level (',unit2,')'])
    
function y = std0(x)
    NN = size(x);
    if NN(1) == 1;
        N1 = NN(2);
    else
        N1 = NN(1);
    end
    y = std(x)/sqrt(N1);
        
    
function y = Hill(beta,x)
 y = beta(3)*x.^beta(1)./(x.^beta(1)+beta(2)^beta(1));
