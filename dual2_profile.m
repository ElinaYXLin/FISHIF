function dual2_profile(nucleus_protein_profile,nucleus_protein2_profile,nucleus_RNA_profile,foci_RNA_profile,image_folder,N_cycle,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Protein - RNA regulation curve plot: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameter setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_pro = 50;
w_pro = 0.1;
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
    pro1 = varargin{2}{1};
    pro2 = varargin{2}{2};
else
    pro1 = 'Protein 1';
    pro2 = 'Protein 2';
end
if length(varargin) > 2
    unit1 = varargin{3}{1};
    unit2 = varargin{3}{2};
else
    unit1 = 'A.U.';
    unit2 = 'A.U.';
end
if length(varargin) > 3
    figure_n = varargin{4};
else
    figure_n = [36,37,38];
end

if isempty(subtp)
    pmin = 0;
    p2min = 0;
else
    pmin = min(nucleus_protein_profile(reg_I,2));
    p2min = min(nucleus_protein2_profile(reg_I,2));
end
pmax = max(nucleus_protein_profile(reg_I,2));
pmax = pmax-pmin;
pmin = 0;
pbin = (pmax-pmin)/N_pro;
pwindow = (pmax-pmin)*w_pro;
pro_center = (pmin+pbin/2):pbin:(pmax-pbin/2);

p2max = max(nucleus_protein2_profile(reg_I,2));
p2max = p2max-p2min;
p2min = 0;
p2bin = (p2max-p2min)/N_pro;
p2window = (p2max-p2min)*w_pro;
pro2_center = (p2min+p2bin/2):p2bin:(p2max-p2bin/2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Background subtraction: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pro_center = pro_center-min(nucleus_protein_profile(:,2));
p_profile = nucleus_protein_profile(:,2)-min(nucleus_protein_profile(:,2));

pro2_center = pro2_center-min(nucleus_protein2_profile(:,2));
p2_profile = nucleus_protein2_profile(:,2)-min(nucleus_protein2_profile(:,2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Regulation curve plot: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(figure_n(1))
clf
RNA_mean = zeros(length(pro_center),length(pro2_center));
RNA_std = zeros(length(pro_center),length(pro2_center));
    for Ip = 1:length(pro_center)
        for Ip2 = 1:length(pro2_center)
            RNA_mean(Ip,Ip2) = mean(nucleus_RNA_profile(reg_I&(p_profile >= (pro_center(Ip)-pwindow/2))&(p_profile <= (pro_center(Ip)+pwindow/2)) & (p2_profile >= (pro2_center(Ip2)-p2window/2))&(p2_profile <= (pro2_center(Ip2)+p2window/2)),4));
            RNA_std(Ip,Ip2) = std(nucleus_RNA_profile(reg_I&(p_profile >= (pro_center(Ip)-pwindow/2))&(p_profile <= (pro_center(Ip)+pwindow/2)) & (p2_profile >= (pro2_center(Ip2)-p2window/2))&(p2_profile <= (pro2_center(Ip2)+p2window/2)),4));
        end
    end
    subplot(1,2,1)
    imagesc(pro2_center,pro_center,RNA_mean)
    axis xy
    title(['Nucleus dual regulation (RNA level) curve: ',image_folder,', cycle = ',num2str(N_cycle)],'Interpreter','none')
    ylabel([pro1,' level (',unit1,')'])
    xlabel([pro2,' level (',unit1,')'])
    hc = colorbar;
    ylabel(hc,['RNA level (',unit2,')'])
    
    subplot(1,2,2)
    imagesc(pro2_center,pro_center,RNA_std)
    axis xy
    title(['Nucleus dual regulation (RNA level) std: ',image_folder,', cycle = ',num2str(N_cycle)],'Interpreter','none')
    ylabel([pro1,' level (',unit1,')'])
    xlabel([pro2,' level (',unit1,')'])
    hc = colorbar;
    ylabel(hc,['RNA std (',unit2,')'])

    
    
figure(figure_n(2))
clf
RNA_mean = zeros(length(pro_center),length(pro2_center));
RNA_std = zeros(length(pro_center),length(pro2_center));
    for Ip = 1:length(pro_center)
        for Ip2 = 1:length(pro2_center)
            RNA_mean(Ip,Ip2) = mean(nucleus_RNA_profile(reg_I&(p_profile >= (pro_center(Ip)-pwindow/2))&(p_profile <= (pro_center(Ip)+pwindow/2)) & (p2_profile >= (pro2_center(Ip2)-p2window/2))&(p2_profile <= (pro2_center(Ip2)+p2window/2)),3));
            RNA_std(Ip,Ip2) = std(nucleus_RNA_profile(reg_I&(p_profile >= (pro_center(Ip)-pwindow/2))&(p_profile <= (pro_center(Ip)+pwindow/2)) & (p2_profile >= (pro2_center(Ip2)-p2window/2))&(p2_profile <= (pro2_center(Ip2)+p2window/2)),3));
        end
    end
    subplot(1,2,1)
    imagesc(pro2_center,pro_center,RNA_mean)
    axis xy
    title(['Nucleus dual regulation (foci #) curve: ',image_folder,', cycle = ',num2str(N_cycle)],'Interpreter','none')
    ylabel([pro1,' level (',unit1,')'])
    xlabel([pro2,' level (',unit1,')'])
    hc = colorbar;
    ylabel(hc,'# of foci per nucleus')
    
    subplot(1,2,2)
    axis xy
    imagesc(pro2_center,pro_center,RNA_std)
    title(['Nucleus dual regulation (foci #) std: ',image_folder,', cycle = ',num2str(N_cycle)],'Interpreter','none')
    ylabel([pro1,' level (',unit1,')'])
    xlabel([pro2,' level (',unit1,')'])
    hc = colorbar;
    ylabel(hc,'std of # of foci per nucleus')
    
        
figure(figure_n(3))
clf
RNA_mean = zeros(length(pro_center),length(pro2_center));
RNA_std = zeros(length(pro_center),length(pro2_center));
    I_foci = find(foci_RNA_profile(:,2) > 0);
    I_foci = I_foci((nucleus_protein_profile(foci_RNA_profile(I_foci,2),1) >= reg_min)&(nucleus_protein_profile(foci_RNA_profile(I_foci,2),1) <= reg_max) & (nucleus_protein2_profile(foci_RNA_profile(I_foci,2),1) >= reg_min)&(nucleus_protein2_profile(foci_RNA_profile(I_foci,2),1) <= reg_max));
    for Ip = 1:length(pro_center)
        for Ip2 = 1:length(pro2_center)
            RNA_mean(Ip,Ip2) = mean(foci_RNA_profile(I_foci((p_profile(foci_RNA_profile(I_foci,2)) >= (pro_center(Ip)-pwindow/2))&(p_profile(foci_RNA_profile(I_foci,2)) <= (pro_center(Ip)+pwindow/2)) & (p2_profile(foci_RNA_profile(I_foci,2)) >= (pro2_center(Ip2)-p2window/2))&(p2_profile(foci_RNA_profile(I_foci,2)) <= (pro2_center(Ip2)+p2window/2))),3));
            RNA_std(Ip,Ip2) = std(foci_RNA_profile(I_foci((p_profile(foci_RNA_profile(I_foci,2)) >= (pro_center(Ip)-pwindow/2))&(p_profile(foci_RNA_profile(I_foci,2)) <= (pro_center(Ip)+pwindow/2)) & (p2_profile(foci_RNA_profile(I_foci,2)) >= (pro2_center(Ip2)-p2window/2))&(p2_profile(foci_RNA_profile(I_foci,2)) <= (pro2_center(Ip2)+p2window/2))),3));
        end
    end
    subplot(1,2,1)
    imagesc(pro2_center,pro_center,RNA_mean)
    axis xy
    title(['Foci dual regulation (RNA level) curve: ',image_folder,', cycle = ',num2str(N_cycle)],'Interpreter','none')
    ylabel([pro1,' level (',unit1,')'])
    xlabel([pro2,' level (',unit1,')'])
    hc = colorbar;
    ylabel(hc,['RNA level (',unit2,')'])
    
    subplot(1,2,2)
    imagesc(pro2_center,pro_center,RNA_std)
    axis xy
    title(['Foci dual regulation (RNA level) std: ',image_folder,', cycle = ',num2str(N_cycle)],'Interpreter','none')
    ylabel([pro1,' level (',unit1,')'])
    xlabel([pro2,' level (',unit1,')'])
    hc = colorbar;
    ylabel(hc,['RNA std (',unit2,')'])

    
    