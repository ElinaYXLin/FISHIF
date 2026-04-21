function GFPIF_profile(nucleus_protein_profile,nucleus_GFP_profile,image_folder,N_cycle)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Protein - GFP comparison curve plot: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameter setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_pro = 100;
pmax = max(nucleus_protein_profile(:,2));
pmin = min(nucleus_protein_profile(:,2));
pbin = (pmax-pmin)/N_pro;
pro_x = pmin:pbin:pmax;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% GFP-IF fitting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[p,S] = polyfit(nucleus_protein_profile(:,2),nucleus_GFP_profile(:,2),1);
[pro_y,delta] = polyval(p,pro_x,S);
R = corrcoef(nucleus_protein_profile(:,2),nucleus_GFP_profile(:,2));
if numel(R) == 1
    R(1,2) = nan;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Comparison curve plot: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(25)
clf
    plot(nucleus_protein_profile(:,2),nucleus_GFP_profile(:,2),'bo')
    hold on
    plot(pro_x,pro_y,'r')
    
    title(['Nucleus GFP-IF comparison curve: ',image_folder,', cycle = ',num2str(N_cycle)])
    xlabel('Antibody level (A.U.)')
    ylabel('GFP level (A.U.)')
    legend('Individual nuclei',['Linear fitting',', I(GFP) = ',num2str(p(1)),'*I(IF) + ',num2str(p(2)),', r = ',num2str(R(1,2)),', dev/L = ',num2str(mean(delta)/(pmax-pmin)/p(1))])
    plot(pro_x,pro_y-2*delta,'g--',pro_x,pro_y+2*delta,'g--')
    hold off

