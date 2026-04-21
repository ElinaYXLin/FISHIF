function [pro_center0,pro2_center0,RNAI_out,RNAN_out,pro_center,pro2_center,null_rate] = dual2_profile_show(nucleus_protein_profile_ab,nucleus_protein2_profile_ab,nucleus_RNA_profile,foci_RNA_profile,image_folder,N_cycle,sub_pos,varargin)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A function to re-plot transcription regulation curve and fit it to a hill function %%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Parameter setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
reg_min = 0.25;%0.25;
reg_max = 0.75;%0.75;
single_add = '_new';

p_profile = nucleus_protein_profile_ab(:,2);
p2_profile = nucleus_protein2_profile_ab(:,2);
reg_I = (nucleus_protein_profile_ab(:,1) >= reg_min) & (nucleus_protein_profile_ab(:,1) <= reg_max) & nucleus_protein_profile_ab(:,3);

N_pro0 = 20;
w_pro0 = 0.1;
pmin0 = 0;
pmax0 = 5e-8;%max(nucleus_protein_profile_ab(reg_I,2));
pbin0 = (pmax0-pmin0)/N_pro0;
pwindow0 = (pmax0-pmin0)*w_pro0;
pro_center0 = (pmin0+pbin0/2):pbin0:(pmax0-pbin0/2);
p2min0 = 0;
p2max0 = 5e-8;%max(nucleus_protein2_profile_ab(reg_I,2));
p2bin0 = (p2max0-p2min0)/N_pro0;
p2window0 = (p2max0-p2min0)*w_pro0;
pro2_center0 = (p2min0+p2bin0/2):p2bin0:(p2max0-p2bin0/2);

N_pro = 50;
w_pro = 0.1;
pmin = 0;
pmax = 5e-8;
pshow_max = 5e-8;
pbin = (pmax-pmin)/N_pro;
pwindow = (pmax-pmin)*w_pro;
pro_center = (pmin+pbin/2):pbin:(pmax-pbin/2);
p2min = 0;
p2max = 5e-8;
p2show_max = 5e-8;
p2bin = (p2max-p2min)/N_pro;
p2window = (p2max-p2min)*w_pro;
pro2_center = (p2min+p2bin/2):p2bin:(p2max-p2bin/2);

if ~isempty(varargin)
    p_name = varargin{1}(1);
    p2_name = varargin{1}(2);
else
    p_name = 'Protein 1';
    p2_name = 'Protein 2';
end

if length(varargin) > 1
    fn = varargin{2};
else
    fn = [42,420,45,450,47,470];
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Re-plot the transcription regulaton curve using total nuclear foci intensity %%
RNA_mean = zeros(length(pro_center0),length(pro2_center0));
RNA_std = zeros(length(pro_center0),length(pro2_center0));

figure(fn(1))
subaxis(sub_pos(1),sub_pos(2),sub_pos(3), 'Spacing', 0, 'PaddingRight',0.01, 'PaddingLeft',0.02, 'PaddingTop',0.045, 'PaddingBottom',0.04, 'Margin', 0.01);
    for Ip = 1:length(pro_center0)
        for Ip2 = 1:length(pro2_center0)
            RNA_mean(Ip,Ip2) = mean(nucleus_RNA_profile(reg_I & (p_profile >= (pro_center0(Ip)-pwindow0/2))&(p_profile <= (pro_center0(Ip)+pwindow0/2)) & (p2_profile >= (pro2_center0(Ip2)-p2window0/2))&(p2_profile <= (pro2_center0(Ip2)+p2window0/2)),4));
            RNA_std(Ip,Ip2) = std0(nucleus_RNA_profile(reg_I & (p_profile >= (pro_center0(Ip)-pwindow0/2))&(p_profile <= (pro_center0(Ip)+pwindow0/2)) & (p2_profile >= (pro2_center0(Ip2)-p2window0/2))&(p2_profile <= (pro2_center0(Ip2)+p2window0/2)),4));
        end
    end
    RNAI_out = RNA_mean;
    imagesc(pro2_center0,pro_center0,RNA_mean)
    axis xy
    title(['Nucleus dual regulation (RNA level) curve: ',image_folder,', cycle = ',num2str(N_cycle)],'Interpreter','none')
    ylabel([p_name,' level (M)'])
    xlabel([p2_name,' level (M)'])
    hc = colorbar;
    ylabel(hc,['RNA level (#)'])

    
figure(fn(2))
    clf
    subplot(1,2,1)
    imagesc(pro2_center0,pro_center0,RNA_mean)
    axis xy
    title(['Nucleus dual regulation (RNA level) curve: ',image_folder,', cycle = ',num2str(N_cycle)],'Interpreter','none')
    ylabel([p_name,' level (M)'])
    xlabel([p2_name,' level (M)'])
    hc = colorbar;
    ylabel(hc,['RNA level (#)'])
    
    subplot(1,2,2)
    imagesc(pro2_center0,pro_center0,RNA_std)
    axis xy
    title(['Nucleus dual regulation (RNA level) error: ',image_folder,', cycle = ',num2str(N_cycle)],'Interpreter','none')
    ylabel([p_name,' level (M)'])
    xlabel([p2_name,' level (M)'])
    hc = colorbar;
    ylabel(hc,['RNA std (#)'])
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Re-plot the transcription regulaton curve using nuclear foci # %%%%%%%%%%%%%%%%%
RNA_mean = zeros(length(pro_center0),length(pro2_center0));
RNA_std = zeros(length(pro_center0),length(pro2_center0));

figure(fn(3))
subaxis(sub_pos(1),sub_pos(2),sub_pos(3), 'Spacing', 0, 'PaddingRight',0.01, 'PaddingLeft',0.02, 'PaddingTop',0.045, 'PaddingBottom',0.04, 'Margin', 0.01);
    for Ip = 1:length(pro_center0)
        for Ip2 = 1:length(pro2_center0)
            RNA_mean(Ip,Ip2) = mean(nucleus_RNA_profile(reg_I & (p_profile >= (pro_center0(Ip)-pwindow0/2))&(p_profile <= (pro_center0(Ip)+pwindow0/2)) & (p2_profile >= (pro2_center0(Ip2)-p2window0/2))&(p2_profile <= (pro2_center0(Ip2)+p2window0/2)),3));
            RNA_std(Ip,Ip2) = std0(nucleus_RNA_profile(reg_I & (p_profile >= (pro_center0(Ip)-pwindow0/2))&(p_profile <= (pro_center0(Ip)+pwindow0/2)) & (p2_profile >= (pro2_center0(Ip2)-p2window0/2))&(p2_profile <= (pro2_center0(Ip2)+p2window0/2)),3));
        end
    end
    RNAN_out = RNA_mean;
    imagesc(pro2_center0,pro_center0,RNA_mean)
    axis xy
    title(['Nucleus dual regulation (foci #) curve: ',image_folder,', cycle = ',num2str(N_cycle)],'Interpreter','none')
    ylabel([p_name,' level (M)'])
    xlabel([p2_name,' level (M)'])
    hc = colorbar;
    ylabel(hc,['# of foci per nucleus'])

    
figure(fn(4))
    clf
    subplot(1,2,1)
    imagesc(pro2_center0,pro_center0,RNA_mean)
    axis xy
    title(['Nucleus dual regulation (foci #) curve: ',image_folder,', cycle = ',num2str(N_cycle)],'Interpreter','none')
    ylabel([p_name,' level (M)'])
    xlabel([p2_name,' level (M)'])
    hc = colorbar;
    ylabel(hc,['# of foci per nucleus'])
    
    subplot(1,2,2)
    imagesc(pro2_center0,pro_center0,RNA_std)
    axis xy
    title(['Nucleus dual regulation (foci #) error: ',image_folder,', cycle = ',num2str(N_cycle)],'Interpreter','none')
    ylabel([p_name,' level (M)'])
    xlabel([p2_name,' level (M)'])
    hc = colorbar;
    ylabel(hc,['std of # of foci per nucleus'])
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Plot the silent foci % vs protein concentration %%%%%%%%%%%%%%%%%%%%%%%%
null_rate = zeros(length(pro_center),length(pro2_center));

figure(fn(5))
subaxis(sub_pos(1),sub_pos(2),sub_pos(3), 'Spacing', 0, 'PaddingRight',0.01, 'PaddingLeft',0.02, 'PaddingTop',0.045, 'PaddingBottom',0.04, 'Margin', 0.01);
    for Ip = 1:length(pro_center)
        for Ip2 = 1:length(pro2_center)
            N_temp = nucleus_RNA_profile(reg_I & (p_profile >= (pro_center(Ip)-pwindow/2))&(p_profile <= (pro_center(Ip)+pwindow/2)) & (p2_profile >= (pro2_center(Ip2)-p2window/2))&(p2_profile <= (pro2_center(Ip2)+p2window/2)),3);
            if length(N_temp) >= 5% && nnz(nucleus_RNA_profile(:,3) > 2) <= 10
                null_rate(Ip,Ip2) = (nnz(N_temp == 0)*2+nnz(N_temp == 1))/length(N_temp)/2;
    %         elseif length(N_temp) >= 5 && nnz(nucleus_RNA_profile(:,3) > 2) >10
    %             null_rate(Ip) = (nnz(N_temp == 0)*4+nnz(N_temp == 1)*3+nnz(N_temp == 2)*2+nnz(N_temp == 3)*1)/length(N_temp)/4;
            else
                null_rate(Ip,Ip2) = nan;
            end
        end
    end
    null_rate = 1-null_rate;
    imagesc(pro2_center,pro_center,null_rate)
    axis xy
    title(['Foci active rate: ',image_folder,', cycle = ',num2str(N_cycle)],'Interpreter','none')
    ylabel([p_name,' level (M)'])
    xlabel([p2_name,' level (M)'])
    hc = colorbar;
    ylabel(hc,['Active rate'])

    
figure(fn(6))
    clf
    imagesc(pro2_center,pro_center,null_rate)
    axis xy
    title(['Foci active rate: ',image_folder,', cycle = ',num2str(N_cycle)],'Interpreter','none')
    ylabel([p_name,' level (M)'])
    xlabel([p2_name,' level (M)'])
    hc = colorbar;
    ylabel(hc,['Active rate'])
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
    

