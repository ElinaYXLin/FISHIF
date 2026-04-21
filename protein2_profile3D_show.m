function [lambda_I,lambda_C,Hill2_I,Hill2_C,pro1_out,pro2_out] = protein2_profile3D_show(nucleus_protein_profile,cytoplasmic_protein_profile,quanti_p,nucleus_protein_profile_ab,nucleus_protein2_profile,cytoplasmic_protein2_profile,quanti_p2,nucleus_protein2_profile_ab,mask_stack,image_folder,N_cycle,sub_pos,flip_EL,pro_name)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A function to re-plot two protein profiles and fit them to given curves 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Parameter setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ccode = [0.7,0.7,0.7 ; 0,0,0.3 ; 0,0.3,0 ; 0,0,0; 0,0,0.7 ; 0,0.7,0.7 ; 1,0,0];
ccode2 = [0.5,0.5,0.5 ; 0,0,0.5 ; 0,0.5,0 ; 0,0,0; 0,0,0.5 ; 0,0.5,0.5 ; 0,1,0];
p_name = pro_name{1};
p2_name = pro_name{2};
Lmin = 0;
Lnbin = 0.05;
Lmax = 1;
Lnwindow = 0.05;
nu_L = Lmin:Lnbin:Lmax;
nu_mean = zeros(size(nu_L));
nu_std = zeros(size(nu_L));
nu2_mean = zeros(size(nu_L));
nu2_std = zeros(size(nu_L));
standard_dorsal = 'Standard protein/Bcd_dorsal.xls';
standard_ventral = 'Standard protein/Bcd_ventral.xls';
bcd_dxy = xlsread(standard_dorsal);
bcd_vxy = xlsread(standard_ventral);

index_true = logical(nucleus_protein_profile(:,4));
if flip_EL
    nucleus_protein_profile(:,1) = 1-nucleus_protein_profile(:,1);
    cytoplasmic_protein_profile(:,1) = 1-cytoplasmic_protein_profile(:,1);
    nucleus_protein_profile_ab(:,1) = 1-nucleus_protein_profile_ab(:,1);
    nucleus_protein2_profile(:,1) = 1-nucleus_protein2_profile(:,1);
    cytoplasmic_protein2_profile(:,1) = 1-cytoplasmic_protein2_profile(:,1);
    nucleus_protein2_profile_ab(:,1) = 1-nucleus_protein2_profile_ab(:,1);
end
exp_C = @(beta,x) beta(1)*exp(-beta(2)*x)+beta(3);
b2 = 5;
hill_C = @(beta,x) beta(3)*beta(2).^beta(1)./(x.^beta(1)+beta(2)^beta(1))+beta(4);
% hill_C = @(beta,x) beta(1)*exp(-beta(2)*x)+beta(3);
h0 = 4;
c0 = 0.5;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Raw profile plotting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nucleus_profile = nucleus_protein_profile;
nucleus2_profile = nucleus_protein2_profile;
for Lcenter = 1:length(nu_L)
    nu_map = (nucleus_profile(index_true,1) >= nu_L(Lcenter)-Lnwindow) & (nucleus_profile(index_true,1) <= nu_L(Lcenter)+Lnwindow);
    nu_mean(Lcenter) = mean(nucleus_profile(nu_map,2));
    nu_std(Lcenter) = std0(nucleus_profile(nu_map,2));
    nu_map = (nucleus2_profile(index_true,1) >= nu_L(Lcenter)-Lnwindow) & (nucleus2_profile(index_true,1) <= nu_L(Lcenter)+Lnwindow);
    nu2_mean(Lcenter) = mean(nucleus2_profile(nu_map,2));
    nu2_std(Lcenter) = std0(nucleus2_profile(nu_map,2));
end

mean_index = (~isnan(nu_mean));
[max_mean,mean_start] = max(nu_mean);
min_mean = min(nu_mean);
if mean_start > 1
    mean_index(1:(mean_start-1)) = false;
end
try
beta1 = nlinfit(nu_L(mean_index),nu_mean(mean_index),exp_C,[max_mean-min_mean,b2,min_mean]);
catch
beta1 = nan(1,3);
end

d_max = mean(bcd_dxy(1:5,2));
d_min = mean(bcd_dxy(end-4:end,2));
d_x0 = mean(bcd_dxy(1:5,1));
d_y = (bcd_dxy(:,2)-d_min)/(d_max-d_min)*exp_C(beta1,d_x0)+min_mean;

v_max = mean(bcd_vxy(1:5,2));
v_min = mean(bcd_vxy(end-4:end,2));
v_x0 = mean(bcd_vxy(1:5,1));
v_y = (bcd_vxy(:,2)-v_min)/(v_max-v_min)*exp_C(beta1,v_x0)+min_mean;

mean2_index = (~isnan(nu2_mean));
[max_mean2,mean2_start] = max(nu2_mean);
[min_mean2,mean2_end] = min(nu2_mean);
if mean2_start > 1
    mean2_index(1:(mean2_start-1)) = false;
end
if mean2_end < length(mean2_index)
    mean2_index((mean2_end+1):end) = false;
end
try
beta2 = nlinfit(nu_L(mean2_index),nu2_mean(mean2_index),hill_C,[h0,c0,max_mean2-min_mean2,min_mean2]);
% beta2 = nlinfit(nu_L(mean2_index),nu2_mean(mean2_index),hill_C,[max_mean2-min_mean2,b2,min_mean2]);
catch
beta2 = nan(1,4);
% beta2 = nan(1,3);
end

figure(2)
subaxis(sub_pos(1),sub_pos(2),sub_pos(3), 'Spacing', 0, 'PaddingRight',0.01, 'PaddingLeft',0.02, 'PaddingTop',0.045, 'PaddingBottom',0.04, 'Margin', 0.01);
    plot(nucleus_profile(:,1),nucleus_profile(:,2),'Marker','.','Color',ccode(1,:),'LineStyle','none');
    hold on
    plot(nucleus_profile(index_true,1),nucleus_profile(index_true,2),'Marker','.','Color',ccode(2,:),'LineStyle','none');
    plot(cytoplasmic_protein_profile(:,1),cytoplasmic_protein_profile(:,2),'Marker','*','Color',ccode(3,:),'LineStyle','none');
    errorbar(nu_L,nu_mean,nu_std,'Marker','.','Color',ccode(4,:),'LineStyle','none');
    plot(bcd_dxy(:,1),d_y,'Marker','X','Color',ccode(5,:));
    plot(bcd_vxy(:,1),v_y,'Marker','X','Color',ccode(6,:));
    plot(nu_L,exp_C(beta1,nu_L),'Color',ccode(7,:))

    plot(nucleus2_profile(:,1),nucleus2_profile(:,2),'Marker','.','Color',ccode2(1,:),'LineStyle','none');
    hold on
    plot(nucleus2_profile(index_true,1),nucleus2_profile(index_true,2),'Marker','.','Color',ccode2(2,:),'LineStyle','none');
    plot(cytoplasmic_protein2_profile(:,1),cytoplasmic_protein2_profile(:,2),'Marker','*','Color',ccode2(3,:),'LineStyle','none');
    errorbar(nu_L,nu2_mean,nu2_std,'Marker','.','Color',ccode2(4,:),'LineStyle','none');
%     plot(bcd_dxy(:,1),d_y,'Marker','X','Color',ccode(5,:));
%     plot(bcd_vxy(:,1),v_y,'Marker','X','Color',ccode(6,:));
    plot(nu_L,hill_C(beta2,nu_L),'Color',ccode2(7,:))

    hold off
    title([image_folder,', cycle = ',num2str(N_cycle)],'Interpreter','none')
    xlabel('AP axis (normalized)')
    ylabel('Intensity (A.U.)')
    legend_text = {[p_name,' nucleus concentration'],[p_name,' good nucleus concentration'],[p_name,' cytoplasmic concentration'],[p_name,' averaged nuclear profile'],[p_name,' standard dorsal profile'],[p_name,' standard ventral profile'],[p_name,' fitting result: lambda = ',num2str(1./beta1(2))]};
    legend_text2 = {[p2_name,' nucleus concentration'],[p2_name,' good nucleus concentration'],[p2_name,' cytoplasmic concentration'],[p2_name,' averaged nuclear profile'],[p2_name,' fitting result: h = ',num2str(beta2(1))]};
    legend_text = cat(2,legend_text,legend_text2);
    legend(legend_text)
    legend('hide')
    xlim([0,1]);
cmax0 = max(nucleus_protein_profile(:,2));
lambda_I = [beta1(1),1./beta1(2),beta1(3),cmax0];
c2max0 = max(nucleus_protein2_profile(:,2));
% Hill2_I = [beta2(1),1./beta2(2),beta2(3),c2max0];
Hill2_I = [beta2,c2max0];

figure(20)
clf
    plot(nucleus_profile(:,1),nucleus_profile(:,2),'Marker','.','Color',ccode(1,:),'LineStyle','none');
    hold on
    plot(nucleus_profile(index_true,1),nucleus_profile(index_true,2),'Marker','.','Color',ccode(2,:),'LineStyle','none');
    plot(cytoplasmic_protein_profile(:,1),cytoplasmic_protein_profile(:,2),'Marker','*','Color',ccode(3,:),'LineStyle','none');
    errorbar(nu_L,nu_mean,nu_std,'Marker','.','Color',ccode(4,:),'LineStyle','none');
    plot(bcd_dxy(:,1),d_y,'Marker','X','Color',ccode(5,:));
    plot(bcd_vxy(:,1),v_y,'Marker','X','Color',ccode(6,:));
    plot(nu_L,exp_C(beta1,nu_L),'Color',ccode(7,:))

    plot(nucleus2_profile(:,1),nucleus2_profile(:,2),'Marker','.','Color',ccode2(1,:),'LineStyle','none');
%     hold on
    plot(nucleus2_profile(index_true,1),nucleus2_profile(index_true,2),'Marker','.','Color',ccode2(2,:),'LineStyle','none');
    plot(cytoplasmic_protein2_profile(:,1),cytoplasmic_protein2_profile(:,2),'Marker','*','Color',ccode2(3,:),'LineStyle','none');
    errorbar(nu_L,nu2_mean,nu2_std,'Marker','.','Color',ccode2(4,:),'LineStyle','none');
%     plot(bcd_dxy(:,1),d_y,'Marker','X','Color',ccode(5,:));
%     plot(bcd_vxy(:,1),v_y,'Marker','X','Color',ccode(6,:));
    plot(nu_L,hill_C(beta2,nu_L),'Color',ccode2(7,:))

    hold off
    title([image_folder,', cycle = ',num2str(N_cycle)],'Interpreter','none')
    xlabel('AP axis (normalized)')
    ylabel('Intensity (A.U.)')
    legend_text = {[p_name,' nucleus concentration'],[p_name,' good nucleus concentration'],[p_name,' cytoplasmic concentration'],[p_name,' averaged nuclear profile'],[p_name,' standard dorsal profile'],[p_name,' standard ventral profile'],[p_name,' fitting result: lambda = ',num2str(1./beta1(2))]};
    legend_text2 = {[p2_name,' nucleus concentration'],[p2_name,' good nucleus concentration'],[p2_name,' cytoplasmic concentration'],[p2_name,' averaged nuclear profile'],[p2_name,' fitting result: h = ',num2str(beta2(1))]};
    legend_text = cat(2,legend_text,legend_text2);
    legend(legend_text)
    legend('hide')
    xlim([0,1]);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Scaled profile plotting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nucleus_profile = nucleus_protein_profile_ab;
nucleus2_profile = nucleus_protein2_profile_ab;
for Lcenter = 1:length(nu_L)
    nu_map = (nucleus_profile(index_true,1) >= nu_L(Lcenter)-Lnwindow) & (nucleus_profile(index_true,1) <= nu_L(Lcenter)+Lnwindow);
    nu_mean(Lcenter) = mean(nucleus_profile(nu_map,2));
    nu_std(Lcenter) = std0(nucleus_profile(nu_map,2));
    nu_map = (nucleus2_profile(index_true,1) >= nu_L(Lcenter)-Lnwindow) & (nucleus2_profile(index_true,1) <= nu_L(Lcenter)+Lnwindow);
    nu2_mean(Lcenter) = mean(nucleus2_profile(nu_map,2));
    nu2_std(Lcenter) = std0(nucleus2_profile(nu_map,2));
end
pro1_out = [nu_L',nu_mean',nu_std'];
pro2_out = [nu_L',nu2_mean',nu2_std'];

mean_index = (~isnan(nu_mean));
[max_mean,mean_start] = max(nu_mean);
min_mean = min(nu_mean);
if mean_start > 1
    mean_index(1:(mean_start-1)) = false;
end
try
beta1 = nlinfit(nu_L(mean_index),nu_mean(mean_index),exp_C,[max_mean-min_mean,b2,min_mean]);
catch
beta1 = nan(1,3);
end

d_max = mean(bcd_dxy(1:5,2));
d_min = mean(bcd_dxy(end-4:end,2));
d_x0 = mean(bcd_dxy(1:5,1));
d_y = (bcd_dxy(:,2)-d_min)/(d_max-d_min)*exp_C(beta1,d_x0)+min_mean;

v_max = mean(bcd_vxy(1:5,2));
v_min = mean(bcd_vxy(end-4:end,2));
v_x0 = mean(bcd_vxy(1:5,1));
v_y = (bcd_vxy(:,2)-v_min)/(v_max-v_min)*exp_C(beta1,v_x0)+min_mean;

mean2_index = (~isnan(nu2_mean));
[max_mean2,mean2_start] = max(nu2_mean);
[min_mean2,mean2_end] = min(nu2_mean);
if mean2_start > 1
    mean2_index(1:(mean2_start-1)) = false;
end
if mean2_end < length(mean2_index)
    mean2_index((mean2_end+1):end) = false;
end
try
beta2 = nlinfit(nu_L(mean2_index),nu2_mean(mean2_index),hill_C,[h0,c0,max_mean2-min_mean2,min_mean2]);
% beta2 = nlinfit(nu_L(mean2_index),nu2_mean(mean2_index),hill_C,[max_mean2-min_mean2,b2,min_mean2]);
catch
beta2 = nan(1,4);
% beta2 = nan(1,3);
end

figure(22)
subaxis(sub_pos(1),sub_pos(2),sub_pos(3), 'Spacing', 0, 'PaddingRight',0.01, 'PaddingLeft',0.02, 'PaddingTop',0.045, 'PaddingBottom',0.04, 'Margin', 0.01);
    plot(nucleus_profile(:,1),nucleus_profile(:,2),'Marker','.','Color',ccode(1,:),'LineStyle','none');
    hold on
    plot(nucleus_profile(index_true,1),nucleus_profile(index_true,2),'Marker','.','Color',ccode(2,:),'LineStyle','none');
    errorbar(nu_L,nu_mean,nu_std,'Marker','.','Color',ccode(4,:),'LineStyle','none');
    plot(bcd_dxy(:,1),d_y,'Marker','X','Color',ccode(5,:));
    plot(bcd_vxy(:,1),v_y,'Marker','X','Color',ccode(6,:));
    plot(nu_L,exp_C(beta1,nu_L),'Color',ccode(7,:))

    plot(nucleus2_profile(:,1),nucleus2_profile(:,2),'Marker','.','Color',ccode2(1,:),'LineStyle','none');
%     hold on
    plot(nucleus2_profile(index_true,1),nucleus2_profile(index_true,2),'Marker','.','Color',ccode2(2,:),'LineStyle','none');
    errorbar(nu_L,nu2_mean,nu_std,'Marker','.','Color',ccode2(4,:),'LineStyle','none');
%     plot(bcd_dxy(:,1),d_y,'Marker','X','Color',ccode(5,:));
%     plot(bcd_vxy(:,1),v_y,'Marker','X','Color',ccode(6,:));
    plot(nu_L,hill_C(beta2,nu_L),'Color',ccode2(7,:))

    hold off
    title([image_folder,', cycle = ',num2str(N_cycle)],'Interpreter','none')
    xlabel('AP axis (normalized)')
    ylabel('Concentration (M)')
    legend_text = {[p_name,' nucleus concentration'],[p_name,' good nucleus concentration'],[p_name,' cytoplasmic concentration'],[p_name,' averaged nuclear profile'],[p_name,' standard dorsal profile'],[p_name,' standard ventral profile'],[p_name,' fitting result: lambda = ',num2str(1./beta1(2))]};
    legend_text2 = {[p2_name,' nucleus concentration'],[p2_name,' good nucleus concentration'],[p2_name,' cytoplasmic concentration'],[p2_name,' averaged nuclear profile'],[p2_name,' fitting result: h = ',num2str(beta2(1))]};
    legend_text = cat(2,legend_text,legend_text2);
    legend(legend_text)
    legend('hide')
    xlim([0,1]);
cmax = max(nucleus_protein_profile_ab(:,2));
lambda_C = [beta1(1),1./beta1(2),beta1(3),cmax];
c2max = max(nucleus_protein2_profile_ab(:,2));
% Hill2_C = [beta2(1),1./beta2(2),beta2(3),c2max];
Hill2_C = [beta2,c2max];

figure(220)
clf
    plot(nucleus_profile(:,1),nucleus_profile(:,2),'Marker','.','Color',ccode(1,:),'LineStyle','none');
    hold on
    plot(nucleus_profile(index_true,1),nucleus_profile(index_true,2),'Marker','.','Color',ccode(2,:),'LineStyle','none');
    errorbar(nu_L,nu_mean,nu_std,'Marker','.','Color',ccode(4,:),'LineStyle','none');
    plot(bcd_dxy(:,1),d_y,'Marker','X','Color',ccode(5,:));
    plot(bcd_vxy(:,1),v_y,'Marker','X','Color',ccode(6,:));
    plot(nu_L,exp_C(beta1,nu_L),'Color',ccode(7,:))

    plot(nucleus2_profile(:,1),nucleus2_profile(:,2),'Marker','.','Color',ccode2(1,:),'LineStyle','none');
%     hold on
    plot(nucleus2_profile(index_true,1),nucleus2_profile(index_true,2),'Marker','.','Color',ccode2(2,:),'LineStyle','none');
    errorbar(nu_L,nu2_mean,nu_std,'Marker','.','Color',ccode2(4,:),'LineStyle','none');
%     plot(bcd_dxy(:,1),d_y,'Marker','X','Color',ccode(5,:));
%     plot(bcd_vxy(:,1),v_y,'Marker','X','Color',ccode(6,:));
    plot(nu_L,hill_C(beta2,nu_L),'Color',ccode2(7,:))

    hold off
    title([image_folder,', cycle = ',num2str(N_cycle)],'Interpreter','none')
    xlabel('AP axis (normalized)')
    ylabel('Concentration (M)')
    legend_text = {[p_name,' nucleus concentration'],[p_name,' good nucleus concentration'],[p_name,' cytoplasmic concentration'],[p_name,' averaged nuclear profile'],[p_name,' standard dorsal profile'],[p_name,' standard ventral profile'],[p_name,' fitting result: lambda = ',num2str(1./beta1(2))]};
    legend_text2 = {[p2_name,' nucleus concentration'],[p2_name,' good nucleus concentration'],[p2_name,' cytoplasmic concentration'],[p2_name,' averaged nuclear profile'],[p2_name,' fitting result: h = ',num2str(beta2(1))]};
    legend_text = cat(2,legend_text,legend_text2);
    legend(legend_text)
    legend('hide')
    xlim([0,1]);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Concentration heat map plotting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(14)
subaxis(sub_pos(1),sub_pos(2),sub_pos(3), 'Spacing', 0, 'PaddingRight',0.01, 'PaddingLeft',0.02, 'PaddingTop',0.045, 'PaddingBottom',0.04, 'Margin', 0.01);
    cmap = [[0,0,0];colormap];
    cmax = max(nucleus_protein_profile_ab(:,2));
    cmin = 2e-9;
    cind_profile = ceil((size(cmap,1)-2)*(log10(nucleus_protein_profile_ab(:,2))-log10(cmin))/(log10(cmax)-log10(cmin)))+2;
    cind_profile(cind_profile < 2) = 2;
    cind_profile(cind_profile > size(cmap,1)) = size(cmap,1);
    cind_profile = [1;cind_profile];
    
    nu_cind = cind_profile(mask_stack+1);
    nu_cind2D = round(sum(nu_cind-1,3)./(sum(mask_stack > 0,3)+(max(mask_stack,[],3) == 0)))+1;
    
    nu_image = zeros([size(nu_cind2D),3]);
    for i_channel = 1:3
        c_temp = cmap(:,i_channel);
        nu_image(:,:,i_channel) = c_temp(nu_cind2D);
    end
    imshow(nu_image)
    
    ctick0 = [2,5,10,20,50,100,200]*1e-9;
    ctick_label0 = {'2nM','5nM','10nM','20nM','50nM','100nM','200nM'};
    ctick = [ctick0(ctick0 < cmax),cmax];
    ctick_position = (size(cmap,1)-2)*((log10(ctick)-log10(cmin))/(log10(cmax)-log10(cmin)))+2;
    ctick_label = cat(2,ctick_label0(ctick0 < cmax),{[num2str(cmax/1e-9),'nM']});
    colorbar('YTick',ctick_position,'YTickLabel',ctick_label)
    title([image_folder,' (',p_name,'), cycle = ',num2str(N_cycle)],'Interpreter','none')

    
figure(140)
    clf
    imshow(nu_image)
    colorbar('YTick',ctick_position,'YTickLabel',ctick_label)
    title([image_folder,' (',p_name,'), cycle = ',num2str(N_cycle)],'Interpreter','none')

%%%---------------------------------------------------------------------%%%
%%%---------------------------------------------------------------------%%%
    
figure(16)
subaxis(sub_pos(1),sub_pos(2),sub_pos(3), 'Spacing', 0, 'PaddingRight',0.01, 'PaddingLeft',0.02, 'PaddingTop',0.045, 'PaddingBottom',0.04, 'Margin', 0.01);
    cmap = [[0,0,0];colormap];
    cmax = max(nucleus_protein2_profile_ab(:,2));
    cmin = 2e-9;
    cind_profile = ceil((size(cmap,1)-2)*(log10(nucleus_protein2_profile_ab(:,2))-log10(cmin))/(log10(cmax)-log10(cmin)))+2;
    cind_profile(cind_profile < 2) = 2;
    cind_profile(cind_profile > size(cmap,1)) = size(cmap,1);
    cind_profile = [1;cind_profile];
    
    nu_cind = cind_profile(mask_stack+1);
    nu_cind2D = round(sum(nu_cind-1,3)./(sum(mask_stack > 0,3)+(max(mask_stack,[],3) == 0)))+1;
    
    nu_image = zeros([size(nu_cind2D),3]);
    for i_channel = 1:3
        c_temp = cmap(:,i_channel);
        nu_image(:,:,i_channel) = c_temp(nu_cind2D);
    end
    imshow(nu_image)
    
    ctick0 = [2,5,10,20,50,100,200]*1e-9;
    ctick_label0 = {'2nM','5nM','10nM','20nM','50nM','100nM','200nM'};
    ctick = [ctick0(ctick0 < cmax),cmax];
    ctick_position = (size(cmap,1)-2)*((log10(ctick)-log10(cmin))/(log10(cmax)-log10(cmin)))+2;
    ctick_label = cat(2,ctick_label0(ctick0 < cmax),{[num2str(cmax/1e-9),'nM']});
    colorbar('YTick',ctick_position,'YTickLabel',ctick_label)
    title([image_folder,' (',p2_name,'), cycle = ',num2str(N_cycle)],'Interpreter','none')

    
figure(160)
    clf
    imshow(nu_image)
    colorbar('YTick',ctick_position,'YTickLabel',ctick_label)
    title([image_folder,' (',p2_name,'), cycle = ',num2str(N_cycle)],'Interpreter','none')
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Fluctuation plot: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = polyfit(nucleus_protein_profile((nucleus_protein_profile(:,1) < 0.5) & index_true,2),nucleus_protein_profile((nucleus_protein_profile(:,1) < 0.5) & index_true,3),1);
p_post = polyfit(nucleus_protein_profile((nucleus_protein_profile(:,1) > 0.8) & index_true,2),nucleus_protein_profile((nucleus_protein_profile(:,1) > 0.8) & index_true,3),1);
p2 = polyfit(nucleus_protein2_profile((nucleus_protein2_profile(:,1) < 0.5) & index_true,2),nucleus_protein2_profile((nucleus_protein2_profile(:,1) < 0.5) & index_true,3),1);
p2_post = polyfit(nucleus_protein2_profile((nucleus_protein2_profile(:,1) > 0.8) & index_true,2),nucleus_protein2_profile((nucleus_protein2_profile(:,1) > 0.8) & index_true,3),1);

figure(62)
subaxis(sub_pos(1),sub_pos(2),sub_pos(3), 'Spacing', 0, 'PaddingRight',0.01, 'PaddingLeft',0.02, 'PaddingTop',0.045, 'PaddingBottom',0.04, 'Margin', 0.01);
    plot(nucleus_protein_profile(:,2),nucleus_protein_profile(:,3),'b.')
    hold on
    plot(nucleus_protein_profile(index_true,2),nucleus_protein_profile(index_true,3),'g.')
    x0 = xlim;
    plot(x0,p(1)*(x0+p_post(2)/p_post(1)),'r')
    
    plot(nucleus_protein2_profile(:,2),nucleus_protein2_profile(:,3),'y.')
    plot(nucleus_protein2_profile(index_true,2),nucleus_protein2_profile(index_true,3),'c.')
    plot(x0,p2(1)*(x0+p2_post(2)/p2_post(1)),'m')
    hold off
    xlabel('Mean Intensity (A.U.)');
    ylabel('Intensity varience (A.U.)');
    Cmax = max(nucleus_protein_profile_ab(:,2));
    C2max = max(nucleus_protein2_profile_ab(:,2));
    title([image_folder,char(10),', ',p_name,'_max = ',num2str(Cmax,4),' M, ',p2_name,'_max = ',num2str(C2max,4),' M'],'Interpreter','none')
    legend([p_name,' data points'],[p_name,' data points cleaned'],[p_name,' linear fit'],[p2_name,' data points'],[p2_name,' data points cleaned'],[p2_name,' linear fit'])
    legend('hide')

figure(620)
clf
    plot(nucleus_protein_profile(:,2),nucleus_protein_profile(:,3),'b.')
    hold on
    plot(nucleus_protein_profile(index_true,2),nucleus_protein_profile(index_true,3),'g.')
    x0 = xlim;
    plot(x0,p(1)*(x0+p_post(2)/p_post(1)),'r')
    
    plot(nucleus_protein2_profile(:,2),nucleus_protein2_profile(:,3),'y.')
    plot(nucleus_protein2_profile(index_true,2),nucleus_protein2_profile(index_true,3),'c.')
    plot(x0,p2(1)*(x0+p2_post(2)/p2_post(1)),'m')
    hold off
    xlabel('Mean Intensity (A.U.)');
    ylabel('Intensity varience (A.U.)');
    Cmax = max(nucleus_protein_profile_ab(:,2));
    C2max = max(nucleus_protein2_profile_ab(:,2));
    title([image_folder,char(10),', ',p_name,'_max = ',num2str(Cmax,4),' M, ',p2_name,'_max = ',num2str(C2max,4),' M'],'Interpreter','none')
    legend([p_name,' data points'],[p_name,' data points cleaned'],[p_name,' linear fit'],[p2_name,' data points'],[p2_name,' data points cleaned'],[p2_name,' linear fit'])
    legend('hide')
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%










