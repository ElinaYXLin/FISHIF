function [lambda_I,lambda_C,Hill2_I,Hill2_C,pro1_out,pro2_out] = protein2GFP_profile3D_show(nucleus_protein_profile,cytoplasmic_protein_profile,quanti_p,nucleus_protein_profile_ab,foci_data,fake_data,nucleus_protein2_profile,cytoplasmic_protein2_profile,quanti_p2,nucleus_protein2_profile_ab,foci_data2,fake_data2,mask_stack,image_folder,N_cycle,sub_pos,flip_EL,pro_name,ir,h,resolution)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A function to re-plot two protein profiles and fit them to given curves 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Parameter setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fLmin = 0;
fLmax = 1;
z_size = size(mask_stack,3);
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
% hill_C = @(beta,x) beta(3)*beta(2).^beta(1)./(x.^beta(1)+beta(2)^beta(1))+beta(4);
hill_C = @(beta,x) beta(1)*exp(-beta(2)*x)+beta(3);
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
% beta2 = nlinfit(nu_L(mean2_index),nu2_mean(mean2_index),hill_C,[h0,c0,max_mean2-min_mean2,min_mean2]);
beta2 = nlinfit(nu_L(mean2_index),nu2_mean(mean2_index),hill_C,[max_mean2-min_mean2,b2,min_mean2]);
catch
% beta2 = nan(1,4);
beta2 = nan(1,3);
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
Hill2_I = [beta2(1),1./beta2(2),beta2(3),c2max0];
% Hill2_I = [beta2,c2max0];

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
% beta2 = nlinfit(nu_L(mean2_index),nu2_mean(mean2_index),hill_C,[h0,c0,max_mean2-min_mean2,min_mean2]);
beta2 = nlinfit(nu_L(mean2_index),nu2_mean(mean2_index),hill_C,[max_mean2-min_mean2,b2,min_mean2]);
catch
% beta2 = nan(1,4);
beta2 = nan(1,3);
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
Hill2_C = [beta2(1),1./beta2(2),beta2(3),c2max];
% Hill2_C = [beta2,c2max];

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




%% Scaled GFP vs Bcd: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nucleus_profile0 = nucleus_protein_profile;
nucleus2_profile0 = nucleus_protein2_profile;
pn0 = polyfit(nucleus_profile0(index_true,2),nucleus2_profile0(index_true,2),1);
Rn0temp = corrcoef([nucleus_profile0(index_true,2),nucleus2_profile0(index_true,2)]);
Rn02 = Rn0temp(1,2)^2;
n0x = [min(nucleus_profile0(:,2)),max(nucleus_profile0(:,2))];
n0y = pn0(1)*n0x+pn0(2);

cytoplasmic_profile = cytoplasmic_protein_profile(:,2);
cytoplasmic2_profile = cytoplasmic_protein2_profile(:,2);
pc = polyfit(cytoplasmic_profile,cytoplasmic2_profile,1);
Rctemp = corrcoef([cytoplasmic_profile,cytoplasmic2_profile]);
Rc2 = Rctemp(1,2)^2;
cx = [min(cytoplasmic_profile),max(cytoplasmic_profile)];
cy = pc(1)*cx+pc(2);

h_area = sum(h{ir}(:));
foci_data{ir}(:,3) = h_area*foci_data{ir}(:,9).*foci_data{ir}(:,13);
foci_data2{ir}(:,3) = h_area*foci_data2{ir}(:,9).*foci_data2{ir}(:,13);
fake_data{ir}(:,3) = h_area*fake_data{ir}(:,9).*fake_data{ir}(:,13);
fake_data2{ir}(:,3) = h_area*fake_data2{ir}(:,9).*fake_data2{ir}(:,13);
f_index = (foci_data{ir}(:,7)>=fLmin) & (foci_data{ir}(:,7)<=fLmax) & (~isnan(foci_data{ir}(:,4))) & (~isnan(foci_data{ir}(:,3))) & (~isnan(fake_data{ir}(:,4))) & (~isnan(fake_data{ir}(:,3))) & (foci_data{ir}(:,8) > 1) & (foci_data{ir}(:,8) < z_size);
f2_index = (foci_data2{ir}(:,7)>=fLmin) & (foci_data2{ir}(:,7)<=fLmax) & (~isnan(foci_data2{ir}(:,4))) & (~isnan(foci_data2{ir}(:,3))) & (~isnan(fake_data2{ir}(:,4))) & (~isnan(fake_data2{ir}(:,3))) & (foci_data2{ir}(:,8) > 1) & (foci_data2{ir}(:,8) < z_size);

nucleus_profile = foci_data{ir}(f_index & f2_index,12)*quanti_p(1)-quanti_p(2);
nucleus2_profile = foci_data2{ir}(f_index & f2_index,12)*quanti_p2(1)-quanti_p2(2);
pn = polyfit(nucleus_profile,nucleus2_profile,1);
Rntemp = corrcoef([nucleus_profile,nucleus2_profile]);
Rn2 = Rntemp(1,2)^2;
nx = [min(nucleus_profile),max(nucleus_profile)];
ny = pn(1)*nx+pn(2);

foci_profile = foci_data{ir}(f_index & f2_index,4)./h_area./foci_data{ir}(f_index & f2_index,9)*quanti_p(1)-quanti_p(2);
foci2_profile = foci_data2{ir}(f_index & f2_index,4)./h_area./foci_data2{ir}(f_index & f2_index,9)*quanti_p2(1)-quanti_p2(2);
pf = polyfit(foci_profile,foci2_profile,1);
Rftemp = corrcoef([foci_profile,foci2_profile]);
Rf2 = Rftemp(1,2)^2;
fx = [min(foci_profile),max(foci_profile)];
fy = pf(1)*fx+pf(2);

enrich_profile = (foci_data{ir}(f_index & f2_index,4)-foci_data{ir}(f_index & f2_index,3))./h_area./foci_data{ir}(f_index & f2_index,9)./foci_data{ir}(f_index & f2_index,16)*quanti_p(1);
enrich2_profile = (foci_data2{ir}(f_index & f2_index,4)-foci_data2{ir}(f_index & f2_index,3))./h_area./foci_data2{ir}(f_index & f2_index,9)./foci_data{ir}(f_index & f2_index,16)*quanti_p2(1);
pe = polyfit(enrich_profile,enrich2_profile,1);
Retemp = corrcoef([enrich_profile,enrich2_profile]);
Re2 = Retemp(1,2)^2;
ex = [min(enrich_profile),max(enrich_profile)];
ey = pe(1)*ex+pe(2);

fake_profile = fake_data{ir}(f_index & f2_index,4)./h_area./fake_data{ir}(f_index & f2_index,9)*quanti_p(1)-quanti_p(2);
fake2_profile = fake_data2{ir}(f_index & f2_index,4)./h_area./fake_data2{ir}(f_index & f2_index,9)*quanti_p2(1)-quanti_p2(2);
pfa = polyfit(fake_profile,fake2_profile,1);
Rfatemp = corrcoef([fake_profile,fake2_profile]);
Rfa2 = Rfatemp(1,2)^2;
fax = [min(fake_profile),max(fake_profile)];
fay = pfa(1)*fax+pfa(2);

enrichf_profile = (fake_data{ir}(f_index & f2_index,4)-fake_data{ir}(f_index & f2_index,3))./h_area./fake_data{ir}(f_index & f2_index,9)./fake_data{ir}(f_index & f2_index,16)*quanti_p(1);
enrichf2_profile = (fake_data2{ir}(f_index & f2_index,4)-fake_data2{ir}(f_index & f2_index,3))./h_area./fake_data2{ir}(f_index & f2_index,9)./fake_data{ir}(f_index & f2_index,16)*quanti_p2(1);
pef = polyfit(enrichf_profile,enrichf2_profile,1);
Reftemp = corrcoef([enrichf_profile,enrichf2_profile]);
Ref2 = Reftemp(1,2)^2;
efx = [min(enrichf_profile),max(enrichf_profile)];
efy = pef(1)*efx+pef(2);

figure(23)
subaxis(sub_pos(1),sub_pos(2),sub_pos(3), 'Spacing', 0, 'PaddingRight',0.01, 'PaddingLeft',0.02, 'PaddingTop',0.045, 'PaddingBottom',0.04, 'Margin', 0.01);
    plot(nucleus_profile0(index_true,2),nucleus2_profile0(index_true,2),'Marker','.','Color',ccode(1,:),'LineStyle','none','DisplayName','Nucleus mean intensity0');
    hold on
    plot(n0x,n0y,'Color',ccode2(1,:),'LineStyle','--','DisplayName',['Nuclear0 linear fit: y = ',num2str(pn0(1),'%5.2f'),' * x + ',num2str(pn0(2),'%5.2f'),', R^2 = ',num2str(Rn02,'%5.2f')]);

    plot(nucleus_profile,nucleus2_profile,'Marker','+','Color',ccode(2,:),'LineStyle','none','DisplayName','Nucleus mean intensity');
    plot(nx,ny,'Color',ccode2(2,:),'DisplayName',['Nuclear linear fit: y = ',num2str(pn(1),'%5.2f'),' * x + ',num2str(pn(2),'%5.2f'),', R^2 = ',num2str(Rn2,'%5.2f')]);

    plot(cytoplasmic_profile,cytoplasmic2_profile,'Marker','x','Color',ccode(3,:),'LineStyle','none','DisplayName','Cytoplasmic mean intensity');
    plot(cx,cy,'Color',ccode2(3,:),'LineStyle','-.','DisplayName',['Nuclear linear fit: y = ',num2str(pc(1),'%5.2f'),' * x + ',num2str(pc(2),'%5.2f'),', R^2 = ',num2str(Rc2,'%5.2f')]);

    plot(foci_profile,foci2_profile,'Marker','*','Color',ccode(4,:),'LineStyle','none','DisplayName','Foci mean intensity');
    plot(fx,fy,'Color',ccode2(4,:),'DisplayName',['Foci linear fit: y = ',num2str(pf(1),'%5.2f'),' * x + ',num2str(pf(2),'%5.2f'),', R^2 = ',num2str(Rf2,'%5.2f')]);
    
    plot(enrich_profile,enrich2_profile,'Marker','o','Color',ccode(5,:),'LineStyle','none','DisplayName','Enrichment mean intensity');
    plot(ex,ey,'Color',ccode2(5,:),'DisplayName',['Enrich linear fit: y = ',num2str(pe(1),'%5.2f'),' * x + ',num2str(pe(2),'%5.2f'),', R^2 = ',num2str(Re2,'%5.2f')]);

    plot(fake_profile,fake2_profile,'Marker','*','Color',ccode2(6,:),'LineStyle','none','DisplayName','Fake mean intensity');
    plot(fax,fay,'Color',ccode2(6,:),'DisplayName',['Fake linear fit: y = ',num2str(pfa(1),'%5.2f'),' * x + ',num2str(pfa(2),'%5.2f'),', R^2 = ',num2str(Rfa2,'%5.2f')]);
    
    plot(enrichf_profile,enrichf2_profile,'Marker','o','Color',ccode(7,:),'LineStyle','none','DisplayName','Fake enrichment mean intensity');
    plot(efx,efy,'Color',ccode2(7,:),'DisplayName',['Fake enrich linear fit: y = ',num2str(pef(1),'%5.2f'),' * x + ',num2str(pef(2),'%5.2f'),', R^2 = ',num2str(Ref2,'%5.2f')]);

    hold off
    title([image_folder,', cycle = ',num2str(N_cycle)],'Interpreter','none')
    xlabel([p_name,' Intensity'])
    ylabel([p2_name,' Intensity'])
%     legend('show')
    
figure(230)
clf
    plot(nucleus_profile0(index_true,2),nucleus2_profile0(index_true,2),'Marker','.','Color',ccode(1,:),'LineStyle','none','DisplayName','Nucleus mean intensity0');
    hold on
    plot(n0x,n0y,'Color',ccode2(1,:),'LineStyle','--','DisplayName',['Nuclear0 linear fit: y = ',num2str(pn0(1),'%5.2f'),' * x + ',num2str(pn0(2),'%5.2f'),', R^2 = ',num2str(Rn02,'%5.2f')]);

    plot(nucleus_profile,nucleus2_profile,'Marker','+','Color',ccode(2,:),'LineStyle','none','DisplayName','Nucleus mean intensity');
    plot(nx,ny,'Color',ccode2(2,:),'DisplayName',['Nuclear linear fit: y = ',num2str(pn(1),'%5.2f'),' * x + ',num2str(pn(2),'%5.2f'),', R^2 = ',num2str(Rn2,'%5.2f')]);

    plot(cytoplasmic_profile,cytoplasmic2_profile,'Marker','x','Color',ccode(3,:),'LineStyle','none','DisplayName','Cytoplasmic mean intensity');
    plot(cx,cy,'Color',ccode2(3,:),'LineStyle','-.','DisplayName',['Nuclear linear fit: y = ',num2str(pc(1),'%5.2f'),' * x + ',num2str(pc(2),'%5.2f'),', R^2 = ',num2str(Rc2,'%5.2f')]);

    plot(foci_profile,foci2_profile,'Marker','*','Color',ccode(4,:),'LineStyle','none','DisplayName','Foci mean intensity');
    plot(fx,fy,'Color',ccode2(4,:),'DisplayName',['Foci linear fit: y = ',num2str(pf(1),'%5.2f'),' * x + ',num2str(pf(2),'%5.2f'),', R^2 = ',num2str(Rf2,'%5.2f')]);
    
    plot(enrich_profile,enrich2_profile,'Marker','o','Color',ccode(5,:),'LineStyle','none','DisplayName','Enrichment mean intensity');
    plot(ex,ey,'Color',ccode2(5,:),'DisplayName',['Enrich linear fit: y = ',num2str(pe(1),'%5.2f'),' * x + ',num2str(pe(2),'%5.2f'),', R^2 = ',num2str(Re2,'%5.2f')]);

    plot(fake_profile,fake2_profile,'Marker','*','Color',ccode2(6,:),'LineStyle','none','DisplayName','Fake mean intensity');
    plot(fax,fay,'Color',ccode2(6,:),'DisplayName',['Fake linear fit: y = ',num2str(pfa(1),'%5.2f'),' * x + ',num2str(pfa(2),'%5.2f'),', R^2 = ',num2str(Rfa2,'%5.2f')]);
    
    plot(enrichf_profile,enrichf2_profile,'Marker','o','Color',ccode(7,:),'LineStyle','none','DisplayName','Fake enrichment mean intensity');
    plot(efx,efy,'Color',ccode2(7,:),'DisplayName',['Fake enrich linear fit: y = ',num2str(pef(1),'%5.2f'),' * x + ',num2str(pef(2),'%5.2f'),', R^2 = ',num2str(Ref2,'%5.2f')]);

    hold off
    title([image_folder,', cycle = ',num2str(N_cycle)],'Interpreter','none')
    xlabel([p_name,' Intensity'])
    ylabel([p2_name,' Intensity'])
    legend('show')
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








