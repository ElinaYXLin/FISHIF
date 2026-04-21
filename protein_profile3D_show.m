function [lambda_I,lambda_C,pro_out] = protein_profile3D_show(nucleus_protein_profile,cytoplasmic_protein_profile,quanti_p,nucleus_protein_profile_ab,mask_stack,image_folder,N_cycle,sub_pos,flip_EL)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A function to re-plot protein profile and fit it to a given curve %%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Parameter setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ccode = [0.7,0.7,0.7 ; 0,0,0.3 ; 0,0.3,0 ; 0,0,0; 0,0,0.7 ; 0,0.7,0.7 ; 1,0,0];
Lmin = 0;
Lnbin = 0.05;
Lmax = 1;
Lnwindow = 0.05;
nu_L = Lmin:Lnbin:Lmax;
nu_mean = zeros(size(nu_L));
nu_std = zeros(size(nu_L));
standard_dorsal = 'Standard protein/Bcd_dorsal.xls';
standard_ventral = 'Standard protein/Bcd_ventral.xls';
bcd_dxy = xlsread(standard_dorsal);
bcd_vxy = xlsread(standard_ventral);

index_true = logical(nucleus_protein_profile(:,4));
if flip_EL
    nucleus_protein_profile(:,1) = 1-nucleus_protein_profile(:,1);
    cytoplasmic_protein_profile(:,1) = 1-cytoplasmic_protein_profile(:,1);
    nucleus_protein_profile_ab(:,1) = 1-nucleus_protein_profile_ab(:,1);
end
exp_C = @(beta,x) beta(1)*exp(-beta(2)*x)+beta(3);
b2 = 5;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Raw profile plotting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nucleus_profile = nucleus_protein_profile;
for Lcenter = 1:length(nu_L)
    nu_map = (nucleus_profile(index_true,1) >= nu_L(Lcenter)-Lnwindow) & (nucleus_profile(index_true,1) <= nu_L(Lcenter)+Lnwindow);
    nu_mean(Lcenter) = mean(nucleus_profile(nu_map,2));
    nu_std(Lcenter) = std0(nucleus_profile(nu_map,2));
end

mean_index = (~isnan(nu_mean));
[max_mean,mean_start] = max(nu_mean);
min_mean = min(nu_mean);
if mean_start > 1
    mean_index(1:(mean_start-1)) = false;
end
beta1 = nlinfit(nu_L(mean_index),nu_mean(mean_index),exp_C,[max_mean-min_mean,b2,min_mean]);

d_max = mean(bcd_dxy(1:5,2));
d_min = mean(bcd_dxy(end-4:end,2));
d_x0 = mean(bcd_dxy(1:5,1));
d_y = (bcd_dxy(:,2)-d_min)/(d_max-d_min)*exp_C(beta1,d_x0)+min_mean;

v_max = mean(bcd_vxy(1:5,2));
v_min = mean(bcd_vxy(end-4:end,2));
v_x0 = mean(bcd_vxy(1:5,1));
v_y = (bcd_vxy(:,2)-v_min)/(v_max-v_min)*exp_C(beta1,v_x0)+min_mean;

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

    hold off
    title([image_folder,', cycle = ',num2str(N_cycle)],'Interpreter','none')
    xlabel('AP axis (normalized)')
    ylabel('Intensity (A.U.)')
    legend('Nucleus concentration','Good nucleus concentration','Cytoplasmic concentration','Averaged nuclear profile','Standard dorsal profile','Standard ventral profile',['Fitting result: lambda = ',num2str(1./beta1(2))])
    legend('hide')
    xlim([0,1]);
cmax0 = max(nucleus_protein_profile(:,2));
lambda_I = [beta1(1),1./beta1(2),beta1(3),cmax0];

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

    hold off
    title([image_folder,', cycle = ',num2str(N_cycle)],'Interpreter','none')
    xlabel('AP axis (normalized)')
    ylabel('Intensity (A.U.)')
    legend('Nucleus concentration','Good nucleus concentration','Cytoplasmic concentration','Averaged nuclear profile','Standard dorsal profile','Standard ventral profile',['Fitting result: lambda = ',num2str(1./beta1(2))])
    legend('hide')
    xlim([0,1]);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Scaled profile plotting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nucleus_profile = nucleus_protein_profile_ab;
for Lcenter = 1:length(nu_L)
    nu_map = (nucleus_profile(index_true,1) >= nu_L(Lcenter)-Lnwindow) & (nucleus_profile(index_true,1) <= nu_L(Lcenter)+Lnwindow);
    nu_mean(Lcenter) = mean(nucleus_profile(nu_map,2));
    nu_std(Lcenter) = std0(nucleus_profile(nu_map,2));
end
pro_out = [nu_L',nu_mean',nu_std'];

mean_index = (~isnan(nu_mean));
[max_mean,mean_start] = max(nu_mean);
min_mean = min(nu_mean);
if mean_start > 1
    mean_index(1:(mean_start-1)) = false;
end
beta1 = nlinfit(nu_L(mean_index),nu_mean(mean_index),exp_C,[max_mean-min_mean,b2,min_mean]);

d_max = mean(bcd_dxy(1:5,2));
d_min = mean(bcd_dxy(end-4:end,2));
d_x0 = mean(bcd_dxy(1:5,1));
d_y = (bcd_dxy(:,2)-d_min)/(d_max-d_min)*exp_C(beta1,d_x0)+min_mean;

v_max = mean(bcd_vxy(1:5,2));
v_min = mean(bcd_vxy(end-4:end,2));
v_x0 = mean(bcd_vxy(1:5,1));
v_y = (bcd_vxy(:,2)-v_min)/(v_max-v_min)*exp_C(beta1,v_x0)+min_mean;

figure(22)
subaxis(sub_pos(1),sub_pos(2),sub_pos(3), 'Spacing', 0, 'PaddingRight',0.01, 'PaddingLeft',0.02, 'PaddingTop',0.045, 'PaddingBottom',0.04, 'Margin', 0.01);
    plot(nucleus_profile(:,1),nucleus_profile(:,2),'Marker','.','Color',ccode(1,:),'LineStyle','none');
    hold on
    plot(nucleus_profile(index_true,1),nucleus_profile(index_true,2),'Marker','.','Color',ccode(2,:),'LineStyle','none');
    errorbar(nu_L,nu_mean,nu_std,'Marker','.','Color',ccode(4,:),'LineStyle','none');
    plot(bcd_dxy(:,1),d_y,'Marker','X','Color',ccode(5,:));
    plot(bcd_vxy(:,1),v_y,'Marker','X','Color',ccode(6,:));
    plot(nu_L,exp_C(beta1,nu_L),'Color',ccode(7,:))

    hold off
    title([image_folder,', cycle = ',num2str(N_cycle)],'Interpreter','none')
    xlabel('AP axis (normalized)')
    ylabel('Concentration (M)')
    legend('Nucleus concentration','Good nucleus concentration','Averaged nuclear profile','Standard dorsal profile','Standard ventral profile',['Fitting result: lambda = ',num2str(1./beta1(2))])
    legend('hide')
    xlim([0,1]);
cmax = max(nucleus_protein_profile_ab(:,2));
lambda_C = [beta1(1),1./beta1(2),beta1(3),cmax];

figure(220)
clf
    plot(nucleus_profile(:,1),nucleus_profile(:,2),'Marker','.','Color',ccode(1,:),'LineStyle','none');
    hold on
    plot(nucleus_profile(index_true,1),nucleus_profile(index_true,2),'Marker','.','Color',ccode(2,:),'LineStyle','none');
    errorbar(nu_L,nu_mean,nu_std,'Marker','.','Color',ccode(4,:),'LineStyle','none');
    plot(bcd_dxy(:,1),d_y,'Marker','X','Color',ccode(5,:));
    plot(bcd_vxy(:,1),v_y,'Marker','X','Color',ccode(6,:));
    plot(nu_L,exp_C(beta1,nu_L),'Color',ccode(7,:))

    hold off
    title([image_folder,', cycle = ',num2str(N_cycle)],'Interpreter','none')
    xlabel('AP axis (normalized)')
    ylabel('Concentration (M)')
    legend('Nucleus concentration','Good nucleus concentration','Averaged nuclear profile','Standard dorsal profile','Standard ventral profile',['Fitting result: lambda = ',num2str(1./beta1(2))])
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
    title([image_folder,', cycle = ',num2str(N_cycle)],'Interpreter','none')

    
figure(140)
    clf
    imshow(nu_image)
    colorbar('YTick',ctick_position,'YTickLabel',ctick_label)
    title([image_folder,', cycle = ',num2str(N_cycle)],'Interpreter','none')
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Fluctuation plot: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = polyfit(nucleus_protein_profile((nucleus_protein_profile(:,1) < 0.5) & index_true,2),nucleus_protein_profile((nucleus_protein_profile(:,1) < 0.5) & index_true,3),1);
p_post = polyfit(nucleus_protein_profile((nucleus_protein_profile(:,1) > 0.8) & index_true,2),nucleus_protein_profile((nucleus_protein_profile(:,1) > 0.8) & index_true,3),1);
mean_post = mean(nucleus_protein_profile((nucleus_protein_profile(:,1) >= 0.8) & index_true,2));

figure(62)
subaxis(sub_pos(1),sub_pos(2),sub_pos(3), 'Spacing', 0, 'PaddingRight',0.01, 'PaddingLeft',0.02, 'PaddingTop',0.045, 'PaddingBottom',0.04, 'Margin', 0.01);
    plot(nucleus_protein_profile(:,2),nucleus_protein_profile(:,3),'b.')
    hold on
    plot(nucleus_protein_profile(index_true,2),nucleus_protein_profile(index_true,3),'g.')
    x0 = xlim;
    plot(x0,p(1)*(x0+p_post(2)/p_post(1)),'r')
    plot(mean_post*[1,1],ylim*0.9,'k--')
    hold off
    xlabel('Mean Intensity (A.U.)');
    ylabel('Intensity varience (A.U.)');
    Cmax = max(nucleus_protein_profile_ab(:,2));
    title([image_folder,char(10),', Cmax = ',num2str(Cmax,4),' M'],'Interpreter','none')
    legend('Data points','Data points cleaned','Linear fit','Lower limit')
    legend('hide')

figure(620)
clf
    plot(nucleus_protein_profile(:,2),nucleus_protein_profile(:,3),'b.')
    hold on
    plot(nucleus_protein_profile(index_true,2),nucleus_protein_profile(index_true,3),'g.')
    x0 = xlim;
    plot(x0,p(1)*(x0+p_post(2)/p_post(1)),'r')
    plot(mean_post*[1,1],ylim*0.9,'k--')
    hold off
    xlabel('Mean Intensity (A.U.)');
    ylabel('Intensity varience (A.U.)');
    Cmax = max(nucleus_protein_profile_ab(:,2));
    title([image_folder,char(10),', Cmax = ',num2str(Cmax,4),' M'],'Interpreter','none')
    legend('Data points','Data points cleaned','Linear fit','Lower limit')
    legend('hide')
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%










