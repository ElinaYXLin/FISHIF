function [nucleus_profile,cytoplasmic_profile] = GFP_profile(seg_bw,cyto_bw,max_image,signal_channel,DAPI_channel,WGA_channel,image_folder,N_cycle)

%% Nucleus/cytoplasmic protein profile calculation: %%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Parameter setting and initiation: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
standard_dorsal = 'Standard protein/Bcd_dorsal.xls';
standard_ventral = 'Standard protein/Bcd_ventral.xls';
subtp = 'background subtraction';
bin2 = 1;
Nb2 = 8;
bcd_dxy = xlsread(standard_dorsal);
bcd_vxy = xlsread(standard_ventral);
%nucleus_profile = zeros(0);
%cytoplasmic_profile = zeros(0);
Lmin = 0;
Lcbin = 0.02;
Lmax = 1;
Lcwindow = 0.02;
cyto_L = Lmin:Lcbin:Lmax;
cyto_mean = zeros(size(cyto_L));

Lmin = 0;
Lnbin = 0.05;
Lmax = 1;
Lnwindow = 0.05;
nu_L = Lmin:Lnbin:Lmax;
nu_mean = zeros(size(nu_L));
nu_std = zeros(size(nu_L));

max_image = double(max_image);
start_limit = 0.3;
end_limit = 1-start_limit;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Nucleus region recognition/processing: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
signal_prop = regionprops(imerode(seg_bw, strel('disk',4)),max_image(:,:,signal_channel),'MeanIntensity','Centroid');

%%% Normalized coordinate calculation: %%%=================================
if size(signal_prop,1)>0
    xy = [signal_prop.Centroid];
    nucleus_xy = [xy(1:2:length(xy)-1);xy(2:2:length(xy))]';
    
    distXY = pdist2(nucleus_xy,nucleus_xy);
    axis_length = max(max(distXY));
    [extreme0,extreme1] = find(distXY >= axis_length,1); %%% find the anterior and posterior pole nucleus
    L2_extreme = dot((nucleus_xy(extreme1,:)-nucleus_xy(extreme0,:)),(nucleus_xy(extreme1,:)-nucleus_xy(extreme0,:)),2);
    nucleus_distance = dot((nucleus_xy-repmat(nucleus_xy(extreme0,:),size(nucleus_xy,1),1)),repmat((nucleus_xy(extreme1,:)-nucleus_xy(extreme0,:)),size(nucleus_xy,1),1),2)/L2_extreme;
    
    nucleus_profile = [nucleus_distance,[signal_prop.MeanIntensity]'];
end
%%% =======================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Cytoplasmic region recognition/processing: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[ylim,xlim] = size(cyto_bw);
xmap = repmat([1:xlim],ylim,1)-nucleus_xy(extreme0,1);
ymap = repmat([1:ylim]',1,xlim)-nucleus_xy(extreme0,2);
Lmap = (xmap*(nucleus_xy(extreme1,1)-nucleus_xy(extreme0,1))+ymap*(nucleus_xy(extreme1,2)-nucleus_xy(extreme0,2)))/L2_extreme;

for Lcenter = 1:length(cyto_L)
    cyto_map = (Lmap >= cyto_L(Lcenter)-Lcwindow) & (Lmap <= cyto_L(Lcenter)+Lcwindow) & cyto_bw;
    cyto_mean(Lcenter) = sum(sum(cyto_map.*max_image(:,:,signal_channel)))/sum(sum(cyto_map));
end
cytoplasmic_profile = [cyto_L',cyto_mean'];

for Lcenter = 1:length(nu_L)
    nu_map = (nucleus_profile(:,1) >= nu_L(Lcenter)-Lnwindow) & (nucleus_profile(:,1) <= nu_L(Lcenter)+Lnwindow);
    nu_mean(Lcenter) = mean(nucleus_profile(nu_map,2));
    nu_std(Lcenter) = std(nucleus_profile(nu_map,2));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Standard profile reshaping: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b2 = -3;
nu_x = nucleus_profile(:,1);
standard_x = [bcd_dxy(:,1);bcd_vxy(:,1)];
if mean(nucleus_profile(nucleus_profile(:,1) <= start_limit,2)) < mean(nucleus_profile(nucleus_profile(:,1) >= end_limit,2))
    nu_x = 1-nucleus_profile(:,1);
    bcd_dxy(:,1) = 1-bcd_dxy(:,1);
    bcd_vxy(:,1) = 1-bcd_vxy(:,1);
end
beta1 = nlinfit(nu_x,nucleus_profile(:,2),@exp_C,[(max(nucleus_profile(:,2))-min(nucleus_profile(:,2))),b2,min(nucleus_profile(:,2))]);
%protein_min = beta1(3);
%protein_max = beta1(1)+beta1(3);
protein_min = min(nu_mean);
protein_max = max(nu_mean);

beta2 = nlinfit(standard_x,[bcd_dxy(:,2);bcd_vxy(:,2)],@exp_C,[(max([bcd_dxy(:,2);bcd_vxy(:,2)])-min([bcd_dxy(:,2);bcd_vxy(:,2)])),b2,min([bcd_dxy(:,2);bcd_vxy(:,2)])]);
%standard_min = beta2(3);
%standard_max = beta2(1)+beta2(3);
standard_min = min(min(bcd_dxy(:,2)),min(bcd_vxy(:,2)));
standard_max = max(max(bcd_dxy(:,2)),min(bcd_vxy(:,2)));

bcd_dxy(:,2) = (bcd_dxy(:,2)-standard_min)*(protein_max-protein_min)/(standard_max-standard_min)+protein_min;
bcd_vxy(:,2) = (bcd_vxy(:,2)-standard_min)*(protein_max-protein_min)/(standard_max-standard_min)+protein_min;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Image output: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(22)
    clf
    if isempty(subtp);
        sub_nu = 0;
        sub_cyto = 0;
    else
        sub_nu = min(nucleus_profile(:,2));
        sub_cyto = min(cytoplasmic_profile(:,2));
    end
    plot(nucleus_profile(:,1),nucleus_profile(:,2)-sub_nu,'o',cytoplasmic_profile(:,1),cytoplasmic_profile(:,2)-sub_cyto,'*',bcd_dxy(:,1),bcd_dxy(:,2)-sub_nu,'X',bcd_vxy(:,1),bcd_vxy(:,2)-sub_nu,'X')
    hold on
    errorbar(nu_L,nu_mean-sub_nu,nu_std,'k');
    hold off
    title(['GFP concentration profile: ',image_folder,', cycle = ',num2str(N_cycle),' (',subtp,')'])
    xlabel('AP axis (normalized)')
    ylabel('Concentration (A.U.)')
    legend('Nucleus concentration','Cytoplasmic concentration','Standard dorsal profile','Standard ventral profile','Averaged nuclear profile')

figure(24)
    clf
    nulabel = bwlabel(seg_bw)+1;
    nuintensity = [-inf;log2(nucleus_profile(:,2)-sub_nu)];
    imagesc(nuintensity(nulabel))%,[min(nucleus_profile(:,2)),max(nucleus_profile(:,2))])
    %axis off
    prolim = max(nuintensity);
    %colormap(jet(1024));
    hcb = colorbar;
    set(hcb,'YTickMode','manual') 
    set(hcb,'YTick',prolim+[-(Nb2-1):1:0]*bin2) 
    nylabel = cell(1,Nb2);
    for ib2 = 1:Nb2
        nylabel{ib2} = ['1/',num2str(2^((Nb2-ib2)*bin2))];
    end
    set(hcb,'YTickLabel',nylabel);
    title(['Nucleus GFP concentration plot: ',image_folder,', cycle = ',num2str(N_cycle),' (',subtp,')'])