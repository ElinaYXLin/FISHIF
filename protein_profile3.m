function [nucleus_profile,cytoplasmic_profile,quanti_p] = protein_profile3(seg_bw,cyto_bw,max_image,signal_channel,mask_stack,signal_stack,image_folder,N_cycle,resolution,varargin)


%% Nucleus/cytoplasmic protein profile calculation: %%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Parameter setting and initiation: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
standard_dorsal = 'Standard protein/Bcd_dorsal.xls';
standard_ventral = 'Standard protein/Bcd_ventral.xls';
subtp = 'background subtraction';
bin2 = 1;
Nb2 = 8;
r_edge = 10;
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

if isempty(varargin)
    typetxt = 'Protein';
    n1 = 2;
    n2 = 14;
else
    typetxt = varargin{1};
    n1 = varargin{2};
    n2 = varargin{3};
end

sigma = 1.2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Nucleus region recognition/processing: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
temp = zeros(0);
Inten = zeros(0);
Inten2 = zeros(0);
mask2 = bwlabel(seg_bw);

for I_layer = 1:size(mask_stack,3)
    sg_prop = regionprops((imerode(mask_stack(:,:,I_layer), strel('disk',r_edge))).*mask2,signal_stack(:,:,I_layer),'MeanIntensity','PixelValues');
    sg2_prop = regionprops((imerode(mask_stack(:,:,I_layer), strel('disk',r_edge))).*mask2,signal_stack(:,:,I_layer).^2,'MeanIntensity');
    if I_layer == 1
        temp0 = cell(size(mask_stack,3),max(max(mask2)));
        temp = zeros(size(mask_stack,3),max(max(mask2)));
        temp2 = zeros(size(mask_stack,3),max(max(mask2)));
    end
    if ~isempty(sg_prop)
        temp0(I_layer,1:length(sg2_prop)) = {sg_prop.PixelValues};  %%%%% Binomial partition preparation
        temp(I_layer,1:length(sg2_prop)) = [sg_prop.MeanIntensity];
        temp2(I_layer,1:length(sg2_prop)) = [sg2_prop.MeanIntensity]-[sg_prop.MeanIntensity].^2;
    end
end
% sg_prop = regionprops((imerode(seg_bw, strel('disk',r_edge))).*mask2,'Centroid');
sg_prop = regionprops((bwmorph(seg_bw, 'thin',r_edge)).*mask2,'Centroid');
centr = cell2mat({sg_prop.Centroid}');
temp(isnan(temp)) = 0;
[Inten,max_index] = max(temp,[],1);
Inten2 = temp2(max_index+size(temp,1)*[0:(size(temp,2)-1)]);
temp00 = temp0(max_index+size(temp,1)*[0:(size(temp,2)-1)]);

%%%%% Binomial partition analysis: %%%%% ----------------------------------
dF = zeros(length(temp00),1);
FF = zeros(length(temp00),1);
SF = zeros(length(temp00),1);
for I_region = 1:length(temp00)
    Inten_value = temp00{I_region};
    size_region = size(Inten_value,1);
    if floor(size_region/2) < ceil(size_region/2)
        size_region = size_region-1;
    end
    dF(I_region) = ((sum(Inten_value(1:(size_region/2)))-sum(Inten_value((size_region/2+1):size_region)))/2).^2;
    FF(I_region) = sum(Inten_value(1:size_region));
    SF(I_region) = size_region;
end
%%%%% ---------------------------------------------------------------------


%%% Normalized coordinate calculation: %%%=================================
if ~isempty(centr)
    xy = centr;
    nucleus_xy = xy(:,1:2);
    if size(seg_bw,2)/size(seg_bw,1) > 1.4
        distXY = pdist2(nucleus_xy,nucleus_xy);
        axis_length = max(max(distXY));
        [extreme0,extreme1] = find(distXY >= axis_length,1); %%% find the anterior and posterior pole nucleus
        x0 = nucleus_xy(extreme0,1);
        y0 = nucleus_xy(extreme0,2);
        x1 = nucleus_xy(extreme1,1);
        y1 = nucleus_xy(extreme1,2);
        L2_extreme = (x1-x0)^2+(y1-y0)^2;
    else
        I1 = nucleus_xy(:,1) < size(seg_bw,2)/2;
        I2 = nucleus_xy(:,1) >= size(seg_bw,2)/2;
        if std(nucleus_xy(I1,2)) <= std(nucleus_xy(I2,2))
            [x0,I0] = min(nucleus_xy(:,1));
            y0 = nucleus_xy(I0,2);
            x1 = 2*size(seg_bw,2)-x0;
            y1 = y0;
            axis_length = x1-x0;
            L2_extreme = axis_length^2;
        else
            [x0,I0] = max(nucleus_xy(:,1));
            y0 = nucleus_xy(I0,2);
            x1 = 2*1-x0;
            y1 = y0;
            axis_length = x1-x0;
            L2_extreme = axis_length^2;
        end
    end
    nucleus_distance = dot((nucleus_xy-repmat([x0,y0],size(nucleus_xy,1),1)),repmat([(x1-x0),(y1-y0)],size(nucleus_xy,1),1),2)/L2_extreme;
    
    nucleus_profile = [nucleus_distance,Inten'];
    fluc_profile = [nucleus_distance,Inten2'];   %%% intensity variance profile
end
%%% =======================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Cytoplasmic region recognition/processing: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[ylim0,xlim0] = size(cyto_bw);
xmap = repmat([1:xlim0],ylim0,1)-x0;
ymap = repmat([1:ylim0]',1,xlim0)-y0;
clear xlim0 ylim0
Lmap = (xmap*(x1-x0)+ymap*(y1-y0))/L2_extreme;

for Lcenter = 1:length(cyto_L)
    cyto_map = (Lmap >= cyto_L(Lcenter)-Lcwindow) & (Lmap <= cyto_L(Lcenter)+Lcwindow) & cyto_bw;
    cyto_mean(Lcenter) = sum(sum(cyto_map.*max_image(:,:,signal_channel)))/sum(sum(cyto_map));
end
cytoplasmic_profile = [cyto_L',cyto_mean'];

for Lcenter = 1:length(nu_L)
    nu_map = (nucleus_profile(:,1) >= nu_L(Lcenter)-Lnwindow) & (nucleus_profile(:,1) <= nu_L(Lcenter)+Lnwindow);
    nu_mean(Lcenter) = mean(nucleus_profile(nu_map,2));
    nu_std(Lcenter) = std0(nucleus_profile(nu_map,2));
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
%beta1 = nlinfit(nu_x,nucleus_profile(:,2),@exp_C,[(max(nucleus_profile(:,2))-min(nucleus_profile(:,2))),b2,min(nucleus_profile(:,2))]);
%protein_min = beta1(3);
%protein_max = beta1(1)+beta1(3);
protein_min = min(nu_mean);
protein_max = max(nu_mean);

%beta2 = nlinfit(standard_x,[bcd_dxy(:,2);bcd_vxy(:,2)],@exp_C,[(max([bcd_dxy(:,2);bcd_vxy(:,2)])-min([bcd_dxy(:,2);bcd_vxy(:,2)])),b2,min([bcd_dxy(:,2);bcd_vxy(:,2)])]);
%standard_min = beta2(3);
%standard_max = beta2(1)+beta2(3);
standard_min = min(min(bcd_dxy(:,2)),min(bcd_vxy(:,2)));
standard_max = max(max(bcd_dxy(:,2)),min(bcd_vxy(:,2)));

bcd_dxy(:,2) = (bcd_dxy(:,2)-standard_min)*(protein_max-protein_min)/(standard_max-standard_min)+protein_min;
bcd_vxy(:,2) = (bcd_vxy(:,2)-standard_min)*(protein_max-protein_min)/(standard_max-standard_min)+protein_min;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Image output: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(n1)
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
    title([typetxt,' concentration profile: ',image_folder,', cycle = ',num2str(N_cycle),' (',subtp,')'],'Interpreter','none')
    xlabel('AP axis (normalized)')
    ylabel('Concentration (A.U.)')
    legend('Nucleus concentration','Cytoplasmic concentration','Standard dorsal profile','Standard ventral profile','Averaged nuclear profile')

figure(n2)
    clf
    nulabel = bwlabel(seg_bw)+1;
    if max(nulabel(:)) > 20
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
        title(['Nucleus ',typetxt,' concentration plot: ',image_folder,', cycle = ',num2str(N_cycle),' (',subtp,')'],'Interpreter','none')
    end
    

%%% Noise output: %%%======================================================
if size(nucleus_profile,1) >= 2
    p = polyfit(nucleus_profile(:,2),fluc_profile(:,2),1);
else
    p = [0,0];
end
figure(62)
    clf
    plot(nucleus_profile(:,2),fluc_profile(:,2),'b.')
    hold on
    x0 = xlim;
    plot(x0,p(1)*x0+p(2),'r')
    hold off
    xlabel('Mean Intensity');
    ylabel('Intensity varience');
    Cmax = (max(nucleus_profile(:,2))-min(nucleus_profile(:,2)))./p(1)./resolution./resolution./6.02e8./4./pi./sigma.^2;
    title(['Protein Intensity noise plot (',image_folder,'): Var = ',num2str(p(1)),' * I + ',num2str(p(2)),', Cmax = ',num2str(Cmax),' M'],'Interpreter','none')
    legend('Data points','Linear fit')
    quanti_p = [p(1).*resolution.*resolution.*6.02e8.*4.*pi.*sigma.^2,p(2)/p(1)];
%%% =======================================================================

%%% Binomial plot: %%% ====================================================
[FF,IX] = sort(FF);
dF = dF(IX);
bin_size = 20;
if length(FF) >= bin_size
    xbin = zeros(1,floor(length(FF)/bin_size)-2);
    ybin = zeros(1,floor(length(FF)/bin_size)-2);
    sbin = zeros(1,floor(length(FF)/bin_size)-2);
    for Ibin = 0:(floor(length(FF)/bin_size)-3)
        ybin(Ibin+1) = mean(dF((Ibin*bin_size+1):((Ibin+1)*bin_size)));
        xbin(Ibin+1) = mean(FF((Ibin*bin_size+1):((Ibin+1)*bin_size)));
        sbin(Ibin+1) = std0(dF((Ibin*bin_size+1):((Ibin+1)*bin_size)));
    end

    pF = polyfit(xbin,ybin,1);
    %pF(1) = (mean(xbin.*ybin)-mean(xbin).*mean(ybin))/(std(xbin,1)^2);
    %pF(2) = 0;
else
    xbin = FF;
    ybin = dF;
    sbin = zeros(size(FF));
    pF = [0,0];
end

figure(63)
    clf
    plot(FF,dF,'b.')
    hold on
    x0 = xlim;
    plot(x0,pF(1)*x0+pF(2),'r')
    
    errorbar(xbin,ybin,sbin,'g.')
    
    hold off
    xlabel('Total cell intensity (FF)');
    ylabel('Binomial partition variation (dF^2)');
    Cmax = max((FF-pF(2)/pF(1))./SF)./7.26./pF(1)./resolution./resolution./6.02e8./sigma.^2;
    title(['Protein binomial partition noise plot (',image_folder,'): dF^2 = ',num2str(pF(1)),' * FF + ',num2str(pF(2)),', Cmax = ',num2str(Cmax),' M'],'Interpreter','none')
    legend('Binomial partition noise','Linear fit of binomial partition noise')
    

%%% =======================================================================


    nucleus_profile = [nucleus_profile,fluc_profile(:,2)];
    
    
    function y = std0(x)
        y = std(x)/sqrt(length(x));
        