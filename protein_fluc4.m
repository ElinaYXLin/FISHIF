function [nucleus_profile,fluc_profile] = protein_fluc4(seg_bw,cyto_bw,max_image,signal_channel,mask_stack,signal_stack,image_folder,resolution,varargin)

%% Nucleus/cytoplasmic protein fluctuation/concentration calculation: %%%%%
%%% seg_bw: 2D nucleus mask;
%%% cyto_bw: 2D cytoplasm mask (not used in this program);
%%% max_image: maximum projection image (not used in this program);
%%% signal_channel: protein channel # (not used in this program);
%%% mask_stack: quasi-3D nucleus mask;
%%% signal_stack: 3D protein signal stack;
%%% image_folder: input image folder (name of data);
%%% resolution: xy resolution of the image;
%%% varargin: for result output
%%%            h1: alternative axes handle for multinormial fluctuation plot (variance vs mean).
%%%            h2: alternative axes handle for binormial fluctuation plot (variance vs mean).
%%%            t1: alternative text handle to put title for h1.
%%%            t2: alternative text handle to put title for h2.
%%%            h3: alternative axes handle for multinormial sub-partition sliding window plot (slope vs window size).
%%%            c3: color code for h3.
%%%            h4: alternative axes handle for multinormial sub-partition sliding window plot (variance vs mean for different window sizes).
%%%            c4: color code for h4.
%%%            window_show: the window size range for h4.
%%%            h5: alternative axes handle for correlation plot.

%%% nucleus_profile: protein intensity list for each nuclei.
%%% fluc_profile: protein fluctuation list for each nuclei.
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

%max_image = double(max_image);
start_limit = 0.3;
end_limit = 1-start_limit;
resolution0 = 0.097;
sigma = 1.16*resolution0./resolution;
sigma_x = 0.985;
sigma_y = 1.37;
%sigma = 1.23;
cz = sqrt(2);
conf95 = erfinv(0.95)*sqrt(2);
Gb = Con_calculator(10000,10000,sigma_x,sigma_y);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Nucleus region recognition/processing: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Multinormial partition analysis: %%% ==================================
temp = zeros(0);
Inten = zeros(0);
Inten2 = zeros(0);
mask2 = bwlabel(imclearborder(seg_bw));
for I_layer = 1:size(mask_stack,3)
    sg_prop = regionprops((imerode(mask_stack(:,:,I_layer), strel('disk',5))).*mask2,signal_stack(:,:,I_layer),'MeanIntensity','PixelValues');
    sg2_prop = regionprops((imerode(mask_stack(:,:,I_layer), strel('disk',5))).*mask2,signal_stack(:,:,I_layer).^2,'MeanIntensity');
    if I_layer == 1
        temp0 = cell(size(mask_stack,3),max(max(mask2)));
        temp = zeros(size(mask_stack,3),max(max(mask2)));
        temp2 = zeros(size(mask_stack,3),max(max(mask2)));
    end
    if ~isempty(sg_prop)
        temp0(I_layer,1:length(sg2_prop)) = {sg_prop.PixelValues}; 
        temp(I_layer,1:length(sg2_prop)) = [sg_prop.MeanIntensity];
        temp2(I_layer,1:length(sg2_prop)) = [sg2_prop.MeanIntensity]-[sg_prop.MeanIntensity].^2;
    end
end
sg_prop = regionprops(mask2,'Centroid');
centr = [sg_prop.Centroid];
temp(isnan(temp)) = 0;
[Inten,max_index] = max(temp,[],1);
Inten2 = temp2(max_index+size(temp,1)*[0:(size(temp,2)-1)]);

refine_signal = zeros(size(seg_bw));
refine_mask = imerode(seg_bw, strel('disk',5));
for I_layer = 1:size(mask_stack,3)
    refine_signal = refine_signal+ismember(mask2,find(max_index == I_layer)).*signal_stack(:,:,I_layer);
end

refine_signal_cyto = signal_stack(:,:,ceil(median(max_index)));
refine_temp = imerode((~seg_bw), strel('disk',5));
refine_mask_cyto = false(size(seg_bw));
refine_mask_cyto(2:(end-1),2:(end-1)) = refine_temp(2:(end-1),2:(end-1));



%%%%% Correlation calculation: 
refine_mask = imclearborder(refine_mask);
r_prop = regionprops(refine_mask,refine_signal,'Image','PixelValues','PixelList');
Crx = 5;
Cry = 5;
mean_corr = zeros(2*Crx+1,2*Cry+1);
%refine2 = bwlabel(refine_mask);

for I_region = 1:length(r_prop)
    region1 = r_prop(I_region).Image;
    im_size = size(r_prop(I_region).Image);
    area1 = zeros(im_size);
    area1(region1) = r_prop(I_region).PixelValues;
    area2 = (area1 - sum(sum(area1))./sum(sum(region1))).*region1;
    corrm = xcorr2(area2)./xcorr2(double(region1));
    try
        corrm2 = corrm((im_size(1)-Crx):(im_size(1)+Crx),(im_size(2)-Cry):(im_size(2)+Cry));
        corrm2 = corrm2./corrm2(Crx+1,Cry+1);
        mean_corr = mean_corr+corrm2;
    catch
    end
end
[XX,YY] = meshgrid(-Cry:Cry,-Crx:Crx);
mean_corr = mean_corr./length(r_prop);
%%% =======================================================================


%%% Binomial partition analysis: %%% ======================================
dF = zeros(length(sg_prop),1);
FF = zeros(length(sg_prop),1);
SF = zeros(length(sg_prop),1);
for I_region = 1:length(sg_prop)
    Inten_value = temp0{max_index+size(temp,1)*(I_region-1)};
    size_region = size(Inten_value,1);
    if floor(size_region/2) < ceil(size_region/2)
        size_region = size_region-1;
    end
    dF(I_region) = ((sum(Inten_value(1:(size_region/2)))-sum(Inten_value((size_region/2+1):size_region)))/2).^2;
    FF(I_region) = sum(Inten_value(1:size_region));
    SF(I_region) = size_region;
end
%%% =======================================================================


%%% Multinormial sub-partition analysis: %%% ==============================
seg = 3:1:61;
k_seg = zeros(size(seg));
kse_seg = zeros(size(seg));
k_seg_cyto = zeros(size(seg));
kse_seg_cyto = zeros(size(seg));
k_cy1 = zeros(size(seg));
kse_cy1 = zeros(size(seg));

h4 = zeros(0);
window_show = zeros(0);

if length(varargin) > 6
    h4 = varargin{7};
    mark_use = varargin{8};
    if mark_use == 'r'
        mark_use_cyto = 'm';
    elseif mark_use == 'b'
        mark_use_cyto = 'c';
    else
        mark_use_cyto = 'g';
    end
    window_show = varargin{9};
else
    figure(65)
    clf
    window_show = 11:10:41;
    mh = length(window_show);
    h4 = zeros(1,mh);
    for i_h = 1:mh
        h4(i_h) = subplot(2,ceil(mh./2),i_h);
    end

    if image_folder((end-1)) == 'A'
        mark_use = 'r';
        mark_use_cyto = 'm';
    elseif image_folder((end-1)) == 'B'
        mark_use = 'b';
        mark_use_cyto = 'c';
    else
        mark_use = 'k';
        mark_use_cyto = 'g';
    end
end


for I_seg = 1:length(seg)
    sub_div = strel('disk', seg(I_seg));
    sub_center = 1;
    %sub_center = getnhood(strel('disk', 2));
    sub_hood = getnhood(sub_div);
    xy_div = size(sub_hood);
    xy_center = size(sub_center);
    
    sub_hood(((xy_div(1)-xy_center(1))/2+1):(xy_div(1)+xy_center(1))/2,((xy_div(2)-xy_center(2))/2+1):(xy_div(2)+xy_center(2))/2) = sub_hood(((xy_div(1)-xy_center(1))/2+1):(xy_div(1)+xy_center(1))/2,((xy_div(2)-xy_center(2))/2+1):(xy_div(2)+xy_center(2))/2)-sub_center;
    sub_hood = sub_hood./sum(sum(sub_hood));
    
    seg_mask = imerode(imclearborder(refine_mask), imfill(logical(sub_hood),'holes'));
    mean_seg = conv2(refine_signal,sub_hood,'same');
    var_seg = (refine_signal-mean_seg).^2;
    %var_seg = conv2(refine_signal.^2,sub_hood,'same')-mean_seg.^2;
    
    seg_mask_cyto = imerode(refine_mask_cyto, imfill(logical(sub_hood),'holes'));
    mean_seg_cyto = conv2(refine_signal_cyto,sub_hood,'same');
    %var_seg = (refine_signal-mean_seg).^2;
    var_seg_cyto = conv2(refine_signal_cyto.^2,sub_hood,'same')-mean_seg_cyto.^2;
    
    
    
    if length(find(seg_mask)) > 10
        p = polyfit(mean_seg(seg_mask),var_seg(seg_mask),1);
        %k_seg(I_seg) = p(1)./(1-Con_calculator(seg(I_seg),seg(I_seg),sigma_x,sigma_y,sub_hood);
        k_seg(I_seg) = p(1);
        cov_xy = cov(mean_seg(seg_mask),var_seg(seg_mask));
        kse_seg(I_seg) = conf95*sqrt((var(mean_seg(seg_mask)).*var(var_seg(seg_mask))-cov_xy(1,2).^2)./(length(mean_seg(seg_mask))-2))./var(mean_seg(seg_mask));
        %kse_seg(I_seg) = sqrt((var(mean_seg(seg_mask)).*var(var_seg(seg_mask))-cov_xy(1,2).^2)./50)./var(mean_seg(seg_mask));
    else
        k_seg(I_seg) = NaN;
        kse_seg(I_seg) = NaN;
    end
    
    if length(find(seg_mask_cyto)) > 10
        pp = polyfit(mean_seg_cyto(seg_mask_cyto),var_seg_cyto(seg_mask_cyto),1);
        k_seg_cyto(I_seg) = pp(1);
        cov_xy_cyto = cov(mean_seg_cyto(seg_mask_cyto),var_seg_cyto(seg_mask_cyto));
        kse_seg_cyto(I_seg) = conf95*sqrt((var(mean_seg_cyto(seg_mask_cyto)).*var(var_seg_cyto(seg_mask_cyto))-cov_xy_cyto(1,2).^2)./(length(mean_seg_cyto(seg_mask_cyto))-2))./var(mean_seg_cyto(seg_mask_cyto));
        %kse_seg_cyto(I_seg) = sqrt((var(mean_seg_cyto(seg_mask_cyto)).*var(var_seg_cyto(seg_mask_cyto))-cov_xy_cyto(1,2).^2)./50)./var(mean_seg_cyto(seg_mask_cyto));
    else
        k_seg_cyto(I_seg) = NaN;
        kse_seg_cyto(I_seg) = NaN;
    end
    
    
    mm_seg_cyto = mean_seg_cyto(seg_mask_cyto);
    dd_seg_cyto = var_seg_cyto(seg_mask_cyto);
    x00 = [min(mm_seg_cyto),max(mm_seg_cyto)];
    if ~isempty(x00)
        Nbin = 40;
        wbin =0.05*(x00(2)-x00(1));
        xbin = x00(1)+[0:1/(Nbin-1):1]*(x00(2)-x00(1));
        ybin = zeros(1,Nbin);
        sbin = zeros(1,Nbin);
        for Ibin = 1:Nbin
            ybin(Ibin) = mean(dd_seg_cyto((mm_seg_cyto >= (xbin(Ibin)-wbin)) & (mm_seg_cyto <= (xbin(Ibin)+wbin))));
            sbin(Ibin) = std0(dd_seg_cyto((mm_seg_cyto >= (xbin(Ibin)-wbin)) & (mm_seg_cyto <= (xbin(Ibin)+wbin))));
        end
    end
    if all(~isnan(ybin(1:6)))
        pcy1 = polyfit(xbin(1:6),ybin(1:6),1);
        k_cy1(I_seg) = pcy1(1);
        cov_cy1 = cov(xbin(1:6),ybin(1:6));
        kse_cy1(I_seg) = conf95*sqrt((var(xbin(1:6)).*var(ybin(1:6))-cov_cy1(1,2).^2)./(length(mean_seg(seg_mask))-2))./var(xbin(1:6));
    else
        k_cy1(I_seg) = NaN;
        kse_cy1(I_seg) = NaN;
    end
    
    %%%%% Sample image show
    if any(ismember(window_show,size(sub_hood,1)))
        axes(h4(ismember(window_show,size(sub_hood,1))));
        hold on
        title(['Subpartition estimation',10,image_folder((strfind(image_folder,'Desktop')+8):(end-3)),10,'(red: anterior, blue: posterior, magenta: anterior cyto, cyan: posterior cyto, window size = ',num2str(window_show(ismember(window_show,size(sub_hood,1)))),')'])
        plot(mean_seg(seg_mask),var_seg(seg_mask),['.',mark_use],mean_seg_cyto(seg_mask_cyto),var_seg_cyto(seg_mask_cyto),['.',mark_use_cyto])
        xlabel('Mean intensity');
        ylabel('Intensity variance');
        
        mm_seg = mean_seg(seg_mask);
        dd_seg = var_seg(seg_mask);
        x00 = [min(mm_seg),max(mm_seg)];
        if ~isempty(x00)
            Nbin = 40;
            wbin =0.05*(x00(2)-x00(1));
            xbin = x00(1)+[0:1/(Nbin-1):1]*(x00(2)-x00(1));
            ybin = zeros(1,Nbin);
            sbin = zeros(1,Nbin);
            for Ibin = 1:Nbin
                ybin(Ibin) = mean(dd_seg((mm_seg >= (xbin(Ibin)-wbin)) & (mm_seg <= (xbin(Ibin)+wbin))));
                sbin(Ibin) = std0(dd_seg((mm_seg >= (xbin(Ibin)-wbin)) & (mm_seg <= (xbin(Ibin)+wbin))));
            end
            errorbar(xbin,ybin,sbin,'o','MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63],'MarkerSize',10);
            if mark_use == 'r'
                mark_use2 = 'r';
            elseif mark_use == 'b'
                mark_use2 = 'b';
            else
                mark_use2 = 'k';

            end
            plot(xlim,p(1)*xlim+p(2),mark_use2, 'LineWidth',2);
        end
        
        
        mm_seg_cyto = mean_seg_cyto(seg_mask_cyto);
        dd_seg_cyto = var_seg_cyto(seg_mask_cyto);
        x00 = [min(mm_seg_cyto),max(mm_seg_cyto)];
        if ~isempty(x00)
            Nbin = 40;
            wbin =0.05*(x00(2)-x00(1));
            xbin = x00(1)+[0:1/(Nbin-1):1]*(x00(2)-x00(1));
            ybin = zeros(1,Nbin);
            sbin = zeros(1,Nbin);
            for Ibin = 1:Nbin
                ybin(Ibin) = mean(dd_seg_cyto((mm_seg_cyto >= (xbin(Ibin)-wbin)) & (mm_seg_cyto <= (xbin(Ibin)+wbin))));
                sbin(Ibin) = std0(dd_seg_cyto((mm_seg_cyto >= (xbin(Ibin)-wbin)) & (mm_seg_cyto <= (xbin(Ibin)+wbin))));
            end
            errorbar(xbin,ybin,sbin,'o','MarkerEdgeColor','k','MarkerFaceColor',[.63 1 .49],'MarkerSize',10);
            if mark_use_cyto == 'm'
                mark_use2_cyto = 'm';
            elseif mark_use_cyto == 'c'
                mark_use2_cyto = 'c';
            else
                mark_use2_cyto = 'k';

            end
            plot(xlim,pp(1)*xlim+pp(2),mark_use2_cyto, 'LineWidth',2);
        end
    end
    
end
%%% =======================================================================


%%% Multinormial sub-clustering analysis: %%% =============================
refine_mask = imclearborder(refine_mask);
r_prop = regionprops(refine_mask,refine_signal,'Image','PixelValues','PixelList','PixelIdxList');
rr = 2;
n_cluster = 2:20;
k_cluster = zeros(size(n_cluster));
kse_cluster = zeros(size(n_cluster));
a_cluster = zeros(size(n_cluster));

h6 = zeros(0);
window_show = zeros(0);

if length(varargin) > 10
    h6 = varargin{11};
    mark_use = varargin{12};
    k_show = varargin{13};
else
    figure(67)
    clf
    k_show = 2:5:20;
    mh = length(k_show);
    h6 = zeros(1,mh);
    for i_h = 1:mh
        h6(i_h) = subplot(2,ceil(mh./2),i_h);
    end
    if image_folder((end-1)) == 'A'
        mark_use = 'r';
    elseif image_folder((end-1)) == 'B'
        mark_use = 'b';
    else
        mark_use = 'k';
    end
end

for I_cluster = 1:length(n_cluster)
    %%%%% Clustering mask generation:
    new_mask = zeros(size(refine_mask));
    for I_region = 1:length(r_prop)
        ccxy = [r_prop(I_region).PixelList,r_prop(I_region).PixelValues];
        mmd = (max(ccxy(:,3))-min(ccxy(:,3)));
        mmd = mmd+(mmd == 0);
        ccxy_norm = [(ccxy(:,1)-min(ccxy(:,1)))./(max(ccxy(:,1))-min(ccxy(:,1))),(ccxy(:,2)-min(ccxy(:,2)))./(max(ccxy(:,2))-min(ccxy(:,2))),rr*(ccxy(:,3)-min(ccxy(:,3)))./mmd];
        IDX = kmeans(ccxy_norm,n_cluster(I_cluster));
        new_mask([r_prop(I_region).PixelIdxList]) = IDX+n_cluster(I_cluster)*(I_region-1);
    end
    %%%%% Mean/variance calculation:
    ccg_prop = regionprops(new_mask,refine_signal,'MeanIntensity','Area');
    ccg2_prop = regionprops(new_mask,refine_signal.^2,'MeanIntensity');
    ccmean = [ccg_prop.MeanIntensity];
    ccvar = [ccg2_prop.MeanIntensity]-[ccg_prop.MeanIntensity].^2;
    ccg = polyfit(ccmean,ccvar,1);
    k_cluster(I_cluster) = ccg(1);
    cov_cxy = cov(ccmean,ccvar);
    kse_cluster(I_cluster) = conf95*sqrt((var(ccmean).*var(ccvar)-cov_cxy(1,2).^2)./(length(ccmean)-2))./var(ccmean);
    a_cluster(I_cluster) = mean([ccg_prop.Area]);
    %%%%% Sample show:
    if any(ismember(k_show,n_cluster(I_cluster)))
        axes(h6(ismember(k_show,n_cluster(I_cluster))));
        hold on
        title(['Subclustering estimation',10,image_folder((strfind(image_folder,'Desktop')+8):(end-3)),10,'(red: anterior, blue: posterior, magenta: anterior cyto, cyan: posterior cyto, cluster # = ',num2str(n_cluster(I_cluster)),', mean area = ',num2str(a_cluster(I_cluster)),')'])
        plot(ccmean,ccvar,['.',mark_use])
        xlabel('Mean intensity');
        ylabel('Intensity variance');
        
        x00 = [min(ccmean),max(ccmean)];
        if ~isempty(x00)
            Nbin = 20;
            wbin =0.05*(x00(2)-x00(1));
            xbin = x00(1)+[0:1/(Nbin-1):1]*(x00(2)-x00(1));
            ybin = zeros(1,Nbin);
            sbin = zeros(1,Nbin);
            for Ibin = 1:Nbin
                ybin(Ibin) = mean(ccvar((ccmean >= (xbin(Ibin)-wbin)) & (ccmean <= (xbin(Ibin)+wbin))));
                sbin(Ibin) = std0(ccvar((ccmean >= (xbin(Ibin)-wbin)) & (ccmean <= (xbin(Ibin)+wbin))));
            end
            errorbar(xbin,ybin,sbin,'o','MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63],'MarkerSize',10);
            if mark_use == 'r'
                mark_use2 = 'r';
            elseif mark_use == 'b'
                mark_use2 = 'b';
            else
                mark_use2 = 'k';

            end
            plot(xlim,ccg(1)*xlim+ccg(2),mark_use2, 'LineWidth',2);
        end
    end
end

%%% =======================================================================



%%% Down-sampling method: %%% =============================================
% refine_mask = imclearborder(refine_mask);
% mask20 = bwlabel(refine_mask);
% size_down = 1:10;   %%%%% Down sampling size
% num_th = 10;
% k_down = zeros(size(size_down));
% for I_down = 1:length(size_down)
%     mean_down = zeros(0);
%     var_down = zeros(0);
%     for I_start = 1:size_down(I_down)
%         for J_start = 1:size_down(I_down)
%             temp_mask = zeros(size(refine_mask));
%             temp_mask(I_start:size_down(I_down):end,J_start:size_down(I_down):end) = refine_mask(I_start:size_down(I_down):end,J_start:size_down(I_down):end);
%             
%             down1_prop = regionprops(mask20,refine_signal.*temp_mask,'MeanIntensity');
%             down2_prop = regionprops(mask20,refine_signal.^2.*temp_mask,'MeanIntensity');
%             down0_prop = regionprops(mask20,temp_mask,'MeanIntensity','Area');
%             
%             temp_mean = [down1_prop.MeanIntensity]./[down0_prop.MeanIntensity];
%             temp_var = [down2_prop.MeanIntensity]./[down0_prop.MeanIntensity]-temp_mean.^2;
%             temp_num = [down0_prop.MeanIntensity].*[down0_prop.Area];
%             temp_var = temp_var./(1-1./temp_num);
%             
%             mean_down = cat(2,mean_down,temp_mean(logical((~isnan(temp_mean)).*(temp_num >= num_th))));
%             var_down = cat(2,var_down,temp_var(logical((~isnan(temp_mean)).*(temp_num >= num_th))));
%         end
%     end
%     
%     if length(mean_down) > 10
%         p = polyfit(mean_down,var_down,1);
%         k_down(I_down) = p(1);
%     else
%         k_down(I_down) = NaN;
%     end
% end
%%% =======================================================================


%%% Normalized coordinate calculation: %%%=================================
if ~isempty(centr)
    xy = centr;
    nucleus_xy = [xy(1:2:length(xy)-1);xy(2:2:length(xy))]';
    
    distXY = pdist2(nucleus_xy,nucleus_xy);
    axis_length = max(max(distXY));
    [extreme0,extreme1] = find(distXY >= axis_length,1); %%% find the anterior and posterior pole nucleus
    L2_extreme = dot((nucleus_xy(extreme1,:)-nucleus_xy(extreme0,:)),(nucleus_xy(extreme1,:)-nucleus_xy(extreme0,:)),2);
    nucleus_distance = dot((nucleus_xy-repmat(nucleus_xy(extreme0,:),size(nucleus_xy,1),1)),repmat((nucleus_xy(extreme1,:)-nucleus_xy(extreme0,:)),size(nucleus_xy,1),1),2)/L2_extreme;
    
    nucleus_profile = [nucleus_distance,Inten'];   %%% mean intensity profile
    fluc_profile = [nucleus_distance,Inten2'];   %%% intensity variance profile
end
%%% =======================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%% Noise fit: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = polyfit(nucleus_profile(:,2),fluc_profile(:,2),1);
pF = polyfit(FF,dF,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






%% Image output: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(varargin)
    figure(62)
    clf
else
    axes(varargin{1});
    cla
end
    plot(nucleus_profile(:,2),fluc_profile(:,2),'.')
    hold on
    x0 = xlim;
    plot(x0,p(1)*x0+p(2))
    hold off
    xlabel('Mean Intensity','FontName','Arial','FontSize',12,'FontWeight','bold');
    ylabel('Intensity varience','FontName','Arial','FontSize',12,'FontWeight','bold');
    Cmax = (max(nucleus_profile(:,2))-p(2)./p(1))./p(1)./resolution./resolution./6.02e8./4./pi./sigma.^2;
    if length(varargin) <= 2 || isempty(varargin{3})
        title(['Protein Intensity noise plot (',image_folder,'): Var = ',num2str(p(1)),' * I + ',num2str(p(2)),', Cmax = ',num2str(Cmax),' M, I0 = ',num2str(p(1)*4*pi*sigma.^2.*cz)],'FontName','Arial','FontSize',14,'FontWeight','bold')
    else
        set(varargin{3},'String',['Multinormial plot: Var = ',num2str(p(1)),' * I + ',num2str(p(2)),', Cmax = ',num2str(Cmax),' M, I0 = ',num2str(p(1)*4*pi*sigma.^2.*cz)])
    end
    legend('Data points','Linear fit')
    
    
if isempty(varargin)
    figure(63)
    clf
else
    axes(varargin{2});
    cla
end
    plot(FF,dF,'.')
    hold on
    x0 = xlim;
    plot(x0,pF(1)*x0+pF(2))
    
    x00 = [min(FF),max(FF)];
    Nbin = 5;
    wbin =0.05*(x00(2)-x00(1));
    xbin = x00(1)+[0:1/(Nbin-1):1]*(x00(2)-x00(1));
    ybin = zeros(1,Nbin);
    sbin = zeros(1,Nbin);
    for Ibin = 1:Nbin
        ybin(Ibin) = mean(dF((FF >= (xbin(Ibin)-wbin)) & (FF <= (xbin(Ibin)+wbin))));
        sbin(Ibin) = std0(dF((FF >= (xbin(Ibin)-wbin)) & (FF <= (xbin(Ibin)+wbin))));
    end
    errorbar(xbin,ybin,sbin,'r.')
    
    hold off
    xlabel('Total cell intensity (FF)','FontName','Arial','FontSize',12,'FontWeight','bold');
    ylabel('Binomial partition variation (dF^2)','FontName','Arial','FontSize',12,'FontWeight','bold');
    Cmax = max((FF-pF(2)/pF(1))./SF)./4./pF(1)./resolution./resolution./6.02e8;
    if length(varargin) <= 2 || isempty(varargin{4})
        title(['Protein binomial partition noise plot (',image_folder,'): dF^2 = ',num2str(pF(1)),' * FF + ',num2str(pF(2)),', Cmax = ',num2str(Cmax),' M, I0 = ',num2str(pF(1)*8*pi*sigma.^2./Gb.*cz)],'FontName','Arial','FontSize',14,'FontWeight','bold')
    else
        set(varargin{4},'String',['Binormial plot: dF^2 = ',num2str(pF(1)),' * FF + ',num2str(pF(2)),', Cmax = ',num2str(Cmax),' M, I0 = ',num2str(pF(1)*8*pi*sigma.^2./Gb.*cz)])
    end
    legend('Binomial partition noise','Linear fit of binomial partition noise')

if isempty(varargin) || length(varargin) <= 4
    figure(64)
    clf
    hold on
    title(['Particle intensity estimation vs. segmentation window size',image_folder])
    errorbar(seg(~isnan(k_seg)),k_seg(~isnan(k_seg))*4*pi*sigma.^2.*cz,kse_seg(~isnan(k_seg))*4*pi*sigma.^2.*cz,'o')
    errorbar(seg(~isnan(k_seg_cyto)),k_seg_cyto(~isnan(k_seg_cyto))*4*pi*sigma.^2.*cz,kse_seg_cyto(~isnan(k_seg_cyto))*4*pi*sigma.^2.*cz,'^')
    plot(xlim,p(1)*4*pi*sigma.^2.*cz*[1,1],'--')
else
    axes(varargin{5});
    mark_use = varargin{6};
    if mark_use == 'r'
        mark_use_cyto = 'm';
    elseif mark_use == 'b'
        mark_use_cyto = 'c';
    else
        mark_use_cyto = 'g';
    end
    hold on
    title(['Particle intensity estimation vs. segmentation window size',10,image_folder(1:(end-3)),'(red: anterior, blue: posterior, magenta: anterior cyto, cyan: posterior cyto)'])
    errorbar(seg(~isnan(k_seg)),k_seg(~isnan(k_seg))*4*pi*sigma.^2.*cz,kse_seg(~isnan(k_seg))*4*pi*sigma.^2.*cz,['o',mark_use])
    errorbar(seg(~isnan(k_seg_cyto)),k_seg_cyto(~isnan(k_seg_cyto))*4*pi*sigma.^2.*cz,kse_seg_cyto(~isnan(k_seg_cyto))*4*pi*sigma.^2.*cz,['o',mark_use_cyto])
    errorbar(seg(~isnan(k_cy1)),k_cy1(~isnan(k_cy1))*4*pi*sigma.^2.*cz,kse_cy1(~isnan(k_cy1))*4*pi*sigma.^2.*cz,['^',mark_use_cyto])
    plot([1,seg(end)],p(1)*4*pi*sigma.^2.*cz*[1,1],['--',mark_use])
end
    ylim0 = ylim;
    ylim([0,ylim0(2)]);
    xlabel('Segmentation window size');
    ylabel('Predicted Single particle intensity');
    
    
% if isempty(varargin) || length(varargin) <= 9
%     figure(66)
%     clf
%     title(['Particle intensity estimation vs. down sampling size',image_folder])
%     plot(size_down(~isnan(k_down)),k_down(~isnan(k_down))*4*pi*sigma.^2.*cz,'o',xlim,p(1)*4*pi*sigma.^2.*cz*[1,1],'--')
% else
%     axes(varargin{10});
%     hold on
%     title(['Particle intensity estimation vs. down sampling size',10,image_folder(1:(end-3)),'(red: anterior, blue: posterior)'])
%     plot(size_down(~isnan(k_down)),k_down(~isnan(k_down))*4*pi*sigma.^2.*cz,['o',varargin{11}],[1,size_down(end)],p(1)*4*pi*sigma.^2.*cz*[1,1],['--',varargin{11}])
% end
%     ylim0 = ylim;
%     ylim([0,ylim0(2)]);
%     xlabel('Down sampling size');
%     ylabel('Predicted Single particle intensity');
    
if isempty(varargin) || length(varargin) <= 9
    figure(66)
    clf
    surf(XX,YY,mean_corr);
    title(['Mean nuclei signal correlation',image_folder])
else
    axes(varargin{10});
    %hold on
    surf(XX,YY,mean_corr);
    title(['Mean nuclei signal correlation',image_folder])
end


if isempty(varargin) || length(varargin) <= 13
    figure(68)
    clf
    hold on
    title(['Particle intensity estimation vs. sub-clustering area size: ',image_folder])
    errorbar(a_cluster,k_cluster*4*pi*sigma.^2.*cz,kse_cluster*4*pi*sigma.^2.*cz,'o')
    plot(xlim,p(1)*4*pi*sigma.^2.*cz*[1,1],'--')
else
    axes(varargin{14});
    mark_use = varargin{12};
    hold on
    title(['Particle intensity estimation vs. sub-clustering area size:',10,image_folder(1:(end-3)),'(red: anterior, blue: posterior)'])
    errorbar(a_cluster,k_cluster*4*pi*sigma.^2.*cz,kse_cluster*4*pi*sigma.^2.*cz,['o',mark_use])
    plot(xlim,p(1)*4*pi*sigma.^2.*cz*[1,1],['--',mark_use])
end
    ylim0 = ylim;
    ylim([0,ylim0(2)]);
    xlabel('sub-clustering area size');
    ylabel('Predicted Single particle intensity');
    
    
    
    
    
    
    
function y = std0(x)
    NN = size(x);
    if NN(1) == 1;
        N1 = NN(2);
    else
        N1 = NN(1);
    end
    y = std(x)/sqrt(N1);
