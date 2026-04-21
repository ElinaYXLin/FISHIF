function [nucleus_profile,cytoplasmic_profile,quanti_p,varargout] = protein_profile3D2(nucleus_DAPI_profile,max_image,EL_info,signal_channel,mask_stack,signal_stack,image_folder,N_cycle,resolution,resolutionz,flip_EL,varargin)


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
resolution0 = 0.0915;
L_ratio = (resolution/resolution0);
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
c_radius = 100;

Lmin = 0;
Lnbin = 0.05;
Lmax = 1;
Lnwindow = 0.05;
nu_L = Lmin:Lnbin:Lmax;
nu_mean = zeros(size(nu_L));
nu_std = zeros(size(nu_L));

% max_image = double(max_image);
start_limit = 0.3;
end_limit = 1-start_limit;

typetxt = 'Protein';
n1 = 2;
n2 = 14;
n3 = 62;
n4 = 63;
sigma = 1.35;
sigmaz = 0.8;
cz0 = 7;

if ~isempty(varargin)
    if ~isempty(varargin{1})
        typetxt = varargin{1};
    end
    if ~isempty(varargin{2})
        n1 = varargin{2};
    end
    if ~isempty(varargin{3})
        n2 = varargin{3};
    end
    if ~isempty(varargin{4})
        n3 = varargin{4};
    end
    if ~isempty(varargin{5})
        n4 = varargin{5};
    end
    if length(varargin) > 5 && ~isempty(varargin{6})
        sigma = varargin{6}(1);
        sigmaz = varargin{6}(2);
    end
    if length(varargin) > 6
        flip_axis = varargin{7};
    else
        flip_axis = false;
    end
    if length(varargin) > 7 && ~isempty(varargin{8})
        EL_range = varargin{8};
        if length(EL_range) < 3
            EL_range2 = [0.8,1,0.7,0.8];
            if flip_axis
                EL_range2 = 1-EL_range2;
            end
            EL_range = [EL_range,EL_range2];
        elseif length(EL_range) < 5
            EL_range2 = [0.7,0.8];
            if flip_axis
                EL_range2 = 1-EL_range2;
            end
            EL_range = [EL_range,EL_range2];
        end
    else
        EL_range = [0,0.5,0.8,1,0.7,0.8];
        if flip_axis
            EL_range = 1-EL_range;
        end
    end
    
    if length(varargin) > 8 && ~isempty(varargin{9})
        Inten_protein = varargin{9};
    else
        Inten_protein = [];
    end
    
    if length(varargin) > 9 && ~isempty(varargin{10})
        tbk0 = varargin{10};
    else
        tbk0 = false;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Nucleus region recognition/processing: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
temp = zeros(0);
Inten = zeros(0);
Inten2 = zeros(0);
%mask2 = max(mask_stack,[],3);
seg_bw = max(mask_stack,[],3);
embryo_region = bwconvhull(seg_bw);
cyto_bw = embryo_region & (~seg_bw);
cyto_bw = imerode(cyto_bw, strel('disk',floor(5/L_ratio)));


for I_layer = 1:size(mask_stack,3)
    sg_prop = regionprops(double(imerode(logical(mask_stack(:,:,I_layer)), strel('disk',r_edge))).*mask_stack(:,:,I_layer),double(signal_stack(:,:,I_layer)),'MeanIntensity','PixelValues');
    sg2_prop = regionprops(double(imerode(logical(mask_stack(:,:,I_layer)), strel('disk',r_edge))).*mask_stack(:,:,I_layer),double(signal_stack(:,:,I_layer)).^2,'MeanIntensity');
    if I_layer == 1
        temp0 = cell(size(mask_stack,3),max(mask_stack(:)));
        temp = zeros(size(mask_stack,3),max(mask_stack(:)));
        temp2 = zeros(size(mask_stack,3),max(mask_stack(:)));
    end
    if ~isempty(sg_prop)
        temp0(I_layer,1:length(sg2_prop)) = {sg_prop.PixelValues};  %%%%% Binomial partition preparation
        temp(I_layer,1:length(sg2_prop)) = [sg_prop.MeanIntensity];
        temp2(I_layer,1:length(sg2_prop)) = [sg2_prop.MeanIntensity]-[sg_prop.MeanIntensity].^2;
    end
end
sg_prop = regionprops(mask_stack,'Centroid');
centr = cell2mat({sg_prop.Centroid}');
temp(isnan(temp)) = 0;
% [Inten,max_index] = max(temp,[],1);
max_index = nucleus_DAPI_profile(:,3)';
Inten = temp(max_index+size(temp,1)*[0:(size(temp,2)-1)]);
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

% flip_EL = false;
%%% Normalized coordinate calculation: %%%=================================
if ~isempty(centr)
    xy = centr;
    nucleus_xy = xy(:,1:2);
    x0 = EL_info(1);
    y0 = EL_info(2);
    x1 = EL_info(3);
    y1 = EL_info(4);
    L2_extreme = EL_info(5);
% % %     if size(seg_bw,2)/size(seg_bw,1) > 1.4
% % %         distXY = pdist2(nucleus_xy,nucleus_xy);
% % %         axis_length = max(max(distXY));
% % %         [extreme0,extreme1] = find(distXY >= axis_length,1); %%% find the anterior and posterior pole nucleus
% % %         x0 = nucleus_xy(extreme0,1);
% % %         y0 = nucleus_xy(extreme0,2);
% % %         x1 = nucleus_xy(extreme1,1);
% % %         y1 = nucleus_xy(extreme1,2);
% % %         L2_extreme = (x1-x0)^2+(y1-y0)^2;
% % %     else
% % %         I1 = nucleus_xy(:,1) < size(seg_bw,2)/2;
% % %         I2 = nucleus_xy(:,1) >= size(seg_bw,2)/2;
% % %         if std(nucleus_xy(I1,2)) <= std(nucleus_xy(I2,2))
% % %             [x0,I0] = min(nucleus_xy(:,1));
% % %             y0 = nucleus_xy(I0,2);
% % %             x1 = 2*size(seg_bw,2)-x0;
% % %             y1 = y0;
% % %             axis_length = x1-x0;
% % %             L2_extreme = axis_length^2;
% % %         else
% % %             [x0,I0] = max(nucleus_xy(:,1));
% % %             y0 = nucleus_xy(I0,2);
% % %             x1 = 2*1-x0;
% % %             y1 = y0;
% % %             axis_length = x1-x0;
% % %             L2_extreme = axis_length^2;
% % %         end
% % %     end
    nucleus_distance = 1-dot((nucleus_xy-repmat([x0,y0],size(nucleus_xy,1),1)),repmat([(x1-x0),(y1-y0)],size(nucleus_xy,1),1),2)/L2_extreme;
    if flip_EL == true
        nucleus_distance = 1-nucleus_distance;
    end
    
    nucleus_profile = [nucleus_distance,Inten'];
    fluc_profile = [nucleus_distance,Inten2'];   %%% intensity variance profile
end
index_true = (max_index' > 1) & (max_index' < size(mask_stack,3));
%%% =======================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Cytoplasmic region recognition/processing: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % [ylim0,xlim0] = size(cyto_bw);
% % xmap = repmat([1:xlim0],ylim0,1)-x0;
% % ymap = repmat([1:ylim0]',1,xlim0)-y0;
% % clear xlim0 ylim0
% % Lmap = (xmap*(x1-x0)+ymap*(y1-y0))/L2_extreme;
% % 
% % for Lcenter = 1:length(cyto_L)
% %     cyto_map = (Lmap >= cyto_L(Lcenter)-Lcwindow) & (Lmap <= cyto_L(Lcenter)+Lcwindow) & cyto_bw;
% %     cyto_mean(Lcenter) = sum(sum(cyto_map.*max_image(:,:,signal_channel)))/sum(sum(cyto_map));
% % end
% % cytoplasmic_profile = [cyto_L',cyto_mean'];
% % 
% % for Lcenter = 1:length(nu_L)
% %     nu_map = (nucleus_profile(:,1) >= nu_L(Lcenter)-Lnwindow) & (nucleus_profile(:,1) <= nu_L(Lcenter)+Lnwindow);
% %     nu_mean(Lcenter) = mean(nucleus_profile(nu_map,2));
% %     nu_std(Lcenter) = std0(nucleus_profile(nu_map,2));
% %     nu_map2 = (nucleus_profile(:,1) >= nu_L(Lcenter)-Lnwindow) & (nucleus_profile(:,1) <= nu_L(Lcenter)+Lnwindow) & index_true;
% %     nu_mean2(Lcenter) = mean(nucleus_profile(nu_map2,2));
% %     nu_std2(Lcenter) = std0(nucleus_profile(nu_map2,2));
% % end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Cytoplasmic region recognition/processing: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % exp_mask = imerode(logical(mask_stack), strel('disk',r_edge));
exp_mask = imdilate(logical(mask_stack), strel('disk',r_edge/2));
ijk = round(xy(:,[2,1,3]));
cytoplasmic_distance  = nucleus_distance;
cytoplasmic_Inten = zeros(size(cytoplasmic_distance));
cytoplasmic_var = zeros(size(cytoplasmic_distance));
cmask0 = getnhood(strel('disk',c_radius));
cid = (size(cmask0,1)-1)/2;
cjd = (size(cmask0,2)-1)/2;
ckd = cz0;
imax = size(mask_stack,1);
jmax = size(mask_stack,2);
kmax = size(mask_stack,3);
for Ic = 1:size(ijk,1)
    dimin = min(ijk(Ic,1)-1,cid);
    dimax = min(imax-ijk(Ic,1),cid);
    djmin = min(ijk(Ic,2)-1,cjd);
    djmax = min(jmax-ijk(Ic,2),cjd);
    dkmin = min(ijk(Ic,3)-1,ckd);
    dkmax = min(kmax-ijk(Ic,3),ckd);
    cmask = cmask0(cid+1-dimin:cid+1+dimax,cjd+1-djmin:cjd+1+djmax);
    try
    nmask = exp_mask(ijk(Ic,1)-dimin:ijk(Ic,1)+dimax,ijk(Ic,2)-djmin:ijk(Ic,2)+djmax,ijk(Ic,3)-dkmin:ijk(Ic,3)+dkmax);
    sig0 = signal_stack(ijk(Ic,1)-dimin:ijk(Ic,1)+dimax,ijk(Ic,2)-djmin:ijk(Ic,2)+djmax,ijk(Ic,3)-dkmin:ijk(Ic,3)+dkmax);
    cytoplasmic_Inten(Ic) = mean(sig0(repmat(cmask,[1,1,size(nmask,3)]) & (~nmask)));
    cytoplasmic_var(Ic) = var(double(sig0(repmat(cmask,[1,1,size(nmask,3)]) & (~nmask))));
    end
end
cytoplasmic_profile = [cytoplasmic_distance,cytoplasmic_Inten,cytoplasmic_var];

for Lcenter = 1:length(nu_L)
    nu_map = (nucleus_profile(:,1) >= nu_L(Lcenter)-Lnwindow) & (nucleus_profile(:,1) <= nu_L(Lcenter)+Lnwindow);
    nu_mean(Lcenter) = mean(nucleus_profile(nu_map,2));
    nu_std(Lcenter) = std0(nucleus_profile(nu_map,2));
    nu_map2 = (nucleus_profile(:,1) >= nu_L(Lcenter)-Lnwindow) & (nucleus_profile(:,1) <= nu_L(Lcenter)+Lnwindow) & index_true;
    nu_mean2(Lcenter) = mean(nucleus_profile(nu_map2,2));
    nu_std2(Lcenter) = std0(nucleus_profile(nu_map2,2));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%% Standard profile reshaping: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b2 = 3*flip_axis;

nu_x = nucleus_profile(:,1);
standard_x = [bcd_dxy(:,1);bcd_vxy(:,1)];
if mean(nucleus_profile(nucleus_profile(:,1) <= start_limit,2)) < mean(nucleus_profile(nucleus_profile(:,1) >= end_limit,2))
    nu_x = 1-nucleus_profile(:,1);
    bcd_dxy(:,1) = 1-bcd_dxy(:,1);
    bcd_vxy(:,1) = 1-bcd_vxy(:,1);
end
% beta1 = nlinfit(nu_x,nucleus_profile(:,2),@exp_C,[(max(nucleus_profile(:,2))-min(nucleus_profile(:,2))),b2,min(nucleus_profile(:,2))]);
%protein_min = beta1(3);
%protein_max = beta1(1)+beta1(3);
protein_min = min(nu_mean);
protein_max = max(nu_mean);

% beta2 = nlinfit(standard_x,[bcd_dxy(:,2);bcd_vxy(:,2)],@exp_C,[(max([bcd_dxy(:,2);bcd_vxy(:,2)])-min([bcd_dxy(:,2);bcd_vxy(:,2)])),b2,min([bcd_dxy(:,2);bcd_vxy(:,2)])]);
%standard_min = beta2(3);
%standard_max = beta2(1)+beta2(3);
standard_min = min(min(bcd_dxy(:,2)),min(bcd_vxy(:,2)));
standard_max = max(max(bcd_dxy(:,2)),min(bcd_vxy(:,2)));

bcd_dxy(:,2) = (bcd_dxy(:,2)-standard_min)*(protein_max-protein_min)/(standard_max-standard_min)+protein_min;
bcd_vxy(:,2) = (bcd_vxy(:,2)-standard_min)*(protein_max-protein_min)/(standard_max-standard_min)+protein_min;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear signal_stack

%% Image output: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(n1)
    clf
    if isempty(subtp)
        sub_nu = 0;
        sub_cyto = 0;
    else
        sub_nu = min(nucleus_profile(:,2));
        sub_cyto = min(cytoplasmic_profile(:,2));
    end
    plot(nucleus_profile(:,1),nucleus_profile(:,2)-sub_nu,'bo',cytoplasmic_profile(:,1),cytoplasmic_profile(:,2)-sub_cyto,'g*',bcd_dxy(:,1),bcd_dxy(:,2)-sub_nu,'X',bcd_vxy(:,1),bcd_vxy(:,2)-sub_nu,'X')
    hold on
    errorbar(nu_L,nu_mean-sub_nu,nu_std,'c');

    plot(nucleus_profile(index_true,1),nucleus_profile(index_true,2)-sub_nu,'ro')
    errorbar(nu_L,nu_mean-sub_nu,nu_std,'k');
    
    hold off
    title([typetxt,' concentration profile: ',image_folder,', cycle = ',num2str(N_cycle),' (',subtp,')'],'Interpreter','none')
    xlabel('AP axis (normalized)')
    ylabel('Concentration (A.U.)')
    legend('Nucleus concentration','Cytoplasmic concentration','Standard dorsal profile','Standard ventral profile','Averaged nuclear profile','Nucleus concentration cleaned','Averaged nuclear profile cleaned')

% figure(n2)
%     clf
%     nulabel = seg_bw+1;
%     nuintensity = [-inf;log2(nucleus_profile(:,2)-sub_nu)];
%     imagesc(nuintensity(nulabel))%,[min(nucleus_profile(:,2)),max(nucleus_profile(:,2))])
%     %axis off
%     prolim = max(nuintensity);
%     %colormap(jet(1024));
%     hcb = colorbar;
%     set(hcb,'YTickMode','manual') 
%     set(hcb,'YTick',prolim+[-(Nb2-1):1:0]*bin2) 
%     nylabel = cell(1,Nb2);
%     for ib2 = 1:Nb2
%         nylabel{ib2} = ['1/',num2str(2^((Nb2-ib2)*bin2))];
%     end
%     set(hcb,'YTickLabel',nylabel);
%     title(['Nucleus ',typetxt,' concentration plot: ',image_folder,', cycle = ',num2str(N_cycle),' (',subtp,')'],'Interpreter','none')
    

%%% Noise output: %%%======================================================
%p = polyfit(nucleus_profile(:,2),fluc_profile(:,2),1);
% p = polyfit(nucleus_profile(nucleus_profile(:,1) < 0.5,2),fluc_profile(nucleus_profile(:,1) < 0.5,2),1);
% p_post = polyfit(nucleus_profile(nucleus_profile(:,1) > 0.8,2),fluc_profile(nucleus_profile(:,1) > 0.8,2),1);
% % if ~flip_axis
    Itrue = (nucleus_profile(:,1) < EL_range(2)) & (nucleus_profile(:,1) >= EL_range(1)) & index_true;
    Itrue2 = (nucleus_profile(:,1) <= EL_range(4)) & (nucleus_profile(:,1) > EL_range(3)) & index_true;
    Itrue3 = (nucleus_profile(:,1) <= EL_range(6)) & (nucleus_profile(:,1) > EL_range(5)) & index_true;
% % else
% %     Itrue = (nucleus_profile(:,1) <= 1-EL_range(1)) & (nucleus_profile(:,1) > 1-EL_range(2)) & index_true;
% %     Itrue2 = (nucleus_profile(:,1) < 1-EL_range(3)) & (nucleus_profile(:,1) >= 1-EL_range(4)) & index_true;
% %     Itrue3 = (nucleus_profile(:,1) < 1-EL_range(5)) & (nucleus_profile(:,1) >= 1-EL_range(6)) & index_true;
% % end
    
p = polyfit(nucleus_profile(Itrue,2),fluc_profile(Itrue,2),1);
if any(Itrue2) && tbk0
    p_post = polyfit(nucleus_profile(Itrue2,2),fluc_profile(Itrue2,2),1);
    mean_post = 0;
elseif any(Itrue2)
    p_post = polyfit(nucleus_profile(Itrue2,2),fluc_profile(Itrue2,2),1);
    mean_post = min(mean(nucleus_profile(Itrue2,2)),mean(nucleus_profile(Itrue3,2)));
else
    p_post = p;
    mean_post = 0;
end
    
pc = polyfit(cytoplasmic_profile(index_true,2),cytoplasmic_profile(index_true,3),1);

% quanti_p = [p(1).*resolution.*resolution.*resolutionz.*6.02e8.*4.*pi.*sigma.^2.*2.*sqrt(pi).*sigmaz,p(2)/p(1)];
quanti_p = [p(1).*resolution.*resolution.*resolutionz.*6.02e8.*4.*pi.*sigma.^2.*2.*sqrt(pi).*sigmaz,p_post(2)/p_post(1)];
% quanti_p = [p(1).*resolution.*resolution.*resolutionz.*6.02e8.*4.*pi.*sigma.^2.*2.*sqrt(pi).*sigmaz,-mean_post];
quanti_pc = [pc(1).*resolution.*resolution.*resolutionz.*6.02e8.*4.*pi.*sigma.^2.*2.*sqrt(pi).*sigmaz,pc(2)/pc(1)];
quanti_p = [quanti_p,quanti_pc,resolution^2*resolutionz*6.02e8*sqrt(2*pi)*sigmaz*Inten_protein];

nucleus_profile = [nucleus_profile,fluc_profile(:,2),max_index'];

% nucleus_profile_ab = ((nucleus_profile(:,2)+p_post(2)/p_post(1))*p_post(1)-nucleus_profile(:,3))/(p_post(1)-p(1))/quanti_p(1);
% nucleus_profile_ab = (nucleus_profile(:,2)-mean(nucleus_profile(nucleus_profile(:,1) > 0.8,2)))/quanti_p(1);
nucleus_profile_ab = (nucleus_profile(:,2)-mean_post)/quanti_p(1);
nucleus_profile_ab(nucleus_profile_ab < 0) = 0;
varargout{1} = nucleus_profile_ab;


figure(n3)
    clf
    plot(nucleus_profile(:,2),fluc_profile(:,2),'b.','DisplayName','Nuclear data points')
    hold on
    plot(nucleus_profile(index_true,2),fluc_profile(index_true,2),'g.','DisplayName','Nuclear data points cleaned')
    x0 = xlim;
    plot(x0,p(1)*(x0+p_post(2)/p_post(1)),'r','DisplayName','Nuclear linear fit')
    if ~isempty(Inten_protein)
        plot(x0,p(1)*quanti_p(5)/quanti_p(1)*(x0+p_post(2)/p_post(1)),'r:','DisplayName','Nuclear linear fit2')
    end
    
    plot(cytoplasmic_profile(:,2),cytoplasmic_profile(:,3),'k.','DisplayName','Cytoplasmic data points')
    plot(cytoplasmic_profile(index_true,2),cytoplasmic_profile(index_true,3),'c.','DisplayName','Cytoplasmic data points cleaned')
    plot(x0,pc(1)*(x0+pc(2)/pc(1)),'m','DisplayName','Cytoplasmic linear fit')
    hold off
    xlabel('Mean Intensity');
    ylabel('Intensity varience');
    Cmax = max(nucleus_profile_ab);
    if ~isempty(Inten_protein)
        Cmax2 = Cmax*quanti_p(1)/quanti_p(5);
    else
        Cmax2 = nan;
    end
    title(['Protein Intensity noise plot (',image_folder,'): Var = ',num2str(p(1)),' * I + ',num2str(p(2)),', Cmax = ',num2str(Cmax),' M, Cmax2 = ',num2str(Cmax2),' M'],'Interpreter','none')
    legend('show')
%%% =======================================================================

%%% Binomial plot: %%% ====================================================
[FF,IX] = sort(FF);
dF = dF(IX);
bin_size = 20;
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

figure(n4)
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



    
figure(n2)
    clf
%     nucleus_profile_ab = (nucleus_profile(:,2)+quanti_p(2))/quanti_p(1);
%     nucleus_profile_ab(nucleus_profile_ab < 0) = 0;
    cmap = [[0,0,0];colormap];
    cmax = max(nucleus_profile_ab);
    cmin = 2e-9;
    cind_profile = ceil((size(cmap,1)-2)*(log10(nucleus_profile_ab)-log10(cmin))/(log10(cmax)-log10(cmin)))+2;
    cind_profile(cind_profile < 2) = 2;
    cind_profile(cind_profile > size(cmap,1)) = size(cmap,1);
    cind_profile = [1;cind_profile];
    
    nu_cind = cind_profile(mask_stack+1);
    nu_cind2D = double(round(sum(nu_cind-1,3)./(sum(mask_stack > 0,3)+(max(mask_stack,[],3) == 0)))+1);
%    nu_cind2D(max(mask_stack,[],3) == 0) = 1;
    
    nu_image = zeros([size(nu_cind2D),3]);
    for i_channel = 1:3
        c_temp = cmap(:,i_channel);
        try
        nu_image(:,:,i_channel) = c_temp(nu_cind2D);
        end
    end
    
    try
        imshow(nu_image)
    catch
        imshow_lxt(nu_image)
    end
    
    ctick0 = [2,5,10,20,50,100,200]*1e-9;
    ctick_label0 = {'2nM','5nM','10nM','20nM','50nM','100nM','200nM'};
    ctick = [ctick0(ctick0 < cmax),cmax];
    ctick_position = (size(cmap,1)-2)*((log10(ctick)-log10(cmin))/(log10(cmax)-log10(cmin)))+2;
    ctick_label = cat(2,ctick_label0(ctick0 < cmax),{[num2str(cmax/1e-9),'nM']});
%     colorbar('YTick',ctick_position,'YTickLabel',ctick_label)
    colorbar('Ticks',ctick_position/max(ctick_position),'TickLabels',ctick_label)    
    
    title(['Nucleus ',typetxt,' concentration plot: ',image_folder,', cycle = ',num2str(N_cycle),' (',subtp,')'],'Interpreter','none')
    
    
    
    function y = std0(x)
        y = std(x)/sqrt(length(x));
        