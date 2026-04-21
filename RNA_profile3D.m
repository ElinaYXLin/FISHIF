function [nucleus_profile,foci_profile,cytoplasmic_profile] = RNA_profile3D(mask_stack,nucleus_DAPI_profile,cyto_bw,foci_bw3D,max_image,EL_info,SS,RNA_stack,resolution,image_folder,N_cycle,varargin)

%% Nucleus/cytoplasmic RNA profile calculation: %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Parameter setting and initialization: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%nucleus_profile = zeros(0);
%cytoplasmic_profile = zeros(0);
resolution0 = 0.09;
L_ratio = (resolution/resolution0);
f_radius = max(ceil(6/L_ratio),3);
Lmin = 0;
Lbin = 0.02;
Lmax = 1;
Lwindow = 0.02;
cyto_L = Lmin:Lbin:Lmax;
nucleus_bin = 0:0.01:1;%0.025:0.05:0.975;
average_radius = 0.03;
intensity_Nbin = 50;
cyto_mean = zeros(size(cyto_L));
% max_image = double(max_image);
% RNA_stack = double(RNA_stack);
if length(varargin) > 1
    flip_EL = varargin{2};
else
    flip_EL = false;
end

if length(varargin) < 3 || isempty(varargin{3})
    use_linear = false;
else
    use_linear = varargin{3};
end
%%% generate PSF matrix: ==================================================
% p_z = 4;
% p_r = 4;
% PSF = false([size(getnhood(strel('disk',p_r))),2*p_z+1]);
% for iz = 1:2*p_z+1
%     r_xy = round(sqrt(p_r^2-(iz-p_z-1)^2));
%     temp = getnhood(strel('disk',r_xy));
%     PSF(((size(PSF,1)-size(temp,1))/2+1):(size(PSF,1)+size(temp,1))/2,((size(PSF,2)-size(temp,2))/2+1):(size(PSF,2)+size(temp,2))/2,iz) = temp;
% end
% PSF = double(PSF);

p_z = 2;
p_r = 3;
PSF = logical(repmat(getnhood(strel('disk',2*p_r)),[1,1,p_z*2+1]));
PSF2 = getnhood(strel('disk',2*p_r));
r_integrate = 100;
cz0 = 7;
r_edge = 10;
c_radius = 100;

%%% =======================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Nucleus region recognition/processing: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nucleus_prop = regionprops(mask_stack,RNA_stack,'MeanIntensity','Centroid'); %%% Nucleus mean intensity calculation
nucleusI_prop = regionprops(mask_stack.*double(~imdilate(logical(foci_bw3D),PSF)),RNA_stack,'MeanIntensity'); %%% Nucleus background total intensity calculation
% nucleusI_prop = regionprops(mask_stack,RNA_stack.*(~imfilter(foci_bw3D,PSF,'same','corr')),'MeanIntensity'); %%% Nucleus background total intensity calculation
% nucleusA_prop = regionprops(mask_stack,(~imfilter(foci_bw3D,PSF,'same','corr')),'MeanIntensity'); %%% Nucleus background area calculation
% nucleus_background = [nucleusI_prop.MeanIntensity]./([nucleusA_prop.MeanIntensity]+([nucleusA_prop.MeanIntensity] == 0)); %%% Nucleus mean background calculation
nucleus_background = [nucleusI_prop.MeanIntensity]; %%% Nucleus mean background calculation



%%% Calculate raw foci intensity: %%%======================================
if size(SS,2) > 1
    Rfoci_prop = regionprops(foci_bw3D,max_image,'MaxIntensity','MeanIntensity','Area','Centroid'); %%%%% Raw foci intensity calculation
    Sfoci_prop = regionprops(foci_bw3D,SS,'MaxIntensity'); %%%%% equivalent foci area calculation
else
    [yind,xind,zind] = ind2sub(size(foci_bw3D),find(foci_bw3D));
    Rfoci_prop = struct('MaxIntensity',num2cell(max_image'),'MeanIntensity',num2cell(max_image'),'Centroid',mat2cell([xind,yind,zind],ones(1,size(xind,1)))');
    Sfoci_prop = struct('MaxIntensity',num2cell(SS'));
end
SSfoci = [Sfoci_prop.MaxIntensity];
% if ~isempty(varargin)
%     Ifoci0 = varargin{1};
% else
    Ifoci0 = [1,1];
% end

%%% =======================================================================

%%% Calculate the foci properties of nuclei: %%%===========================
Nfoci_nucleus = zeros(length(nucleus_prop),1); %%%%% Initialization of nucleus foci # profile
Ifoci_nucleus = zeros(length(nucleus_prop),1); %%%%% Initialization of nucleus foci intensity profile
if size(SS,2) > 1
    N_prop = regionprops(foci_bw3D,mask_stack,'MaxIntensity'); %%%%% Link foci to nuclei
else
    N_prop = struct('MaxIntensity',num2cell(mask_stack(foci_bw3D)'));
end
background0 = [0,nucleus_background];
background0(:) = 0;
% % foci_intensity2 = round(([Rfoci_prop.MaxIntensity]-background0(([N_prop.MaxIntensity]+1)))./Ifoci0(2)); %%%%% Refined foci maximal intensity
% % foci_intensity = round(([Rfoci_prop.MeanIntensity]-background0(([N_prop.MaxIntensity]+1))).*SSfoci./Ifoci0(1)); %%%%% Refined foci total intensity
foci_intensity2 = ([Rfoci_prop.MaxIntensity]-background0(([N_prop.MaxIntensity]+1)))./Ifoci0(2); %%%%% Refined foci maximal intensity
foci_intensity = ([Rfoci_prop.MeanIntensity]-background0(([N_prop.MaxIntensity]+1))).*SSfoci./Ifoci0(1); %%%%% Refined foci total intensity

for I_nucleus = 1:length(nucleus_prop)
    Nfoci_nucleus(I_nucleus) = sum([N_prop.MaxIntensity] == I_nucleus); %%%%% Nucleus foci # profile
    Ifoci_nucleus(I_nucleus) = sum(([N_prop.MaxIntensity] == I_nucleus).*foci_intensity); %%%%% Nucleus foci intensity profile
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%  Distance calculation and recording for nucleus: %%%===================
nucleus_distance = zeros(0);
nucleus_profile = zeros(0);
if size(nucleus_prop,1)>0
    xy = cell2mat({nucleus_prop.Centroid}');
    nucleus_xy = xy(:,1:2);
    x0 = EL_info(1);
    y0 = EL_info(2);
    x1 = EL_info(3);
    y1 = EL_info(4);
    L2_extreme = EL_info(5);
% % %     %%%%% Normalized coordinate calculation:
% % %     if size(mask_stack,2)/size(mask_stack,1) > 1.4
% % %         distXY = pdist2(nucleus_xy,nucleus_xy);
% % %         axis_length = max(max(distXY));
% % %         [extreme0,extreme1] = find(distXY >= axis_length,1); %%% find the anterior and posterior pole nucleus
% % %         x0 = nucleus_xy(extreme0,1);
% % %         y0 = nucleus_xy(extreme0,2);
% % %         x1 = nucleus_xy(extreme1,1);
% % %         y1 = nucleus_xy(extreme1,2);
% % %         L2_extreme = (x1-x0)^2+(y1-y0)^2;
% % %     else
% % %         I1 = nucleus_xy(:,1) < size(mask_stack,2)/2;
% % %         I2 = nucleus_xy(:,1) >= size(mask_stack,2)/2;
% % %         if std(nucleus_xy(I1,2)) <= std(nucleus_xy(I2,2))
% % %             [x0,I0] = min(nucleus_xy(:,1));
% % %             y0 = nucleus_xy(I0,2);
% % %             x1 = 2*size(mask_stack,2)-x0;
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
    if flip_EL
        nucleus_distance = 1-nucleus_distance;
    end
    
    foci_bw2D = max(foci_bw3D,[],3);
    bg_bw2D = ~imdilate(logical(foci_bw2D),PSF2);
    nucleus_background2 = get_bg(mask_stack,[nucleus_xy(:,[2,1]),nucleus_DAPI_profile(:,3)],RNA_stack,bg_bw2D,r_integrate);
    
    %%%%% nucleus data recording:    
    nucleus_profile = [nucleus_distance,[nucleus_prop.MeanIntensity]',Nfoci_nucleus,Ifoci_nucleus,nucleus_background',nucleus_background2'];
end
%%% =======================================================================

%%%  Distance calculation and recording for foci: %%%======================
foci_distance = zeros(0);
foci_profile = zeros(0);
if size(Rfoci_prop,1)>0
    xy = cell2mat({Rfoci_prop.Centroid}');
    foci_xy = xy(:,1:2);
    %%%%% Normalized coordinate calculation:
    foci_distance = 1-dot((foci_xy-repmat([x0,y0],size(foci_xy,1),1)),repmat([(x1-x0),(y1-y0)],size(foci_xy,1),1),2)/L2_extreme;
    if flip_EL
        foci_distance = 1-foci_distance;
    end
    %%%%% nucleus data recording:    
    foci_profile = [foci_distance,double([N_prop.MaxIntensity]'),foci_intensity',foci_intensity2',find(foci_bw3D)];
end
%%% =======================================================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Cytoplasmic region recognition/processing: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
exp_mask = imdilate(logical(mask_stack), strel('disk',r_edge/2));
xy = cell2mat({nucleus_prop.Centroid}');
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
    nmask = exp_mask(ijk(Ic,1)-dimin:ijk(Ic,1)+dimax,ijk(Ic,2)-djmin:ijk(Ic,2)+djmax,ijk(Ic,3)-dkmin:ijk(Ic,3)+dkmax);
    sig0 = RNA_stack(ijk(Ic,1)-dimin:ijk(Ic,1)+dimax,ijk(Ic,2)-djmin:ijk(Ic,2)+djmax,ijk(Ic,3)-dkmin:ijk(Ic,3)+dkmax);
    cytoplasmic_Inten(Ic) = mean(sig0(repmat(cmask,[1,1,size(nmask,3)]) & (~nmask)));
    cytoplasmic_var(Ic) = var(double(sig0(repmat(cmask,[1,1,size(nmask,3)]) & (~nmask))));
end
cytoplasmic_profile = [cytoplasmic_distance,cytoplasmic_Inten,cytoplasmic_var];


% % % [ylim,xlim] = size(cyto_bw);
% % % xmap = repmat([1:xlim],ylim,1)-x0;
% % % ymap = repmat([1:ylim]',1,xlim)-y0;
% % % Lmap = 1-(xmap*(x1-x0)+ymap*(y1-y0))/L2_extreme;
% % % 
% % % for Lcenter = 1:length(cyto_L)
% % %     cyto_map = (Lmap >= cyto_L(Lcenter)-Lwindow) & (Lmap <= cyto_L(Lcenter)+Lwindow) & cyto_bw & (~imdilate(logical(max(foci_bw3D,[],3)),strel('disk',f_radius)));
% % %     cyto_mean(Lcenter) = sum(sum(cyto_map.*max(RNA_stack,[],3)))/sum(sum(cyto_map));
% % % end
% % % if flip_EL
% % %     cyto_L = 1-cyto_L;
% % % end
% % % cytoplasmic_profile = [cyto_L',cyto_mean'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Image output: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(varargin) || isempty(varargin{1})
    figure(5)
else
    figure(varargin{1}(1))
end
    bin_max = min(nucleus_bin+average_radius,1);
    bin_min = max(nucleus_bin-average_radius,0);
    fi0 = zeros(size(nucleus_bin));
    fi1 = zeros(size(nucleus_bin));
    for I_bin = 1:length(nucleus_bin)
        fi0(I_bin) = mean(nucleus_profile((nucleus_distance >= bin_min(I_bin))&(nucleus_distance <= bin_max(I_bin)),4));
        fi1(I_bin) = std(nucleus_profile((nucleus_distance >= bin_min(I_bin))&(nucleus_distance <= bin_max(I_bin)),4));
    end
    
    plot(nucleus_profile(:,1),nucleus_profile(:,2),'o',nucleus_profile(:,1),nucleus_profile(:,4),'+',nucleus_profile(:,1),nucleus_profile(:,5),'^',cytoplasmic_profile(:,1),cytoplasmic_profile(:,2),'*',nucleus_bin',fi0,nucleus_bin',fi1)
    title(['FISH profile: ',image_folder,', cycle = ',num2str(N_cycle)],'Interpreter','none')
    xlabel('AP axis (normalized)')
    ylabel('Intensity (A.U.)')
    legend('Nucleus mean intensity','Nucleus foci intensity','Nucleus background mean intensity','Cytoplasmic mean intensity','Mean nucleus foci intensity','Std of nucleus foci intensity')

if isempty(varargin) || isempty(varargin{1})
    figure(6)
else
    figure(varargin{1}(2))
end
    bin_max = min(nucleus_bin+average_radius,1);
    bin_min = max(nucleus_bin-average_radius,0);
    fn0 = zeros(size(nucleus_bin));
    fn1 = zeros(size(nucleus_bin));
    fn2 = zeros(size(nucleus_bin));
    fn3 = zeros(size(nucleus_bin));
    for I_bin = 1:length(nucleus_bin)
        fn0(I_bin) = sum((nucleus_distance(Nfoci_nucleus == 0) >= bin_min(I_bin)).*(nucleus_distance(Nfoci_nucleus == 0) <= bin_max(I_bin)));
        fn1(I_bin) = sum((nucleus_distance(Nfoci_nucleus == 1) >= bin_min(I_bin)).*(nucleus_distance(Nfoci_nucleus == 1) <= bin_max(I_bin)));
        fn2(I_bin) = sum((nucleus_distance(Nfoci_nucleus == 2) >= bin_min(I_bin)).*(nucleus_distance(Nfoci_nucleus == 2) <= bin_max(I_bin)));
        fn3(I_bin) = sum((nucleus_distance(Nfoci_nucleus  > 2) >= bin_min(I_bin)).*(nucleus_distance(Nfoci_nucleus  > 2) <= bin_max(I_bin)));
    end
    %fn0 = hist(nucleus_distance(Nfoci_nucleus == 0), nucleus_bin);
    %fn1 = hist(nucleus_distance(Nfoci_nucleus == 1), nucleus_bin);
    %fn2 = hist(nucleus_distance(Nfoci_nucleus == 2), nucleus_bin);
    %fn3 = hist(nucleus_distance(Nfoci_nucleus > 2), nucleus_bin);
    fn_all = fn0+fn1+fn2+fn3+(fn0+fn1+fn2+fn3 == 0);
    plot(nucleus_bin',(1-fn0./fn_all)*100,'m',nucleus_bin',fn1./fn_all*100,'b',nucleus_bin',fn2./fn_all*100,'g',nucleus_bin',fn3./fn_all*100,'r')
    title(['Nucleus foci distribution: ',image_folder,', cycle = ',num2str(N_cycle)],'Interpreter','none')
    xlabel('AP axis (normalized)')
    ylabel('%')
    legend('nuclei with active foci','nuclei with 1 foci','nuclei with 2 foci','nuclei with >3 foci')

if isempty(varargin) || isempty(varargin{1})
    figure(7)
else
    figure(varargin{1}(3))
end
    foci_n = hist(foci_distance,nucleus_bin);
    plot(nucleus_bin',foci_n/(sum(foci_n)+(sum(foci_n) == 0))*100)
    title(['Foci position distribution: ',image_folder,', cycle = ',num2str(N_cycle)],'Interpreter','none')
    xlabel('AP axis (normalized)')
    ylabel('%')

if isempty(varargin) || isempty(varargin{1})
    figure(8)
else
    figure(varargin{1}(4))
end
    [Ifoci_n,xout] = hist(foci_intensity,intensity_Nbin);
    bar(xout,Ifoci_n/(sum(Ifoci_n)+(sum(Ifoci_n) == 0))*100)
    title(['Foci intensity distribution: ',image_folder,', cycle = ',num2str(N_cycle)],'Interpreter','none')
    xlabel('Intensity (A.U.)')
    ylabel('%')

if isempty(varargin) || isempty(varargin{1})
    figure(9)
else
    figure(varargin{1}(5))
end
    plot(foci_intensity,foci_intensity2,'.')
    title(['Foci maximal intensity vs total intensity: ',image_folder,', cycle = ',num2str(N_cycle)],'Interpreter','none')
    xlabel('Maximal intensity (A.U.)')
    ylabel('Total intensity (A.U.)')

clear nucleus_prop nucleusI_prop nucleusA_prop nucleus_background Rfoci_prop

if isempty(varargin) || isempty(varargin{1})
    figure(11)
else
    figure(varargin{1}(6))
end
    im0 = logical(ismember(mask_stack,find(Nfoci_nucleus == 0)));
    im1 = logical(ismember(mask_stack,find(Nfoci_nucleus == 1)));
    imfate(:,:,3) = max(double(im0 | im1),[],3);
    clear im1
    im2 = logical(ismember(mask_stack,find(Nfoci_nucleus == 2)));
    imfate(:,:,2) = max(double(im0 | im2),[],3);
    clear im2
    im3 = logical(ismember(mask_stack,find(Nfoci_nucleus > 2)));
    imfate(:,:,1) = max(double(im0 | im3),[],3);
    clear im3
    
    try
        imshow(imfate);
    catch
        imshow_lxt(imfate);
    end
    title(['Nuclei transcription status: ',image_folder,' (cycle = ',num2str(N_cycle),', white: 0 foci, blue: 1 foci, green: 2 foci, red: > 2 foci)'],'Interpreter','none');

    
    
    
    
    
    
function nucleus_background2 = get_bg(mask_stack,nucleus_xyz,RNA_stack,bg_bw2D,r_integrate)

embryo_mask = bwconvhull(logical(max(mask_stack,[],3)));
nucleus_xyz = round(nucleus_xyz);
xmin = max(1,nucleus_xyz(:,1)-r_integrate);
xmax = min(size(mask_stack,1),nucleus_xyz(:,1)+r_integrate);
ymin = max(1,nucleus_xyz(:,2)-r_integrate);
ymax = min(size(mask_stack,2),nucleus_xyz(:,2)+r_integrate);

nucleus_background2 = zeros(1,length(xmin));
for ii = 1:length(xmin)
    bg_mask = embryo_mask(xmin(ii):xmax(ii),ymin(ii):ymax(ii)) & bg_bw2D(xmin(ii):xmax(ii),ymin(ii):ymax(ii));
    bg_I = RNA_stack(xmin(ii):xmax(ii),ymin(ii):ymax(ii),nucleus_xyz(ii,3));
    nucleus_background2(ii) = sum(bg_I(bg_mask))/sum(bg_mask(:));
end


