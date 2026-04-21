function [nucleus_profile,foci_profile,cytoplasmic_profile] = RNA_profile(seg_bw,cyto_bw,foci_bw,max_image,signal_channel,resolution,image_folder,N_cycle,varargin)

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
max_image = double(max_image);

if length(varargin) > 2 && varargin{3} > 0 && varargin{3} < 1
    scale0 = varargin{3};
else
    scale0 = 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Nucleus region recognition/processing: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nucleus_prop = regionprops(seg_bw,max_image(:,:,signal_channel),'MeanIntensity','Centroid'); %%% Nucleus mean intensity calculation
nucleusI_prop = regionprops(seg_bw,max_image(:,:,signal_channel).*(~imdilate(foci_bw,strel('disk',f_radius))),'MeanIntensity'); %%% Nucleus background total intensity calculation
nucleusA_prop = regionprops(seg_bw,(~imdilate(foci_bw,strel('disk',f_radius))),'MeanIntensity'); %%% Nucleus background area calculation
nucleus_background = [nucleusI_prop.MeanIntensity]./([nucleusA_prop.MeanIntensity]+([nucleusA_prop.MeanIntensity] == 0)); %%% Nucleus mean background calculation

%%% Calculate raw foci intensity: %%%======================================
Rfoci_prop = regionprops(foci_bw,max_image(:,:,signal_channel),'MaxIntensity','MeanIntensity','Area','Centroid'); %%%%% Raw foci intensity calculation
if ~isempty(varargin) && ~isempty(varargin{1})
    Sfoci_prop = regionprops(foci_bw,varargin{1},'MaxIntensity'); %%%%% equivalent foci area calculation
    SSfoci = [Sfoci_prop.MaxIntensity];
    if length(varargin) > 1 && ~isempty(varargin{2})
        Ifoci0 = varargin{2};
    else
        Ifoci0 = [1,1];
    end
else
    SSfoci = [Rfoci_prop.Area];
    Ifoci0 = [1,1];
end
%%% =======================================================================

%%% Calculate the foci properties of nuclei: %%%===========================
Nfoci_nucleus = zeros(length(nucleus_prop),1); %%%%% Initialization of nucleus foci # profile
Ifoci_nucleus = zeros(length(nucleus_prop),1); %%%%% Initialization of nucleus foci intensity profile
N_prop = regionprops(foci_bw,bwlabel(seg_bw),'MaxIntensity'); %%%%% Link foci to nuclei
background0 = [0,nucleus_background];
foci_intensity2 = round(([Rfoci_prop.MaxIntensity]-background0(([N_prop.MaxIntensity]+1)))./Ifoci0(2)); %%%%% Refined foci maximal intensity
foci_intensity = round(([Rfoci_prop.MeanIntensity]-background0(([N_prop.MaxIntensity]+1))).*SSfoci./Ifoci0(1)); %%%%% Refined foci total intensity

for I_nucleus = 1:length(nucleus_prop)
    Nfoci_nucleus(I_nucleus) = sum([N_prop.MaxIntensity] == I_nucleus); %%%%% Nucleus foci # profile
    Ifoci_nucleus(I_nucleus) = sum(([N_prop.MaxIntensity] == I_nucleus).*foci_intensity); %%%%% Nucleus foci intensity profile
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%  Distance calculation and recording for nucleus: %%%===================
nucleus_distance = zeros(0);
nucleus_profile = zeros(0);
if size(nucleus_prop,1)>0
    xy = [nucleus_prop.Centroid];
    nucleus_xy = [xy(1:2:length(xy)-1);xy(2:2:length(xy))]';
    %%%%% Normalized coordinate calculation:
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

    %%%%% nucleus data recording:    
    nucleus_profile = [nucleus_distance,[nucleus_prop.MeanIntensity]',Nfoci_nucleus,Ifoci_nucleus,nucleus_background'];
end
%%% =======================================================================

%%%  Distance calculation and recording for foci: %%%======================
foci_distance = zeros(0);
foci_profile = zeros(0);
if size(Rfoci_prop,1)>0
    xy = [Rfoci_prop.Centroid];
    foci_xy = [xy(1:2:length(xy)-1);xy(2:2:length(xy))]';
    %%%%% Normalized coordinate calculation:
    foci_distance = dot((foci_xy-repmat([x0,y0],size(foci_xy,1),1)),repmat([(x1-x0),(y1-y0)],size(foci_xy,1),1),2)/L2_extreme;
    %%%%% nucleus data recording:    
    foci_profile = [foci_distance,[N_prop.MaxIntensity]',foci_intensity',foci_intensity2'];
end
%%% =======================================================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Cytoplasmic region recognition/processing: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[ylim,xlim] = size(cyto_bw);
xmap = repmat([1:xlim],ylim,1)-x0;
ymap = repmat([1:ylim]',1,xlim)-y0;
Lmap = (xmap*(x1-x0)+ymap*(y1-y0))/L2_extreme;

for Lcenter = 1:length(cyto_L)
    cyto_map = (Lmap >= cyto_L(Lcenter)-Lwindow) & (Lmap <= cyto_L(Lcenter)+Lwindow) & cyto_bw & (~imdilate(foci_bw,strel('disk',f_radius)));
    cyto_mean(Lcenter) = sum(sum(cyto_map.*max_image(:,:,signal_channel)))/sum(sum(cyto_map));
end
cytoplasmic_profile = [cyto_L',cyto_mean'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Image output: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(5)
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

figure(6)
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

figure(7)
    foci_n = hist(foci_distance,nucleus_bin);
    plot(nucleus_bin',foci_n/(sum(foci_n)+(sum(foci_n) == 0))*100)
    title(['Foci position distribution: ',image_folder,', cycle = ',num2str(N_cycle)],'Interpreter','none')
    xlabel('AP axis (normalized)')
    ylabel('%')

figure(8)
    [Ifoci_n,xout] = hist(foci_intensity,intensity_Nbin);
    bar(xout,Ifoci_n/(sum(Ifoci_n)+(sum(Ifoci_n) == 0))*100)
    title(['Foci intensity distribution: ',image_folder,', cycle = ',num2str(N_cycle)],'Interpreter','none')
    xlabel('Intensity (A.U.)')
    ylabel('%')

figure(9)
    plot(foci_intensity,foci_intensity2,'.')
    title(['Foci maximal intensity vs total intensity: ',image_folder,', cycle = ',num2str(N_cycle)],'Interpreter','none')
    xlabel('Maximal intensity (A.U.)')
    ylabel('Total intensity (A.U.)')

clear nucleus_prop nucleusI_prop nucleusA_prop nucleus_background Rfoci_prop

figure(11)
    im0 = logical(ismember(bwlabel(seg_bw),find(Nfoci_nucleus == 0)));
    im1 = logical(ismember(bwlabel(seg_bw),find(Nfoci_nucleus == 1)));
    imfate(:,:,3) = double(im0 | im1);
    clear im1
    im2 = logical(ismember(bwlabel(seg_bw),find(Nfoci_nucleus == 2)));
    imfate(:,:,2) = double(im0 | im2);
    clear im2
    im3 = logical(ismember(bwlabel(seg_bw),find(Nfoci_nucleus > 2)));
    imfate(:,:,1) = double(im0 | im3);
    clear im3
    
    if scale0 < 1
        imfate1 = imresize(imfate,scale0);
    else
        imfate1 = imfate;
    end

    try
        imshow(imfate1);
    catch
        imshow_lxt(imfate1);
    end
    title(['Nuclei transcription status: ',image_folder,' (cycle = ',num2str(N_cycle),', white: 0 foci, blue: 1 foci, green: 2 foci, red: > 2 foci)'],'Interpreter','none');
