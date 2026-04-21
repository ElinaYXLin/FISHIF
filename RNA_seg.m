function [foci_bw,psize0,g_threshold0] = RNA_seg(seg_bw,max_image,signal_channel,WGA_channel,resolution,image_folder,N_cycle,varargin)

%% Transcription foci recognition: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Parameter setting and initialization: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
resolution0 = 0.13;
% L_ratio = (resolution/resolution0);
% psize = 0.8:0.1:2;
L_ratio = 1;
psize = 2;
Vpmax = zeros(size(psize));
Ngmax = zeros(size(psize));
Ng_th = 1e5;
max_image0 = double(max_image(:,:,signal_channel))/double(max(max(max_image(:,:,signal_channel))));
fit_radius = 3; %%% Radius of fit window for local gradient calculation (for point recognition threshold selection)
test_threshold = [0:0.002:0.5]; %%% Point recognition threshold scaning range

if ~isempty(varargin) && varargin{1} > 0 && varargin{1} < 1
    scale0 = varargin{1};
else
    scale0 = 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Point width optimization: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
for Ip = 1:length(psize)
    H = -fspecial('log',15,psize(Ip)/L_ratio);
    Ng = zeros(1,length(test_threshold)); %%% Initialization of Ng (# of recognized regions)
    Tg = zeros(1,length(test_threshold)); %%% Initialization of Tg (Threshold value)
    fit_th = zeros(1,length(test_threshold)-2*fit_radius); %%% Initialization of threshold range (for point recognition threshold profile fit)
    fit_g = zeros(1,length(test_threshold)-2*fit_radius); %%% Initialization of # of recognized regions under certain threshold (for point recognition threshold profile fit)

%%% Application of point recognition filter: %%%===========================
    g = imfilter(max_image0,H,'replicate');
    g(g<0) = 0;
%%% =======================================================================

%%% Profile (# of recognized regions under certain threshold) collection: %
    for Ig = 1:length(test_threshold)
        bw_g = im2bw(g,test_threshold(Ig));
        g_prop = regionprops(bw_g,'Area');
        Ng(Ig) = length(g_prop);
        Tg(Ig) = test_threshold(Ig);
    end
%%% =======================================================================

%%% Gradient of the profile using linear fit: %%%==========================
    for Ig = (1+fit_radius):(length(test_threshold)-fit_radius)
        fit_th(Ig-fit_radius) = test_threshold(Ig);
        temp_fit = polyfit(test_threshold((Ig-fit_radius):(Ig+fit_radius)),Ng((Ig-fit_radius):(Ig+fit_radius)),1);
        fit_g(Ig-fit_radius) = temp_fit(1);
    end
%%% =======================================================================

    %%% Selection of the threshold as the first local maxima of the gradient: %
    g_max = imregionalmin(abs(fit_g));
    g_max(1) = 0;
    II_max = intersect(find(g_max),find(Ng < Ng_th)-fit_radius);
    if ~isempty(II_max)
        Vpmax(Ip) = fit_g(II_max(1));
        Ngmax(Ip) = Ng(II_max(1)+fit_radius);
    else
        Vpmax(Ip) = -inf;
        Ngmax(Ip) = 0;
    end
end
%%% =======================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Selection of point spread width: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Find the maxima of point recognition region and reinitialization: %%%==
figure(10)
[AX,H1,H2] = plotyy(psize,Vpmax,psize,Ngmax);
title(['Vmax, Ng vs pixel size (',image_folder,') cycle = ',num2str(N_cycle)],'Interpreter','none')
xlabel('Point size (pixel)')
set(get(AX(1),'Ylabel'),'String','Ng gradient') 
set(get(AX(2),'Ylabel'),'String','# of recognized regions') 

%[~,Ip0] = min(abs(Vpmax));
[~,Ip0] = max(Vpmax);
psize0 = psize(Ip0);
H = -fspecial('log',15,psize(Ip)/L_ratio);
Ng = zeros(1,length(test_threshold)); %%% Initialization of Ng (# of recognized regions)
Tg = zeros(1,length(test_threshold)); %%% Initialization of Tg (Threshold value)
fit_th = zeros(1,length(test_threshold)-2*fit_radius); %%% Initialization of threshold range (for point recognition threshold profile fit)
fit_g = zeros(1,length(test_threshold)-2*fit_radius); %%% Initialization of # of recognized regions under certain threshold (for point recognition threshold profile fit)
g=imfilter(max_image0,H,'replicate');
g(g<0) = 0;

for Ig = 1:length(test_threshold)
    bw_g = im2bw(g,test_threshold(Ig));
    g_prop = regionprops(bw_g,'Area');
    Ng(Ig) = length(g_prop);
    Tg(Ig) = test_threshold(Ig);
end
%%% =======================================================================

%%% Gradient of the profile using linear fit: %%%==========================
for Ig = (1+fit_radius):(length(test_threshold)-fit_radius)
    fit_th(Ig-fit_radius) = test_threshold(Ig);
    temp_fit = polyfit(test_threshold((Ig-fit_radius):(Ig+fit_radius)),Ng((Ig-fit_radius):(Ig+fit_radius)),1);
    fit_g(Ig-fit_radius) = temp_fit(1);
end
%%% =======================================================================

%%% Selection of the threshold as the first local maxima of the gradient: %
g_max = imregionalmin(abs(fit_g));
g_max(1) = 0;
II_max = intersect(find(g_max),find(Ng < Ng_th)-fit_radius);
g_threshold0 = fit_th(II_max(1));
%%% =======================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Plot of the profile with the selected threshold: %%%%%%%%%%%%%%%%%%%%%%%
figure(3)
plot(Tg,Ng,g_threshold0*[1,1],[min(Ng),max(Ng)],'--')
xlabel('Threshold value (normalized)')
ylabel('# of recognized regions')
title(['Threshold - # profile and the selected threshold value (',image_folder,') ','cycle = ',num2str(N_cycle),', threshold = ',num2str(g_threshold0),', pixel size = ',num2str(psize(Ip0))],'Interpreter','none')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Transcription foci recognition: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
foci_bw = im2bw(g,g_threshold0);
foci_bw = imfill(foci_bw,'holes');

%%% Image output: %%%======================================================
foci_bw0 = imdilate(foci_bw,strel('disk',4));
bw_perim_g = bwperim(foci_bw0);
bw_perim_WGA = bwperim(seg_bw);
new_image(:,:,1:2) = max_image(:,:,[signal_channel,WGA_channel]);
new_image(:,:,2) = 0;
new_image(:,:,3) = 0;
overlay = imoverlay(new_image, bw_perim_g, [1,1,1]);
overlay = imoverlay(overlay,bw_perim_WGA,[0,0,1]);

if scale0 < 1
    overlay1 = imresize(overlay,scale0);
else
    overlay1 = overlay;
end

figure(4)
try
    imshow(overlay1)
catch
    imshow_lxt(overlay1)
end
title([image_folder,' (cycle = ',num2str(N_cycle),', white: transcription foci recognition, blue: nucleus recognition)'],'Interpreter','none')
%%% =======================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


