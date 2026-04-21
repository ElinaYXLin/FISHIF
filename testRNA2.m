psize = 0.5:0.1:2.0;
Vpmax = zeros(size(psize));

for Ip = 1:length(psize)
    H = -fspecial('log',15,psize(Ip)/(resolution/0.13));
    fit_radius = 10;
    test_threshold = [0.01:0.002:1];
    fit_t = zeros(1,length(test_threshold)-2*fit_radius);
    fit_g = zeros(length(test_threshold)-2*fit_radius,2);

    g=imfilter(max_image(:,:,1),H,'replicate');
    g(g<0) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

    for Ig = 1:length(test_threshold)
        bw_g = im2bw(g,test_threshold(Ig));
        g_prop = regionprops(bw_g,g,'MaxIntensity');
        Ng(Ig) = length(g_prop);
        Tg(Ig) = test_threshold(Ig);
    end

    for Ig = (1+fit_radius):(length(test_threshold)-fit_radius)
        fit_t(Ig-fit_radius) = test_threshold(Ig);
        fit_g(Ig-fit_radius,:) = polyfit(test_threshold((Ig-fit_radius):(Ig+fit_radius)),Ng((Ig-fit_radius):(Ig+fit_radius)),1);
    end
    fit_gg = fit_g(:,1)';
    g_min = imregionalmin(fit_gg);
    g_max = imregionalmax(fit_gg);
    g_min(1) = 0;
    %g_max(1) = 0;
    [~,Imin] = min(g_min.*fit_gg);
    I1 = find(g_max,2);
    I2 = find(g_max(1:Imin),1,'last');
    I_threshold = floor(I1(2));
    %g_threshold0 = fit_t(I_threshold);
    Vpmax(Ip) = fit_gg(I_threshold);%-fit_gg(I2);
end

figure(10)
plot(psize,Vpmax)
title('Vmax vs pixel size')

%%%%%----------------------------------------------------------------------
[~,Ip0] = max(Vpmax);
Ip0 = Ip0;
H = -fspecial('log',15,psize(Ip0)/(resolution/0.13));
fit_radius = 10;
test_threshold = [0.01:0.001:1];
fit_t = zeros(1,length(test_threshold)-2*fit_radius);
fit_g = zeros(length(test_threshold)-2*fit_radius,2);

g=imfilter(max_image(:,:,1),H,'replicate');
g(g<0) = 0;

for Ig = 1:length(test_threshold)
    bw_g = im2bw(g,test_threshold(Ig));
    g_prop = regionprops(bw_g,g,'MaxIntensity');
    Ng(Ig) = length(g_prop);
    Tg(Ig) = test_threshold(Ig);
end

for Ig = (1+fit_radius):(length(test_threshold)-fit_radius)
    fit_t(Ig-fit_radius) = test_threshold(Ig);
    fit_g(Ig-fit_radius,:) = polyfit(test_threshold((Ig-fit_radius):(Ig+fit_radius)),Ng((Ig-fit_radius):(Ig+fit_radius)),1);
end
fit_gg = fit_g(:,1)';
g_min = imregionalmin(fit_gg);
g_max = imregionalmax(fit_gg);
g_min(1) = 0;
g_max(1) = 0;
[~,Imin] = min(g_min.*fit_gg);
I1 = find(g_max,2);
I2 = find(g_max(1:Imin),1,'last');
I_threshold = floor(I1(2));
g_threshold0 = fit_t(I_threshold);
%%%%%----------------------------------------------------------------------

figure(3)
plot(Tg,Ng,g_threshold0*[1,1],[min(Ng),max(Ng)],'--')
title(['Threshold = ',num2str(g_threshold0),', Pixel size = ',num2str(psize(Ip0))])
%[x,y,button] = ginput(1);
%if button == 1
%    g_threshold0 = x;
%end
%[AX,H1,H2] = plotyy(Tg,Ng,fit_t,fit_gg);
%set(H1,'LineStyle','none','Marker','o')
%set(H2,'LineStyle','none','Marker','*')

%%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
%reform_all = imdilate(cumulated_all,strel('disk',1));
%reform_D = bwdist(~reform_all);
%reform_L = watershed(-reform_D);
%reform_all(reform_L == 0) = 0;
%reform_all = reform_all | cumulated_all;
reform_all = seg_bw;

bw_g = im2bw(g,g_threshold0);
%bw_g = imdilate(bw_g,strel('disk',1));
bw_g = imfill(bw_g,'holes');
g_prop = regionprops(bw_g,reform_all,'Area','MeanIntensity','MaxIntensity','Eccentricity');
g_threshold = 0;%max([g_prop.MaxIntensity])/2;
keepIdx = find([g_prop.MaxIntensity] >= g_threshold);
bw_g0 = logical(ismember(bwlabel(bw_g),keepIdx));
bw_g0 = imdilate(bw_g0,strel('disk',4));
bw_perim_g = bwperim(bw_g0);
bw_perim_WGA = bwperim(reform_all);
new_image(:,:,1:2) = max_image(:,:,1:2);
new_image(:,:,3) = 0;
overlay = imoverlay(new_image, bw_perim_g, [1,1,1]);
overlay = imoverlay(overlay,bw_perim_WGA,[0,0,1]);
figure(1)
imshow(overlay)
%%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
%max_signal = zeros(0);
%signal_channel = 1;
%for I_layer = 1:length(layer_name)
%    temp_image = imread([image_folder,layer_name(I_layer).name]);   %%% get images from every layer
%    gI_prop = regionprops(bw_g,temp_image(:,:,signal_channel),'MaxIntensity');
%    max_signal = cat(1,max_signal,[gI_prop.MaxIntensity]);
%end

%for Iregion = 1:size(max_signal,2)
%    figure(2)
%    plot(max_signal(:,Iregion))
%    waitforbuttonpress
%end
%%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%