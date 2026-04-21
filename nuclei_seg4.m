function [seg_bw,cyto_bw,max_image,N_cycle,WGA_th0,DAPI_th0] = nuclei_seg4(image_folder,image_type,WGA_channel,DAPI_channel,resolution)
%clear
%% Segmentation of nuclei image stacks: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This function recognize/segment fly embryo nuclei of a given stack of 
%%% images using WGA and DAPI signals, and output the nucleus segmentation 
%%% and cytoplasmic region segmentation in a binary image, as well as a 
%%% maximum projection image of the original stack for further processing.
%%% The results are also plotted on an overlaied image.

%%% seg_bw (image array): a binary image representing the nuclei recognition result.
%%% cyto_bw (image array): a binary image representing the cytoplasm recognition result.
%%% max_image (multi-channel image array): the maximum projection image of the original image stack.
%%% image_folder (string): name of the folder that contains all images in the stack.
%%% image_type (string): type of the image files in the folder (for example, '*.tif').
%%% WGA_channel (integer: 1,2,3): channel # of the WGA signal.
%%% DAPI_channel (integer: 1,2,3): channel # of the DAPI signal.
%%% resolution (real number): pixel size in um.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Parameter setting and initialization: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%image_folder = '10052010/stacks/EXPD2_20101005_60X_Z1_4_006/';
%image_folder = '11152010/stacks/11152010_4_FISHhb_001/';
%image_folder = '10022010/stacks/expd2_10022010_4_005/';
%image_type = '*.tif';
%WGA_channel = 2;
%DAPI_channel = 1;
%resolution = 0.09;
resolution0 = 0.09;
L_ratio = (resolution/resolution0);
%WGA_area = floor(100/L_ratio^2);
%WGA_max = floor(5000/L_ratio^2);
%DAPI_max = floor(5000/L_ratio^2);
H = -fspecial('log',15,5);
bound_threshold = 50;
ex_win = 200;
WGA_threshold0 = 0.5;
threshold_limit = 1.5;
threshold_step = 0.1;
DAPI_threshold = 1;
convex_threshold = 0.3;
convex_threshold2 = 0.1;
r_exp = 0;
area_ratio = 0.25;
area_ratio2 = 2;
se = strel('disk',floor(20/L_ratio));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Image loading: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
layer_name = dir([image_folder,image_type]);
if length(layer_name) < 1
    error(['Folder error: ',image_folder,'has no image!'])
end

for I_layer = 1:length(layer_name)
    temp_image = imread([image_folder,layer_name(I_layer).name]);   %%% get images from every layer
    if ~exist('max_image')
        max_image = temp_image;
    else
        max_image = max(max_image,temp_image);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    



new_WGA = imtophat(max_image(:,:,WGA_channel), se); %%% tophat filtering for WGA channel
new_DAPI = imtophat(max_image(:,:,DAPI_channel), se); %%% tophat filtering for DAPI channel
%for i = 1:2 %from leo's code, gaussian filtration step
%    new_WGA = imfilter(new_WGA,fspecial('gaussian',5,1),'symmetric','conv');
%    new_DAPI = imfilter(new_DAPI,fspecial('gaussian',5,1),'symmetric','conv');
%end
%new_WGA = imclose(imopen(new_WGA,strel('disk',floor(2/L_ratio))),strel('disk',floor(2/L_ratio)));
%new_DAPI = imclose(imopen(new_DAPI,strel('disk',floor(2/L_ratio))),strel('disk',floor(2/L_ratio)));

DAPI_th0 = DAPI_threshold*graythresh(new_DAPI);
bw_DAPI = im2bw(new_DAPI,DAPI_th0);
bw_DAPI = imfill(bw_DAPI,'holes');
bw_DAPI = imopen(bw_DAPI, strel('disk',floor(5/L_ratio)));
bw_DAPI = imerode(bw_DAPI, strel('disk',floor(3/L_ratio)));
imsize = size(bw_DAPI);

dmax_WGA = double(new_WGA);
fWGA = imfilter(dmax_WGA ,H,'replicate','conv');
fWGA(fWGA < 0) = 0;
%fWGA = fWGA.*rangefilt(fWGA);
fWGA = fWGA/max(max(fWGA));
WGA_th0 = WGA_threshold0*graythresh(fWGA);
bw_WGA = im2bw(fWGA,WGA_th0);
bw_WGA = bwareaopen(bw_WGA, bound_threshold);
bw_WGA = imclose(bw_WGA, strel('disk',floor(8/L_ratio)));
bw_WGA = ~bw_WGA;

%bw_WGA = imopen(bw_WGA, strel('disk',floor(8/L_ratio)));
WGA_prop = regionprops(logical(bw_WGA),bw_DAPI,'MaxIntensity','Area','Image','ConvexArea','ConvexImage','BoundingBox');
if size(WGA_prop(:),1)>0
    keepIdx = find(([WGA_prop.ConvexArea]-[WGA_prop.Area])./[WGA_prop.ConvexArea] < convex_threshold);
    keepIdx2 = find([WGA_prop.MaxIntensity] > 0);
    keepIdx3 = find(([WGA_prop.Area] > 0.25*median([WGA_prop.Area]))&([WGA_prop.Area] < 2*median([WGA_prop.Area])));
    keepIdx = intersect(keepIdx,keepIdx2);
    keepIdx = intersect(keepIdx,keepIdx3);
    WGA_all = logical(ismember(bwlabel(bw_WGA),keepIdx)); %stores objects that fit in the keepIdx
    for I_temp = keepIdx
        x1 = uint16(WGA_prop(I_temp).BoundingBox(2));
        x2 = uint16(x1+WGA_prop(I_temp).BoundingBox(4)-1);
        y1 = uint16(WGA_prop(I_temp).BoundingBox(1));
        y2 = uint16(y1+WGA_prop(I_temp).BoundingBox(3)-1);
        WGA_all(x1:x2,y1:y2) = WGA_all(x1:x2,y1:y2) | WGA_prop(I_temp).ConvexImage;
    end
end
bw_WGA = WGA_all;
bw_WGA = imopen(bw_WGA, strel('disk',floor(5/L_ratio)));
%bw_WGA = bwareaopen(bw_WGA, ceil(0.25*median([WGA_prop.Area])));
%bw_WGA = bwareaclose(bw_WGA, floor(2*median([WGA_prop.Area])));
seg_bw = bwmorph(bw_WGA,'thicken',ceil(5/L_ratio));
all_prop = regionprops(seg_bw,'Area');
N_cycle = round(log2(size(all_prop(:),1)))+2;   %%% Calculate the nuclear cycle number

%%% Cytoplasm region recognition: %%%======================================
embryo_region = false(size(seg_bw));
if any(any(seg_bw))
    temp_label = 2*seg_bw;
    temp_label(1,1) = 1;
    embryo_prop = regionprops(temp_label,'ConvexImage','BoundingBox');
    x1 = uint16(embryo_prop(2).BoundingBox(2));
    x2 = uint16(x1+embryo_prop(2).BoundingBox(4)-1);
    y1 = uint16(embryo_prop(2).BoundingBox(1));
    y2 = uint16(y1+embryo_prop(2).BoundingBox(3)-1);
    embryo_region(x1:x2,y1:y2) = embryo_prop(2).ConvexImage;
end

DAPI_prop = regionprops(logical(bw_DAPI),seg_bw,'MaxIntensity','Centroid');
r_more = ceil(sqrt(area_ratio2*median([all_prop.Area])/pi));
more_bw = zeros(imsize);
if size(DAPI_prop(:),1)>0
    keepIdx = find([DAPI_prop.MaxIntensity] == 0);
    for I_more = 1:length(keepIdx)
        more_bw(ceil(DAPI_prop(keepIdx(I_more)).Centroid(2)),ceil(DAPI_prop(keepIdx(I_more)).Centroid(1))) = 1;
    end
end
more_bw = logical(conv2(more_bw,double(getnhood(strel('disk',r_more))),'same'));
cyto_bw = embryo_region & (~seg_bw) & (~more_bw);
cyto_bw = imerode(cyto_bw, strel('disk',floor(5/L_ratio)));
%bw_DAPI = imopen(bw_DAPI, strel('disk',floor(10/L_ratio)));
%%% =======================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Output: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(51)
bw_perim_WGA = bwperim(seg_bw);
bw_perim_DAPI = bwperim(bw_DAPI);
overlay = imoverlay(adapthisteq(max_image(:,:,WGA_channel)), bw_perim_WGA, [0,1,0]);
overlay = imoverlay(overlay, bw_perim_DAPI, [1,0,0]);
imshow(overlay)   %%% show the embryo boundary recognition
title([image_folder,' (green: nucleus recognition, red: DAPI signal), cycle = ',num2str(N_cycle)]);
%%%waitforbuttonpress
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%imshow(bw_WGA)


