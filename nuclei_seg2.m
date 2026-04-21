function [seg_bw,cyto_bw,max_image] = nuclei_seg(image_folder,image_type,WGA_channel,DAPI_channel,resolution)
%clear
%% Segmentation of nuclei image stacks: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This function recognize/segment fly embryo nuclei of a given stack of 
%%% images using WGA and DAPI signals, and output the nucleus segmentation 
%%% and cytoplasmic region segmentation in a binary image, as well as a 
%%% maximum projection image of the original stack for further processing.

%%% seg_bw (image array): a binary image representing the nuclei recognition result.
%%% cyto_bw (image array): a binary image representing the cytoplasm recognition result.
%%% max_image (multi-channel image array): the maximum projection image of the original image stack.
%%% image_folder (string): name of the folder that contains all images in the stack.
%%% image_type (string): type of the image files in the folder (for example, '*.tif').
%%% WGA_channel (integer: 1,2,3): channel # of the WGA signal.
%%% DAPI_channel (integer: 1,2,3): channel # of the DAPI signal.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Parameter setting and initiation: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%image_folder = '10022010/stacks/expd2_10022010_4_005/';
%image_folder = '10272010/stacks/10262010_1_bcd_40x_15/';
%image_type = '*.tif';
%WGA_channel = 2;
%DAPI_channel = 3;
%resolution = 0.45;
L_ratio = (resolution/0.09);
WGA_area = floor(400/L_ratio^2);
DAPI_area = floor(400/L_ratio^2);
DAPI_area2 = WGA_area;
WGA_max = floor(20000/L_ratio^2);
DAPI_max = floor(20000/L_ratio^2);
threshold_DAPI = 1;
threshold_WGA = 1;
convex_threshold = 0.1;
se = strel('disk',floor(20/L_ratio));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Image loading: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
layer_name = dir([image_folder,image_type]);
if length(layer_name) < 1
    error(['Folder error: ',image_folder,'has no image!'])
end

cumulated_WGA = false(0);
cumulated_DAPI = false(0);

for I_layer = 1:length(layer_name)
    temp_image = imread([image_folder,layer_name(I_layer).name]);   %%% get images from every layer
    temp_DAPI = temp_image(:,:,DAPI_channel);
    %max_DAPI(I_layer) = max(max(temp_DAPI));
    temp_WGA = temp_image(:,:,WGA_channel);
    %max_WGA(I_layer) = max(max(temp_WGA));
    if ~exist('max_image')
        max_image = temp_image;
    else
        max_image = max(max_image,temp_image);
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
end

temp_DAPI = max_image(:,:,DAPI_channel);
temp_WGA = max_image(:,:,WGA_channel);

%% Image analysis/segmentation: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    new_WGA = imtophat(temp_WGA, se); %%% tophat filtering for WGA channel
    new_DAPI = imtophat(temp_DAPI, se); %%% tophat filtering for DAPI channel
    for i = 1:2 %from leo's code, gaussian filtration step
        new_WGA = imfilter(new_WGA,fspecial('gaussian',2,1),'symmetric','conv');
        new_DAPI = imfilter(new_DAPI,fspecial('gaussian',2,1),'symmetric','conv');
    end
    bw_WGA = im2bw(new_WGA,graythresh(new_WGA)*threshold_WGA);
    bw_WGA = ~bw_WGA;
    %bw_WGA = imfill(bw_WGA,'holes');
    bw_WGA = imopen(bw_WGA, strel('disk',floor(5/L_ratio)));
    bw_WGA = bwareaopen(bw_WGA, WGA_area);
    bw_WGA = bwareaclose(bw_WGA, WGA_max);
    bw_WGA = imclearborder(bw_WGA, 4); %border nuclei removed
    bw_WGA0 = bw_WGA;
    bw_WGA = imfill(bw_WGA,'holes');
    
    bw_DAPI = im2bw(new_DAPI,graythresh(new_DAPI)*threshold_DAPI);
    bw_DAPI = bwareaopen(bw_DAPI, DAPI_area);
    bw_DAPI = imopen(bw_DAPI, strel('disk',floor(5/L_ratio)));
    bw_DAPI = imerode(bw_DAPI, strel('disk',floor(5/L_ratio)));
    bw_DAPI = imfill(bw_DAPI,'holes');
    bw_DAPI = bwareaclose(bw_DAPI, DAPI_max);
    bw_DAPI0 = bw_DAPI;
    
%%% Segments refinement: %%%===============================================
    bw_WGA = imdilate(bw_WGA, strel('disk',floor(20/L_ratio)));
    label_WGA = bwlabel(bw_WGA); %%% convert bw_WGA to label matrix
    DAPI_prop = regionprops(logical(bw_DAPI),bw_WGA,'ConvexImage','BoundingBox','MaxIntensity'); %%% calculate the covex image and the bouding box of DAPI signal
    WGA_prop = regionprops(logical(bw_WGA),bw_DAPI,'Area','MaxIntensity','ConvexArea','Centroid','ConvexImage','BoundingBox'); %%% calculate the area, overlaps with DAPI area, convex area, center coordinates, convex image and the bounding box of each recognized region

    if size(DAPI_prop(:),1)>0
        for I_DAPI = 1:length(DAPI_prop)
            x1 = uint16(DAPI_prop(I_DAPI).BoundingBox(2));
            x2 = uint16(x1+DAPI_prop(I_DAPI).BoundingBox(4)-1);
            y1 = uint16(DAPI_prop(I_DAPI).BoundingBox(1));
            y2 = uint16(y1+DAPI_prop(I_DAPI).BoundingBox(3)-1);
            bw_DAPI(x1:x2,y1:y2) = bw_DAPI(x1:x2,y1:y2)|DAPI_prop(I_DAPI).ConvexImage;
        end
        keepIdx = find([DAPI_prop.MaxIntensity] == 0);
        bw_DAPI_only = logical(ismember(bwlabel(bw_DAPI),keepIdx)); %stores objects that fit in the keepIdx
    end
    bw_DAPI_only = bwareaopen(bw_DAPI_only, DAPI_area2);
    
    if size(WGA_prop(:),1)>0
        for I_WGA = 1:length(WGA_prop)
            x1 = uint16(WGA_prop(I_WGA).BoundingBox(2));
            x2 = uint16(x1+WGA_prop(I_WGA).BoundingBox(4)-1);
            y1 = uint16(WGA_prop(I_WGA).BoundingBox(1));
            y2 = uint16(y1+WGA_prop(I_WGA).BoundingBox(3)-1);
            label_WGA(x1:x2,y1:y2) = max(label_WGA(x1:x2,y1:y2),I_WGA*WGA_prop(I_WGA).ConvexImage);
        end
        keepIdx = find([WGA_prop.MaxIntensity] >= 0);
        keepIdx2 = find(([WGA_prop.ConvexArea]-[WGA_prop.Area])./[WGA_prop.Area] < convex_threshold);
        keepIdx = intersect(keepIdx,keepIdx2);
        bw_WGA = logical(ismember(label_WGA,keepIdx)); %stores objects that fit in the keepIdx
        bw_WGA_negative = logical(label_WGA)-bw_WGA;
    end
    %bw_WGA = imdilate(bw_WGA,strel('disk',floor(5/L_ratio)));
    bw_WGA = bwareaclose(bw_WGA, WGA_max); 
    
        
    if ~exist('cumulated_all')
        %cumulated_WGA = bw_WGA;
        cumulated_DAPI = bw_DAPI;
        cumulated_all = bw_WGA | bw_DAPI_only;
        cumulated_negative = bw_WGA_negative;
    else
        %cumulated_WGA = cumulated_WGA | bw_WGA;
        cumulated_DAPI = cumulated_DAPI | bw_DAPI;
        cumulated_all = cumulated_all | bw_WGA | bw_DAPI_only;
        cumulated_negative = cumulated_negative | bw_WGA_negative;
    end
%%%========================================================================


seg_bw = cumulated_all; %%% nucleus recognition result output

%%% Cytoplasm region recognition: %%%======================================
embryo_region = false(size(seg_bw));
temp_label = 2*seg_bw;
temp_label(1,1) = 1;
embryo_prop = regionprops(temp_label,'ConvexImage','BoundingBox');
x1 = uint16(embryo_prop(2).BoundingBox(2));
x2 = uint16(x1+embryo_prop(2).BoundingBox(4)-1);
y1 = uint16(embryo_prop(2).BoundingBox(1));
y2 = uint16(y1+embryo_prop(2).BoundingBox(3)-1);
embryo_region(x1:x2,y1:y2) = embryo_prop(2).ConvexImage;
cyto_bw = embryo_region & (~seg_bw);
%%% =======================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Output: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
bw_perim_WGA = bwperim(cumulated_all);
bw_perim_DAPI = bwperim(cumulated_DAPI);
bw_perim_negative = bwperim(cumulated_negative);
overlay = imoverlay(adapthisteq(temp_WGA), bw_perim_WGA, [0,1,0]);
overlay = imoverlay(overlay, bw_perim_DAPI, [1,0,0]);
overlay = imoverlay(overlay, bw_perim_negative, [0,0,1]);
imshow(overlay)   %%% show the embryo boundary recognition
title([image_folder,' (green: nucleus recognition, red: DAPI signal)']);
%waitforbuttonpress
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%