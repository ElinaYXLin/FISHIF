%function [seg_bw,cyto_bw,max_image] = nuclei_seg(image_folder,image_type,WGA_channel,DAPI_channel,resolution)
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
WGA_area = floor(200/L_ratio^2);
DAPI_area = floor(200/L_ratio^2);
DAPI_area2 = WGA_area;
WGA_max = floor(5000/L_ratio^2);
DAPI_max = floor(5000/L_ratio^2);
threshold_DAPI = 1;
threshold_WGA = 1;
convex_threshold = 0.1;
convex_threshold2 = 0.1;
r_exp = 2;
area_ratio = 0.3;
se = strel('disk',floor(20/L_ratio));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Image loading: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
layer_name = dir([image_folder,image_type]);
if length(layer_name) < 1
    error(['Folder error: ',image_folder,'has no image!'])
end

cumulated_WGA = false(0);
cumulated_DAPI = false(0);
cumulated_all = false(0);

%%% Maximal intensity/threshold comparison: %%% ===========================
for I_layer = 1:length(layer_name)
    temp_image = imread([image_folder,layer_name(I_layer).name]);   %%% get images from every layer
    temp_DAPI = temp_image(:,:,DAPI_channel);
    max_DAPI(I_layer) = max(max(temp_DAPI));
    temp_WGA = temp_image(:,:,WGA_channel);
    max_WGA(I_layer) = max(max(temp_WGA));
    new_WGA = imtophat(temp_WGA, se); %%% tophat filtering for WGA channel
    new_DAPI = imtophat(temp_DAPI, se); %%% tophat filtering for DAPI channel
    for i = 1:2 %from leo's code, gaussian filtration step
        new_WGA = imfilter(new_WGA,fspecial('gaussian',2,1),'symmetric','conv');
        new_DAPI = imfilter(new_DAPI,fspecial('gaussian',2,1),'symmetric','conv');
    end
    %new_WGA = imclose(imopen(new_WGA,strel('disk',floor(5/L_ratio))),strel('disk',floor(5/L_ratio)));
    %new_DAPI = imclose(imopen(new_DAPI,strel('disk',floor(5/L_ratio))),strel('disk',floor(5/L_ratio)));
    
    th_WGA(I_layer) = graythresh(new_WGA);
    th_DAPI(I_layer) = graythresh(new_DAPI);
end

max_DAPI0 = max(max_DAPI);
max_WGA0 = max(max_WGA);
th_DAPI0 = (max(th_DAPI)+min(th_DAPI))/2;
th_WGA0 = max(th_WGA);
%%% =======================================================================

for I_layer = 1:length(layer_name)
    temp_image = imread([image_folder,layer_name(I_layer).name]);   %%% get images from every layer
    temp_DAPI = temp_image(:,:,DAPI_channel);
    max_DAPI(I_layer) = max(max(temp_DAPI));
    temp_WGA = temp_image(:,:,WGA_channel);
    max_WGA(I_layer) = max(max(temp_WGA));
    if ~exist('max_image')
        max_image = temp_image;
    else
        max_image = max(max_image,temp_image);
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

%% Image analysis/segmentation: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    new_WGA = imtophat(temp_WGA, se); %%% tophat filtering for WGA channel
    new_DAPI = imtophat(temp_DAPI, se); %%% tophat filtering for DAPI channel
    for i = 1:2 %from leo's code, gaussian filtration step
        new_WGA = imfilter(new_WGA,fspecial('gaussian',2,1),'symmetric','conv');
        new_DAPI = imfilter(new_DAPI,fspecial('gaussian',2,1),'symmetric','conv');
    end
    %new_WGA = imclose(imopen(new_WGA,strel('disk',floor(5/L_ratio))),strel('disk',floor(5/L_ratio)));
    %new_DAPI = imclose(imopen(new_DAPI,strel('disk',floor(5/L_ratio))),strel('disk',floor(5/L_ratio)));
    
    bw_WGA = im2bw(new_WGA,th_WGA0*threshold_WGA);
    bw_WGA = ~bw_WGA;
    %bw_WGA = imfill(bw_WGA,'holes');
    bw_WGA = imerode(bw_WGA, strel('disk',floor(2/L_ratio)));
    bw_WGA = imclose(bw_WGA, strel('disk',floor(5/L_ratio)));
    bw_WGA = imdilate(bw_WGA, strel('disk',floor(2/L_ratio)));
    bw_WGA = imopen(bw_WGA, strel('disk',floor(5/L_ratio)));
    bw_WGA = bwareaopen(bw_WGA, WGA_area);
    bw_WGA = bwareaclose(bw_WGA, WGA_max);
    %bw_WGA = imclearborder(bw_WGA, 4); %border nuclei removed
    bw_WGA0 = bw_WGA;
    
    bw_DAPI = im2bw(new_DAPI,th_DAPI0*threshold_DAPI);
    bw_DAPI = bwareaopen(bw_DAPI, DAPI_area);
    bw_DAPI = imerode(bw_DAPI, strel('disk',floor(2/L_ratio)));
    bw_DAPI = imclose(bw_DAPI, strel('disk',floor(5/L_ratio)));
    bw_DAPI = imdilate(bw_DAPI, strel('disk',floor(2/L_ratio)));
    bw_DAPI = imopen(bw_DAPI, strel('disk',floor(5/L_ratio)));
    bw_DAPI = imerode(bw_DAPI, strel('disk',floor(5/L_ratio)));
    bw_DAPI = imfill(bw_DAPI,'holes');
    bw_DAPI = bwareaclose(bw_DAPI, DAPI_max);
    bw_DAPI0 = bw_DAPI;
    
%%% Segments refinement: %%%===============================================
%%%%% DAPI refinement: %%%%%-----------------------------------------------
    label_WGA = bwlabel(bw_WGA); %%% convert bw_WGA to label matrix
    %bw_DAPI = imdilate(bw_DAPI, strel('disk',floor(5/L_ratio)));
    DAPI_prop = regionprops(logical(bw_DAPI),bw_WGA,'ConvexImage','BoundingBox','MaxIntensity'); %%% calculate the covex image and the bouding box of DAPI signal
    if size(DAPI_prop(:),1)>0
        %for I_DAPI = 1:length(DAPI_prop)
        %    x1 = uint16(DAPI_prop(I_DAPI).BoundingBox(2));
        %    x2 = uint16(x1+DAPI_prop(I_DAPI).BoundingBox(4)-1);
        %    y1 = uint16(DAPI_prop(I_DAPI).BoundingBox(1));
        %    y2 = uint16(y1+DAPI_prop(I_DAPI).BoundingBox(3)-1);
        %    bw_DAPI(x1:x2,y1:y2) = bw_DAPI(x1:x2,y1:y2)|DAPI_prop(I_DAPI).ConvexImage;
        %end
        keepIdx = find([DAPI_prop.MaxIntensity] == 0);
        bw_DAPI_only = logical(ismember(bwlabel(bw_DAPI),keepIdx)); %stores objects that fit in the keepIdx
    else
        bw_DAPI_only = zeros(size(bw_DAPI));
    end
    bw_DAPI_only = bwareaopen(bw_DAPI_only, DAPI_area2);
    
    if isempty(cumulated_DAPI)
        cumulated_DAPI = bw_DAPI;
    else
        cumulated_DAPI = cumulated_DAPI | bw_DAPI;
    end
%%%%% ---------------------------------------------------------------------

%%%%% WGA refinement: %%%%%------------------------------------------------    
    WGA_prop = regionprops(logical(bw_WGA),cumulated_DAPI,'Area','MaxIntensity','ConvexArea','Centroid','ConvexImage','Image','BoundingBox'); %%% calculate the area, overlaps with DAPI area, convex area, center coordinates, convex image and the bounding box of each recognized region
    
    if size(WGA_prop(:),1)>0
        keepIdx = find([WGA_prop.MaxIntensity] > 0);
        keepIdx2 = find(([WGA_prop.ConvexArea]-[WGA_prop.Area])./[WGA_prop.Area] < convex_threshold);
        keepIdx3 = setdiff(keepIdx,keepIdx2);
        keepIdx = intersect(keepIdx,keepIdx2);
        bw_WGA = logical(ismember(label_WGA,keepIdx)); %stores objects that fit in the keepIdx
        for I_WGA = keepIdx
            x1 = uint16(WGA_prop(I_WGA).BoundingBox(2));
            x2 = uint16(x1+WGA_prop(I_WGA).BoundingBox(4)-1);
            y1 = uint16(WGA_prop(I_WGA).BoundingBox(1));
            y2 = uint16(y1+WGA_prop(I_WGA).BoundingBox(3)-1);
            bw_WGA(x1:x2,y1:y2) = bw_WGA(x1:x2,y1:y2) | WGA_prop(I_WGA).ConvexImage;
        end
    end
    %bw_WGA = imdilate(bw_WGA,strel('disk',floor(5/L_ratio)));
    bw_WGA = bwareaclose(bw_WGA, WGA_max); 
%%%%% ---------------------------------------------------------------------

%%%%% Mask cumulation/refinement: %%%%%------------------------------------
    if isempty(cumulated_all)
        %cumulated_WGA = bw_WGA;
        cumulated_all = bw_WGA;% | bw_DAPI_only;
    else
        temp_all = cumulated_all | bw_WGA;% | bw_DAPI_only;
        temp_prop = regionprops(temp_all,'Area','ConvexArea','ConvexImage','BoundingBox');
        if size(WGA_prop(:),1)>0
            keepIdx = find(([temp_prop.ConvexArea]-[temp_prop.Area])./[temp_prop.Area] < convex_threshold);
            temp_all = logical(ismember(bwlabel(temp_all),keepIdx)); %stores objects that fit in the keepIdx
            for I_temp = keepIdx
                x1 = uint16(temp_prop(I_temp).BoundingBox(2));
                x2 = uint16(x1+temp_prop(I_temp).BoundingBox(4)-1);
                y1 = uint16(temp_prop(I_temp).BoundingBox(1));
                y2 = uint16(y1+temp_prop(I_temp).BoundingBox(3)-1);
                temp_all(x1:x2,y1:y2) = temp_all(x1:x2,y1:y2) | temp_prop(I_temp).ConvexImage;
            end
        end
        
        region_center = zeros(size(temp_all));
        all_prop = regionprops(cumulated_all,'Area','Centroid');
        if size(all_prop(:),1)>0
            for I_all = 1:length(all_prop)
                xcenter = ceil(all_prop(I_all).Centroid(1));
                ycenter = ceil(all_prop(I_all).Centroid(2));
                region_center(ycenter,xcenter) = 1;
            end
        end
        
        temp_prop = regionprops(temp_all,region_center,'Area','MeanIntensity');
        if size(temp_prop(:),1)>0
            keepIdx = find([temp_prop.Area].*[temp_prop.MeanIntensity] <= 1);
            cumulated_all = logical(ismember(temp_all,keepIdx)) | cumulated_all;
        end
    end
%%%%% ---------------------------------------------------------------------

    all_prop = regionprops(cumulated_all,'Area');
    if size(all_prop(:),1)>0
        WGA_area = floor(area_ratio*mean([all_prop.Area]));
        %keepIdx = find([all_prop.Area] >= 0.3*mean([all_prop.Area]));
        %cumulated_all = logical(ismember(cumulated_all,keepIdx)); %stores objects that fit in the keepIdx
    end
    cumulated_all = bwareaopen(cumulated_all, WGA_area);
%%%========================================================================


figure(1)
bw_perim_WGA = bwperim(bw_WGA0);
bw_perim_DAPI = bwperim(cumulated_DAPI);
overlay = imoverlay(adapthisteq(temp_WGA), bw_perim_WGA, [0,1,0]);
overlay = imoverlay(overlay, bw_perim_DAPI, [1,0,0]);
imshow(overlay)   %%% show the embryo boundary recognition
title([image_folder,' (green: nucleus recognition, red: DAPI signal), ',num2str(I_layer),'/',num2str(length(layer_name))]);
waitforbuttonpress
end

%%% Segment expansion: %%%=================================================
reform_all = imdilate(cumulated_all,strel('disk',r_exp));
all_prop = regionprops(reform_all,'Area','ConvexArea');
if size(all_prop(:),1)>0
    keepIdx = find(([all_prop.ConvexArea]-[all_prop.Area])./[all_prop.Area] >= convex_threshold2);
    con_all = logical(ismember(bwlabel(reform_all),keepIdx));
    %rem_all = reform_all & (~con_all);
    raw_shed = (~watershed(-bwdist(~con_all))).*con_all;
    no_shed = bwlabel(raw_shed).*imdilate(imerode(con_all, strel('disk',floor(20/L_ratio))),strel('disk',floor(10/L_ratio)));
    raw_shed(logical(ismember(bwlabel(raw_shed),setdiff(unique(no_shed),0)))) = 0;
    reform_all = reform_all&(~raw_shed);

end
cumulated_all = reform_all;


all_prop = regionprops(cumulated_all,'Area');
if size(all_prop(:),1)>0
    WGA_area = floor(area_ratio*mean([all_prop.Area]));
end
cumulated_all = bwareaopen(cumulated_all, WGA_area);
%%% =======================================================================

seg_bw = cumulated_all; %%% nucleus recognition result output

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
cyto_bw = embryo_region & (~seg_bw);
cyto_bw = imerode(cyto_bw, strel('disk',floor(5/L_ratio)));
cumulated_DAPI = imopen(cumulated_DAPI, strel('disk',floor(10/L_ratio)));
%%% =======================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Output: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
bw_perim_WGA = bwperim(cumulated_all);
bw_perim_DAPI = bwperim(cumulated_DAPI);
overlay = imoverlay(adapthisteq(max_image(:,:,WGA_channel)), bw_perim_WGA, [0,1,0]);
overlay = imoverlay(overlay, bw_perim_DAPI, [1,0,0]);
imshow(overlay)   %%% show the embryo boundary recognition
title([image_folder,' (green: nucleus recognition, red: DAPI signal)']);
%waitforbuttonpress
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%