%function [seg_bw,cyto_bw,max_image] = nuclei_seg3(image_folder,image_type,WGA_channel,DAPI_channel,resolution)
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
ex_win = 200;
WGA_threshold0 = 1;
threshold_limit = 1.5;
threshold_step = 0.1;
DAPI_threshold = 1.2;
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

bw_DAPI = im2bw(new_DAPI,DAPI_threshold*graythresh(new_DAPI));
bw_DAPI = imfill(bw_DAPI,'holes');
bw_DAPI = imclose(bw_DAPI, strel('disk',floor(5/L_ratio)));
bw_DAPI = imerode(bw_DAPI, strel('disk',floor(3/L_ratio)));
imsize = size(bw_DAPI);

dmax_WGA = double(new_WGA);
fWGA = imfilter(dmax_WGA ,H,'replicate','conv');
fWGA(fWGA < 0) = 0;
%fWGA = fWGA.*rangefilt(fWGA);
fWGA = fWGA/max(max(fWGA));
bw_WGA = im2bw(fWGA,0.5*graythresh(fWGA));
bw_WGA = imclose(bw_WGA, strel('disk',floor(5/L_ratio)));
bw_WGA = ~bw_WGA;
bw_WGA = imopen(bw_WGA, strel('disk',floor(8/L_ratio)));
WGA_prop = regionprops(logical(bw_WGA),bw_DAPI,'MaxIntensity','Eccentricity','Area','Image','ConvexArea','ConvexImage','BoundingBox');
if size(WGA_prop(:),1)>0
    keepIdx = find(([WGA_prop.ConvexArea]-[WGA_prop.Area])./[WGA_prop.ConvexArea] < convex_threshold);
    keepIdx2 = find([WGA_prop.MaxIntensity] > 0);
    keepIdx = intersect(keepIdx,keepIdx2);
    WGA_all = logical(ismember(bwlabel(bw_WGA),keepIdx)); %stores objects that fit in the keepIdx
end
bw_WGA = WGA_all;

DAPI_prop = regionprops(logical(bw_DAPI),bw_WGA,'MaxIntensity','Image','BoundingBox'); %%% calculate the covex image and the bouding box of DAPI signal
if size(DAPI_prop(:),1)>0
    keepIdxD = find([DAPI_prop.MaxIntensity] == 0);
    for I_DAPI = keepIdxD
        WGA_threshold = WGA_threshold0;
        N_WGA = 0;
        x10 = uint16(DAPI_prop(I_DAPI).BoundingBox(2));
        x20 = uint16(x10+DAPI_prop(I_DAPI).BoundingBox(4)-1);
        y10 = uint16(DAPI_prop(I_DAPI).BoundingBox(1));
        y20 = uint16(y10+DAPI_prop(I_DAPI).BoundingBox(3)-1);
        
        x1w = max((x10-ex_win),1);
        x2w = min((x20+ex_win),imsize(1));
        y1w = max((y10-ex_win),1);
        y2w = min((y20+ex_win),imsize(2));
        
        x1d = x10-x1w+1;
        x2d = x20-x1w+1;
        y1d = y10-y1w+1;
        y2d = y20-y1w+1;
        bw_DAPI00 = zeros(x2w-x1w+1,y2w-y1w+1);
        bw_DAPI00(x1d:x2d,y1d:y2d) = DAPI_prop(I_DAPI).Image;
        
        WGA00 = fWGA(x1w:x2w,y1w:y2w);
        while (~N_WGA)&&(WGA_threshold <= threshold_limit)
            bw_WGA00 = im2bw(WGA00,WGA_threshold*graythresh(WGA00));
            bw_WGA00 = imclose(bw_WGA00, strel('disk',floor(7/L_ratio)));
            bw_WGA00 = ~bw_WGA00;
            bw_WGA00 = imopen(bw_WGA00, strel('disk',floor(5/L_ratio)));
            bw_WGA00 = imclearborder(bw_WGA00, 4); %border nuclei removed
            %bw_WGA00 = bwareaopen(bw_WGA00, WGA_area);
            %bw_WGA00 = bwareaclose(bw_WGA00, WGA_max);
            WGA_prop = regionprops(logical(bw_WGA00),bw_DAPI00,'MaxIntensity');
            N_WGA = size(WGA_prop(:),1);
            if size(WGA_prop(:),1)>0
                keepIdx = find([WGA_prop.MaxIntensity] > 0);
                N_WGA = N_WGA*length(keepIdx);
            end
            WGA_threshold = WGA_threshold+threshold_step;
        end
        
        if (size(WGA_prop(:),1) > 0)&&(~isempty(keepIdx))
            bw_WGA00 = logical(ismember(bwlabel(bw_WGA00),keepIdx));
            bw_WGA(x1w:x2w,y1w:y2w) = bw_WGA(x1w:x2w,y1w:y2w) | bw_WGA00;
        else
            WGA_threshold = WGA_threshold0;
            bw_WGA00 = im2bw(WGA00,WGA_threshold*graythresh(WGA00));
            bw_WGA00 = imclose(bw_WGA00, strel('disk',floor(5/L_ratio)));
            bw_WGA11 = imclose(bw_WGA00, strel('disk',floor(10/L_ratio)));
            %bw_WGA11 = imfill(imclose(bw_WGA00, strel('disk',floor(10/L_ratio))),'holes');
            re_WGA00 = zeros(size(bw_WGA00));
            %re_WGA00 = bw_WGA11;
            re_prop = regionprops(logical(bw_WGA00),'ConvexImage','ConvexArea','BoundingBox');
            if size(re_prop(:),1)>0
                for I_re = 1:length(re_prop)
                    x1 = uint16(re_prop(I_re).BoundingBox(2));
                    x2 = uint16(x1+re_prop(I_re).BoundingBox(4)-1);
                    y1 = uint16(re_prop(I_re).BoundingBox(1));
                    y2 = uint16(y1+re_prop(I_re).BoundingBox(3)-1);
                    if any(any(re_prop(I_re).ConvexImage.*bw_DAPI00(x1:x2,y1:y2)))%&&(re_prop(I_re).ConvexArea < WGA_max)
                        re_WGA00(x1:x2,y1:y2) = re_WGA00(x1:x2,y1:y2) | re_prop(I_re).ConvexImage;
                    end
                end
                bw_WGA00 = re_WGA00 & (~bw_WGA11);
                %bw_WGA00 = imerode(bw_WGA00, strel('disk',floor(5/L_ratio)));
            end
            WGA2_prop = regionprops(logical(bw_WGA00),bw_DAPI00,'MaxIntensity');
            if (size(WGA2_prop(:),1) > 0) && any([WGA2_prop.MaxIntensity] > 0)
                keepIdx0 = find([WGA2_prop.MaxIntensity] > 0);
                %bw_WGA00 = imdilate(logical(ismember(bwlabel(bw_WGA00),keepIdx0)),strel('disk',floor(5/L_ratio))) | bw_DAPI00;
                bw_WGA00 = logical(ismember(bwlabel(bw_WGA00),keepIdx0)) | bw_DAPI00;
                bw_WGA(x1w:x2w,y1w:y2w) = bw_WGA(x1w:x2w,y1w:y2w) | bw_WGA00;
            
            %else
            %figure(1)
            %bw_perim_WGA = bwperim(bw_WGA(x1w:x2w,y1w:y2w));
            %bw_perim_DAPI = bwperim(bw_DAPI00);
            %overlay = imoverlay(adapthisteq(WGA00), bw_perim_DAPI, [1,0,0]);
            %overlay = imoverlay(overlay, bw_perim_WGA, [0,1,0]);
            %imshow(overlay)   %%% show the embryo boundary recognition
            %title([' (green: nucleus recognition, red: DAPI signal), ',num2str(I_DAPI),'/',num2str(length(DAPI_prop))]);
            %waitforbuttonpress
            end
        end
    end
end


WGA_prop = regionprops(logical(bw_WGA),bw_DAPI,'MaxIntensity','Eccentricity','Area','Image','ConvexArea','ConvexImage','BoundingBox');
if size(WGA_prop(:),1)>0
    keepIdx = find(([WGA_prop.ConvexArea]-[WGA_prop.Area])./[WGA_prop.ConvexArea] < convex_threshold);
    keepIdx3 = keepIdx;
    keepIdx2 = find([WGA_prop.MaxIntensity] > 0);
    keepIdx = intersect(keepIdx,keepIdx2);
    WGA_all = logical(ismember(bwlabel(bw_WGA),keepIdx)); %stores objects that fit in the keepIdx

    if ~isempty(keepIdx3)
        for I_temp = keepIdx3
            x1 = uint16(WGA_prop(I_temp).BoundingBox(2));
            x2 = uint16(x1+WGA_prop(I_temp).BoundingBox(4)-1);
            y1 = uint16(WGA_prop(I_temp).BoundingBox(1));
            y2 = uint16(y1+WGA_prop(I_temp).BoundingBox(3)-1);
            WGA_temp = imerode(WGA_prop(I_temp).Image, strel('disk',floor(5/L_ratio)));
            WGA_prop00 = regionprops(logical(WGA_temp),bw_DAPI(x1:x2,y1:y2),'MaxIntensity','ConvexImage','BoundingBox');
            keepIdx4 = find([WGA_prop00.MaxIntensity] > 0);
            WGA_temp = logical(ismember(bwlabel(WGA_temp),keepIdx4));
            if ~isempty(keepIdx4)
                for I_temp2 = keepIdx4
                    x11 = uint16(WGA_prop00(I_temp2).BoundingBox(2));
                    x21 = uint16(x11+WGA_prop00(I_temp2).BoundingBox(4)-1);
                    y11 = uint16(WGA_prop00(I_temp2).BoundingBox(1));
                    y21 = uint16(y11+WGA_prop00(I_temp2).BoundingBox(3)-1);
                    WGA_temp(x11:x21,y11:y21) = WGA_temp(x11:x21,y11:y21) | WGA_prop00(I_temp2).ConvexImage;
                end
            end
            WGA_temp = imdilate(WGA_prop(I_temp).Image, strel('disk',floor(3/L_ratio)));
            WGA_all(x1:x2,y1:y2) = WGA_all(x1:x2,y1:y2) | WGA_temp;
        end
    end
end

bw_WGA = WGA_all;
WGA_prop = regionprops(logical(bw_WGA),imopen(bwlabel(bw_DAPI), strel('disk',floor(5/L_ratio))),'MaxIntensity','Area','ConvexImage','BoundingBox');
if size(WGA_prop(:),1)>0
    keepIdx = find(([WGA_prop.Area] > area_ratio*mean([WGA_prop.Area])) & ([WGA_prop.Area] < area_ratio2*mean([WGA_prop.Area])));
    keepIdx2 = find([WGA_prop.MaxIntensity] > 0);
    keepIdx = intersect(keepIdx,keepIdx2);
    WGA_all = logical(ismember(bwlabel(bw_WGA),keepIdx)); %stores objects that fit in the keepIdx
    if ~isempty(keepIdx)
        for I_temp = keepIdx
            x1 = uint16(WGA_prop(I_temp).BoundingBox(2));
            x2 = uint16(x1+WGA_prop(I_temp).BoundingBox(4)-1);
            y1 = uint16(WGA_prop(I_temp).BoundingBox(1));
            y2 = uint16(y1+WGA_prop(I_temp).BoundingBox(3)-1);
            WGA_all(x1:x2,y1:y2) = WGA_all(x1:x2,y1:y2) | WGA_prop(I_temp).ConvexImage;
        end
    end
end
%bw_WGA = WGA_all;
%bw_WGA = imopen(bw_WGA, strel('disk',floor(5/L_ratio)));
%bw_WGA = bwareaopen(bw_WGA, ceil(0.25*mean([WGA_prop.Area])));
%bw_WGA = bwareaclose(bw_WGA, floor(2*mean([WGA_prop.Area])));
cumulated_all = WGA_all;


%%% Segment expansion: %%%=================================================
%reform_all = imdilate(cumulated_all,strel('disk',r_exp));
%all_prop = regionprops(reform_all,'Area','ConvexArea');
%if size(all_prop(:),1)>0
%    keepIdx = find(([all_prop.ConvexArea]-[all_prop.Area])./[all_prop.Area] >= convex_threshold2);
%    r_more = ceil(sqrt(mean([all_prop(keepIdx).Area])/pi));
%    con_all = logical(ismember(bwlabel(reform_all),keepIdx));
%    %rem_all = reform_all & (~con_all);
%    raw_shed = (~watershed(-bwdist(~con_all))).*con_all;
%    no_shed = bwlabel(raw_shed).*imdilate(imerode(con_all, strel('disk',floor(r_more/2))),strel('disk',floor(r_more/4)));
%    raw_shed(logical(ismember(bwlabel(raw_shed),setdiff(unique(no_shed),0)))) = 0;
%    reform_all = reform_all&(~raw_shed);
%end
%cumulated_all = reform_all;

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

DAPI_prop = regionprops(logical(bw_DAPI),seg_bw,'MaxIntensity','Centroid');
r_more = ceil(sqrt(area_ratio2*mean([all_prop.Area])/pi));
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
figure(1)
bw_perim_WGA = bwperim(seg_bw);
bw_perim_DAPI = bwperim(bw_DAPI);
overlay = imoverlay(adapthisteq(max_image(:,:,WGA_channel)), bw_perim_WGA, [0,1,0]);
overlay = imoverlay(overlay, bw_perim_DAPI, [1,0,0]);
imshow(overlay)   %%% show the embryo boundary recognition
title([image_folder,' (green: nucleus recognition, red: DAPI signal)']);
%waitforbuttonpress
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%imshow(bw_WGA)


