function [seg_bw,cyto_bw,max_image,N_cycle,WGA_th0,DAPI_th0] = nuclei_segDAPI(image_folder,image_type,DAPI_channel,resolution,varargin)
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
%%% varargin (real number): {1} scaling of the output image
%%%                         {2} minimum area of nuclei
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Parameter setting and initialization: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
resolution0 = 0.0915;
% L_ratio = (resolution/resolution0);
L_ratio = 1;
%se = strel('disk',floor(20/L_ratio));
if ~isempty(varargin) && varargin{1} > 0 && varargin{1} < 1
    scale0 = varargin{1};
else
    scale0 = 1;
end

if length(varargin) >= 2 && varargin{2} > 0
    Smin = varargin{2};
else
    Smin = 0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Image loading: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
layer_name = dir([image_folder,image_type]);
if length(layer_name) < 1
    error(['Folder error: ',image_folder,'has no image!'])
end

for I_layer = 1:length(layer_name)
    temp_image = imread([image_folder,layer_name(I_layer).name]);   %%% get images from every layer
    %imstack(:,:,I_layer) = imfilter(double(temp_image(:,:,DAPI_channel)),fspecial('gaussian',10,3/L_ratio),'symmetric','conv');
    if ~exist('max_image')
        max_image = temp_image;
    else
        max_image = max(max_image,temp_image);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

% imstack = double(imstack)./max(max(max(double(imstack))));   %%% Normalization
% immask = zeros(size(imstack));
% DAPI_th0 = graythresh(imstack)*1.2;
max_DAPI = max_image(:,:,DAPI_channel);
%max_DAPI = imtophat(max_image(:,:,DAPI_channel),strel('disk',floor(10/L_ratio)));
max_DAPI = imfilter(max_DAPI,fspecial('gaussian',10,3/L_ratio),'symmetric','conv');
DAPI_th0 = graythresh(max_DAPI);
WGA_th0 = DAPI_th0;
%% Raw mask calculation: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for I_layer = 1:size(imstack,3)
%     immask2 = im2bw(imstack(:,:,I_layer),DAPI_th0);
%     immask2 = imclose(immask2,strel('disk',round(3/L_ratio)));
%     immask2 = imopen(immask2,strel('disk',round(20/L_ratio)));
%     immask2 = bwareaopen(immask2,round(1000/L_ratio/L_ratio));
%     immask(:,:,I_layer) = immask2;
% end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 2D Sup mask fine segmentation: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% maskall = logical(max(immask,[],3));
maskall = im2bw(max_DAPI,DAPI_th0*1);
maskall = imopen(maskall,strel('disk',round(3/L_ratio)));
maskall = imclose(maskall,strel('disk',round(3/L_ratio)));
maskall = imopen(maskall,strel('disk',round(20/L_ratio)));
maskall = bwareaopen(maskall,round(1000/L_ratio/L_ratio));

se_mask = regionprops(maskall,'Area');
se_area = [se_mask.Area];
area1 = geomean(se_area(se_area > 1000/L_ratio/L_ratio));   %%% calculate the mean area for a single nucleus
seg_bw = reseg(maskall,area1);
if Smin > 0
    seg_bw = bwareaopen(seg_bw,Smin);
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Cytoplasm region recognition: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
embryo_region = bwconvhull(seg_bw);
cyto_bw = embryo_region & (~seg_bw);
cyto_bw = imerode(cyto_bw, strel('disk',floor(5/L_ratio)));
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Nuclei cycle estimation: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
all_prop = regionprops(seg_bw,'Area');
N_cycle = round(log2(size(all_prop(:),1)))+2;   %%% Calculate the nuclear cycle number
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Output: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(51)
% bw_perim_WGA = bwperim(seg_bw);
% bw_perim_DAPI = bwperim(bw_DAPI);
% overlay = imoverlay(adapthisteq(max_image(:,:,DAPI_channel)), bw_perim_WGA, [0,1,0]);
% overlay = imoverlay(overlay, bw_perim_DAPI, [1,0,0]);
% imshow(overlay)   %%% show the embryo boundary recognition
% title([image_folder,' (green: nucleus recognition, red: DAPI signal), cycle = ',num2str(N_cycle)]);

bw_perim_DAPI = bwperim(seg_bw);
overlay = zeros([size(seg_bw),3]);
overlay(:,:,1) = 0;
overlay(:,:,2) = bw_perim_DAPI;
overlay(:,:,3) = double(max_image(:,:,DAPI_channel))./max(max(double(max_image(:,:,DAPI_channel))));

if scale0 < 1
    overlay1 = imresize(overlay,scale0);
else
    overlay1 = overlay;
end

try
    imshow(overlay1)
catch
    imshow_lxt(overlay1)
end
title([image_folder,' (green: nucleus recognition, blue: DAPI signal), cycle = ',num2str(N_cycle)],'Interpreter','none');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear imstack

