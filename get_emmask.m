function em_mask = get_emmask(image_folder,varargin)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A function to extract embryo area mask %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% image_folder (input): embryo image folder;
%% em_mask (output): embryo area mask;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Data loading: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
im_tail = 'stack*.tif';
fname = dir([image_folder,im_tail]);
im0 = imread([image_folder,fname(end).name]);
if isempty(varargin) || isempty(varargin{1})
    th0 = 0.6;
else
    th0 = varargin{1};
end
if length(varargin) < 2 || isempty(varargin{2})
    ch0 = 1;
else
    ch0 = varargin{2};
end
% % % im0 = max_image;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Image processing: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % mask0 = im2bw(im0(:,:,1),0.1*graythresh(im0(:,:,1)));
im000 = mean(double(im0(:,:,ch0)),3);
im00 = im000/max(im000(:));
im00(im2bw(im00,graythresh(im00))) = graythresh(im00);
mask0 = im2bw(im00,th0*graythresh(im00));
% mask0 = imopen(mask0,strel('disk',10));
mask0 = imclose(mask0,strel('disk',50));
mask1 = bwareaopen(mask0,1000);
mask1 = imfill(mask1,'holes');
mask1 = imopen(mask1,strel('disk',5));
prop0 = regionprops(mask1,'Area');
[~,Imax] = max([prop0.Area]);
em_mask = bwlabel(mask1) == Imax;
% em_mask = bwconvhull(em_mask);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



