max_WGA = max_image(:,:,2);
max_DAPI = max_image(:,:,1);

max_WGA = imtophat(max_WGA, strel('disk',floor(40/L_ratio))); %%% tophat filtering for WGA channel
max_DAPI = imtophat(max_DAPI, strel('disk',floor(40/L_ratio))); %%% tophat filtering for DAPI channel

%for i = 1:2 %from leo's code, gaussian filtration step
%    max_WGA = imfilter(max_WGA,fspecial('gaussian',8,1),'symmetric','conv');
%    max_DAPI = imfilter(max_DAPI,fspecial('gaussian',8,1),'symmetric','conv');
%end
for i = 1:2
    max_WGA = imclose(imopen(max_WGA,strel('disk',floor(5/L_ratio))),strel('disk',floor(5/L_ratio)));
    max_DAPI = imclose(imopen(max_DAPI,strel('disk',floor(5/L_ratio))),strel('disk',floor(5/L_ratio)));
end

bw_max = im2bw(max_WGA,1.2*graythresh(max_WGA));
bw_max = imopen(bw_max, strel('disk',floor(5/L_ratio)));
bw_max = imclose(bw_max, strel('disk',floor(5/L_ratio)));

bw_max2 = im2bw(max_DAPI,1.05*graythresh(max_DAPI));
bw_max2 = imclose(bw_max2, strel('disk',floor(8/L_ratio)));
bw_max2 = imopen(bw_max2, strel('disk',floor(8/L_ratio)));
bw_max2 = imfill(bw_max2,'holes');
bw_max2 = bwareaopen(bw_max2, 300);

%DAPI_prop = regionprops(logical(bw_max2),'ConvexImage','BoundingBox');
%if size(DAPI_prop(:),1)>0
%    for I_DAPI = 1:size(DAPI_prop(:),1)
%        x1 = uint16(DAPI_prop(I_DAPI).BoundingBox(2));
%        x2 = uint16(x1+DAPI_prop(I_DAPI).BoundingBox(4)-1);
%        y1 = uint16(DAPI_prop(I_DAPI).BoundingBox(1));
%        y2 = uint16(y1+DAPI_prop(I_DAPI).BoundingBox(3)-1);
%        bw_max2(x1:x2,y1:y2) = bw_max2(x1:x2,y1:y2) | DAPI_prop(I_DAPI).ConvexImage;
%    end
%end
%bw_max2 = imerode(bw_max2, strel('disk',floor(5/L_ratio)));

bw_max0 = ~bw_max;
bw_max0 = bwareaopen(bw_max0, 150);
%bw_max0 = bwareaclose(bw_max0, WGA_max);
bw_max00 = bwareaopen(bw_max0, WGA_max);
max_prop = regionprops(logical(bw_max0),bw_max2,'MaxIntensity');
keepIdx = find([max_prop.MaxIntensity] == 0);
bw_max0 = logical(ismember(bwlabel(bw_max0),keepIdx));
bw_max0 = bw_max0 | bw_max00;
bw_max0 = ~bw_max0;
bw_max0 = bwareaopen(bw_max0, 300);
imshow(bw_max0)
