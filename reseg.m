function bwout = reseg(immask2,area1)    



    se_mask = regionprops(immask2,'Area');
    se_area = [se_mask.Area];
    bwout = false(size(immask2));   %%% Output mask initialization

    r1 = sqrt(area1/pi);   %%% estimated radius of a single nucleus

%% Partition: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bwraw = immask2;
%bwraw = bwareaopen(bwraw,round(area1*0.5));   %%% remove small areas

%%% Estimation of the number of nuclei included in each recognized regions
se_mask = regionprops(bwraw,'Area');   %%% area calculation
nmask = round([se_mask.Area]./area1);   %%% calculate the estimated # of nuclei for each recognized area
lmask = bwlabel(bwraw);   %%% mask labeling

%%% Filter out the merged nuclei
bwout(ismember(lmask,find(nmask <= 1))) = true;   %%% Add the single nuclei areas to the output mask
bwraw(ismember(lmask,find(nmask <= 1))) = false;   %%% Remove the single nuclei areas from the raw mask

if any(any(bwraw))
    %%% Watersheding to chop merged nuclei
    D= bwdist(~bwraw);
    g = imfilter(double(bwraw),double(getnhood(strel('disk',round(r1*0.75)))),'symmetric','corr');
    g2 = imimposemin(-D,imdilate(imregionalmax(g),strel('disk',5)));
    L2 = watershed(g2);
    bwnew = bwraw & (~(L2 == 0));
    bwout = bwout | bwnew;
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Convex the recognized regions: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% se_mask = regionprops(bwout,'Area');
% se_area = [se_mask.Area];
% bwout = ismember(bwlabel(bwout),find(se_area > 1000));

%bwout0 = bwout;
se_mask = regionprops(bwout,'ConvexImage','BoundingBox');
bwout = zeros(size(bwout));
bwp = zeros(size(bwout));
if size(se_mask(:),1)>0
    for I_temp = 1:length(se_mask)
        x1 = uint16(se_mask(I_temp).BoundingBox(2));
        x2 = uint16(x1+se_mask(I_temp).BoundingBox(4)-1);
        y1 = uint16(se_mask(I_temp).BoundingBox(1));
        y2 = uint16(y1+se_mask(I_temp).BoundingBox(3)-1);
        bwout(x1:x2,y1:y2) = bwout(x1:x2,y1:y2) + double(se_mask(I_temp).ConvexImage);
        bwp(x1:x2,y1:y2) = bwp(x1:x2,y1:y2) | bwperim(se_mask(I_temp).ConvexImage);
    end
end

bwout = logical(bwout == 1) & (~bwp);
bwout = imerode(bwout,strel('disk',5));
bwout = bwconvhull(bwout,'objects');
bwout = bwmorph(bwout,'thicken',5);

