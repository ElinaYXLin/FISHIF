function [overlay,g] = foci_show(max_image0,psize0,L_ratio,g_threshold0,seg_bw,peak_on)

H = -fspecial('log',15,psize0/L_ratio);
g = imfilter(max_image0,H,'replicate');
g(g<0) = 0;
foci_bw = im2bw(g,g_threshold0);
foci_bw = imfill(foci_bw,'holes');
%foci_bw0 = imdilate(foci_bw,strel('disk',4));
foci_bw0 = bwmorph(foci_bw,'thicken',4);
bw_perim_g = bwperim(foci_bw0);
bw_perim_WGA = bwperim(seg_bw);
if peak_on
    bw_peak = imregionalmax(g) & foci_bw;
    overlay(:,:,1) = max_image0+double(bw_peak);
    overlay(:,:,2) = double(bw_perim_g)+double(bw_peak);
    overlay(:,:,3) = double(bw_perim_WGA)+double(bw_peak);
else
    overlay(:,:,1) = max_image0;
    overlay(:,:,2) = double(bw_perim_g);
    overlay(:,:,3) = double(bw_perim_WGA);
end
