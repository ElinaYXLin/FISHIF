function [overlay,g] = cf_show(max_image0,vlayer,psize0,L_ratio,g_threshold0,peak_on,vcontrast,vcontrmin,vmask,varargin)

g = zeros(size(max_image0),'uint16');
H = -fspecial('log',15,psize0/L_ratio);
for I_layer = 1:size(max_image0,3)
    g(:,:,I_layer) = imfilter(max_image0(:,:,I_layer),H,'replicate');
end
% g(g<0) = 0;

if ~isempty(varargin)
    temp = varargin{1};
else
    temp = false;
end
overlay = cf_show1(max_image0,g,vlayer,g_threshold0,peak_on,vcontrast,vcontrmin,vmask,temp);




% foci_bw = im2bw(g(:,:,vlayer),g_threshold0);
% foci_bw = imfill(foci_bw,'holes');
% %foci_bw0 = imdilate(foci_bw,strel('disk',4));
% foci_bw0 = bwmorph(foci_bw,'thicken',4);
% bw_perim_g = bwperim(foci_bw0);
% 
% if peak_on
%     bw_peak = imregionalmax(g(:,:,vlayer)) & foci_bw;
%     overlay(:,:,1) = (max_image0(:,:,vlayer)-vcontrmin)/(vcontrast-vcontrmin)+double(bw_peak);
%     overlay(:,:,2) = vmask*double(bw_perim_g)+double(bw_peak);
%     overlay(:,:,3) = 0;
% else
%     overlay(:,:,1) = (max_image0(:,:,vlayer)-vcontrmin)/(vcontrast-vcontrmin);
%     overlay(:,:,2) = vmask*double(bw_perim_g);
%     overlay(:,:,3) = 0;
% end
