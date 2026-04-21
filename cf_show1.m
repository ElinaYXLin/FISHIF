function overlay = cf_show1(max_image0,g,vlayer,g_threshold0,peak_on,vcontrast,vcontrmin,vmask,varargin)

global mask_out RNA_mask

foci_bw = im2bw(g(:,:,vlayer),g_threshold0);
foci_bw = imfill(foci_bw,'holes');
%%foci_bw0 = imdilate(foci_bw,strel('disk',4));
%foci_bw0 = bwmorph(foci_bw,'thicken',4);
%bw_perim_g = bwperim(foci_bw0);
bw_perim_g = foci_bw;
%peak_on = true;
if peak_on
    %bw_peak = imregionalmax(max_image0(:,:,vlayer));% & foci_bw;
    %bw_peak = imregionalmax(g);
    %bw_peak = bw_peak(:,:,vlayer);
    if ~isempty(mask_out)
        bw_peak = imdilate(mask_out(:,:,vlayer),strel('disk',1));
    else
        bw_peak = false(size(max_image0,1),size(max_image0,2));
    end
    if ~isempty(RNA_mask) && ~isempty(varargin) && varargin{1}
        bw_peak2 = bwperim(imdilate(RNA_mask(:,:,vlayer),strel('disk',5)));
    else
        bw_peak2 = false(size(max_image0,1),size(max_image0,2));
    end
    
    overlay(:,:,1) = max((max_image0(:,:,vlayer)-vcontrmin)/(vcontrast-vcontrmin),bw_peak);
    overlay(:,:,2) = max(vmask*double(bw_perim_g),bw_peak | bw_peak2);
    overlay(:,:,3) = double(bw_peak);
else
    overlay(:,:,1) = (max_image0(:,:,vlayer)-vcontrmin)/(vcontrast-vcontrmin);
    overlay(:,:,2) = vmask*double(bw_perim_g);
    overlay(:,:,3) = 0;
end

