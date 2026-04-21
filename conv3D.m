function im_out = conv3D(im_in)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A function to generate 3D convex hull without merging %%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

im_out = false(size(im_in));
for I_layer = 1:size(im_in,3)
    temp_im = logical(im_in(:,:,I_layer));
    if any(temp_im(:))
        im_prop = regionprops(temp_im, 'Centroid');
        temp_conv = bwconvhull(temp_im,'objects');
        centroid_xy = round(cell2mat({im_prop.Centroid}'));
        temp_centroid = false(size(temp_im));
        temp_centroid(sub2ind(size(temp_im),centroid_xy(:,2),centroid_xy(:,1))) = true;
        D= bwdist(~temp_conv);
        g = temp_centroid;
        g2 = imimposemin(-D,g);
        L2 = watershed(g2);
        temp_conv = temp_conv & (L2 ~= 0);
        label_conv = bwlabel(temp_conv);
        raw_list = label_conv(temp_centroid);
        im_out(:,:,I_layer) = ismember(label_conv,raw_list);
    end
end