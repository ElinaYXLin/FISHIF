function im_out = dilate3D(im_in,r)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A function to dilate the 3D mask in xy without merging %%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% im_out (Output): Dilated 3D image (with label);                       %%
%% im_in (Input): Original 3D image (with label);                        %%
%% r (Input): Dilate radius (>0: dilation, <0: erosion).                 %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

im_out = zeros(size(im_in));

if r > 0
    for I_layer = 1:size(im_in,3)
        temp_im = bwlabel(bwmorph(logical(im_in(:,:,I_layer)),'thicken',r));
        im_prop = regionprops(temp_im,im_in(:,:,I_layer),'MaxIntensity');
        label_index = [0,[im_prop.MaxIntensity]];
        im_out(:,:,I_layer) = label_index(temp_im+1);
    end
elseif r < 0
    for I_layer = 1:size(im_in,3)
        temp_im = imerode(logical(im_in(:,:,I_layer)),strel('disk',-r));
        im_out(:,:,I_layer) = im_in(:,:,I_layer).*temp_im;
    end
else
    im_out = im_in;
end


