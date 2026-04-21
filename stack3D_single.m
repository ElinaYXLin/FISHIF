function stack0 = stack3D_single(image_folder,channel0)

image_type = '*.tif';
imlist = dir([image_folder,image_type]); %%% get the image list from image folder
N_layer = length(imlist);
stack0 = zeros(0,'uint16');

%%% Image stack loading: %=================================================
for I_layer = 1:N_layer
    temp = imread([image_folder,imlist(I_layer).name]);
    stack0 = cat(3,stack0,temp(:,:,channel0));
end
%%% =======================================================================
