function [mask_stack,signal_stack,RNA_stack,varargout] = mask3D(seg_bw,signal_channel,DAPI_channel,RNA_channel,image_folder)

image_type = '*.tif';
imlist = dir([image_folder,image_type]); %%% get the image list from image folder
thresh3 = 0.5;
N_layer = length(imlist);
DAPI_stack = zeros([size(seg_bw),N_layer]);
signal_stack = zeros([size(seg_bw),N_layer]);
mask_stack = zeros([size(seg_bw),N_layer]);
RNA_stack = zeros([size(seg_bw),N_layer]);
z_DAPI = zeros(0);


%%% Image stack loading: %=================================================
for I_layer = 1:N_layer
    temp = imread([image_folder,imlist(I_layer).name]);
    DAPI_stack(:,:,I_layer) = temp(:,:,DAPI_channel);
    signal_stack(:,:,I_layer) = temp(:,:,signal_channel);
    RNA_stack(:,:,I_layer) = temp(:,:,RNA_channel);
    STATS = regionprops(seg_bw,DAPI_stack(:,:,I_layer),'MeanIntensity');
    z_DAPI(I_layer,:) = double([STATS.MeanIntensity]);
end
%%% =======================================================================

%%% 3D mask generation: ===================================================
z_DAPI = z_DAPI./repmat(max(z_DAPI,[],1),size(z_DAPI,1),1);
region3 = z_DAPI >= thresh3;
region3 = cat(2,zeros(N_layer,1),region3);
bwn = bwlabel(seg_bw);

for I_layer = 1:N_layer
    temp_DAPI = region3(I_layer,:);
    mask_stack(:,:,I_layer) = temp_DAPI(bwn+1);
end
%%% =======================================================================

varargout = {DAPI_stack};