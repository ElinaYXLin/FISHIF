function max_spot = manual_foci3D(input_folder,RNA_channel,mask_out,foci_list,varargin)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A function to fit added foci with 3D coordinates %%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% new_list: parameters of new foci.
%% input_folder: folder for FISH images.
%% channel_value: FISH channel.
%% add_id: locations of added foci candidates.
%% foci_list: existing foci information.
%% varargin: {1} [zmin,zmax] ==> z position range
%%           {2} xrange
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Data loading: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
image_tail = '*.tif';
rxy = 10;

lsm_stack = dir([input_folder,image_tail]); %%% image file list loading
immax = length(lsm_stack);
imstack0 = zeros(size(mask_out));
for I_layer = 1:immax
    temp0 = imread([input_folder,lsm_stack(I_layer).name]);
    imstack0(:,:,I_layer) = temp0(:,:,RNA_channel);
end

if ~isempty(varargin) && ~isempty(varargin{1})
    zmin = varargin{1}(1);
    zmax = varargin{1}(2);
else
    zmin = 1;
    zmax = immax;
end    

if length(varargin) >= 2 && ~isempty(varargin{2})
    xrange = varargin{2};
else
    xrange = 5;
end


mask_out(:,:,1:(zmin-1)) = false;
mask_out(:,:,(zmax+1):end) = false;

foci_xy = max(round(foci_list(:,6:7)),1);
id_list = sub2ind(size(temp0(:,:,1)),foci_xy(:,1),foci_xy(:,2));
foci_2D = false(size(temp0(:,:,1)));
foci_2D(id_list) = true;
mask_out2D = repmat(foci_2D,[1,1,immax]);
mask_out2D = mask_out2D & (~mask_out);

peak_parameter = ptrack(imstack0,imstack0,mask_out,mask_out2D,false,xrange);   %%% Track/fit mRNA spots
max_spot = peak_parameter;