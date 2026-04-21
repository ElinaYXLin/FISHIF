function max_spot = manual_foci(input_folder,RNA_channel,add_id,foci_list,varargin)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A function to fit added foci %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% new_list: parameters of new foci.
%% input_folder: folder for FISH images.
%% channel_value: FISH channel.
%% add_id: locations of added foci candidates.
%% foci_list: existing foci information.
%% varargin: {1} [zmin,zmax] ==> z position range
%%           {2} xrange
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Data loading: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
image_tail = '*.tif';
rxy = 10;

lsm_stack = dir([input_folder,image_tail]); %%% image file list loading
immax = length(lsm_stack);
imstack0 = zeros([size(add_id),immax],'uint16');
for I_layer = 1:immax
    temp0 = imread([input_folder,lsm_stack(I_layer).name]);
    imstack0(:,:,I_layer) = temp0(:,:,RNA_channel);
end
imstack = imfilter(imstack0,fspecial('gaussian',8,2),'symmetric','conv');

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


spot_id0 = find(add_id)';
% spot_id = repmat(spot_id0,[immax,1])+repmat([0:immax-1]',[1,numel(spot_id0)])*numel(add_id);
spot_id = repmat(spot_id0,[(zmax-zmin+1),1])+repmat([(zmin-1):(zmax-1)]',[1,numel(spot_id0)])*numel(add_id);
[~,Iz] = max(imstack(spot_id));
idz = spot_id(Iz+[0:(numel(spot_id0)-1)]*(zmax-zmin+1));
mask_out = false(size(imstack0));
mask_out(idz) = true;

foci_xy = max(round(foci_list(:,6:7)),1);
id_list = sub2ind(size(add_id),foci_xy(:,1),foci_xy(:,2));
foci_2D = false(size(add_id));
foci_2D(id_list) = true;
mask_out2D = repmat(foci_2D | add_id,[1,1,immax]);
mask_out2D = mask_out2D & (~mask_out);

peak_parameter = ptrack(imstack0,imstack,mask_out,mask_out2D,false,xrange);   %%% Track/fit mRNA spots
max_spot = spfilter(peak_parameter,rxy,rxy,1);   %%% Fine filter to get rid of noise
