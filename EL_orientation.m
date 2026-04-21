function flip_EL = EL_orientation(nucleus_DAPI_profile,EL_info,mask_stack,signal_stack,flip_axis,varargin)


%% Determine the orientation of EL axis: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Parameter setting and initiation: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(varargin) && ~isempty(varargin{1})
    EL_range = varargin{1};
else
    EL_range = [0.15,0.4,0.6,0.85];
end

resolution0 = 0.083;
if length(varargin) < 2 || isempty(varargin{2})
    L_ratio = 1;
else
    L_ratio = varargin{2}/resolution0;
end
r_edge = round(10/L_ratio);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Nucleus region recognition/processing: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
temp = zeros(0);
Inten = zeros(0);


for I_layer = 1:size(mask_stack,3)
    sg_prop = regionprops(double(imerode(logical(mask_stack(:,:,I_layer)), strel('disk',r_edge))).*mask_stack(:,:,I_layer),double(signal_stack(:,:,I_layer)),'MeanIntensity','PixelValues');
    sg2_prop = regionprops(double(imerode(logical(mask_stack(:,:,I_layer)), strel('disk',r_edge))).*mask_stack(:,:,I_layer),double(signal_stack(:,:,I_layer)).^2,'MeanIntensity');
    if I_layer == 1
        temp = zeros(size(mask_stack,3),max(mask_stack(:)));
    end
    if ~isempty(sg_prop)
        temp(I_layer,1:length(sg2_prop)) = [sg_prop.MeanIntensity];
    end
end
sg_prop = regionprops(mask_stack,'Centroid');
centr = cell2mat({sg_prop.Centroid}');
temp(isnan(temp)) = 0;
% [Inten,max_index] = max(temp,[],1);
max_index = nucleus_DAPI_profile(:,3)';
Inten = temp(max_index+size(temp,1)*[0:(size(temp,2)-1)]);


flip_EL = false;
%%% Normalized coordinate calculation: %%%=================================
if ~isempty(centr)
    xy = centr;
    nucleus_xy = xy(:,1:2);
    x0 = EL_info(1);
    y0 = EL_info(2);
    x1 = EL_info(3);
    y1 = EL_info(4);
    L2_extreme = EL_info(5);
    nucleus_distance = 1-dot((nucleus_xy-repmat([x0,y0],size(nucleus_xy,1),1)),repmat([(x1-x0),(y1-y0)],size(nucleus_xy,1),1),2)/L2_extreme;
    
    IA = (nucleus_distance >= EL_range(1)) & (nucleus_distance <= EL_range(2));
    IP = (nucleus_distance >= EL_range(3)) & (nucleus_distance <= EL_range(4));
    if xor(mean(Inten(IP)) > mean(Inten(IA)),flip_axis)
%         nucleus_distance = 1-nucleus_distance;
        flip_EL = true;
    end
end
%%% =======================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


     