function [nucleus_expand_profile,r_edge] = DAPI_profile3D_expand(mask_stack,signal_stack,varargin)


%% Nucleus DAPI expand profile calculation: %%%%%%%%%%%%%%%%%%%%%%
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Parameter setting and initiation: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(varargin) || isempty(varargin{1})
    r_edge00 = 15;
    r_edge0 = (-r_edge00-1)*2:2:r_edge00*2;
else
    r_edge0 = varargin{1};
end
r_edge0 = sort(r_edge0);
r_edge = r_edge0(2:end);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Nucleus region recognition/processing: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_nucleus = max(mask_stack(:));
nucleus_expand_profile = zeros(length(r_edge),N_nucleus,size(mask_stack,3));

Npid = feature('numcores');
if Npid > 30
    Npid = 30;
end
pid = parpool(Npid);

for I_layer = 1:size(mask_stack,3)
    temp_mean = zeros(length(r_edge0),N_nucleus);
    temp_area = zeros(length(r_edge0),N_nucleus);
    
    for ii = 1:length(r_edge0)
        if r_edge0(ii) < 0
            sg_prop = regionprops(double(imerode(logical(mask_stack(:,:,I_layer)), strel('disk',-r_edge0(ii)))).*double(mask_stack(:,:,I_layer)),double(signal_stack(:,:,I_layer)),'MeanIntensity','Area');
        elseif r_edge0(ii) == 0
            sg_prop = regionprops(double(mask_stack(:,:,I_layer)),double(signal_stack(:,:,I_layer)),'MeanIntensity','Area');
        else
            sg_prop = regionprops(double(dilate3D(mask_stack(:,:,I_layer),r_edge0(ii))),double(signal_stack(:,:,I_layer)),'MeanIntensity','Area');
        end
        
        temp_mean(ii,1:length(sg_prop)) = [sg_prop.MeanIntensity];
        temp_area(ii,1:length(sg_prop)) = [sg_prop.Area];
    end
    
    temp_mean(isnan(temp_mean)) = 0;
    temp_area(isnan(temp_area)) = 0;
    
    nucleus_expand_profile(:,:,I_layer) = diff(temp_mean.*temp_area)./diff(temp_area);
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
delete(pid)


