function mask1 = mask_adjust(mask0,nucleus_expand_profile,r_edge,expand_th)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A function to adjust the size of each nucleus (erosion/dilation) based on the intensity profile of boundary pixels %%
%% mask1 (output): Adjusted mask for output
%% mask0 (input): original mask for adjustment
%% nucleus_expand_profile (input): nucleus boundary pixel intensity profile for different radius adjustments
%% r_edge (input): radius adjustment values
%% expand_th (input): threshold for relative boundary pixel intensity
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Calculate the relative intensity of boundary pixels for different radius adjustments: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nucleus_expand_max = repmat(max(nucleus_expand_profile),size(nucleus_expand_profile,1),1,1);
nucleus_expand_min = repmat(min(nucleus_expand_profile),size(nucleus_expand_profile,1),1,1);
nucleus_expand_ratio = (nucleus_expand_profile-nucleus_expand_min)./(nucleus_expand_max-nucleus_expand_min);

mask1 = zeros(size(mask0));
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Determine the radius expansion/shrinking: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Npid = feature('numcores');
if Npid > 30
    Npid = 30;
end
pid = parpool(Npid);

parfor I_layer = 1:size(mask0,3)
    ratio0 = nucleus_expand_ratio(:,:,I_layer);
    mask00 = mask0(:,:,I_layer);
    mask11 = zeros(size(mask00));
    
    [~,I0] = max(flipud(ratio0));
    I_true0 = zeros(size(ratio0));
    I_true0(I0+(0:(size(ratio0,2)-1))*size(ratio0,1)) = 1;
    I_true0 = logical(cumsum(flipud(I_true0)));
    
    ratio1 = ratio0.*I_true0; ratio1(~I_true0) = 1;
    [~,I0] = min(ratio1);
    I_true1 = zeros(size(ratio0));
    I_true1(I0+(0:(size(ratio0,2)-1))*size(ratio0,1)) = 1;
    I_true1 = ~logical(cumsum(I_true1));
    
    ratio2 = ratio1.*I_true1;
    [~,I_th] = max(cumsum(ratio2 > expand_th));
    I_nan = all(isnan(ratio0)); 
    I_th(I_nan) = nan;
    N_th = sort(unique(I_th(~I_nan)),'descend');
    
    %%% Radius expansion/shrinking:
    for ii = 1:length(N_th)
        r0 = r_edge(N_th(ii));
        mask_temp = dilate3D(mask00,r0);
        
        I_adjust = ismember(mask_temp,find(I_th == N_th(ii)));
        mask11(I_adjust) = mask_temp(I_adjust);
    end
    
    mask1(:,:,I_layer) = mask11;
end
delete(pid)



