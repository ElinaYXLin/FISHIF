function flip_EL = EL_orientation2(EL_info,mask_stack,file_name,protein_RNA_mismatch,flip_axis,varargin)


%% Determine the orientation of EL axis using foci intensity: %%%%
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Parameter setting and initiation: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(varargin) && ~isempty(varargin{1})
    EL_range = varargin{1};
else
    EL_range = [0.15,0.4,0.6,0.85];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Nucleus and foci processing: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[foci_list,~,~] = xlsread(file_name);
Inten_foci = prod(foci_list(:,1:3),2);
i_foci = min(max(round(foci_list(:,6)),1),size(mask_stack,1));
j_foci = min(max(round(foci_list(:,7)),1),size(mask_stack,2));
k_foci = min(max(round(foci_list(:,8))+protein_RNA_mismatch,1),size(mask_stack,3));
ind_foci = sub2ind(size(mask_stack),i_foci,j_foci,k_foci);
N_foci = mask_stack(ind_foci);

sg_prop = regionprops(mask_stack,'Centroid');
centr = cell2mat({sg_prop.Centroid}');

Inten = zeros(size(centr,1),1);
for ii = 1:size(centr,1)
    Inten(ii) = sum(Inten_foci(N_foci == ii));
end
Inten(isnan(Inten)) = 0;

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


     