function mask_label = label3D_GPU(raw3D)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A function to label the stitched 2D segmentation mask %%%%%%%%%%%%%%%%%%
%% raw3D: stitched 2D mask.                              %%%%%%%%%%%%%%%%%%
%% mask_label: labeled 3D mask.                          %%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% raw3D = gpuArray(raw3D);
mask_label = zeros(size(raw3D),'uint16');
% max_label = 0;   %%% maximal label #
if size(raw3D,3) <= 1
   thick_thresh = 1;   %%% thickness threshold
else
    thick_thresh = 3;   %%% thickness threshold
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% For the first layer: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mask_label(:,:,1) = bwlabel(raw3D(:,:,1));
max_label = max(mask_label(:));
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% 2D mask stitching: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for I_layer = 1:(size(raw3D,3)-1)
    temp_label1 = gpuArray(mask_label(:,:,I_layer));
    temp_label2 = gpuArray(zeros(size(raw3D(:,:,1)),'uint16'));
    temp1 = gpuArray(raw3D(:,:,I_layer));
    temp2 = gpuArray(raw3D(:,:,I_layer+1));
    bwlabel1 = uint16(bwlabel(temp1));
    bwlabel2 = uint16(bwlabel(temp2));
    prop1 = regionprops(temp1,'Centroid');
    prop2 = regionprops(temp2,'Centroid');
    
    if ~isempty(prop1)
        xy_center0 = round(cell2mat({prop1.Centroid}'));   %%% center coordinates matrix for layer I
        ind_center0 = sub2ind(size(bwlabel1),xy_center0(:,2),xy_center0(:,1));   %%% center linear indices matrix for layer I
        bw2temp_prop = regionprops(bwlabel1,temp_label1,'MaxIntensity');   %%% convert bwlabel1 labels to temp_label1 labels
        bw2temp = [bw2temp_prop.MaxIntensity];
        
        for I_area = 1:length(prop2)

            x_center = round(prop2(I_area).Centroid(1));
            y_center = round(prop2(I_area).Centroid(2));
            if temp_label1(y_center,x_center)
                temp_label2(bwlabel2 == I_area) = temp_label1(y_center,x_center);
            end

            ind_list = find(bwlabel2(ind_center0) == I_area);
            current_ind = temp_label2(find(bwlabel2 == I_area,1));
            if (~current_ind) && (~isempty(ind_list))
                current_ind = min(bw2temp(ind_list));
                temp_label2(bwlabel2 == I_area) = current_ind;
            elseif (~current_ind) && isempty(ind_list)
                max_label = max_label+1;
                current_ind = max_label;
                temp_label2(bwlabel2 == I_area) = current_ind;
            end
            ind_label = setdiff(bw2temp(ind_list),current_ind);
            if ~isempty(ind_label)
                temp_label1(ismember(temp_label1,ind_label)) = current_ind;
                mask_label(ismember(mask_label,ind_label)) = current_ind;
                bw2temp(ind_list) = current_ind;
            end
        end
        mask_label(:,:,I_layer+1) = gather(temp_label2);
    else
        mask_label(:,:,I_layer+1) = gather((bwlabel2+max_label)).*uint16(raw3D(:,:,I_layer+1));
        max_label = max(mask_label(:));
    end
%     figure; imshow(mask_label(:,:,1));
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Remove the fake regions that are too thin: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
thickness_region = zeros(1,max(mask_label(:)));
for I_layer = 1:size(mask_label,3)
    region_plus = setdiff(unique(mask_label(:,:,I_layer)),[0]);
    thickness_region(region_plus) = thickness_region(region_plus)+1;
end
thin_region = find(thickness_region < thick_thresh);
mask_label(ismember(mask_label,thin_region)) = 0;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Resign the label #: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
new_label = zeros(1,max_label+1,'uint16');
new_label(sort(unique(mask_label))+1) = [0:(length(unique(mask_label))-1)];
mask_label = gather(new_label(mask_label+1));
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
