function temp_seg3D(lf_name)
%clear all
close all
tic

%% Parameter setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
in_folder = 'stacks/';
out_folder = 'masks/';
input_name = 'matchlist.xls';
mask_name = 'mask.mat';
standard_record = 'Calibration/Results/standard.mat';
image_type = '*.tif';
%out_folder = 'Results/';
I_max = 256;
low_th = 500;
high_th = 50000;
cir_th = 0.7;
th_ratio = 1.0;
thick_thresh = 3;   %%% thickness threshold
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Data loading: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global standard_data
standard_data = load(standard_record);

if isa(lf_name,'char') && strcmp(lf_name(end-3:end),'.xls')
    list_name = lf_name;
    [num_list, folder_list] = xlsread(list_name);
    folder_list = folder_list(strcmpi('T',folder_list(:,6)),:);
elseif isa(lf_name,'cell')
    folder_list = lf_name;
else
    error('Incorrect input!')
end
[N1,N2] = size(folder_list);

for list_I = 1:N1
    [sub_num, sub_list] = xlsread([folder_list{list_I,1},in_folder,input_name]);
    [M1,M2] = size(sub_list);
    channel_name = eval(folder_list{list_I,5});
    for list_J = 1:M1
        image_folder = [folder_list{list_I,1},in_folder,sub_list{list_J,3}];
        imlist = [folder_list{list_I,1},out_folder,sub_list{list_J,3},mask_name]; %%% get the raw mask array
        load(imlist);
%% 3D mask refinement: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         mask_stack = label3D(raw3D);   %%% mask stitching and alignment
    %%% mask rethresholding: %%% ==========================================
% %         mask_stack = zeros(size(raw3D));
% %         parfor I_nuclei = 1:max(mask_label(:))
% %             temp_label = mask_label == I_nuclei;
% %             mask_stack(temp_label) = I_nuclei*im2bw(fil3D(temp_label),max(th_3D(temp_label)));
% %         end
%         label_prop = regionprops(mask_label,th_3D,'MaxIntensity');
%         th_max = [1,[label_prop.MaxIntensity]];
%         mask_stack = mask_label.*(fil3D > th_ratio*65535*th_max(mask_label+1));
%         mask_stack = mask_label;
%         clear mask_label th_max
    %%% ===================================================================
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%toc


%% Remove the fake regions that are too thin: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        thickness_region = zeros(1,max(mask_stack(:)));
        for I_layer = 1:size(mask_stack,3)
            region_plus = setdiff(unique(mask_stack(:,:,I_layer)),[0]);
            thickness_region(region_plus) = thickness_region(region_plus)+1;
        end
        thin_region = find(thickness_region <= thick_thresh);
        mask_stack(ismember(mask_stack,thin_region)) = 0;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Resign the label #: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        new_label = zeros(1,max(mask_stack(:))+1);
        new_label(sort(unique(mask_stack))+1) = [0:(length(unique(mask_stack))-1)];
        mask_stack = new_label(mask_stack+1);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%toc
%% Fine segmentation output: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         if ~exist([folder_list{list_I,1},out_folder,sub_list{list_J,3}],'dir')
%             mkdir([folder_list{list_I,1},out_folder,sub_list{list_J,3}])
%         end
%matlabpool
%         for image_I = 1:length(imlist)
%             mask_stack(:,:,image_I) = imfill(mask_stack(:,:,image_I),'holes');
%             mask_stack(:,:,image_I) = mask_stack(:,:,image_I).*bwareaopen(logical(mask_stack(:,:,image_I)), low_th);
%             raw_im = imread([folder_list{list_I,1},in_folder,sub_list{list_J,3},imlist(image_I).name]);
%             im_out = uint16(zeros([size_im([1:2]),3]));
%             im_out(:,:,1) = raw_im(:,:,WGA_channel);
%             im_out(:,:,2) = uint16(bwperim(mask_stack(:,:,image_I))*65535);
%             im_out(:,:,3) = fil3D(:,:,image_I);
%             imwrite(im_out,[folder_list{list_I,1},out_folder,sub_list{list_J,3},imlist(image_I).name])
%         end
%matlabpool close
        save(imlist,'mask_stack','-v7.3')
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end
end

toc



% function mask_label = label3D(raw3D)
%     
% end
