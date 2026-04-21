function DAPI_seg3D2_single(lf_name)
%clear all
% close all
tic

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameter setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
in_folder = 'stacks/';
out_folder = 'masks/';
input_name = 'matchlist.xls';
mask_name = 'mask.mat';
standard_record = 'Calibration/Results/standard.mat';
image_type = '*.tif';
%out_folder = 'Results/';

% imblur0 = 10;   %%% Initial bluring mask radius
imblur0 = 3;   %%% Initial bluring mask radius
blur_fit = [0.0033,-7.5];   %%% Bluring mask radius recalculation parameters
bad_ratio = 1.5;   %%% Threshold for being a merged (bad) nuclei center mask region
I_max = 256;   %%% threshold scanning step for circular mask algorithm
% % low_th = 1500;    %%% lower area limit for circular mask algorithm
low_th = 500;    %%% lower area limit for circular mask algorithm
% low_th = 100;    %%% lower area limit for circular mask algorithm
high_th = 10000;   %%% higher area limit for circular mask algorithm
cir_th = 0.70;   %%% circularity threshold for circular mask algorithm
% conv_th = 0.9;   %%% convexity threshold for circular mask algorithm
sigma_th = 100;  %%% neighboring mask threshold correlation decay length for supplemental mask recognition
merge_th = 10000;   %%% higher area limit (nuclei merge limit) for supplemental mask recognition
th_ratio = 1.0;   %%% threshold resetting ratio for mask refinement on z direction
% back_th = 1;
back_th = 0.6;
% r_refine = 10;
r_refine = 5;

% Advanced Guassian filtering    
% Ex = fspecial('gaussian',100,20);
% Ix = fspecial('gaussian',100,30);
Ex = fspecial('gaussian',100,10);
Ix = fspecial('gaussian',100,20);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Data loading: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% global standard_data
% standard_data = load(standard_record);

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
    I_cut = strfind(folder_list{list_I,1},in_folder);
    if isempty(I_cut)
        [sub_num, sub_list] = xlsread([folder_list{list_I,1},in_folder,input_name]);
    else
        [sub_num, sub_list] = xlsread([folder_list{list_I,1}(1:I_cut-1),in_folder,input_name]);
        I_temp = strmatch(folder_list{list_I,1}(I_cut+length(in_folder):end-1),sub_list(:,3));
        sub_list = sub_list(I_temp,:);
        folder_list{list_I,1} = folder_list{list_I,1}(1:I_cut-1);
    end
    [M1,M2] = size(sub_list);
%     channel_name = eval(folder_list{list_I,5});
    for list_J = 1%:M1
        image_folder = [folder_list{list_I,1},in_folder,sub_list{list_J,3}];
        Mdim = sub_num(list_J,3);
        Nbin1 = sub_num(list_J,2);
        WGA_channel = sub_num(list_J,6);
        DAPI_channel = sub_num(list_J,7);
        RNA_channel = sub_num(list_J,10);
        protein_channel = sub_num(list_J,8);
        resolution = sub_num(list_J,9);
        imblur = imblur0;
        
        imlist = dir([folder_list{list_I,1},in_folder,sub_list{list_J,3},image_type]); %%% get the image list from image folder1
%         max_temp = zeros(1,length(imlist));
% matlabpool
%         parfor image_I = 1:length(imlist)
%             raw_im = imread([folder_list{list_I,1},in_folder,sub_list{list_J,3},imlist(image_I).name]);
%             max_temp(image_I) = max(max(raw_im(:,:,DAPI_channel)));
%         end
% matlabpool close
%         clear raw_im
%         im_max = max(max_temp);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%% segmentation test and mean nuclear size measurement: %%%%%%%%%%%%%%%%%%%
        raw_im = imread([folder_list{list_I,1},in_folder,sub_list{list_J,3},imlist(round(length(imlist)/2)).name]);
        input_im = raw_im(:,:,DAPI_channel);
        size_im = size(input_im);
        H = fspecial('disk',imblur); % Filter Kernel  
        fil_im = imfilter(input_im,H,'same','conv');
        im_max = max(fil_im(:));
        im_temp = fil_im;
        bw_im = false(size(fil_im));

        for I = 1:(I_max-1)
            th_I = I/I_max*double(im_max)/65535;
            bw_temp = imfill(im2bw(im_temp,th_I),'holes');
            bw_prop = regionprops(bw_temp,'Area','Perimeter');%,'ConvexArea');
            bw_area = [bw_prop.Area];
%             bw_conv_area = [bw_prop.ConvexArea];
            bw_perim = [bw_prop.Perimeter];
            ind_true = find(bw_area >= low_th & bw_area <= high_th & 4*pi*bw_area./bw_perim.^2 >= cir_th);% & bw_area./bw_conv_area >= conv_th);
            bw_true = ismember(bwlabel(bw_temp),ind_true);
            bw_im = bw_im | bw_true;
            im_temp(bw_true) = 0;
        end
        bw_im = imopen(bw_im,strel('disk',10));
        temp_prop = regionprops(bw_im,'Area');
        mean_area = mean([temp_prop.Area]);
        imblur = max(3,round(mean_area*blur_fit(1)+blur_fit(2)));
%         r_mean = max(ceil(sqrt(mean_area/pi)),28);
% 
%         Ex = fspecial('gaussian',100,round(r_mean*0.4));
%         Ix = fspecial('gaussian',100,round(r_mean*0.6));
        
        clear raw_im input_im fil_im im_temp bw_im bw_temp bw_prop bw_area bw_perim bw_true temp_prop
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        


        raw3D = false([size_im([1:2]),length(imlist)]); %%% Raw 3D mask
        th_3D = zeros([size_im([1:2]),length(imlist)]); %%% Raw 3D threshold
        fil3D = zeros([size_im([1:2]),length(imlist)]); %%% Smoothed image
        
% % matlabpool
% %         parfor image_I = 1:length(imlist)
        input_im = int16(zeros(size(raw3D)));
        for image_I = 1:length(imlist)
            raw_im = imread([folder_list{list_I,1},in_folder,sub_list{list_J,3},imlist(image_I).name]);
            input_im(:,:,image_I) = raw_im(:,:,DAPI_channel);
        end
        input_im = max(input_im,[],3);
%         input_im = double(max(input_im,[],3));
%         input_im = int16(input_im/max(input_im(:))*65535);
        

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%% Targeting nuclear centers: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            H = fspecial('disk',imblur); % Filter Kernel       
            I_blur = imfilter(input_im,H,0); %Apply Filter
            I_blur =  adapthisteq(I_blur); %Step 1: enhances the contrast of the grayscale image   

            %%% Gaussian difference filtering and segmentation  
            outE = imfilter(single(I_blur),Ex,'replicate'); 
            outI = imfilter(single(I_blur),Ix,'replicate'); 
            outims = outE - outI;  
            W = watershed(max(outims(:))-outims);

            outims(outims < 0) = 0;
            outims = outims/max(outims(:));

            bw_diff = im2bw(outims,1.2*graythresh(outims)); 
            temp2_props = regionprops(bw_diff,'Area');
            mean_area_diff = mean([temp2_props.Area]);
            bad_nu_I = find([temp2_props.Area] > bad_ratio*mean_area_diff);
            bad_bw0 = ismember(bwlabel(bw_diff),bad_nu_I);
            bad_bw = imerode(bad_bw0 & (W ~= 0),strel('disk',5));
            bw_diff(bad_bw0) = bad_bw(bad_bw0);

            %%% Recognize the center points of nuclei
            diff_props = regionprops(bw_diff,'Centroid');
            center_xy = round(cell2mat({diff_props.Centroid}'));
            bw_center = sub2ind(size(bw_diff),center_xy(:,2),center_xy(:,1));
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%% Circular segmentation: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            H = fspecial('disk',imblur); % Filter Kernel  
            fil_im = imfilter(input_im,H,'same','conv');
            th_2D = zeros(size_im([1:2])); %%% Raw 2D threshold
            im_max = max(fil_im(:));
            im_temp = fil_im;
            bw_im = false(size(fil_im));
%%% Cicle search
            for I = 1:(I_max-1)
                th_I = I/I_max;%*double(im_max)/65535;
                bw_temp = imfill(im2bw(im_temp,th_I),'holes');
                bw_prop = regionprops(bw_temp,'Area','Perimeter');
                bw_area = [bw_prop.Area];
                bw_perim = [bw_prop.Perimeter];
                ind_true = find(bw_area >= low_th & bw_area <= high_th & 4*pi*bw_area./bw_perim.^2 >= cir_th);
                bw_true = ismember(bwlabel(bw_temp),ind_true);
                bw_im = bw_im | bw_true;
                im_temp(bw_true) = 0;
                th_2D(bw_true) = th_I;
            end
%             bw_im = imopen(bw_im,strel('disk',10));
            bw_out = bw_im;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%% Supplement segmentation based on center recognition: %%%%%%%%%%%%%%%%%%%
            bw_props = regionprops(bw_im,'Centroid');
            if ~isempty(bw_props)
                center_xy2 = round(cell2mat({bw_props.Centroid}'));
                bw_center2 = sub2ind(size(bw_im),center_xy2(:,2),center_xy2(:,1));
                th_center2 = th_2D(bw_center2);

                bw_center0 = bw_center(~bw_im(bw_center));
                bw_center0_area = ismember(bwlabel(bw_diff),find((~bw_im(bw_center))));
                [center_y0,center_x0] = ind2sub(size(bw_diff),bw_center0);
                center_xy0 = [center_x0,center_y0];

                for I_center = 1:size(bw_center0,1)
                    center_dist2 = (center_xy2(:,1)-center_xy0(I_center,1)).^2+(center_xy2(:,2)-center_xy0(I_center,2)).^2;
                    th_estimate = sum(th_center2.*exp(-center_dist2/2/sigma_th^2))/sum(exp(-center_dist2/2/sigma_th^2));
                    bw_th = imfill(bwareaclose(im2bw(im_temp,th_estimate),merge_th),'holes');
%                     bw_th = im2bw(im_temp,th_estimate);
                    bw_th_label = bwlabel(bw_th);
                    bw_I = (bw_th_label == bw_th_label(bw_center0(I_center))) & bw_th;
                    bw_out = bw_out | bw_I;
                    th_2D(bw_I) = max(th_2D(bw_I),th_estimate);
                end

                bw_cut = false(size(fil_im));
    %             bw_cut([bw_center2;bw_center0]) = true;
                bw_cut(bw_center2) = true;
                bw_cut = bw_cut | (bw_center0_area & bw_out);

                temp_label = bwlabel(bw_out);
                temp_sel = ismember(temp_label,temp_label(bw_cut));
                D= bwdist(~bw_out);
                g2 = imimposemin(-D,bw_cut);
                bw_cut = imdilate(watershed(g2) == 0,strel('disk',1)) & (temp_sel);
                bw_out = bw_out & (~ bw_cut);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%% Mask area refinement and convexation: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                im_refine = imfilter(input_im,fspecial('disk',5),'same','conv');
                Inten_back = sort(im_refine(~bw_out));
                back_thresh = back_th*Inten_back(round(length(Inten_back)*0.9));
                bw_refine0 = bwareaopen(bw_out & (im_refine >= back_thresh),low_th);
                bw_refine0 = imfill(bw_refine0,'holes');
                bw_refine = imopen(bw_refine0,strel('disk',r_refine));
%                 bw_refine = bwconvhull(bw_refine0,'objects');

                D= bwdist(~bw_refine);
                g2 = imimposemin(-D,bw_refine0);
                bw_separate = imdilate(watershed(g2) == 0,strel('disk',1));
                bw_out = bw_refine & (~ bw_separate);
                bw_out = bwconvhull(bw_out,'objects');
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Raw segmentation output: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for image_I = 1:length(imlist)
            raw3D(:,:,image_I) = bw_out;
            th_3D(:,:,image_I) = th_2D;
            fil3D(:,:,image_I) = fil_im;
        end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %         end
% % % %         parfor image_I = 1:size(raw3D,3)
% % % %             raw3D(:,:,image_I) = imopen(raw3D(:,:,image_I),strel('disk',10));
% % % % %             raw3D(:,:,image_I) = imerode(raw3D(:,:,image_I),strel('disk',2));
% % % %         end
% % matlabpool close
        
        
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3D mask tiling and refinement: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         mask_label = label3D(raw3D);   %%% mask tiling
        mask_stack = label3D(raw3D);   %%% mask tiling
        clear raw3D
    %%% mask rethresholding to refine z distribution: %%% =================
%         mask_stack = zeros(size(raw3D));
%         parfor I_nuclei = 1:max(mask_label(:))
%             temp_label = mask_label == I_nuclei;
%             mask_stack(temp_label) = I_nuclei*im2bw(fil3D(temp_label),max(th_3D(temp_label)));
%         end
%         label_prop = regionprops(mask_label,th_3D,'MaxIntensity');
%         th_max = [1,[label_prop.MaxIntensity]];
%         mask_stack = mask_label.*(fil3D > th_ratio*65535*th_max(mask_label+1));
%         mask_stack = mask_label;
%         clear mask_label th_max
    %%% ===================================================================
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%toc


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Resign the label #: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        new_label = zeros(1,max(mask_stack(:))+1);
        new_label(sort(unique(mask_stack))+1) = [0:(length(unique(mask_stack))-1)];
        mask_stack = new_label(mask_stack+1);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%toc


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Segmentation result output: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ~exist([folder_list{list_I,1},out_folder,sub_list{list_J,3}],'dir')
            mkdir([folder_list{list_I,1},out_folder,sub_list{list_J,3}])
        end
        save([folder_list{list_I,1},out_folder,sub_list{list_J,3},mask_name],'mask_stack','-v7.3')
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end
end

toc

