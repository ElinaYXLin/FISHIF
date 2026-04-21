function DAPI_seg3D2_test(lf_name,varargin)
%clear all
close all
tic

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameter setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
in_folder = 'stacks/';
out_folder = 'masks/';
input_name = 'matchlist.xls';
mask_name = 'mask.mat';
mask_imname = 'mask_image.fig';
% standard_record = 'Calibration/Results/standard.mat';
image_type = '*.tif';
%out_folder = 'Results/';

% imblur0 = 10;   %%% Initial bluring mask radius
imblur0 = 3;   %%% Initial bluring mask radius
blur_fit = [0.0033,-7.5];   %%% Bluring mask radius recalculation parameters
bad_ratio = 1.5;   %%% Threshold for being a merged (bad) nuclei center mask region
I_max = 128;   %%% threshold scanning step for circular mask algorithm
% % low_th = 1500;    %%% lower area limit for circular mask algorithm
% low_th = 1000;    %%% lower area limit for circular mask algorithm
% cir_th = 0.75;   %%% circularity threshold for circular mask algorithm
% cir_th = 0.6;   %%% circularity threshold for circular mask algorithm
% conv_th = 0.9;   %%% convexity threshold for circular mask algorithm
sigma_th = 100;  %%% neighboring mask threshold correlation decay length for supplemental mask recognition
merge_th = 10000;   %%% higher area limit (nuclei merge limit) for supplemental mask recognition
th_ratio = 1.0;   %%% threshold resetting ratio for mask refinement on z direction
% % % back_th = 1.5;
r_refine = 10;
% r_refine = 5;
pcore = 0.5;

% Advanced Guassian filtering    
Ex = fspecial('gaussian',100,20);
Ix = fspecial('gaussian',100,30);
% % Ex = fspecial('gaussian',100,10);
% % Ix = fspecial('gaussian',100,20);
% Ex = fspecial('gaussian',100,10);
% Ix = fspecial('gaussian',100,20);

if ~isempty(varargin) && ~isempty(varargin{1})
    channel0 = varargin{1};
else
    channel0 = [];
end

if length(varargin)>= 2 && ~isempty(varargin{2})
    cir_th = varargin{2};
else
    cir_th = 0.8;   %%% circularity threshold for circular mask algorithm
end

if length(varargin)>= 3 && ~isempty(varargin{3})
    list_I0 = varargin{3}{1};
    list_J0 = varargin{3}{2};
else
    list_I0 = [];
    list_J0 = [];
end

if length(varargin)>= 4 && ~isempty(varargin{4})
    r_erode = varargin{4};
else
    r_erode = 5;   %%% radius of erosion after initial recognition
end

if length(varargin)>= 5 && ~isempty(varargin{5})
    low_th = varargin{5}(1);    %%% lower area limit for circular mask algorithm
    high_th = varargin{5}(2);   %%% higher area limit for circular mask algorithm
else
    low_th = 500;    %%% lower area limit for circular mask algorithm
    high_th = 50000;   %%% higher area limit for circular mask algorithm
end

if length(varargin)>= 6 && ~isempty(varargin{6})
    scale0 = varargin{6};
else
    scale0 = 1;   %%% rescaling factor for image size
end

if length(varargin)>= 7 && ~isempty(varargin{7})
    z_thresh = varargin{7};
else
    z_thresh = 3;   %%% threshold of nuclear span on z
end

if length(varargin)>= 8 && ~isempty(varargin{8})
    back_th = varargin{8};
else
    back_th = 0.6;   %%% background threshold multiplier
end

if length(varargin)>= 9 && ~isempty(varargin{9})
    r_close = varargin{9};
else
    r_close = 0;   %%% radius of close after initial recognition
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
warning off

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

if isempty(list_I0)
    list_I00 = 1:N1;
else
    list_I00 = list_I0;
end

for list_I = list_I00
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

    if isempty(list_J0)
        list_J00 = 1:M1;
    elseif ~iscell(list_J0)
        list_J00 = list_J0;
    else
        list_J00 = list_J0{list_I};
    end
    
    for list_J = list_J00
        image_folder = [folder_list{list_I,1},in_folder,sub_list{list_J,3}];
        Mdim = sub_num(list_J,3);
        Nbin1 = sub_num(list_J,2);
        WGA_channel = sub_num(list_J,6);
        if isempty(channel0)
            DAPI_channel = sub_num(list_J,7);
        else
            DAPI_channel = channel0;
        end
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
        
        size0 = size(raw_im);
        if scale0 < 1
            n1 = mod(size0(1),1/scale0);
            n2 = mod(size0(2),1/scale0);
            raw_im = imresize(raw_im(1:end-n1,1:end-n2,:),scale0);
        else
            n1 = 0;
            n2 = 0;
        end
        
        if numel(DAPI_channel) == 1
            input_im = raw_im(:,:,DAPI_channel);
        else
            input_im = uint16(mean(raw_im(:,:,DAPI_channel),3));
        end
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
        
% matlabpool
% pid = parpool(length(imlist));
% % ncore = feature('numcores');
% % npar = min(length(imlist),round(ncore*pcore));
% % pid = parpool(npar);
Npid = feature('numcores');
if Npid > 6
    Npid = 6;
end
% % % pid = parpool(Npid);
        for image_I = 28:length(imlist)
            warning off
            raw_im = imread([folder_list{list_I,1},in_folder,sub_list{list_J,3},imlist(image_I).name]);
            if scale0 < 1
                raw_im = imresize(raw_im(1:end-n1,1:end-n2,:),scale0);
            end
            
            if numel(DAPI_channel) == 1
                input_im = raw_im(:,:,DAPI_channel);
            else
                input_im = uint16(mean(raw_im(:,:,DAPI_channel),3));
            end
        

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
            if ~isempty(center_xy)
                bw_center = sub2ind(size(bw_diff),center_xy(:,2),center_xy(:,1));
            else
                bw_center = [];
            end
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
                th_I = I/I_max*double(im_max)/65535;
                bw_temp = imfill(im2bw(im_temp,th_I),'holes');
                if r_close > 0
                    bw_temp = imopen(bw_temp,strel('disk',3));
                    bw_temp = imclose(bw_temp,strel('disk',r_close));
                end
%                 bw_temp = imopen(imfill(im2bw(im_temp,th_I),'holes'),strel('disk',10));
                bw_prop = regionprops(bw_temp,'Area','Perimeter');
                bw_area = [bw_prop.Area];
                bw_perim = [bw_prop.Perimeter];
                ind_true = find(bw_area >= low_th & bw_area <= high_th & 4*pi*bw_area./bw_perim.^2 >= cir_th);
                bw_true = ismember(bwlabel(bw_temp),ind_true);
                bw_im = bw_im | bw_true;
                im_temp(bw_true) = 0;
                th_2D(bw_true) = th_I;
            end
            
% % %             if r_close > 0
% % %                 bw_im = imclose(bw_im,strel('disk',r_close));
% % %             end
            bw_im = imopen(bw_im,strel('disk',10));
% % %             if r_erode > 0
% % %                 bw_im = imerode(bw_im,strel('disk',r_erode));
% % %             end
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
                
                bw_im2 = imdilate(bw_im,strel('disk',10));
                bw_center0 = bw_center(~bw_im2(bw_center));
                bw_center0_area = ismember(bwlabel(bw_diff),find((~bw_im2(bw_center))));
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
                
                bw_cut = imerode(bw_out,strel('disk',15));
% % %                 bw_cut = false(size(fil_im));
% % %     %             bw_cut([bw_center2;bw_center0]) = true;
% % %                 bw_cut(bw_center2) = true;
% % %                 bw_cut = bw_cut | (bw_center0_area & bw_out);

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
% % %                 if r_close > 0
% % %                     bw_out = imclose(bw_out,strel('disk',r_close));
% % %                 end
                bw_out = imopen(bw_out,strel('disk',10));
                
                bw_out  = ellipse_fit(bw_out,25);   %%% refine the mask using elliptical shape
                bw_out  = ellipse_fit(bw_out,25);   %%% refine the mask using elliptical shape

                im_refine = imfilter(input_im,fspecial('disk',5),'same','conv');
                Inten_back = sort(im_refine(~bw_out));
                back_thresh = back_th*Inten_back(round(length(Inten_back)*0.9));
                bw_refine0 = bwareaopen(bw_out & (im_refine >= back_thresh),low_th);
                bw_refine0 = imfill(bw_refine0,'holes');
                if r_close > 0
                    bw_refine0 = imclose(bw_refine0,strel('disk',r_close));
                else
                    bw_refine0 = imopen(bw_refine0,strel('disk',r_refine));
                end
                if r_erode > 0
                    bw_refine0 = imerode(bw_refine0,strel('disk',r_erode));
                end
                
                bw_refine = bwconvhull(bw_refine0,'objects');

                D= bwdist(~bw_refine);
                g2 = imimposemin(-D,bw_refine0);
                bw_separate = imdilate(watershed(g2) == 0,strel('disk',1));
                bw_out = bw_refine & (~ bw_separate);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Raw segmentation output: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            raw3D(:,:,image_I) = bw_out;
            th_3D(:,:,image_I) = th_2D;
            fil3D(:,:,image_I) = fil_im;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
% % % %         parfor image_I = 1:size(raw3D,3)
% % % %             raw3D(:,:,image_I) = imopen(raw3D(:,:,image_I),strel('disk',10));
% % % % %             raw3D(:,:,image_I) = imerode(raw3D(:,:,image_I),strel('disk',2));
% % % %         end
% matlabpool close
% % % delete(pid)
% toc

im000 = zeros([size(fil3D,1),size(fil3D,2),3]);
figure
    im000(:,:,3) = max(fil3D,[],3)/65535*3; 
    im000(:,:,2) = double(bwperim(max(raw3D,[],3))); 
    imshow(im000)
title([folder_list{list_I,1},out_folder,sub_list{list_J,3}])

if ~exist([folder_list{list_I,1},out_folder,sub_list{list_J,3}],'dir')
    mkdir([folder_list{list_I,1},out_folder,sub_list{list_J,3}])
end

saveas(gcf,[folder_list{list_I,1},out_folder,sub_list{list_J,3},mask_imname])
clear im000
        
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3D mask tiling and refinement: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         mask_label = label3D(raw3D);   %%% mask tiling
%         mask_stack = label3D_GPU2(raw3D);   %%% mask tiling
        mask_stack = label3D(raw3D,z_thresh);   %%% mask tiling
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
        new_label = zeros(1,max(mask_stack(:))+1,'uint16');
        new_label(sort(unique(mask_stack))+1) = [0:(length(unique(mask_stack))-1)];
        mask_stack = new_label(mask_stack+1);
        
        if scale0 < 1
            mask_stack = imresize(mask_stack,1/scale0,'nearest');
            mask_stack = cat(1,mask_stack,false(n1,size(mask_stack,2),size(mask_stack,3)));
            mask_stack = cat(2,mask_stack,false(size(mask_stack,1),n2,size(mask_stack,3)));
        end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%toc


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Segmentation result output: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %         if ~exist([folder_list{list_I,1},out_folder,sub_list{list_J,3}],'dir')
% % %             mkdir([folder_list{list_I,1},out_folder,sub_list{list_J,3}])
% % %         end
        save([folder_list{list_I,1},out_folder,sub_list{list_J,3},mask_name],'mask_stack','-v7.3')
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end
end
warning on
toc


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A function to refine the 2D mask by assuming that the nucleus has an elliptical shape
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function bw_out  = ellipse_fit(bw_in,r0)

s0 = regionprops(bw_in,'MinorAxisLength');
if ~isempty(s0)
    r1 = round(min(median([s0.MinorAxisLength]/2)-10,r0));

    s = regionprops(imerode(bw_in,strel('disk',r1)), 'Orientation', 'MajorAxisLength','MinorAxisLength', 'Eccentricity', 'Centroid');
    ind0 = find(bw_in);
    [y00,x00] = ind2sub(size(bw_in),ind0);
    Itrue = false(size(y00));

    for k = 1:length(s)
        xbar = s(k).Centroid(1);
        ybar = s(k).Centroid(2);
        a = s(k).MajorAxisLength/2+r1+2;
        b = s(k).MinorAxisLength/2+r1+2;
        theta = pi*s(k).Orientation/180;
        R = [cos(theta),sin(theta);-sin(theta),cos(theta)];

        xy11 = [x00-xbar,y00-ybar]*R;
        Itrue = Itrue | (xy11(:,1)/a).^2+(xy11(:,2)/b).^2 <= 1;
    end
    ind1 = ind0(Itrue);
    bw_out = false(size(bw_in)); 
    bw_out(ind1) = true;
else
    bw_out = bw_in;
end

