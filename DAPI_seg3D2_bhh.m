function DAPI_seg3D2_bhh(lf_name,varargin)
%clear all
close all
tic

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameter setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
in_folder = 'stacks/';
result_folder = 'Results/';
out_folder = 'masks/';
input_name = 'matchlist.xls';
mask_name = 'mask.mat';
mask_imname = 'mask_image.fig';
% standard_record = 'Calibration/Results/standard.mat';
image_type = '*.tif';
%out_folder = 'Results/';
mat_tail = '.mat';

% imblur0 = 10;   %%% Initial bluring mask radius
imblur0 = 5;   %%% Initial bluring mask radius
imblur2 = 12;   %%% Initial bluring mask radius
blur_fit = [0.0033,-7.5];   %%% Bluring mask radius recalculation parameters
bad_ratio = 1.5;   %%% Threshold for being a merged (bad) nuclei center mask region
I_max = 128;   %%% threshold scanning step for circular mask algorithm
I_max0 = 128;   %%% threshold max
% % low_th = 1500;    %%% lower area limit for circular mask algorithm
% low_th = 1000;    %%% lower area limit for circular mask algorithm
% cir_th = 0.75;   %%% circularity threshold for circular mask algorithm
% cir_th = 0.6;   %%% circularity threshold for circular mask algorithm
% conv_th = 0.9;   %%% convexity threshold for circular mask algorithm
sigma_th = 100;  %%% neighboring mask threshold correlation decay length for supplemental mask recognition
th_ratio = 1.0;   %%% threshold resetting ratio for mask refinement on z direction
% % % back_th = 1.5;
pcore = 0.5;

% Advanced Guassian filtering    
Ex = fspecial('gaussian',100,10);
Ix = fspecial('gaussian',100,20);
% % % % Ex = fspecial('gaussian',100,20);
% % % % Ix = fspecial('gaussian',100,30);
% % % Ex = fspecial('gaussian',100,15);
% % % Ix = fspecial('gaussian',100,30);
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
% merge_th = 10000;   %%% higher area limit (nuclei merge limit) for supplemental mask recognition
% merge_th = 100000;   %%% higher area limit (nuclei merge limit) for supplemental mask recognition
merge_th = high_th;   %%% higher area limit (nuclei merge limit) for supplemental mask recognition

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
if length(r_close) < 2
    r_close(2:3) = r_close(1);
elseif length(r_close) < 3
    r_close(3) = r_close(2);
end

if length(varargin)>= 10 && ~isempty(varargin{10})
    r_refine = varargin{10};
else
    r_refine = 10;   %%% radius of open after initial recognition
%     r_refine = 5;
end

%%% Methods to locate nuclei and determine the local threshold:     0. Layer-wise detection using criterion of circles; ...
                                                                 ...1. Location using maximum projected previous mask (with mean intensity threshold) and detection using criterion of circles; ...
                                                                 ...2. Location using maximum projected image (with mean intensity threshold) and detection using criterion of circles; ...
                                                                 ...3. Location using previous mask (without mean intensity threshold) and detection using criterion of circles applied on maximum projected image; ...
                                                                 ...4. Location using maximum projected image (with mean intensity threshold) and detection using criterion of circles applied on maximum projected image; ...
                                                                 ...5. Location using maximum projected image (with mean intensity threshold) and detection using a global threshold from maximum projected image; ...
if length(varargin)>= 11 && ~isempty(varargin{11}) 
    t_reuse = varargin{11};
else
    t_reuse = 0;
end

if t_reuse == 1 || t_reuse == 3
    reuse_file = 'quick_mask.mat';
    reuse_diff = 'bw_applied3D_save';
    mask_name = 'quick_mask_reuse.mat';
% % elseif t_reuse == 2
% %     reuse_file = 'quick_mask.mat';
% %     reuse_diff = 'bw_applied3D_save';
% %     mask_name = 'quick_mask_test.mat';
end

r00 = 0.4;
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
        
        %%% Create an ellipitcal mask for imerode
% % %         if exist([folder_list{list_I,1},result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mat_tail])
% % %             load([folder_list{list_I,1},result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mat_tail],'em_mask');
% % %             EL_info = get_EL(em_mask);   %%% get EL information (extreme points, EL length) from the embryo mask
% % %             theta0 = angle((EL_info(3)-EL_info(1))-1i*(EL_info(4)-EL_info(2)))/pi*180;
% % %         else
% % %             theta0 = 0;
% % %         end
% % %         mask_elp0 =  ellipse_mask(10,3.3,theta0);
% % %         mask_elp =  ellipse_mask(18,6,theta0);

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
%         if exist([folder_list{list_I,1},result_folder,sub_list{list_J,3}(1:end-1),mat_tail])
%             S = load([folder_list{list_I,1},result_folder,sub_list{list_J,3}(1:end-1),mat_tail],'max_image');
%             raw_im = S.max_image; clear S
%         else
            raw_im = imread([folder_list{list_I,1},in_folder,sub_list{list_J,3},imlist(round(length(imlist)/2)).name]);
%         end
        
        size0 = size(raw_im);
        if scale0 > 1
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
        
        fil_im2 = imfilter(input_im,fspecial('disk',imblur2),'same','conv');
        
        im_max = max(fil_im(:));
        im_temp = fil_im;
        bw_im = false(size(fil_im));

        for I = 1:(I_max-1)
            th_I = I/I_max0*double(im_max)/65535;
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
        imblur = max(3,round((mean_area*blur_fit(1)+blur_fit(2))*1.5));
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
% % %         [X0,Y0] = meshgrid(1:size_im(2),1:size_im(1));
        old_diff = raw3D;
        diff0 = old_diff(:,:,1);
            
        th_2D0 = zeros(size_im(1:2));
        bw_im0 = false(size_im(1:2));
        I_blur = zeros(size_im(1:2));
        im_max0 = 0;
        th_ratio = 1.2;
        MidI = ones(size_im(1:2));

        if t_reuse == 1 || t_reuse == 3
            S = load([folder_list{list_I,1},out_folder,sub_list{list_J,3},reuse_file],reuse_diff);
            old_diff = eval(['S.',reuse_diff]); clear S
            if scale0 > 1
                old_diff = imresize(old_diff(1:end-n1,1:end-n2,:),scale0);
            end
            
            diff0 = max(old_diff,[],3);
            
            D= -bwdist(~diff0);
            bw_cut = imextendedmin(D,2);
            g2 = imimposemin(D,bw_cut);
            bw_cut = imdilate(watershed(g2) == 0,strel('disk',1));% & (temp_sel);
            diff0 = diff0 & (~ bw_cut);
        end
            
        if t_reuse == 2 || t_reuse == 3 || t_reuse == 4 || t_reuse == 5
            S = load([folder_list{list_I,1},result_folder,sub_list{list_J,3}(1:end-1),mat_tail],'max_image');
            raw_im = S.max_image; clear S
            if scale0 > 1
                raw_im = imresize(raw_im(1:end-n1,1:end-n2,:),scale0);
            end

            if numel(DAPI_channel) == 1
                input_im = raw_im(:,:,DAPI_channel);
            else
                input_im = uint16(mean(raw_im(:,:,DAPI_channel),3));
            end

            H = fspecial('disk',imblur); % Filter Kernel  
        end
            
        if t_reuse == 2 || t_reuse == 4 || t_reuse == 5
            I_blur = imfilter(input_im,H,0); %Apply Filter
            I_blur =  adapthisteq(I_blur); %Step 1: enhances the contrast of the grayscale image   

            %%% Gaussian difference filtering and segmentation  
            outE = imfilter(single(I_blur),Ex,'replicate'); 
            outI = imfilter(single(I_blur),Ix,'replicate'); 
            outims = outE - outI;  
            W = watershed(max(outims(:))-outims);

            outims(outims < 0) = 0;
            
            MidI = imfilter(single(I_blur),double(getnhood(strel('disk',500))),'replicate'); 
            MidI = MidI/max(MidI(:));
            outims = outims./MidI;
            
            outims = outims/max(outims(:));

            bw_diff = im2bw(outims,1*graythresh(outims)); 
% % % % %             temp2_props = regionprops(bw_diff,'Area');
% % % % %             mean_area_diff = mean([temp2_props.Area]);
% % % % %             bad_nu_I = find([temp2_props.Area] > bad_ratio*mean_area_diff);
% % % % %             bad_bw0 = ismember(bwlabel(bw_diff),bad_nu_I);
% % % % %             bad_bw = imerode(bad_bw0 & (W ~= 0),strel('disk',5));
% % % % %             bw_diff(bad_bw0) = bad_bw(bad_bw0);
% % % % %             bw_diff = imclose(bw_diff,strel('disk',10));
% % % % %             bw_diff = imfill(bw_diff,'holes');
            bw_diff = bwconvhull(bw_diff,'objects');
            bw_diff = imopen(bw_diff,strel('disk',10));
            
            old_diff = repmat(bw_diff,[1,1,length(imlist)]);
            diff0 = bw_diff;
        end
        
        if t_reuse == 3 || t_reuse == 4
            fil_im = uint16(double(imfilter(input_im,H,'same','conv'))./double(MidI));
            im_max0 = max(fil_im(:));
            [th_2D0,bw_im0] = circular_th(fil_im,I_max,I_max0,im_max0,r_close,low_th,high_th,cir_th);
            th_ratio = 0.7;
% %             th_ratio = 1;
% %             th_2D0(:) = 0.06;
        end
        
        if t_reuse == 5
            fil_im = uint16(double(imfilter(input_im,H,'same','conv'))./double(MidI));
            bw_im0 = true(size_im(1:2));
            bw0 = im2bw(fil_im,graythresh(fil_im));
            th_2D0 = double(min(fil_im(bw0)))/65535*ones(size_im(1:2));
            th_ratio = 0.4;
        end
        
        prop0 = regionprops(diff0,I_blur,'MeanIntensity');
        Inten_diff0 = [prop0.MeanIntensity];
        label_diff0 = bwlabel(diff0);
        
% matlabpool
% pid = parpool(length(imlist));
% % ncore = feature('numcores');
% % npar = min(length(imlist),round(ncore*pcore));
% % pid = parpool(npar);
Npid = feature('numcores');
if Npid > 18
    Npid = 18;
end
pid = parpool(Npid);
        parfor image_I = 1:length(imlist)
            warning off
            raw_im = imread([folder_list{list_I,1},in_folder,sub_list{list_J,3},imlist(image_I).name]);
            if scale0 > 1
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

            if t_reuse
                if t_reuse == 1 || t_reuse == 3
                    bw_diff = old_diff(:,:,image_I);
                else
                    bw_diff = diff0;
                    prop0 = regionprops(bw_diff,I_blur,'MeanIntensity');
                    Inten_diff = [prop0.MeanIntensity];
                    I_true = find(Inten_diff >= r00*Inten_diff0);
                    bw_diff = ismember(label_diff0,I_true);
                end
            else
                %%% Gaussian difference filtering and segmentation  
                outE = imfilter(single(I_blur),Ex,'replicate'); 
                outI = imfilter(single(I_blur),Ix,'replicate'); 
                outims = outE - outI;  
                W = watershed(max(outims(:))-outims);

                outims(outims < 0) = 0;
                outims = outims/max(outims(:));

                bw_diff = im2bw(outims,1*graythresh(outims)); 
                temp2_props = regionprops(bw_diff,'Area');
                mean_area_diff = mean([temp2_props.Area]);
                bad_nu_I = find([temp2_props.Area] > bad_ratio*mean_area_diff);
                bad_bw0 = ismember(bwlabel(bw_diff),bad_nu_I);
                bad_bw = imerode(bad_bw0 & (W ~= 0),strel('disk',5));
                bw_diff(bad_bw0) = bad_bw(bad_bw0);
                bw_diff = imclose(bw_diff,strel('disk',10));
                bw_diff = imfill(bw_diff,'holes');
                bw_diff = imopen(bw_diff,strel('disk',10));
                bw_diff = bwconvhull(bw_diff,'objects');
            end

            
            %%% Recognize the center points of nuclei
            bw_diff_label = bwlabel(bw_diff);
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
            fil_im = uint16(double(imfilter(input_im,H,'same','conv'))./double(MidI));
% %             fil_im2 = imfilter(input_im,fspecial('disk',imblur2),'same','conv');

            if t_reuse == 3 || t_reuse == 4 || t_reuse == 5
                im_max = im_max0;
                th_2D = th_2D0;
                bw_im = bw_im0;
                bw_out = false(size(bw_im));
            else
                im_max = max(fil_im(:));
                [th_2D,bw_im] = circular_th(fil_im,I_max,I_max0,im_max,r_close,low_th,high_th,cir_th);
                bw_out = bw_im;
            end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%% Supplement segmentation based on center recognition: %%%%%%%%%%%%%%%%%%%
            bw_props = regionprops(bw_im,'Centroid');
            if ~isempty(bw_props) && t_reuse
                center_xy2 = round(cell2mat({bw_props.Centroid}'));
                bw_center2 = sub2ind(size(bw_im),center_xy2(:,2),center_xy2(:,1));
                label0 = find(bw_diff_label(bw_center2) == 0);
                if t_reuse ~= 5
                    bw_im(ismember(bwlabel(bw_im),label0)) = 0;
                end
                if t_reuse ~= 3 && t_reuse ~= 4 && t_reuse ~= 5
                    bw_out = bw_im;
                end
            end    
                
            bw_props = regionprops(bw_im,'Centroid');
            if ~isempty(bw_props)
                if t_reuse ~= 5
                    center_xy2 = round(cell2mat({bw_props.Centroid}'));
                    bw_center2 = sub2ind(size(bw_im),center_xy2(:,2),center_xy2(:,1));
                    th_center2 = th_2D(bw_center2);

    % %                 bw_im2 = imclose(bw_im,strel('disk',5));
                    bw_im2 = bw_out;
                    bw_center0 = bw_center(~bw_im2(bw_center));
    % %                 bw_center0_area = ismember(bwlabel(bw_diff),find((~bw_im2(bw_center))));
                    [center_y0,center_x0] = ind2sub(size(bw_diff),bw_center0);
                    center_xy0 = [center_x0,center_y0];

                    for I_center = 1:size(bw_center0,1)
                        center_dist2 = (center_xy2(:,1)-center_xy0(I_center,1)).^2+(center_xy2(:,2)-center_xy0(I_center,2)).^2;
                        th_estimate = th_ratio*sum(th_center2.*exp(-center_dist2/2/sigma_th^2))/sum(exp(-center_dist2/2/sigma_th^2));
                        bw_th = imfill(bwareaclose(fil_im >= th_estimate*65535,merge_th),'holes');
    % % %                     bw_th = bwareaclose(im2bw(fil_im,th_estimate),merge_th);
    % %                     bw_th = imfill(bwareaclose(im2bw(fil_im,th_estimate),merge_th),'holes');
    % % % %                     bw_th = im2bw(im_temp,th_estimate);
                        label0 = bw_diff_label(bw_center0(I_center));
                        if label0
                            bw_th_label = bwlabel(bw_th);
                            label_select = setdiff(unique(bw_th_label(bw_diff_label == label0)),0);
                            bw_I = ismember(bw_th_label,label_select);% & bw_th;
        % %                     if nnz(bw_I)
        % %                         aa(:,:,3) = fil_im*6; aa(:,:,2) = uint16(bwperim(bw_out)*65535); aa(:,:,1) = uint16(bwperim(bw_I)*65535); imshow(aa); title(num2str(I_center))
        % %                         pause
        % %                     end
                            if isempty(label_select) && (t_reuse == 1 || t_reuse == 3)
                                bw_I = bw_diff_label == label0;
                            end
                            bw_out = bw_out | bw_I;
                            th_2D(bw_I) = max(th_2D(bw_I),th_estimate);
                        end
                    end
                    
                else
                    bw_th = imclose(fil_im >= th_ratio*th_2D(1)*65535,strel('disk',10));
                    bw_th_label = bwlabel(bw_th);
                    label_select = setdiff(unique(bw_th_label(bw_diff)),0);
                    bw_out = ismember(bw_th_label,label_select);% & bw_th;
                end
                
% %                 bw_th = pixel_th(X0,Y0,center_xy2,th_center2);
% %                 bw_out = imfill(bwareaopen(bwareaclose(fil_im >= 1.0*bw_th*65535,merge_th),low_th),'holes');
                
                bw_out = bwareaopen(imopen(bw_out,strel('disk',3)),low_th);
                if r_close(2) > 0
                    bw_out = imclose(bw_out,strel('disk',r_close(2)));
                    if t_reuse ~= 5
                        bw_out = imfill(bw_out,'holes');
                    else
                        bw_out = bw_out | bwareaclose(~bw_out,300);
                    end
                end
%                 bw_out = imopen(bw_out,mask_elp0);
% % %                 bw_cut = false(size(fil_im));
% % %     %             bw_cut([bw_center2;bw_center0]) = true;
% % %                 bw_cut(bw_center2) = true;
% % %                 bw_cut = bw_cut | (bw_center0_area & bw_out);

% % %                 temp_label = bwlabel(bw_out);
% % %                 temp_sel = ismember(temp_label,temp_label(bw_cut));
% %                 D= -bwdist(~bw_out);
% %                 bw_cut = imextendedmin(D,2);
% % % % %                 bw_cut = ellipse_vote(bw_out,fil_im);
% %                 g2 = imimposemin(D,bw_cut);
% %                 bw_cut = imdilate(watershed(g2) == 0,strel('disk',1));% & (temp_sel);
% %                 bw_out = bw_out & (~ bw_cut);
                
% % % %                 bw_cut =  imfindmax(bw_out,bw_im,fil_im2,10,5,theta0);
% %                 im_erode0 = imerode(bw_out,ellipse_mask(8,3,theta0));
% %                 D = -bwdist(~im_erode0);
% %                 bw_cut = imextendedmin(D,2);
                
% %                 bw_out = imclose(bw_out,strel('disk',12));
                
                D= -bwdist(~bw_out);
                bw_cut = imextendedmin(D,2);
                g2 = imimposemin(D,bw_cut);
                bw_cut = imdilate(watershed(g2) == 0,strel('disk',1));% & (temp_sel);
                bw_out = bw_out & (~ bw_cut);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%% Mask area refinement and convexation: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %                 if r_close > 0
% % %                     bw_out = imclose(bw_out,strel('disk',r_close));
% % %                 end
%                 bw_out = imopen(bw_out,strel('disk',10));
                
                bw_out  = ellipse_fit(bw_out,20);   %%% refine the mask using elliptical shape
                bw_out  = ellipse_fit(bw_out,20);   %%% refine the mask using elliptical shape

                if back_th > 0
                    im_refine = imfilter(input_im,fspecial('disk',5),'same','conv');
                    Inten_back = sort(im_refine(~bw_out));
                    back_thresh = back_th*Inten_back(round(length(Inten_back)*0.9));
                    bw_refine0 = bwareaopen(bw_out & (im_refine >= back_thresh),low_th);
                    bw_refine0 = imfill(bw_refine0,'holes');
                else
                    bw_refine0 = bwareaopen(bw_out,low_th);
                    bw_refine0 = imfill(bw_refine0,'holes');
                end
                if r_close(3) > 0
                    bw_refine0 = imclose(bw_refine0,strel('disk',r_close(3)));
                elseif r_refine > 0
                    bw_refine0 = imopen(bw_refine0,strel('disk',r_refine));
                end
                if r_erode > 0
                    bw_refine0 = imerode(bw_refine0,strel('disk',r_erode));
                end
                
                if r_close(3) > 0
                    D= -bwdist(~bw_refine0);
                    bw_separate = imextendedmin(D,2);
                    g2 = imimposemin(D,bw_separate);
                    bw_separate = imdilate(watershed(g2) == 0,strel('disk',1));
                    bw_refine0 = bw_refine0 & (~ bw_separate);
                end
                
                bw_refine = bwconvhull(bw_refine0,'objects');
                
                if t_reuse == 1 || t_reuse == 3
                    bw_center0 = bw_center(~bw_refine(bw_center));
                    label0 = setdiff(bw_diff_label(bw_center0),0);
                    bw_I = ismember(bw_diff_label,label0);
                    bw_refine = bw_refine | bw_I;
                end

                D= -bwdist(~bw_refine);
                bw_separate = imextendedmin(D,2);
                g2 = imimposemin(D,bw_separate);
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
        delete(pid)
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

        if t_reuse == 1 || t_reuse == 3
            if scale0 > 1
                raw3D = imresize(raw3D,1/scale0,'nearest');
                raw3D = cat(1,raw3D,false(n1,size(raw3D,2),size(raw3D,3)));
                raw3D = cat(2,raw3D,false(size(raw3D,1),n2,size(raw3D,3)));
            end
            eval([reuse_diff,' = raw3D;'])
            save([folder_list{list_I,1},out_folder,sub_list{list_J,3},mask_name],reuse_diff,'-v7.3')
            clear raw3D
        else
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

            if scale0 > 1
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
    r1 = round(max(min(median([s0.MinorAxisLength]/2)-10,r0),median([s0.MinorAxisLength]/4)));
%     r1 = round(min(median([s0.MinorAxisLength]/2)-10,r0));
    
    bw_in1 = imerode(bw_in,strel('disk',r1));
    s = regionprops(bw_in1, 'Orientation', 'MajorAxisLength','MinorAxisLength', 'Eccentricity', 'Centroid');
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
    
    s2 = regionprops(bw_in,bw_in1, 'MaxIntensity');
    ind1 = find(~[s2.MaxIntensity]);
    if ~isempty(ind1)
        bw_out = bw_out | ismember(bwlabel(bw_in),ind1);
    end
else
    bw_out = bw_in;
end




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A function to create marker for overlapping elliptical areas by finding local maximum 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function bw_mask =  imfindmax(bw_out,bw_im,fil_im,Lmajor,Lminor,theta0)

bw_out_label = bwlabel(bw_out);
immax0 = imregionalmax(fil_im) & imerode(bw_out,strel('disk',5));
immax0_label = bw_out_label.*immax0;
bw_mask = (imdilate(immax0_label,ellipse_mask(Lmajor,Lminor,theta0)) == bw_out_label) & bw_out;
bw_mask = imclose(bw_mask,ellipse_mask(3,1.5,theta0));

bw_mask = (bw_mask | bw_im) & bw_out;




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A function to calculate the threshold for every pixel in the image %%%%%%%%%%%%%%%%%%% 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function bw_th = pixel_th(X0,Y0,center_xy2,th_center2,varargin)

if ~isempty(varargin)
    sigma_th = varargin{1};
else
    sigma_th = 100;
end

th_sum = zeros(size(X0));
weight_sum = zeros(size(X0));

for I_center = 1:size(center_xy2,1)
    center_dist2 = (X0-center_xy2(I_center,1)).^2+(Y0-center_xy2(I_center,2)).^2;
    th_sum = th_sum+th_center2(I_center)*exp(-center_dist2/2/sigma_th^2);
    weight_sum = weight_sum+exp(-center_dist2/2/sigma_th^2);
end

bw_th = th_sum./weight_sum;





%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A function to find threshold and boundary of nucleus using criterion of circle  %%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [th_2D,bw_im] = circular_th(fil_im,I_max,I_max0,im_max,r_close,low_th,high_th,cir_th)
            
th_2D = zeros(size(fil_im)); %%% Raw 2D threshold
im_temp = fil_im;
bw_im = false(size(fil_im));
%%% Cicle search
for I = 1:(I_max-1)
    th_I = I/I_max0*double(im_max)/65535;
%     bw_temp = imfill(im2bw(im_temp,th_I),'holes');
    bw_temp = imfill(im_temp >= th_I*65535,'holes');
    if r_close(1) > 0
% %                     bw_temp = imopen(bw_temp,strel('disk',3));
        bw_temp = imclose(bw_temp,strel('disk',r_close(1)));
        bw_temp = imfill(bw_temp,'holes');
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

