function DAPI_seg3D(lf_name)
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
low_th = 1500;
high_th = 50000;
merge_th = 100000;
cir_th = 0.8;
th_ratio = 1.0;
imblur = 20;  
sigma_th = 100;

% Advanced Guassian filtering    
Ex = fspecial('gaussian',100,20);
Ix = fspecial('gaussian',100,30);
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
%     channel_name = eval(folder_list{list_I,5});
    for list_J = 1:M1
        image_folder = [folder_list{list_I,1},in_folder,sub_list{list_J,3}];
        Mdim = sub_num(list_J,3);
        Nbin1 = sub_num(list_J,2);
        WGA_channel = sub_num(list_J,6);
        DAPI_channel = sub_num(list_J,7);
        RNA_channel = sub_num(list_J,10);
        protein_channel = sub_num(list_J,8);
        resolution = sub_num(list_J,9);
        
        imlist = dir([folder_list{list_I,1},in_folder,sub_list{list_J,3},image_type]); %%% get the image list from image folder1
        max_temp = zeros(1,length(imlist));
matlabpool
        parfor image_I = 1:length(imlist)
            raw_im = imread([folder_list{list_I,1},in_folder,sub_list{list_J,3},imlist(image_I).name]);
            max_temp(image_I) = max(max(raw_im(:,:,DAPI_channel)));
        end
matlabpool close
        clear raw_im
        im_max = max(max_temp);
        size_im = size(imread([folder_list{list_I,1},in_folder,sub_list{list_J,3},imlist(1).name]));
        raw3D = false([size_im([1:2]),length(imlist)]); %%% Raw 3D mask
        th_3D = zeros([size_im([1:2]),length(imlist)]); %%% Raw 3D threshold
        fil3D = zeros([size_im([1:2]),length(imlist)]); %%% Smoothed image
        
matlabpool
        parfor image_I = 1:length(imlist)
            raw_im = imread([folder_list{list_I,1},in_folder,sub_list{list_J,3},imlist(image_I).name]);
            th_2D = zeros(size_im([1:2])); %%% Raw 2D threshold
            input_im = raw_im(:,:,DAPI_channel);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        

%% Find region centers: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            H = fspecial('disk',imblur); % Filter Kernel       
            I_blur = imfilter(input_im,H,0); %Apply Filter
            I_blur =  adapthisteq(I_blur); %Step 1: enhances the contrast of the grayscale image   


            % Faster method to apply filter -- use straight Gaussians.  
            outE = imfilter(single(I_blur),Ex,'replicate'); 
            outI = imfilter(single(I_blur),Ix,'replicate'); 
            outims = outE - outI;  

            outims(outims < 0) = 0;
            outims = outims/max(outims(:));

            % Recognize the center points of nuclei
            bw_diff = im2bw(outims,1.2*graythresh(outims)); 
            diff_props = regionprops(bw_diff,'Centroid');
            center_xy = round(cell2mat({diff_props.Centroid}'));
            bw_center = sub2ind(size(bw_diff),center_xy(:,2),center_xy(:,1));
            
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% % % %% Segmentation: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %             %fil_im = imfilter(input_im,fspecial('gaussian',15,5),'same','conv');
% % %             fil_im = imfilter(input_im,fspecial('disk',20),'same','conv');
% % % 
% % %             im_temp = fil_im;
% % %             bw_im = false(size(fil_im));
% % %             for I = 1:(I_max-1)
% % %                 th_I = I/I_max*double(im_max)/65535;
% % %                 bw_temp = imfill(im2bw(im_temp,th_I),'holes');
% % %                 %bw_temp = imerode(bw_temp,strel('disk',5));
% % %                 bw_prop = regionprops(bw_temp,'Area','Perimeter');
% % %                 bw_area = [bw_prop.Area];
% % %                 bw_perim = [bw_prop.Perimeter];
% % %                 ind_true = find(bw_area >= low_th & bw_area <= high_th & 4*pi*bw_area./bw_perim.^2 >= cir_th);
% % %                 bw_true = ismember(bwlabel(bw_temp),ind_true);
% % %                 bw_im = bw_im | bw_true;
% % %                 im_temp(bw_true) = 0;
% % %                 th_2D(bw_true) = th_I;
% % %             end
% % % %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Segmentation: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            fil_im = imfilter(input_im,fspecial('disk',20),'same','conv');
            im_max = max(fil_im(:));
            im_temp = fil_im;
            bw_im = false(size(fil_im));
%%% Cicle search
            for I = 1:(I_max-1)
                th_I = I/I_max*double(im_max)/65535;
                bw_temp = imfill(im2bw(im_temp,th_I),'holes');
                %bw_temp = imerode(bw_temp,strel('disk',5));
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

            bw_props = regionprops(bw_im,'Centroid');
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

%%% Mask shape refinement
            im_refine = imfilter(input_im,fspecial('disk',5),'same','conv');
            Inten_back = sort(im_refine(~bw_out));
            back_thresh = 1.5*Inten_back(round(length(Inten_back)*0.9));
            bw_refine0 = bwareaopen(bw_out & (im_refine >= back_thresh),low_th);
            bw_refine = bwconvhull(bw_refine0,'objects');

            D= bwdist(~bw_refine);
            g2 = imimposemin(-D,bw_refine0);
            bw_separate = imdilate(watershed(g2) == 0,strel('disk',1));
            bw_out = bw_refine & (~ bw_separate);
            
            
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Raw segmentation output: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            raw3D(:,:,image_I) = bw_out;
            th_3D(:,:,image_I) = th_2D;
            fil3D(:,:,image_I) = fil_im;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        parfor image_I = 1:size(raw3D,3)
            raw3D(:,:,image_I) = imopen(raw3D(:,:,image_I),strel('disk',10));
%             raw3D(:,:,image_I) = imerode(raw3D(:,:,image_I),strel('disk',2));
        end
        
% % % % %% Customized imopen: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % %         [X,Y] = meshgrid(1:size(raw3D,2),1:size(raw3D,1));
% % % %         parfor image_I = 1:size(raw3D,3)
% % % %             raw_im = imread([folder_list{list_I,1},in_folder,sub_list{list_J,3},imlist(image_I).name]);
% % % %             input_im = raw_im(:,:,DAPI_channel);
% % % %             nu_props = regionprops(raw3D(:,:,image_I),'Centroid','Eccentricity','BoundingBox','Image');
% % % %             nux2_props = regionprops(raw3D(:,:,image_I),X.^2,'MeanIntensity');
% % % %             nuy2_props = regionprops(raw3D(:,:,image_I),Y.^2,'MeanIntensity');
% % % % 
% % % %             nuI_props = regionprops(raw3D(:,:,image_I),double(input_im),'WeightedCentroid','MeanIntensity');
% % % %             nuIx2_props = regionprops(raw3D(:,:,image_I),double(input_im).*X.^2,'MeanIntensity');
% % % %             nuIy2_props = regionprops(raw3D(:,:,image_I),double(input_im).*Y.^2,'MeanIntensity');
% % % % 
% % % %             epsilon = [nu_props.Eccentricity]';
% % % %             nu_xy = cell2mat({nu_props.Centroid}');
% % % %             nu_x2 = [nux2_props.MeanIntensity]';
% % % %             nu_y2 = [nuy2_props.MeanIntensity]';
% % % % 
% % % %             nuI_xy = cell2mat({nuI_props.WeightedCentroid}');
% % % %             nuI_x2 = [nuIx2_props.MeanIntensity]'./[nuI_props.MeanIntensity]';
% % % %             nuI_y2 = [nuIy2_props.MeanIntensity]'./[nuI_props.MeanIntensity]';
% % % % 
% % % %             nu_r = sqrt((nu_x2+nu_y2-nu_xy(:,1).^2-nu_xy(:,2).^2));%*3*2.*(1-epsilon.^2)./(2-epsilon.^2));
% % % %             nuI_r = sqrt((nuI_x2+nuI_y2-nuI_xy(:,1).^2-nuI_xy(:,2).^2));%*3*2.*(1-epsilon.^2)./(2-epsilon.^2));
% % % %             dr = nuI_r-nu_r;
% % % % 
% % % %             nu_bound = cell2mat({nu_props.BoundingBox}');
% % % %             bw_out_new = false(size(input_im));
% % % %             for I_nu = 1:length(nu_props)
% % % %                 x_low = nu_bound(I_nu,1);
% % % %                 y_low = nu_bound(I_nu,2);
% % % %                 x_high = x_low+nu_bound(I_nu,3)-1;
% % % %                 y_high = y_low+nu_bound(I_nu,4)-1;
% % % %                 nu_temp = false(nu_bound(I_nu,4)+2,nu_bound(I_nu,3)+2);
% % % %                 nu_temp(2:end-1,2:end-1) = nu_props(I_nu).Image;
% % % %                 nu_temp = imerode(nu_temp,strel('disk',min(10,round(max(0,-dr(I_nu))))));
% % % %                 bw_out_new(y_low:y_high,x_low:x_high) = bw_out_new(y_low:y_high,x_low:x_high) | nu_temp(2:end-1,2:end-1);
% % % %             end   
% % % %             raw3D(:,:,image_I) = bw_out_new;
% % % %         end
% % % % %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

matlabpool close
        
        
%% 3D mask refinement: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        mask_label = label3D(raw3D);   %%% mask stitching and alignment
        clear raw3D
    %%% mask rethresholding: %%% ==========================================
%         mask_stack = zeros(size(raw3D));
%         parfor I_nuclei = 1:max(mask_label(:))
%             temp_label = mask_label == I_nuclei;
%             mask_stack(temp_label) = I_nuclei*im2bw(fil3D(temp_label),max(th_3D(temp_label)));
%         end
        label_prop = regionprops(mask_label,th_3D,'MaxIntensity');
        th_max = [1,[label_prop.MaxIntensity]];
        mask_stack = mask_label.*(fil3D > th_ratio*65535*th_max(mask_label+1));
        mask_stack = mask_label;
        clear mask_label th_max
    %%% ===================================================================
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%toc


%% Resign the label #: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        new_label = zeros(1,max(mask_stack(:))+1);
        new_label(sort(unique(mask_stack))+1) = [0:(length(unique(mask_stack))-1)];
        mask_stack = new_label(mask_stack+1);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%toc
%% Fine segmentation output: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ~exist([folder_list{list_I,1},out_folder,sub_list{list_J,3}],'dir')
            mkdir([folder_list{list_I,1},out_folder,sub_list{list_J,3}])
        end
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
        save([folder_list{list_I,1},out_folder,sub_list{list_J,3},mask_name],'mask_stack','-v7.3')
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end
end

toc



% function mask_label = label3D(raw3D)
%     
% end
