clear %%all
close all

% function Dualprocess7_reassign
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A function to add the old version of processed nuclear segmentation/RNA identification results to the new version of unprocessed results
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Parameter setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
list_name = 'Duallist.xls';
in_folder = 'stacks/';
old_add = '_old';
input_name = 'matchlist.xls';

stitch_name = 'stitch_parameter.mat';
image_type = '*.tif';

mask_folder = 'masks/';
mask_name = 'mask.mat';
quick_mask_name = 'quick_mask.mat';

hist_folder = 'Histogram/';
hist_folder2 = 'Histogram_RNA2/';

hist_old_tail = '_raw_old.xls';
hist_tail = '_raw.xls';
mat_tail = '.mat';
seg_add = '_seg';
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Data loading: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[num_list, folder_list] = xlsread(list_name);
folder_list = folder_list(strcmpi('T',folder_list(:,6)),:);
[N1,N2] = size(folder_list);
run_list = {2,1,1};
for list_I = 1:2%N1
    [sub_num, sub_list] = xlsread([folder_list{list_I,1},in_folder,input_name]);
    [M1,M2] = size(sub_list);
    
    for list_J = run_list{list_I}
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% coordinates conversion: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        image_folder0 = [folder_list{list_I,1},in_folder,sub_list{list_J,3}(1:end-1),old_add,'/'];
        N_z0 = length(dir([image_folder0,image_type]));
        S0 = load([image_folder0,stitch_name],'x_start','y_start','z_start','imcorr1');
        imcorr1 = S0.imcorr1;
        x_start0 = S0.x_start; x_max0 = sum(size(imcorr1,2)-x_start0+1,2); dx0 = x_max0-min(x_max0);
        y_start0 = S0.y_start; y_max0 = sum(size(imcorr1,1)-y_start0+1,1); dy0 = y_max0-min(y_max0);
        z_start0 = S0.z_start; 
        x_end0 = size(imcorr1,2)*ones(size(x_start0)); x_end0(:,end) = x_end0(:,end)-dx0;
        y_end0 = size(imcorr1,1)*ones(size(y_start0)); y_end0(end,:) = y_end0(end,:)-dy0;
        z_end0 = z_start0+N_z0-1;
        [x_start_im0,y_start_im0,z_start_im0,x_end_im0,y_end_im0,z_end_im0] = tile_info_extract2(x_start0,y_start0,z_start0,x_end0,y_end0,z_end0);
        
        image_folder1 = [folder_list{list_I,1},in_folder,sub_list{list_J,3}];
        N_z1 = length(dir([image_folder1,image_type]));
        S1 = load([image_folder1,stitch_name],'x_start','y_start','z_start');
        x_start1 = S1.x_start; x_max1 = sum(size(imcorr1,2)-x_start1+1,2); dx1 = x_max1-min(x_max1);
        y_start1 = S1.y_start; y_max1 = sum(size(imcorr1,1)-y_start1+1,1); dy1 = y_max1-min(y_max1);
        z_start1 = S1.z_start; 
        x_end1 = size(imcorr1,2)*ones(size(x_start1)); x_end1(:,end) = x_end1(:,end)-dx1;
        y_end1 = size(imcorr1,1)*ones(size(y_start1)); y_end1(end,:) = y_end1(end,:)-dy1;
        z_end1 = z_start1+N_z1-1;
        [x_start_im1,y_start_im1,z_start_im1,x_end_im1,y_end_im1,z_end_im1] = tile_info_extract2(x_start1,y_start1,z_start1,x_end1,y_end1,z_end1);
        
        x_start = max(x_start0,x_start1); x_end = min(x_end0,x_end1);
        y_start = max(y_start0,y_start1); y_end = min(y_end0,y_end1);
        z_start = max(z_start0,z_start1); z_end = min(z_end0,z_end1);
        
        x_start_im_old = x_start_im0+x_start-x_start0;
        y_start_im_old = y_start_im0+y_start-y_start0;
        z_start_im_old = z_start_im0+z_start-z_start0;
        x_end_im_old = x_end_im0+x_end-x_end0;
        y_end_im_old = y_end_im0+y_end-y_end0;
        z_end_im_old = z_end_im0+z_end-z_end0;
        
        x_start_im_new = x_start_im1+x_start-x_start1;
        y_start_im_new = y_start_im1+y_start-y_start1;
        z_start_im_new = z_start_im1+z_start-z_start1;
        x_end_im_new = x_end_im1+x_end-x_end1;
        y_end_im_new = y_end_im1+y_end-y_end1;
        z_end_im_new = z_end_im1+z_end-z_end1;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Add old nuclear segmentation results to the new ones: %%%%%%%%%%%
        S0 = load([folder_list{list_I,1},mask_folder,sub_list{list_J,3}(1:end-1),old_add,'/',quick_mask_name],'bw_applied3D_save');   %%% load old 3D mask
        mask0 = logical(S0.bw_applied3D_save);
        
        S1 = load([folder_list{list_I,1},mask_folder,sub_list{list_J,3},mask_name],'mask_stack');   %%% load new 3D mask
        bw_applied3D_save = logical(S1.mask_stack);
        
        if size(mask0,2) ~= x_end_im0(end,end) || size(mask0,1) ~= y_end_im0(end,end) || size(bw_applied3D_save,2) ~= x_end_im1(end,end) || size(bw_applied3D_save,1) ~= y_end_im1(end,end)
            error('Dimension mismatch')
        end
        for ii = 1:size(x_start,1)
            for jj = 1:size(x_start,2)
                bw_applied3D_save(y_start_im_new(ii,jj):y_end_im_new(ii,jj),x_start_im_new(ii,jj):x_end_im_new(ii,jj),z_start_im_new(ii,jj):z_end_im_new(ii,jj)) = mask0(y_start_im_old(ii,jj):y_end_im_old(ii,jj),x_start_im_old(ii,jj):x_end_im_old(ii,jj),z_start_im_old(ii,jj):z_end_im_old(ii,jj));
            end
        end
        
        save([folder_list{list_I,1},mask_folder,sub_list{list_J,3},quick_mask_name],'bw_applied3D_save','-v7.3')
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Modify the coordinates of RNA recognition results: %%%%%%%%%%%%%%
        [foci_list0,~,~] = xlsread([folder_list{list_I,1},hist_folder,sub_list{list_J,3}(1:(end-1)),hist_old_tail]);
        foci_list = foci_list0;
        xyz = foci_list0(:,6:8);
        Itrue = false(size(foci_list0,1),1);
        for ii = 1:size(x_start,1)
            for jj = 1:size(x_start,2)
                Itrue0 = xyz(:,1) >= y_start_im_old(ii,jj) & xyz(:,1) <= y_end_im_old(ii,jj) & xyz(:,2) >= x_start_im_old(ii,jj) & xyz(:,2) <= x_end_im_old(ii,jj) & xyz(:,3) >= z_start_im_old(ii,jj) & xyz(:,3) <= z_end_im_old(ii,jj);
                Itrue = Itrue | Itrue0;
                foci_list(Itrue0,6) = foci_list0(Itrue0,6)+y_start_im_new(ii,jj)-y_start_im_old(ii,jj);
                foci_list(Itrue0,7) = foci_list0(Itrue0,7)+x_start_im_new(ii,jj)-x_start_im_old(ii,jj);
                foci_list(Itrue0,8) = foci_list0(Itrue0,8)+z_start_im_new(ii,jj)-z_start_im_old(ii,jj);
            end
        end
        foci_list = foci_list(Itrue,:);
        xlswrite([folder_list{list_I,1},hist_folder,sub_list{list_J,3}(1:(end-1)),hist_tail],foci_list)
        
        [foci_list0,~,~] = xlsread([folder_list{list_I,1},hist_folder2,sub_list{list_J,3}(1:(end-1)),hist_old_tail]);
        foci_list = foci_list0;
        xyz = foci_list0(:,6:8);
        Itrue = false(size(foci_list0,1),1);
        for ii = 1:size(x_start,1)
            for jj = 1:size(x_start,2)
                Itrue0 = xyz(:,1) >= y_start_im_old(ii,jj) & xyz(:,1) <= y_end_im_old(ii,jj) & xyz(:,2) >= x_start_im_old(ii,jj) & xyz(:,2) <= x_end_im_old(ii,jj) & xyz(:,3) >= z_start_im_old(ii,jj) & xyz(:,3) <= z_end_im_old(ii,jj);
                Itrue = Itrue | Itrue0;
                foci_list(Itrue0,6) = foci_list0(Itrue0,6)+y_start_im_new(ii,jj)-y_start_im_old(ii,jj);
                foci_list(Itrue0,7) = foci_list0(Itrue0,7)+x_start_im_new(ii,jj)-x_start_im_old(ii,jj);
                foci_list(Itrue0,8) = foci_list0(Itrue0,8)+z_start_im_new(ii,jj)-z_start_im_old(ii,jj);
            end
        end
        foci_list = foci_list(Itrue,:);
        xlswrite([folder_list{list_I,1},hist_folder2,sub_list{list_J,3}(1:(end-1)),hist_tail],foci_list)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end








