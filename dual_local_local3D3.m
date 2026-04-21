function [foci_data,fake_data,h,t_absolute,varargout] = dual_local_local3D3(foci_list0,max_image00,mask_stack,DNA_mask,signal_stack,RNA_stack,nucleus_protein_profile,image_folder,N_cycle,varargin)


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A function to plot the local transcription factor concentration at transcription sites %%
%% foci_data: foci enrichment output data (cell).
%% fake_data: fake foci enrichment output data (cell).

%% for both foci_data and fake_data, here are the list of the meaning of different rows:
% % 1. Foci intensity
% % 2. Protein mean concentration
% % 3. Enrichment background value
% % 4. Enrichment value
% % 5. FISH signal enrichment value
% % 6. Nucleus label
% % 7. Nucleus EL length
% % 8. Nucleus max z slice (for mean concentration extraction)
% % 9. Enrichment integration area ratio
% % 10.Background-ring integration area ratio
% % 11.Foci linear coordinate
% % 12.Protein mean concentration at foci located slice
% % 13.Protein background mean concentration at foci located slice
% % 14.Protein background concentration variance at foci located slice
% % 15.Integral covered ratio assuming a single point source at the center
% % 16.z dimension Integral covered ratio assuming a single point source at the center
% % 17.Nuclear area at the foci slice
% % (previous) 16.Protein background mean concentration (from low DNA intensity region) at foci located slice
% % (previous) 17.Protein background concentration variance (from low DNA intensity region) at foci located slice
% % 18. RNA2 signal enrichment value

%% h: enrichment mask.
%% t_absolute: a logical variable (true: absolute units, false: relative units)
%% varargout: list of the radius of enrichment mask (cell).

%% foci_list0: transcription foci position indices list.
%% max_image00: transcription foci intensity list.
%% mask_stack: 3D nuclei mask (labeled).
%% DNA_mask: mask of high intensity DNA region (high hoechst signal region) in nuclei (in case the protein signal of these areas is needed).
%% signal_stack: 3D transcription factor image stack.
%% RNA_stack: 3D FISH signal image stack.
%% nucleus_protein_profile: nucleus protein mean concentration profile (EL/mean concentration/...).
%% image_folder: image folder name.
%% N_cycle: embryo cycle.
%% varargin: {resolution, resolutionz, indices list of foci of other species (for control experiment), figure handle list, PSF sigmas, foci reuse data, RNA2_stack}.
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Initialization: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%h = fspecial('gaussian', 7, 2);   %%% Local mask for foci protein concentration
Lmin = -1;%0.15;  %%% EL lower limit for data plot
Lmax = 2;%0.85;  %%% EL higher limit for data plot
N_bin = 15;  %%% Number of bins for enrihment data averaging
rbin = 2;
Nmean = 5;
N_fit = 40;
r_edge = 10;
rwidth = 4;
r_range0 = 250;

N_fb = 20;
w_fb = 0.05;
N_pb = 20;
w_pb = 0.05;
ppbin = [-1:0.05:1];   %%% binning of the ratio histogram
Npp = length(ppbin);
int_z = 0;   %%% z radius for integration

if int_z
    signal_stack3 = imfilter(signal_stack,ones(1,1,int_z*2+1)/(int_z*2+1),'replicate','same','conv');
end

r_size = 0:1:9;   %%% enrichment integration radius list

if isempty(varargin)
    resolution = 1;
    resolutionz = 1/2;
    nM_con = 1;
    unit1 = 'A.U.';
    unit2 = 'A.U.';
    t_absolute = false;
else
    resolution = varargin{1};
    resolutionz = varargin{2};
    nM_con = 6.02e8;
    unit1 = '#';
    unit2 = 'M';
    t_absolute = true;
    resolution0 = 0.083;
    r_size = r_size*max(round(resolution0/resolution),1);
    r_range0 = r_range0*max(round(resolution0/resolution),1);
    
%     if length(varargin) > 2
%         foci_mask = repmat(logical(varargin{3}),[1,1,size(mask_stack,3)]);
%         foci_mask_true = true;
%     else
%         foci_mask = false(size(mask_stack));
%         foci_mask_true = false;
%     end

    if length(varargin) > 2 && ~isempty(varargin{3})
        [foci_mask_x,foci_mask_y] = ind2sub(size(mask_stack),find(logical(varargin{3})));
        foci_mask_true = true;
    else
        foci_mask_x = zeros(0);
        foci_mask_y = zeros(0);
        foci_mask_true = false;
    end

    
    if length(varargin) > 3 && ~isempty(varargin{4})
        figure_list = varargin{4};
    else
        figure_list = [73:77,80,81];
    end
    
    
    if length(varargin) > 4 && ~isempty(varargin{5})
        if length(varargin{5}) == 3
            sigmai = varargin{5}(1);
            sigmaj = varargin{5}(2);
            sigmak = varargin{5}(3);
        else
            sigmai = varargin{5}(1);
            sigmaj = varargin{5}(1);
            sigmak = varargin{5}(2);
        end
    else
        sigmai = 1.35;
        sigmaj = 1.35;
        sigmak = 0.8;
    end
    
    if length(varargin) > 5 && ~isempty(varargin{6})
        foci_data0 = varargin{6}{1};
        fake_data0 = varargin{6}{2};
        r_size = varargin{6}{3};
        protein_protein2_xymismatch = varargin{6}{4};
        foci_reuse = true;
    else
        foci_reuse = false;
    end
    
    if length(varargin) > 6
        RNA2_stack = varargin{7};
    else
        RNA2_stack = [];
    end
end

varargout = {r_size};
h = cell(size(r_size));
% mr_size = 3*ones(size(r_size));
mr_size = r_size+1;   %%% enrichment local background integration radius list
mh = cell(size(r_size));
area_h = zeros(size(h));


foci_data = cell(size(r_size));   %%% foci data matrix initialization
fake_data = cell(size(r_size));   %%% control data matrix initialization
p_ratio0 = cell(size(r_size));
fp_ratio0 = cell(size(r_size));
op_ratio0 = cell(size(r_size));
ofp_ratio0 = cell(size(r_size));
nucleus_protein_profile((nucleus_protein_profile(:,1) < Lmin)|(nucleus_protein_profile(:,1) > Lmax),:) = nan(sum((nucleus_protein_profile(:,1) < Lmin)|(nucleus_protein_profile(:,1) > Lmax)),size(nucleus_protein_profile,2));
nucleus_protein_profile = [nan(1,size(nucleus_protein_profile,2));nucleus_protein_profile];

foci_temp = false(size(mask_stack));
foci_temp(foci_list0) = true;
foci_bw = max(foci_temp,[],3);   %%% 2D foci mask
clear foci_temp

mean_signal2D = nan(size(mask_stack,3),max(mask_stack(:)));
nu_centroid = cell(1,size(mask_stack,3));
% area2D = nan(size(mask_stack,3),max(mask_stack(:)));
% % bwi = logical(mask_stack);
% % bw01 = imerode(bwi,strel('disk',r_edge));% & (~imdilate(foci_mask,strel('disk',r_size(ir))));
bw01 = imerode(logical(mask_stack),strel('disk',r_edge));% & (~imdilate(foci_mask,strel('disk',r_size(ir))));
mask_shrink = double(bw01).*mask_stack;

foci_list = foci_list0(mask_shrink(foci_list0) > 0);
foci_raw = foci_list;
foci_nlabel = mask_stack(foci_list);
[raw_x,raw_y,~] = ind2sub(size(mask_stack),foci_list);
if foci_mask_true
    d_foci = min(pdist2([raw_x,raw_y],[foci_mask_x,foci_mask_y]),[],2);
else
    d_foci = inf(size(raw_x));
end

% matlabpool(4)

for iz = 1:size(mask_stack,3)
    z_start = max(1,iz-int_z);
    z_end = min(size(mask_stack,3),iz+int_z);
    temp_props = regionprops(mask_shrink(:,:,z_start:z_end),signal_stack(:,:,z_start:z_end),'MeanIntensity','Area','Centroid');
    if ~isempty(temp_props)
        mean_signal2D(iz,1:length(temp_props)) = [temp_props.MeanIntensity];
%         area2D(iz,1:length(temp_props)) = [temp_props.Area]*resolution*resolution*2*resolutionz*nM_con;
        nu_centroid{iz} = round(cell2mat({temp_props.Centroid}'));
    end
end

for ir = 1:length(r_size)%;r_size(ir)
    h{ir} = double(getnhood(strel('disk',r_size(ir))));
    mh{ir} = double(getnhood(strel('disk',r_size(ir)+mr_size(ir))));
    hxy = size(h{ir});
    mhxy = size(mh{ir});
    mh{ir}((mhxy(1)-hxy(1))/2+1:(mhxy(1)+hxy(1))/2,(mhxy(2)-hxy(2))/2+1:(mhxy(2)+hxy(2))/2) = mh{ir}((mhxy(1)-hxy(1))/2+1:(mhxy(1)+hxy(1))/2,(mhxy(2)-hxy(2))/2+1:(mhxy(2)+hxy(2))/2)-h{ir};
%    mh{ir} = mh{ir}*sum(h{ir}(:))/sum(mh{ir}(:));
    h{ir} = repmat(h{ir}*resolution*resolution*resolutionz*nM_con,[1,1,2*int_z+1]);
    mh{ir} = repmat(mh{ir}*resolution*resolution*resolutionz*nM_con,[1,1,2*int_z+1]);
    
    %h{ir} = h{ir}/sum(sum(h{ir}));
    %area_h(ir) = sum(sum(h{ir}(:,:,1)));
    area_h(ir) = sum(h{ir}(:));

    %%% Fake spot mask generation: %%% ========================================
% %     foci_bw00 = foci_raw;
%     bwi = logical(mask_stack);
%     mask_shrink = imerode(bwi,strel('disk',r_edge));
% %     bw00_temp = imerode(bwi,strel('disk',r_edge+r_size(ir)+mr_size(ir)));
% % 
% %     mask_shrink = false(size(bw00_temp));
% %     for I_layer = 2:(size(bw00_temp,3)-1)
% %         mask_shrink(:,:,I_layer) = bw00_temp(:,:,I_layer) & bw00_temp(:,:,I_layer-1) & bw00_temp(:,:,I_layer+1);
% %     end
% %     mask_shrink = bwi;
    
%     mask_shrink = bw01.*mask_stack;
    
% %     clear bw00_temp
% %     if foci_mask_true
% %         mask_clean = (~imdilate(foci_mask,strel('disk',2*r_size(ir)))).*mask_shrink;
% %     else
% %         mask_clean = mask_shrink;
% %     end
% %     foci_nlabel_clean = mask_clean(foci_bw00);
% %     foci_bw00 = foci_bw00(logical(foci_nlabel_clean));
% %     foci_nlabel_clean = foci_nlabel_clean(logical(foci_nlabel_clean));
    
    if foci_mask_true
        foci_clean_list = d_foci >= 2*r_size(ir);
        foci_bw00 = foci_raw(foci_clean_list);
        foci_nlabel_clean = foci_nlabel(foci_clean_list);
    else
        foci_bw00 = foci_raw;
        foci_nlabel_clean = foci_nlabel;
    end
    
    if foci_reuse
        foci_bw00 = foci_data0{ir}(:,11);
        foci_nlabel_clean = foci_data0{ir}(:,6);
        [x_raw,y_raw,~] = ind2sub(size(mask_stack),foci_bw00);
        dxy = xyshift([x_raw,y_raw],protein_protein2_xymismatch,foci_bw);
        xy_max = size(foci_bw);
        
%         x0_min = max(1,dxy(:,1)+1);
%         x0_max = min(xy_max(1),dxy(:,1)+xy_max(1));
%         y0_min = max(1,dxy(:,2)+1);
%         y0_max = min(xy_max(2),dxy(:,2)+xy_max(2));
%         x1_min = max(1,-dxy(:,1)+1);
%         x1_max = min(xy_max(1),-dxy(:,1)+xy_max(1));
%         y1_min = max(1,-dxy(:,2)+1);
%         y1_max = min(xy_max(2),-dxy(:,2)+xy_max(2));
    else
        [x_raw,y_raw,~] = ind2sub(size(mask_stack),foci_bw00);
    end


%     foci_bw_neg = (~imdilate(foci_bw,strel('disk',2*r_size(ir))));
    foci_bw_removed = mask_shrink.*repmat(double(~imdilate(foci_bw,strel('disk',2*r_size(ir)))),[1,1,size(mask_stack,3)]);
%     foci_bw_removed2 = mask_shrink.*repmat((~imdilate(foci_bw,strel('disk',r_size(ir)))),[1,1,size(mask_stack,3)]);
%     foci_bw_removed3 = foci_bw_removed2.*(~DNA_mask);

    fake_bw00 = nan(size(foci_bw00));
    fn_mean_clean = nan(size(foci_bw00));
    fn_var_clean = nan(size(foci_bw00));
    nucleus_prop_Area2 = nan(size(foci_bw00));
%     fn_mean_clean3 = nan(size(foci_bw00));
%     fn_var_clean3 = nan(size(foci_bw00));

    xmin0 = max(x_raw-r_range0,1);
    xmax0 = min(x_raw+r_range0,size(foci_bw_removed,1));
    ymin0 = max(y_raw-r_range0,1);
    ymax0 = min(y_raw+r_range0,size(foci_bw_removed,2));
    
    for I_f = 1:length(foci_bw00)
        [~,~,foci_layer] = ind2sub(size(mask_stack),foci_bw00(I_f));
% % %         nucleus_prop = (mask_shrink(:,:,foci_layer) == foci_nlabel_clean(I_f)) & foci_bw_neg;
% %         nucleus_prop = foci_bw_removed(:,:,foci_layer) == foci_nlabel_clean(I_f);
        nucleus_prop = false(size(foci_bw_removed,1),size(foci_bw_removed,2));
        nucleus_prop(xmin0(I_f):xmax0(I_f),ymin0(I_f):ymax0(I_f)) = foci_bw_removed(xmin0(I_f):xmax0(I_f),ymin0(I_f):ymax0(I_f),foci_layer) == foci_nlabel_clean(I_f);
        
        nu_xy = nu_centroid{foci_layer}(foci_nlabel_clean(I_f),[2,1]);
        df_center = round(pdist2(nu_xy,[x_raw(I_f),y_raw(I_f)]));
        hsmall = logical(getnhood(strel('disk',max(0,df_center-rwidth))));
        hbig = logical(getnhood(strel('disk',df_center+rwidth)));
        hsxy = size(hsmall);
        hbxy = size(hbig);
        hbig((hbxy(1)-hsxy(1))/2+1:(hbxy(1)+hsxy(1))/2,(hbxy(2)-hsxy(2))/2+1:(hbxy(2)+hsxy(2))/2) = hbig((hbxy(1)-hsxy(1))/2+1:(hbxy(1)+hsxy(1))/2,(hbxy(2)-hsxy(2))/2+1:(hbxy(2)+hsxy(2))/2)-hsmall;
        rbxy = floor(hbxy/2);

        nu_xmin = max(nu_xy(1)-rbxy(1),1);
        nu_ymin = max(nu_xy(2)-rbxy(2),1);
        nu_xmax = min(nu_xy(1)+rbxy(1),size(nucleus_prop,1));
        nu_ymax = min(nu_xy(2)+rbxy(2),size(nucleus_prop,2));
        
        hb_xmin = 1+rbxy(1)+nu_xmin-nu_xy(1);
        hb_ymin = 1+rbxy(2)+nu_ymin-nu_xy(2);
        hb_xmax = 1+rbxy(1)+nu_xmax-nu_xy(1);
        hb_ymax = 1+rbxy(2)+nu_ymax-nu_xy(2);
        
        nucleus_prop(nu_xmin:nu_xmax,nu_ymin:nu_ymax) = nucleus_prop(nu_xmin:nu_xmax,nu_ymin:nu_ymax) & hbig(hb_xmin:hb_xmax,hb_ymin:hb_ymax);
        
%         if foci_reuse
%             nucleus_prop(x1_min(I_f):x1_max(I_f),y1_min(I_f):y1_max(I_f),foci_layer) = nucleus_prop(x0_min(I_f):x0_max(I_f),y0_min(I_f):y0_max(I_f),foci_layer);
%         end
%         nucleus_prop_Area = nnz(nucleus_prop);
%         nucleus_prop_PixelIdxList = find(nucleus_prop)+size(mask_shrink,1)*size(mask_shrink,2)*(foci_layer-1);
        if foci_reuse
            [nucleus_prop_x,nucleus_prop_y] = find(nucleus_prop);
            nucleus_prop_x = nucleus_prop_x+dxy(1);
            nucleus_prop_y = nucleus_prop_y+dxy(2);
            nucleus_prop_id = (nucleus_prop_x >= 1)&(nucleus_prop_x <= xy_max(1))&(nucleus_prop_y >= 1)&(nucleus_prop_y <= xy_max(2));
            nucleus_prop_PixelIdxList = sub2ind(xy_max,nucleus_prop_x(nucleus_prop_id),nucleus_prop_y(nucleus_prop_id))+size(mask_shrink,1)*size(mask_shrink,2)*(foci_layer-1);
            nucleus_prop_Area = numel(nucleus_prop_PixelIdxList);
        else
            nucleus_prop_PixelIdxList = find(nucleus_prop)+size(mask_shrink,1)*size(mask_shrink,2)*(foci_layer-1);
            nucleus_prop_Area = numel(nucleus_prop_PixelIdxList);
        end
        
        layer_min = max(1,foci_layer-int_z);
        layer_max = min(size(mask_stack,3),foci_layer+int_z);
%         nucleus_prop2 = foci_bw_removed2(:,:,layer_min:layer_max) == foci_nlabel_clean(I_f);
%         nucleus_prop_PixelIdxList2 = find(nucleus_prop2)+size(mask_shrink,1)*size(mask_shrink,2)*(layer_min-1);
        if int_z
            fn_mean_clean(I_f) = mean(signal_stack3(nucleus_prop_PixelIdxList));
            fn_var_clean(I_f) = var(double(signal_stack3(nucleus_prop_PixelIdxList)));
        else
            fn_mean_clean(I_f) = mean(signal_stack(nucleus_prop_PixelIdxList));
            fn_var_clean(I_f) = var(double(signal_stack(nucleus_prop_PixelIdxList)));
        end
        
%         nucleus_prop0 = mask_shrink(:,:,layer_min:layer_max) == foci_nlabel_clean(I_f);
%         nucleus_prop_Area2(I_f) = nnz(nucleus_prop0);
        nucleus_prop_Area2(I_f) = nucleus_prop_Area;
% %         nucleus_prop3 = foci_bw_removed3(:,:,layer_min:layer_max) == foci_nlabel_clean(I_f);
% %         nucleus_prop_PixelIdxList3 = find(nucleus_prop3)+size(mask_shrink,1)*size(mask_shrink,2)*(layer_min-1);
%         nucleus_prop_PixelIdxList3 = nucleus_prop_PixelIdxList2(~DNA_mask(nucleus_prop_PixelIdxList2));
%         fn_mean_clean3(I_f) = mean(signal_stack(nucleus_prop_PixelIdxList3));
%         fn_var_clean3(I_f) = var(signal_stack(nucleus_prop_PixelIdxList3));
        if foci_reuse
            fake_bw00(I_f) = fake_data0{ir}(I_f,11);
        else
            try 
                if nucleus_prop_Area > 1
                    tcontinue = true;
                    while tcontinue
                        random_I = ceil(rand(1)*nucleus_prop_Area);
                        I_temp = nucleus_prop_PixelIdxList(random_I);
                        if ~any(fake_bw00 == I_temp)
                            fake_bw00(I_f) = I_temp;
                            tcontinue = false;
                        end
                    end
               else
                    %random_I = nan;
                    foci_bw00(I_f) = 0;
                end
            catch
                foci_bw00(I_f) = 0;
            end
        end
        
    end
    temp_true = logical(foci_bw00);
    foci_bw00 = foci_bw00(temp_true);
    fake_bw00 = fake_bw00(temp_true);
    fn_mean_clean = fn_mean_clean(temp_true);
    fn_var_clean = fn_var_clean(temp_true);
    nucleus_prop_Area2 = nucleus_prop_Area2(temp_true);
%     fn_mean_clean3 = fn_mean_clean3(temp_true);
%     fn_var_clean3 = fn_var_clean3(temp_true);
    foci_data{ir} = zeros(length(foci_bw00),17);   %%% foci data matrix initialization
    fake_data{ir} = zeros(length(fake_bw00),17);   %%% control data matrix initialization
%     clear temp_ture
    %%% =======================================================================


    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Calculate local protein concentration: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     bw01 = imerode(bwi,strel('disk',r_edge));% & (~imdilate(foci_mask,strel('disk',r_size(ir))));
    x_h = floor(size(h{ir},1)/2);
    y_h = floor(size(h{ir},2)/2);
    z_h = floor(size(h{ir},3)/2);
    mx_h = floor(size(mh{ir},1)/2);
    my_h = floor(size(mh{ir},2)/2);
    mz_h = floor(size(mh{ir},3)/2);
    x_lim = size(mask_stack,1);
    y_lim = size(mask_stack,2);
    z_lim = size(mask_stack,3);

    [x_all,y_all,z_all] = ind2sub(size(mask_stack),foci_bw00);
%%%%%
    x_all_low = max(x_all-x_h,1);
    x_all_high = min(x_all+x_h,x_lim);
    y_all_low = max(y_all-y_h,1);
    y_all_high = min(y_all+y_h,y_lim);
    z_all_low = max(z_all-z_h,1);
    z_all_high = min(z_all+z_h,z_lim);
%%%%%
    x_h_low = 1+x_h+x_all_low-x_all;
    x_h_high = 1+x_h+x_all_high-x_all;
    y_h_low = 1+y_h+y_all_low-y_all;
    y_h_high = 1+y_h+y_all_high-y_all;
    z_h_low = 1+z_h+z_all_low-z_all;
    z_h_high = 1+z_h+z_all_high-z_all;
%%%%%
    mx_all_low = max(x_all-mx_h,1);
    mx_all_high = min(x_all+mx_h,x_lim);
    my_all_low = max(y_all-my_h,1);
    my_all_high = min(y_all+my_h,y_lim);
    mz_all_low = max(z_all-mz_h,1);
    mz_all_high = min(z_all+mz_h,z_lim);
%%%%%
    mx_h_low = 1+mx_h+mx_all_low-x_all;
    mx_h_high = 1+mx_h+mx_all_high-x_all;
    my_h_low = 1+my_h+my_all_low-y_all;
    my_h_high = 1+my_h+my_all_high-y_all;
    mz_h_low = 1+mz_h+mz_all_low-z_all;
    mz_h_high = 1+mz_h+mz_all_high-z_all;
%%%%%

    if foci_reuse
        x0 = -dxy(:,1)+x_all_low;
        y0 = -dxy(:,2)+y_all_low;
        x1_min = max(1,-dxy(:,1)+x_all_low);
        x1_max = min(xy_max(1),-dxy(:,1)+x_all_high);
        y1_min = max(1,-dxy(:,2)+y_all_low);
        y1_max = min(xy_max(2),-dxy(:,2)+y_all_high);
        x0_min = x1_min-x0+1;
        x0_max = x1_max-x0+1;
        y0_min = y1_min-y0+1;
        y0_max = y1_max-y0+1;

        mx0 = -dxy(:,1)+mx_all_low;
        my0 = -dxy(:,2)+my_all_low;
        mx1_min = max(1,-dxy(:,1)+mx_all_low);
        mx1_max = min(xy_max(1),-dxy(:,1)+mx_all_high);
        my1_min = max(1,-dxy(:,2)+my_all_low);
        my1_max = min(xy_max(2),-dxy(:,2)+my_all_high);
        mx0_min = mx1_min-mx0+1;
        mx0_max = mx1_max-mx0+1;
        my0_min = my1_min-my0+1;
        my0_max = my1_max-my0+1;
    end


    [ii0,jj0,kk0] = ind2sub(size(h{ir}),find(h{ir})); 
%     if size(h{ir},3) == 1
%         Integral0 = 2*pi*sigmai*sigmaj;
%     else
        Integral0 = 2*pi*sigmai*sigmaj*sqrt(2*pi)*sigmak;
        Integralz0 = sqrt(2*pi)*sigmak;
%     end
    
    for I_ind = 1:length(x_all)
        if foci_reuse
            signal_stack0 = double(signal_stack(x_all_low(I_ind):x_all_high(I_ind),y_all_low(I_ind):y_all_high(I_ind),z_all_low(I_ind):z_all_high(I_ind)));
            signal_stack0(x0_min(I_ind):x0_max(I_ind),y0_min(I_ind):y0_max(I_ind),:) = signal_stack(x1_min(I_ind):x1_max(I_ind),y1_min(I_ind):y1_max(I_ind),z_all_low(I_ind):z_all_high(I_ind));

            RNA_stack0 = double(RNA_stack(x_all_low(I_ind):x_all_high(I_ind),y_all_low(I_ind):y_all_high(I_ind),z_all_low(I_ind):z_all_high(I_ind)));
            
            signal_stack1 = double(signal_stack(mx_all_low(I_ind):mx_all_high(I_ind),my_all_low(I_ind):my_all_high(I_ind),mz_all_low(I_ind):mz_all_high(I_ind)));
            signal_stack1(mx0_min(I_ind):mx0_max(I_ind),my0_min(I_ind):my0_max(I_ind),:) = signal_stack(mx1_min(I_ind):mx1_max(I_ind),my1_min(I_ind):my1_max(I_ind),mz_all_low(I_ind):mz_all_high(I_ind));
            
            foci_data{ir}(I_ind,4) = sum(sum(sum(signal_stack0.*bw01(x_all_low(I_ind):x_all_high(I_ind),y_all_low(I_ind):y_all_high(I_ind),z_all_low(I_ind):z_all_high(I_ind)).*h{ir}(x_h_low(I_ind):x_h_high(I_ind),y_h_low(I_ind):y_h_high(I_ind),z_h_low(I_ind):z_h_high(I_ind)))));
            foci_data{ir}(I_ind,5) = sum(sum(sum(   RNA_stack0.*bw01(x_all_low(I_ind):x_all_high(I_ind),y_all_low(I_ind):y_all_high(I_ind),z_all_low(I_ind):z_all_high(I_ind)).*h{ir}(x_h_low(I_ind):x_h_high(I_ind),y_h_low(I_ind):y_h_high(I_ind),z_h_low(I_ind):z_h_high(I_ind)))));
            foci_data{ir}(I_ind,9) = sum(sum(sum(bw01(x_all_low(I_ind):x_all_high(I_ind),y_all_low(I_ind):y_all_high(I_ind),z_all_low(I_ind):z_all_high(I_ind)).*h{ir}(x_h_low(I_ind):x_h_high(I_ind),y_h_low(I_ind):y_h_high(I_ind),z_h_low(I_ind):z_h_high(I_ind)))))/sum(sum(sum(h{ir}(x_h_low(I_ind):x_h_high(I_ind),y_h_low(I_ind):y_h_high(I_ind),:))));
            foci_data{ir}(I_ind,10) = sum(sum(sum(bw01(mx_all_low(I_ind):mx_all_high(I_ind),my_all_low(I_ind):my_all_high(I_ind),mz_all_low(I_ind):mz_all_high(I_ind)).*mh{ir}(mx_h_low(I_ind):mx_h_high(I_ind),my_h_low(I_ind):my_h_high(I_ind),mz_h_low(I_ind):mz_h_high(I_ind)))))/sum(sum(sum(mh{ir}(mx_h_low(I_ind):mx_h_high(I_ind),my_h_low(I_ind):my_h_high(I_ind),:))));
            f_ratio = sum(sum(sum(bw01(x_all_low(I_ind):x_all_high(I_ind),y_all_low(I_ind):y_all_high(I_ind),z_all_low(I_ind):z_all_high(I_ind)).*h{ir}(x_h_low(I_ind):x_h_high(I_ind),y_h_low(I_ind):y_h_high(I_ind),z_h_low(I_ind):z_h_high(I_ind)))))/sum(sum(sum(bw01(mx_all_low(I_ind):mx_all_high(I_ind),my_all_low(I_ind):my_all_high(I_ind),mz_all_low(I_ind):mz_all_high(I_ind)).*mh{ir}(mx_h_low(I_ind):mx_h_high(I_ind),my_h_low(I_ind):my_h_high(I_ind),mz_h_low(I_ind):mz_h_high(I_ind)))));
            foci_data{ir}(I_ind,3) = sum(sum(sum(signal_stack1.*bw01(mx_all_low(I_ind):mx_all_high(I_ind),my_all_low(I_ind):my_all_high(I_ind),mz_all_low(I_ind):mz_all_high(I_ind)).*mh{ir}(mx_h_low(I_ind):mx_h_high(I_ind),my_h_low(I_ind):my_h_high(I_ind),mz_h_low(I_ind):mz_h_high(I_ind)))))*f_ratio;
            if ~isempty(RNA2_stack)
                RNA2_stack0 = double(RNA2_stack(x_all_low(I_ind):x_all_high(I_ind),y_all_low(I_ind):y_all_high(I_ind),z_all_low(I_ind):z_all_high(I_ind)));
                foci_data{ir}(I_ind,18) = sum(sum(sum(RNA2_stack0.*bw01(x_all_low(I_ind):x_all_high(I_ind),y_all_low(I_ind):y_all_high(I_ind),z_all_low(I_ind):z_all_high(I_ind)).*h{ir}(x_h_low(I_ind):x_h_high(I_ind),y_h_low(I_ind):y_h_high(I_ind),z_h_low(I_ind):z_h_high(I_ind)))));
            end
        else
            signal_stack0 = double(signal_stack(x_all_low(I_ind):x_all_high(I_ind),y_all_low(I_ind):y_all_high(I_ind),z_all_low(I_ind):z_all_high(I_ind)));
            RNA_stack0 = double(RNA_stack(x_all_low(I_ind):x_all_high(I_ind),y_all_low(I_ind):y_all_high(I_ind),z_all_low(I_ind):z_all_high(I_ind)));
            signal_stack1 = double(signal_stack(mx_all_low(I_ind):mx_all_high(I_ind),my_all_low(I_ind):my_all_high(I_ind),mz_all_low(I_ind):mz_all_high(I_ind)));
            
            foci_data{ir}(I_ind,4) = sum(sum(sum(signal_stack0.*bw01(x_all_low(I_ind):x_all_high(I_ind),y_all_low(I_ind):y_all_high(I_ind),z_all_low(I_ind):z_all_high(I_ind)).*h{ir}(x_h_low(I_ind):x_h_high(I_ind),y_h_low(I_ind):y_h_high(I_ind),z_h_low(I_ind):z_h_high(I_ind)))));
            foci_data{ir}(I_ind,5) = sum(sum(sum(   RNA_stack0.*bw01(x_all_low(I_ind):x_all_high(I_ind),y_all_low(I_ind):y_all_high(I_ind),z_all_low(I_ind):z_all_high(I_ind)).*h{ir}(x_h_low(I_ind):x_h_high(I_ind),y_h_low(I_ind):y_h_high(I_ind),z_h_low(I_ind):z_h_high(I_ind)))));
            foci_data{ir}(I_ind,9) = sum(sum(sum(bw01(x_all_low(I_ind):x_all_high(I_ind),y_all_low(I_ind):y_all_high(I_ind),z_all_low(I_ind):z_all_high(I_ind)).*h{ir}(x_h_low(I_ind):x_h_high(I_ind),y_h_low(I_ind):y_h_high(I_ind),z_h_low(I_ind):z_h_high(I_ind)))))/sum(sum(sum(h{ir}(x_h_low(I_ind):x_h_high(I_ind),y_h_low(I_ind):y_h_high(I_ind),:))));
            foci_data{ir}(I_ind,10) = sum(sum(sum(bw01(mx_all_low(I_ind):mx_all_high(I_ind),my_all_low(I_ind):my_all_high(I_ind),mz_all_low(I_ind):mz_all_high(I_ind)).*mh{ir}(mx_h_low(I_ind):mx_h_high(I_ind),my_h_low(I_ind):my_h_high(I_ind),mz_h_low(I_ind):mz_h_high(I_ind)))))/sum(sum(sum(mh{ir}(mx_h_low(I_ind):mx_h_high(I_ind),my_h_low(I_ind):my_h_high(I_ind),:))));
            f_ratio = sum(sum(sum(bw01(x_all_low(I_ind):x_all_high(I_ind),y_all_low(I_ind):y_all_high(I_ind),z_all_low(I_ind):z_all_high(I_ind)).*h{ir}(x_h_low(I_ind):x_h_high(I_ind),y_h_low(I_ind):y_h_high(I_ind),z_h_low(I_ind):z_h_high(I_ind)))))/sum(sum(sum(bw01(mx_all_low(I_ind):mx_all_high(I_ind),my_all_low(I_ind):my_all_high(I_ind),mz_all_low(I_ind):mz_all_high(I_ind)).*mh{ir}(mx_h_low(I_ind):mx_h_high(I_ind),my_h_low(I_ind):my_h_high(I_ind),mz_h_low(I_ind):mz_h_high(I_ind)))));
            foci_data{ir}(I_ind,3) = sum(sum(sum(signal_stack1.*bw01(mx_all_low(I_ind):mx_all_high(I_ind),my_all_low(I_ind):my_all_high(I_ind),mz_all_low(I_ind):mz_all_high(I_ind)).*mh{ir}(mx_h_low(I_ind):mx_h_high(I_ind),my_h_low(I_ind):my_h_high(I_ind),mz_h_low(I_ind):mz_h_high(I_ind)))))*f_ratio;
            if ~isempty(RNA2_stack)
                RNA2_stack0 = double(RNA2_stack(x_all_low(I_ind):x_all_high(I_ind),y_all_low(I_ind):y_all_high(I_ind),z_all_low(I_ind):z_all_high(I_ind)));
                foci_data{ir}(I_ind,18) = sum(sum(sum(RNA2_stack0.*bw01(x_all_low(I_ind):x_all_high(I_ind),y_all_low(I_ind):y_all_high(I_ind),z_all_low(I_ind):z_all_high(I_ind)).*h{ir}(x_h_low(I_ind):x_h_high(I_ind),y_h_low(I_ind):y_h_high(I_ind),z_h_low(I_ind):z_h_high(I_ind)))));
            end
        end
%         nu_mask_temp = (mask_stack(:,:,z_all_low(I_ind):z_all_high(I_ind)) == mask_stack(foci_bw00(I_ind))) & bw01(:,:,z_all_low(I_ind):z_all_high(I_ind));
%         foci_data{ir}(I_ind,12) = sum(sum(sum(signal_stack(:,:,z_all_low(I_ind):z_all_high(I_ind)).*nu_mask_temp)))./sum(sum(sum(nu_mask_temp)));
%         nu_mask_temp(x_all_low(I_ind):x_all_high(I_ind),y_all_low(I_ind):y_all_high(I_ind),:) = nu_mask_temp(x_all_low(I_ind):x_all_high(I_ind),y_all_low(I_ind):y_all_high(I_ind),:) & (~h{ir}(x_h_low(I_ind):x_h_high(I_ind),y_h_low(I_ind):y_h_high(I_ind),z_h_low(I_ind):z_h_high(I_ind)));
%         foci_data{ir}(I_ind,13) = sum(sum(sum(signal_stack(:,:,z_all_low(I_ind):z_all_high(I_ind)).*nu_mask_temp)))./sum(sum(sum(nu_mask_temp)));
        foci_data{ir}(I_ind,12) = mean_signal2D(z_all(I_ind),mask_stack(foci_bw00(I_ind)));
        foci_data{ir}(I_ind,13) = fn_mean_clean(I_ind);
        foci_data{ir}(I_ind,14) = fn_var_clean(I_ind);
        
        exp_ind = (ii0 >= x_h_low(I_ind)) & (ii0 <= x_h_high(I_ind)) & (jj0 >= y_h_low(I_ind)) & (jj0 <= y_h_high(I_ind)) & (kk0 >= z_h_low(I_ind)) & (kk0 <= z_h_high(I_ind));
        foci_data{ir}(I_ind,15) = sum(exp(-(ii0(exp_ind)-x_h-1).^2/2/sigmai^2-(jj0(exp_ind)-y_h-1).^2/2/sigmaj^2-(kk0(exp_ind)-z_h-1).^2/2/sigmak^2))/Integral0;

        foci_data{ir}(I_ind,16) = sum(exp(-(kk0(exp_ind)-z_h-1).^2/2/sigmak^2))/Integralz0/length(kk0);
        foci_data{ir}(I_ind,17) = nucleus_prop_Area2(I_ind);
%         foci_data{ir}(I_ind,16) = fn_mean_clean3(I_ind);
%         foci_data{ir}(I_ind,17) = fn_var_clean3(I_ind);
    end
    
%%%
%%%

    [x_all,y_all,z_all] = ind2sub(size(mask_stack),fake_bw00);
%%%%%
    x_all_low = max(x_all-x_h,1);
    x_all_high = min(x_all+x_h,x_lim);
    y_all_low = max(y_all-y_h,1);
    y_all_high = min(y_all+y_h,y_lim);
    z_all_low = max(z_all-z_h,1);
    z_all_high = min(z_all+z_h,z_lim);
%%%%%
    x_h_low = 1+x_h+x_all_low-x_all;
    x_h_high = 1+x_h+x_all_high-x_all;
    y_h_low = 1+y_h+y_all_low-y_all;
    y_h_high = 1+y_h+y_all_high-y_all;
    z_h_low = 1+z_h+z_all_low-z_all;
    z_h_high = 1+z_h+z_all_high-z_all;
%%%%%
    mx_all_low = max(x_all-mx_h,1);
    mx_all_high = min(x_all+mx_h,x_lim);
    my_all_low = max(y_all-my_h,1);
    my_all_high = min(y_all+my_h,y_lim);
    mz_all_low = max(z_all-mz_h,1);
    mz_all_high = min(z_all+mz_h,z_lim);
%%%%%
    mx_h_low = 1+mx_h+mx_all_low-x_all;
    mx_h_high = 1+mx_h+mx_all_high-x_all;
    my_h_low = 1+my_h+my_all_low-y_all;
    my_h_high = 1+my_h+my_all_high-y_all;
    mz_h_low = 1+mz_h+mz_all_low-z_all;
    mz_h_high = 1+mz_h+mz_all_high-z_all;
%%%%%

    if foci_reuse
        x0 = -dxy(:,1)+x_all_low;
        y0 = -dxy(:,2)+y_all_low;
        x1_min = max(1,-dxy(:,1)+x_all_low);
        x1_max = min(xy_max(1),-dxy(:,1)+x_all_high);
        y1_min = max(1,-dxy(:,2)+y_all_low);
        y1_max = min(xy_max(2),-dxy(:,2)+y_all_high);
        x0_min = x1_min-x0+1;
        x0_max = x1_max-x0+1;
        y0_min = y1_min-y0+1;
        y0_max = y1_max-y0+1;

        mx0 = -dxy(:,1)+mx_all_low;
        my0 = -dxy(:,2)+my_all_low;
        mx1_min = max(1,-dxy(:,1)+mx_all_low);
        mx1_max = min(xy_max(1),-dxy(:,1)+mx_all_high);
        my1_min = max(1,-dxy(:,2)+my_all_low);
        my1_max = min(xy_max(2),-dxy(:,2)+my_all_high);
        mx0_min = mx1_min-mx0+1;
        mx0_max = mx1_max-mx0+1;
        my0_min = my1_min-my0+1;
        my0_max = my1_max-my0+1;
    end

    [ii0,jj0,kk0] = ind2sub(size(h{ir}),find(h{ir})); 
%     if size(h{ir},3) == 1
%         Integral0 = 2*pi*sigmai*sigmaj;
%     else
        Integral0 = 2*pi*sigmai*sigmaj*sqrt(2*pi)*sigmak;
        Integralz0 = sqrt(2*pi)*sigmak;
%     end
    
    for I_ind = 1:length(x_all)
        if foci_reuse
            signal_stack0 = double(signal_stack(x_all_low(I_ind):x_all_high(I_ind),y_all_low(I_ind):y_all_high(I_ind),z_all_low(I_ind):z_all_high(I_ind)));
            signal_stack0(x0_min(I_ind):x0_max(I_ind),y0_min(I_ind):y0_max(I_ind),:) = signal_stack(x1_min(I_ind):x1_max(I_ind),y1_min(I_ind):y1_max(I_ind),z_all_low(I_ind):z_all_high(I_ind));

            RNA_stack0 = double(RNA_stack(x_all_low(I_ind):x_all_high(I_ind),y_all_low(I_ind):y_all_high(I_ind),z_all_low(I_ind):z_all_high(I_ind)));
            
            signal_stack1 = double(signal_stack(mx_all_low(I_ind):mx_all_high(I_ind),my_all_low(I_ind):my_all_high(I_ind),mz_all_low(I_ind):mz_all_high(I_ind)));
            signal_stack1(mx0_min(I_ind):mx0_max(I_ind),my0_min(I_ind):my0_max(I_ind),:) = signal_stack(mx1_min(I_ind):mx1_max(I_ind),my1_min(I_ind):my1_max(I_ind),mz_all_low(I_ind):mz_all_high(I_ind));
            
            fake_data{ir}(I_ind,4) = sum(sum(sum(signal_stack0.*bw01(x_all_low(I_ind):x_all_high(I_ind),y_all_low(I_ind):y_all_high(I_ind),z_all_low(I_ind):z_all_high(I_ind)).*h{ir}(x_h_low(I_ind):x_h_high(I_ind),y_h_low(I_ind):y_h_high(I_ind),z_h_low(I_ind):z_h_high(I_ind)))));
            fake_data{ir}(I_ind,5) = sum(sum(sum(   RNA_stack0.*bw01(x_all_low(I_ind):x_all_high(I_ind),y_all_low(I_ind):y_all_high(I_ind),z_all_low(I_ind):z_all_high(I_ind)).*h{ir}(x_h_low(I_ind):x_h_high(I_ind),y_h_low(I_ind):y_h_high(I_ind),z_h_low(I_ind):z_h_high(I_ind)))));
            fake_data{ir}(I_ind,9) = sum(sum(sum(bw01(x_all_low(I_ind):x_all_high(I_ind),y_all_low(I_ind):y_all_high(I_ind),z_all_low(I_ind):z_all_high(I_ind)).*h{ir}(x_h_low(I_ind):x_h_high(I_ind),y_h_low(I_ind):y_h_high(I_ind),z_h_low(I_ind):z_h_high(I_ind)))))/sum(sum(sum(h{ir}(x_h_low(I_ind):x_h_high(I_ind),y_h_low(I_ind):y_h_high(I_ind),:))));
            fake_data{ir}(I_ind,10) = sum(sum(sum(bw01(mx_all_low(I_ind):mx_all_high(I_ind),my_all_low(I_ind):my_all_high(I_ind),mz_all_low(I_ind):mz_all_high(I_ind)).*mh{ir}(mx_h_low(I_ind):mx_h_high(I_ind),my_h_low(I_ind):my_h_high(I_ind),mz_h_low(I_ind):mz_h_high(I_ind)))))/sum(sum(sum(mh{ir}(mx_h_low(I_ind):mx_h_high(I_ind),my_h_low(I_ind):my_h_high(I_ind),:))));
            f_ratio = sum(sum(sum(bw01(x_all_low(I_ind):x_all_high(I_ind),y_all_low(I_ind):y_all_high(I_ind),z_all_low(I_ind):z_all_high(I_ind)).*h{ir}(x_h_low(I_ind):x_h_high(I_ind),y_h_low(I_ind):y_h_high(I_ind),z_h_low(I_ind):z_h_high(I_ind)))))/sum(sum(sum(bw01(mx_all_low(I_ind):mx_all_high(I_ind),my_all_low(I_ind):my_all_high(I_ind),mz_all_low(I_ind):mz_all_high(I_ind)).*mh{ir}(mx_h_low(I_ind):mx_h_high(I_ind),my_h_low(I_ind):my_h_high(I_ind),mz_h_low(I_ind):mz_h_high(I_ind)))));
            if ~isempty(RNA2_stack)
                RNA2_stack0 = double(RNA2_stack(x_all_low(I_ind):x_all_high(I_ind),y_all_low(I_ind):y_all_high(I_ind),z_all_low(I_ind):z_all_high(I_ind)));
                fake_data{ir}(I_ind,18) = sum(sum(sum(RNA2_stack0.*bw01(x_all_low(I_ind):x_all_high(I_ind),y_all_low(I_ind):y_all_high(I_ind),z_all_low(I_ind):z_all_high(I_ind)).*h{ir}(x_h_low(I_ind):x_h_high(I_ind),y_h_low(I_ind):y_h_high(I_ind),z_h_low(I_ind):z_h_high(I_ind)))));
            end
            try
            fake_data{ir}(I_ind,3) = sum(sum(sum(signal_stack1.*bw01(mx_all_low(I_ind):mx_all_high(I_ind),my_all_low(I_ind):my_all_high(I_ind),mz_all_low(I_ind):mz_all_high(I_ind)).*mh{ir}(mx_h_low(I_ind):mx_h_high(I_ind),my_h_low(I_ind):my_h_high(I_ind),mz_h_low(I_ind):mz_h_high(I_ind)))))*f_ratio;
            catch
                error(['ir = ',num2str(ir),', I_ind = ',num2str(I_ind),', size1 = ',num2str(size(signal_stack1)),', size2 = ',num2str(size(bw01(mx_all_low(I_ind):mx_all_high(I_ind),my_all_low(I_ind):my_all_high(I_ind),mz_all_low(I_ind):mz_all_high(I_ind)))),', size3 = ',num2str(size(mh{ir}(mx_h_low(I_ind):mx_h_high(I_ind),my_h_low(I_ind):my_h_high(I_ind),mz_h_low(I_ind):mz_h_high(I_ind))))])
            end
        else
            signal_stack0 = double(signal_stack(x_all_low(I_ind):x_all_high(I_ind),y_all_low(I_ind):y_all_high(I_ind),z_all_low(I_ind):z_all_high(I_ind)));
            RNA_stack0 = double(RNA_stack(x_all_low(I_ind):x_all_high(I_ind),y_all_low(I_ind):y_all_high(I_ind),z_all_low(I_ind):z_all_high(I_ind)));
            signal_stack1 = double(signal_stack(mx_all_low(I_ind):mx_all_high(I_ind),my_all_low(I_ind):my_all_high(I_ind),mz_all_low(I_ind):mz_all_high(I_ind)));
            
            fake_data{ir}(I_ind,4) = sum(sum(sum(signal_stack0.*bw01(x_all_low(I_ind):x_all_high(I_ind),y_all_low(I_ind):y_all_high(I_ind),z_all_low(I_ind):z_all_high(I_ind)).*h{ir}(x_h_low(I_ind):x_h_high(I_ind),y_h_low(I_ind):y_h_high(I_ind),z_h_low(I_ind):z_h_high(I_ind)))));
            fake_data{ir}(I_ind,5) = sum(sum(sum(   RNA_stack0.*bw01(x_all_low(I_ind):x_all_high(I_ind),y_all_low(I_ind):y_all_high(I_ind),z_all_low(I_ind):z_all_high(I_ind)).*h{ir}(x_h_low(I_ind):x_h_high(I_ind),y_h_low(I_ind):y_h_high(I_ind),z_h_low(I_ind):z_h_high(I_ind)))));
            fake_data{ir}(I_ind,9) = sum(sum(sum(bw01(x_all_low(I_ind):x_all_high(I_ind),y_all_low(I_ind):y_all_high(I_ind),z_all_low(I_ind):z_all_high(I_ind)).*h{ir}(x_h_low(I_ind):x_h_high(I_ind),y_h_low(I_ind):y_h_high(I_ind),z_h_low(I_ind):z_h_high(I_ind)))))/sum(sum(sum(h{ir}(x_h_low(I_ind):x_h_high(I_ind),y_h_low(I_ind):y_h_high(I_ind),:))));
            fake_data{ir}(I_ind,10) = sum(sum(sum(bw01(mx_all_low(I_ind):mx_all_high(I_ind),my_all_low(I_ind):my_all_high(I_ind),mz_all_low(I_ind):mz_all_high(I_ind)).*mh{ir}(mx_h_low(I_ind):mx_h_high(I_ind),my_h_low(I_ind):my_h_high(I_ind),mz_h_low(I_ind):mz_h_high(I_ind)))))/sum(sum(sum(mh{ir}(mx_h_low(I_ind):mx_h_high(I_ind),my_h_low(I_ind):my_h_high(I_ind),:))));
            f_ratio = sum(sum(sum(bw01(x_all_low(I_ind):x_all_high(I_ind),y_all_low(I_ind):y_all_high(I_ind),z_all_low(I_ind):z_all_high(I_ind)).*h{ir}(x_h_low(I_ind):x_h_high(I_ind),y_h_low(I_ind):y_h_high(I_ind),z_h_low(I_ind):z_h_high(I_ind)))))/sum(sum(sum(bw01(mx_all_low(I_ind):mx_all_high(I_ind),my_all_low(I_ind):my_all_high(I_ind),mz_all_low(I_ind):mz_all_high(I_ind)).*mh{ir}(mx_h_low(I_ind):mx_h_high(I_ind),my_h_low(I_ind):my_h_high(I_ind),mz_h_low(I_ind):mz_h_high(I_ind)))));
            fake_data{ir}(I_ind,3) = sum(sum(sum(signal_stack1.*bw01(mx_all_low(I_ind):mx_all_high(I_ind),my_all_low(I_ind):my_all_high(I_ind),mz_all_low(I_ind):mz_all_high(I_ind)).*mh{ir}(mx_h_low(I_ind):mx_h_high(I_ind),my_h_low(I_ind):my_h_high(I_ind),mz_h_low(I_ind):mz_h_high(I_ind)))))*f_ratio;
            if ~isempty(RNA2_stack)
                RNA2_stack0 = double(RNA2_stack(x_all_low(I_ind):x_all_high(I_ind),y_all_low(I_ind):y_all_high(I_ind),z_all_low(I_ind):z_all_high(I_ind)));
                fake_data{ir}(I_ind,18) = sum(sum(sum(RNA2_stack0.*bw01(x_all_low(I_ind):x_all_high(I_ind),y_all_low(I_ind):y_all_high(I_ind),z_all_low(I_ind):z_all_high(I_ind)).*h{ir}(x_h_low(I_ind):x_h_high(I_ind),y_h_low(I_ind):y_h_high(I_ind),z_h_low(I_ind):z_h_high(I_ind)))));
            end
        end
%         nu_mask_temp = (mask_stack(:,:,z_all_low(I_ind):z_all_high(I_ind)) == mask_stack(fake_bw00(I_ind))) & bw01(:,:,z_all_low(I_ind):z_all_high(I_ind));
%         fake_data{ir}(I_ind,12) = sum(sum(sum(signal_stack(:,:,z_all_low(I_ind):z_all_high(I_ind)).*nu_mask_temp)))./sum(sum(sum(nu_mask_temp)));
%         nu_mask_temp(x_all_low(I_ind):x_all_high(I_ind),y_all_low(I_ind):y_all_high(I_ind),:) = nu_mask_temp(x_all_low(I_ind):x_all_high(I_ind),y_all_low(I_ind):y_all_high(I_ind),:) & (~h{ir}(x_h_low(I_ind):x_h_high(I_ind),y_h_low(I_ind):y_h_high(I_ind),z_h_low(I_ind):z_h_high(I_ind)));
%         fake_data{ir}(I_ind,13) = sum(sum(sum(signal_stack(:,:,z_all_low(I_ind):z_all_high(I_ind)).*nu_mask_temp)))./sum(sum(sum(nu_mask_temp)));
        fake_data{ir}(I_ind,12) = mean_signal2D(z_all(I_ind),mask_stack(fake_bw00(I_ind)));
        fake_data{ir}(I_ind,13) = fn_mean_clean(I_ind);
        fake_data{ir}(I_ind,14) = fn_var_clean(I_ind);
        
        exp_ind = (ii0 >= x_h_low(I_ind)) & (ii0 <= x_h_high(I_ind)) & (jj0 >= y_h_low(I_ind)) & (jj0 <= y_h_high(I_ind)) & (kk0 >= z_h_low(I_ind)) & (kk0 <= z_h_high(I_ind));
        fake_data{ir}(I_ind,15) = sum(exp(-(ii0(exp_ind)-x_h-1).^2/2/sigmai^2-(jj0(exp_ind)-y_h-1).^2/2/sigmaj^2-(kk0(exp_ind)-z_h-1).^2/2/sigmak^2))/Integral0;

        fake_data{ir}(I_ind,16) = sum(exp(-(kk0(exp_ind)-z_h-1).^2/2/sigmak^2))/Integralz0/length(kk0);
        fake_data{ir}(I_ind,17) = nucleus_prop_Area2(I_ind);
%         fake_data{ir}(I_ind,16) = fn_mean_clean3(I_ind);
%         fake_data{ir}(I_ind,17) = fn_var_clean3(I_ind);
    end
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Calculate mean protein intensity at different plains: %%%%%%%%%%%%%%%%%%
%     mean3_nu = zeros(size(foci_bw00));
%     for I_layer = 1:size(signal_stack,3)
%         nu_prop = regionprops(imerode(bwi, strel('disk',r_size(ir)+r_edge)).*bwn,conv2(signal_stack(:,:,I_layer),h{ir},'same'),'MeanIntensity');
%         temp_nu = [0,[nu_prop.MeanIntensity]];
%         mean3_nu((foci_layer == I_layer)|(fake_layer == I_layer)) = temp_nu(bwn((foci_layer == I_layer)|(fake_layer == I_layer))+1);
%     end
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Collect foci data (foci intensity, foci local protein concentration, mean protein concentration)
    if length(foci_bw00) > 0
        foci_data{ir}(:,1) = max_image00(ismember(foci_raw,foci_bw00));   %%% foci intensity
        foci_data{ir}(:,2) = nucleus_protein_profile((mask_stack(foci_bw00)+1),2);   %%% protein mean concentration
        foci_data{ir}(:,6) = mask_stack(foci_bw00);   %%% nucleus label
        foci_data{ir}(:,7) = nucleus_protein_profile((mask_stack(foci_bw00)+1),1);   %%% nucleus EL length;
        foci_data{ir}(:,8) = nucleus_protein_profile((mask_stack(foci_bw00)+1),3);   %%% nucleus max z slice;
        foci_data{ir}(:,11) = foci_bw00;   %%% foci linear coordinate;
    end
    

    if length(fake_bw00) > 0
        fake_data{ir}(:,1) = 0;
        fake_data{ir}(:,2) = nucleus_protein_profile((mask_stack(fake_bw00)+1),2);
        fake_data{ir}(:,6) = mask_stack(fake_bw00);   %%% nucleus label
        fake_data{ir}(:,7) = nucleus_protein_profile((mask_stack(fake_bw00)+1),1);   %%% nucleus EL length;
        fake_data{ir}(:,8) = nucleus_protein_profile((mask_stack(fake_bw00)+1),3);   %%% nucleus max z slice;
        fake_data{ir}(:,11) = fake_bw00;   %%% fake foci linear coordinate;
    end
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    foci_data{ir} = foci_data{ir}(~isnan(foci_data{ir}(:,2)),:);
    fake_data{ir} = fake_data{ir}(~isnan(fake_data{ir}(:,2)),:);
    
    p_ratio0{ir} = foci_data{ir}(:,4)./foci_data{ir}(:,3)-1;
    fp_ratio0{ir} = fake_data{ir}(:,4)./fake_data{ir}(:,3)-1;
    op_ratio0{ir} = foci_data{ir}(:,4)-foci_data{ir}(:,3);
    ofp_ratio0{ir} = fake_data{ir}(:,4)-fake_data{ir}(:,3);

    %h{ir} = h{ir}(:,:,1);
    %mh{ir} = mh{ir}(:,:,1);
end

% matlabpool close

%% Figure plotting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
color_RGB0 = [1,0,0;0,1,0;0,0,1;0,1,1;1,0,1;1,1,0];
color_RGB = [color_RGB0;[0,0,0];color_RGB0/2;[0.5,0.5,0.5]];
if length(r_size) > size(color_RGB0,1)
    LL2 = ceil(length(r_size)/size(color_RGB0,1));
    color_RGB  = repmat(color_RGB,LL2,1);
end    

figure(figure_list(1))
clf
subplot(1,2,1)
legend_text = cell(0);
for ir = 1:length(r_size)
    plot(p_ratio0{ir},foci_data{ir}(:,1),'o','color',color_RGB(ir,:))
    hold on
    [xtemp,IX] = sort(p_ratio0{ir});
    ytemp = foci_data{ir}(IX,1);
    N_bin0 = floor(N_bin/rbin);
    N_temp = floor(length(xtemp)/N_bin0);
    x_f = zeros(1,N_bin0);
    y_f = zeros(1,N_bin0);
    x_ferr = zeros(1,N_bin0);
    y_ferr = zeros(1,N_bin0);
    for I_bin = 1:N_bin0
        x_f(I_bin) = mean(xtemp(((I_bin-1)*N_temp+1):I_bin*N_temp));
        x_ferr(I_bin) = std0(xtemp(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y_f(I_bin) = mean(ytemp(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y_ferr(I_bin) = std0(ytemp(((I_bin-1)*N_temp+1):I_bin*N_temp));
    end
    errorbarxy(x_f,y_f,x_ferr,y_ferr,x_ferr,y_ferr,color_RGB(ir,:),color_RGB(ir,:))
    legend_text = cat(2,legend_text,{['local intensity, r = ',num2str(r_size(ir))],['binned data, r = ',num2str(r_size(ir))]});
end
title(['Foci intensity vs foci local protein enrichment: ',image_folder,', cycle = ',num2str(N_cycle)],'Interpreter','none')
xlabel('Protein enrichment')
ylabel('Foci intensity (A.U.)')
legend(legend_text)

    
subplot(1,2,2)
legend_text = cell(0);
for ir = 1:length(r_size)
    plot(foci_data{ir}(:,2),p_ratio0{ir},'*','color',color_RGB(2*ir-1,:))
    hold on
    plot(fake_data{ir}(:,2),fp_ratio0{ir},'.','color',color_RGB(2*ir,:))
    
    [xtemp1,IX1] = sort(foci_data{ir}(:,2));
    [xtemp2,IX2] = sort(fake_data{ir}(:,2));
    ytemp1 = p_ratio0{ir}(IX1);
    ytemp2 = fp_ratio0{ir}(IX2);
    ytemp3 = foci_data{ir}(IX1,4)./fake_data{ir}(IX2,4)-1;
    temp_ind = (~isnan(ytemp1)) & (~isnan(ytemp2)) & (~isnan(ytemp3));
    xtemp1 = xtemp1(temp_ind);
    xtemp2 = xtemp2(temp_ind);
    ytemp1 = ytemp1(temp_ind);
    ytemp2 = ytemp2(temp_ind);
    ytemp3 = ytemp3(temp_ind);
    N_temp = floor(length(xtemp1)/N_bin);
    x1_f = zeros(1,N_bin);
    x2_f = zeros(1,N_bin);
    y1_f = zeros(1,N_bin);
    y2_f = zeros(1,N_bin);
    y3_f = zeros(1,N_bin);
    x1_ferr = zeros(1,N_bin);
    x2_ferr = zeros(1,N_bin);
    y1_ferr = zeros(1,N_bin);
    y2_ferr = zeros(1,N_bin);
    y3_ferr = zeros(1,N_bin);

    for I_bin = 1:N_bin
        x1_f(I_bin) = mean(xtemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
        x1_ferr(I_bin) = std0(xtemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
        x2_f(I_bin) = mean(xtemp2(((I_bin-1)*N_temp+1):I_bin*N_temp));
        x2_ferr(I_bin) = std0(xtemp2(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y1_f(I_bin) = mean(ytemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y1_ferr(I_bin) = std0(ytemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y2_f(I_bin) = mean(ytemp2(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y2_ferr(I_bin) = std0(ytemp2(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y3_f(I_bin) = mean(ytemp3(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y3_ferr(I_bin) = std0(ytemp3(((I_bin-1)*N_temp+1):I_bin*N_temp));
    end
    errorbarxy(x1_f,y1_f,x1_ferr,y1_ferr,x1_ferr,y1_ferr,color_RGB(2*ir-1,:),color_RGB(2*ir-1,:),'o')
    errorbarxy(x2_f,y2_f,x2_ferr,y2_ferr,x2_ferr,y2_ferr,color_RGB(2*ir,:),color_RGB(2*ir,:),'*')
    errorbarxy(x1_f,y3_f,x1_ferr,y3_ferr,x1_ferr,y3_ferr,mean(color_RGB((2*ir-1):(2*ir),:)),mean(color_RGB((2*ir-1):(2*ir),:)),'--')
%     plot(x1_f,y1_f-y2_f,'--','color',mean(color_RGB((2*ir-1):(2*ir),:)))
    
    x_limit = xlim;
    %%% Nonlinear fitting 1: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [~,I_middle] = min(abs(y1_f-(mean(y1_f((end-Nmean+1):end))+mean(y1_f(1:Nmean)))/2));
    beta0 = [2,x1_f(I_middle),(mean(y1_f((end-Nmean+1):end))-mean(y1_f(1:Nmean))),mean(y1_f(1:Nmean))];
    try
        [beta1,r,~,~,~] = nlinfit(x1_f,y1_f,@Hill,beta0);
        x_fit = x1_f(1):(x1_f(end)-x1_f(1))/N_fit:x1_f(end);
        plot(x_fit,Hill(beta1,x_fit),'--','color',color_RGB(2*ir-1,:))
    catch err
        plot(x_limit,[0,0],'--','color',color_RGB(2*ir-1,:))
        beta1 = [nan,nan,nan,nan];
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%% Nonlinear fitting 2: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    try
        [beta2,r,~,~,~] = nlinfit(xtemp1,ytemp1,@Hill,beta0);
        x_fit = xtemp1(1):(xtemp1(end)-xtemp1(1))/N_fit:xtemp1(end);
        plot(x_fit,Hill(beta2,x_fit),'-','color',color_RGB(2*ir-1,:))
    catch err
        plot(x_limit,[0,0],'-','color',color_RGB(2*ir-1,:))
        beta2 = [nan,nan,nan,nan];
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    legend_text = cat(2,legend_text,{['foci spots, r = ',num2str(r_size(ir))],['control spots, r = ',num2str(r_size(ir))],['binned data, r = ',num2str(r_size(ir))],['binned control data, r = ',num2str(r_size(ir))],['difference of binned data, r = ',num2str(r_size(ir))],['Fitting of binned data, r = ',num2str(r_size(ir)),', E1 = ',num2str(beta1(4)),', E2 = ',num2str(beta1(3)+beta1(4)),', th = ',num2str(beta1(2)),', h = ',num2str(beta1(1))],['Fitting of raw data, r = ',num2str(r_size(ir)),', E1 = ',num2str(beta2(4)),', E2 = ',num2str(beta2(3)+beta2(4)),', th = ',num2str(beta2(2)),', h = ',num2str(beta2(1))]});
end
title(['Foci local protein enrichment vs nuclear mean protein intensity: ',image_folder,', cycle = ',num2str(N_cycle)],'Interpreter','none')
ylabel('Protein enrichment')
xlabel(['P_m_e_a_n (',unit2,')'])
legend(legend_text)


figure(figure_list(2))
clf
for ir = 1:length(r_size)
    subplot(2,ceil(length(r_size)/2),ir)
    n_foci = hist(p_ratio0{ir},ppbin);
    n_fake = hist(fp_ratio0{ir},ppbin);
    plot(ppbin,(n_foci/sum(n_foci)*100),ppbin,(n_fake/sum(n_fake)*100));
    title(['Foci local protein enrichment histogram (r = ',num2str(r_size(ir)),'): ',image_folder,', cycle = ',num2str(N_cycle),', foci mean = ',num2str(mean(p_ratio0{ir}((~isnan(p_ratio0{ir}))&(p_ratio0{ir}>0)))),', foci std = ',num2str(std(p_ratio0{ir}((~isnan(p_ratio0{ir}))&(p_ratio0{ir}>0)))),', fake foci mean = ',num2str(mean(fp_ratio0{ir}((~isnan(fp_ratio0{ir}))&(fp_ratio0{ir}>0)))),', fake foci std = ',num2str(std(fp_ratio0{ir}((~isnan(fp_ratio0{ir}))&(fp_ratio0{ir}>0))))],'Interpreter','none')
    ylabel('%')
    xlabel('Protein enrichment')
    legend('foci spots','control spots')
end


figure(figure_list(3))
clf
legend_text = cell(0);
for ir = 1:length(r_size)
    plot(op_ratio0{ir},foci_data{ir}(:,5),'o','color',color_RGB(ir,:))
    hold on
    
    [xtemp,IX] = sort(op_ratio0{ir});
    ytemp = foci_data{ir}(IX,5);
    N_bin0 = floor(N_bin/rbin);
    N_temp = floor(length(xtemp)/N_bin0);
    x_f = zeros(1,N_bin0);
    y_f = zeros(1,N_bin0);
    x_ferr = zeros(1,N_bin0);
    y_ferr = zeros(1,N_bin0);
    for I_bin = 1:N_bin0
        x_f(I_bin) = mean(xtemp(((I_bin-1)*N_temp+1):I_bin*N_temp));
        x_ferr(I_bin) = std0(xtemp(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y_f(I_bin) = mean(ytemp(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y_ferr(I_bin) = std0(ytemp(((I_bin-1)*N_temp+1):I_bin*N_temp));
    end
    errorbarxy(x_f,y_f,x_ferr,y_ferr,x_ferr,y_ferr,color_RGB(ir,:),color_RGB(ir,:))
    
    legend_text = cat(2,legend_text,{['foci spots, r = ',num2str(r_size(ir))],['binned data, r = ',num2str(r_size(ir))]});
end
title(['Foci local protein extra fluorescence vs RNA fluorescence: ',image_folder,', cycle = ',num2str(N_cycle)],'Interpreter','none')
xlabel('P_e_x_t_r_a')
ylabel('Foci fluorescence (A.U.)')
legend(legend_text)


figure(figure_list(4))
clf
subplot(1,2,1)
legend_text = cell(0);
for ir = 1:length(r_size)
    plot(op_ratio0{ir},foci_data{ir}(:,1),'o','color',color_RGB(ir,:))
    hold on
    
    [xtemp,IX] = sort(op_ratio0{ir});
    ytemp = foci_data{ir}(IX,1);
    N_bin0 = floor(N_bin/rbin);
    N_temp = floor(length(xtemp)/N_bin0);
    x_f = zeros(1,N_bin0);
    y_f = zeros(1,N_bin0);
    x_ferr = zeros(1,N_bin0);
    y_ferr = zeros(1,N_bin0);
    for I_bin = 1:N_bin0
        x_f(I_bin) = mean(xtemp(((I_bin-1)*N_temp+1):I_bin*N_temp));
        x_ferr(I_bin) = std0(xtemp(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y_f(I_bin) = mean(ytemp(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y_ferr(I_bin) = std0(ytemp(((I_bin-1)*N_temp+1):I_bin*N_temp));
    end
    errorbarxy(x_f,y_f,x_ferr,y_ferr,x_ferr,y_ferr,color_RGB(ir,:),color_RGB(ir,:))

    legend_text = cat(2,legend_text,{['local intensity, r = ',num2str(r_size(ir))],['binned data, r = ',num2str(r_size(ir))]});
end
title(['Foci local protein intensity vs foci intensity: ',image_folder,', cycle = ',num2str(N_cycle)],'Interpreter','none')
xlabel(['Protein enrichment (',unit1,')'])
ylabel('Foci intensity (A.U.)')
legend(legend_text)

    
subplot(1,2,2)
legend_text = cell(0);
for ir = 1:length(r_size)
    plot(foci_data{ir}(:,2),op_ratio0{ir},'*','color',color_RGB(2*ir-1,:))
    hold on
    plot(fake_data{ir}(:,2),ofp_ratio0{ir},'.','color',color_RGB(2*ir,:))

    [xtemp1,IX1] = sort(foci_data{ir}(:,2));
    [xtemp2,IX2] = sort(fake_data{ir}(:,2));
    ytemp1 = op_ratio0{ir}(IX1);
    ytemp2 = ofp_ratio0{ir}(IX2);
    ytemp3 = foci_data{ir}(IX1,4)-fake_data{ir}(IX2,4);
    temp_ind = (~isnan(ytemp1)) & (~isnan(ytemp2));
    xtemp1 = xtemp1(temp_ind);
    xtemp2 = xtemp2(temp_ind);
    ytemp1 = ytemp1(temp_ind);
    ytemp2 = ytemp2(temp_ind);
    ytemp3 = ytemp3(temp_ind);
    N_temp = floor(length(xtemp1)/N_bin);
    x1_f = zeros(1,N_bin);
    x2_f = zeros(1,N_bin);
    y1_f = zeros(1,N_bin);
    y2_f = zeros(1,N_bin);
    y3_f = zeros(1,N_bin);
    x1_ferr = zeros(1,N_bin);
    x2_ferr = zeros(1,N_bin);
    y1_ferr = zeros(1,N_bin);
    y2_ferr = zeros(1,N_bin);
    y3_ferr = zeros(1,N_bin);
    
    for I_bin = 1:N_bin
        x1_f(I_bin) = mean(xtemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
        x1_ferr(I_bin) = std0(xtemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
        x2_f(I_bin) = mean(xtemp2(((I_bin-1)*N_temp+1):I_bin*N_temp));
        x2_ferr(I_bin) = std0(xtemp2(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y1_f(I_bin) = mean(ytemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y1_ferr(I_bin) = std0(ytemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y2_f(I_bin) = mean(ytemp2(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y2_ferr(I_bin) = std0(ytemp2(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y3_f(I_bin) = mean(ytemp3(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y3_ferr(I_bin) = std0(ytemp3(((I_bin-1)*N_temp+1):I_bin*N_temp));
    end
    errorbarxy(x1_f,y1_f,x1_ferr,y1_ferr,x1_ferr,y1_ferr,color_RGB(2*ir-1,:),color_RGB(2*ir-1,:),'o')
    errorbarxy(x2_f,y2_f,x2_ferr,y2_ferr,x2_ferr,y2_ferr,color_RGB(2*ir,:),color_RGB(2*ir,:),'*')
    errorbarxy(x1_f,y3_f,x1_ferr,y3_ferr,x1_ferr,y3_ferr,mean(color_RGB((2*ir-1):(2*ir),:)),mean(color_RGB((2*ir-1):(2*ir),:)),'--')
%     plot(x1_f,y1_f-y2_f,'--','color',mean(color_RGB((2*ir-1):(2*ir),:)))

    x_limit = xlim;
    %%% Nonlinear fitting 1: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [~,I_middle] = min(abs(y1_f-(mean(y1_f((end-Nmean+1):end))+mean(y1_f(1:Nmean)))/2));
    beta0 = [2,x1_f(I_middle),(mean(y1_f((end-Nmean+1):end))-mean(y1_f(1:Nmean))),mean(y1_f(1:Nmean))];
    try
        [beta1,r,~,~,~] = nlinfit(x1_f,y1_f,@Hill,beta0);
        x_fit = x1_f(1):(x1_f(end)-x1_f(1))/N_fit:x1_f(end);
        plot(x_fit,Hill(beta1,x_fit),'--','color',color_RGB(2*ir-1,:))
    catch err
        plot(x_limit,[0,0],'--','color',color_RGB(2*ir-1,:))
        beta1 = [nan,nan,nan,nan];
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%% Nonlinear fitting 2: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    try
        [beta2,r,~,~,~] = nlinfit(xtemp1,ytemp1,@Hill,beta0);
        x_fit = xtemp1(1):(xtemp1(end)-xtemp1(1))/N_fit:xtemp1(end);
        plot(x_fit,Hill(beta2,x_fit),'-','color',color_RGB(2*ir-1,:))
    catch err
        plot(x_limit,[0,0],'-','color',color_RGB(2*ir-1,:))
        beta2 = [nan,nan,nan,nan];
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    legend_text = cat(2,legend_text,{['foci spots, r = ',num2str(r_size(ir))],['control spots, r = ',num2str(r_size(ir))],['binned data, r = ',num2str(r_size(ir))],['binned control data, r = ',num2str(r_size(ir))],['difference of binned data, r = ',num2str(r_size(ir))],['Fitting of binned data, r = ',num2str(r_size(ir)),', E1 = ',num2str(beta1(4)),', E2 = ',num2str(beta1(3)+beta1(4)),', th = ',num2str(beta1(2)),', h = ',num2str(beta1(1))],['Fitting of raw data, r = ',num2str(r_size(ir)),', E1 = ',num2str(beta2(4)),', E2 = ',num2str(beta2(3)+beta2(4)),', th = ',num2str(beta2(2)),', h = ',num2str(beta2(1))]});
end
title(['Foci local protein intensity vs nuclear mean protein intensity: ',image_folder,', cycle = ',num2str(N_cycle)],'Interpreter','none')
ylabel(['Protein enrichment (',unit1,')'])
xlabel(['P_m_e_a_n (',unit2,')'])
legend(legend_text)


figure(figure_list(5))
clf
for ir = 1:length(r_size)
    subplot(2,ceil(length(r_size)/2),ir)
    [n_foci,xout_foci] = hist(op_ratio0{ir},Npp);
    [n_fake,xout_fake] = hist(ofp_ratio0{ir},Npp);
    plot(xout_foci,(n_foci/sum(n_foci)*100),xout_fake,(n_fake/sum(n_fake)*100));
    title(['Foci local protein intensity histogram (r = ',num2str(r_size(ir)),'): ',image_folder,', cycle = ',num2str(N_cycle),', foci mean = ',num2str(mean(op_ratio0{ir}((~isnan(op_ratio0{ir}))&(op_ratio0{ir}>0)))),', foci std = ',num2str(std(op_ratio0{ir}((~isnan(op_ratio0{ir}))&(op_ratio0{ir}>0)))),', fake foci mean = ',num2str(mean(ofp_ratio0{ir}((~isnan(ofp_ratio0{ir}))&(ofp_ratio0{ir}>0)))),', fake foci std = ',num2str(std(ofp_ratio0{ir}((~isnan(ofp_ratio0{ir}))&(ofp_ratio0{ir}>0))))],'Interpreter','none')
    ylabel('%')
    xlabel(['Protein enrichment (',unit1,')'])
    legend('foci spots','control spots')
end

figure(figure_list(6))
clf
    mean_p = zeros(size(r_size));
    mean_fp = zeros(size(r_size));
    mean_d = zeros(size(r_size));
    std_p = zeros(size(r_size));
    std_fp = zeros(size(r_size));
    std_d = zeros(size(r_size));
    for ir = 1:length(r_size)
%         [~,IX1] = sort(foci_data{ir}(:,2));
%         [~,IX2] = sort(fake_data{ir}(:,2));
%         ytemp1 = p_ratio0{ir}(IX1);
%         ytemp2 = fp_ratio0{ir}(IX2);
%         temp_ind = (~isnan(ytemp1)) & (~isnan(ytemp2));
%         mean_p(ir) = mean(ytemp1(temp_ind));
%         mean_fp(ir) = mean(ytemp2(temp_ind));
        temp_ind = (~isnan(foci_data{ir}(:,4)./fake_data{ir}(:,4)));
        mean_p(ir) = mean(p_ratio0{ir}(~isnan(p_ratio0{ir})));
        mean_fp(ir) = mean(fp_ratio0{ir}(~isnan(fp_ratio0{ir})));
        mean_d(ir) = mean(foci_data{ir}(temp_ind,4)./fake_data{ir}(temp_ind,4)-1);
        std_p(ir) = std0(p_ratio0{ir}(~isnan(p_ratio0{ir})));
        std_fp(ir) = std0(fp_ratio0{ir}(~isnan(fp_ratio0{ir})));
        std_d(ir) = std0(foci_data{ir}(temp_ind,4)./fake_data{ir}(temp_ind,4)-1);
    end
    hold on
    errorbar(r_size,mean_p,std_p,'go-');
    errorbar(r_size,mean_fp,std_fp,'b*-');
    errorbar(r_size,mean_d,std_d,'r--');
%     plot(r_size,mean_p-mean_fp,'r--')
    title(['Mean relative enrichment vs integration radius: ',image_folder,', cycle = ',num2str(N_cycle)],'Interpreter','none')
    ylabel('Protein enrichment')
    xlabel(['Radius (pixel)'])
    legend('foci spots','control spots','Difference')

figure(figure_list(7))
clf
    mean_op = zeros(size(r_size));
    mean_ofp = zeros(size(r_size));
    mean_od = zeros(size(r_size));
    std_op = zeros(size(r_size));
    std_ofp = zeros(size(r_size));
    std_od = zeros(size(r_size));
    for ir = 1:length(r_size)
%         [~,IX1] = sort(foci_data{ir}(:,2));
%         [~,IX2] = sort(fake_data{ir}(:,2));
%         ytemp1 = op_ratio0{ir}(IX1);
%         ytemp2 = ofp_ratio0{ir}(IX2);
%         temp_ind = (~isnan(ytemp1)) & (~isnan(ytemp2));
%         mean_op(ir) = mean(ytemp1(temp_ind));
%         mean_ofp(ir) = mean(ytemp2(temp_ind));
        temp_ind = (~isnan(foci_data{ir}(:,4)-fake_data{ir}(:,4)));
        mean_op(ir) = mean(op_ratio0{ir}(~isnan(op_ratio0{ir})));
        mean_ofp(ir) = mean(ofp_ratio0{ir}(~isnan(ofp_ratio0{ir})));
        mean_od(ir) = mean(foci_data{ir}(temp_ind,4)-fake_data{ir}(temp_ind,4));
        std_op(ir) = std0(op_ratio0{ir}(~isnan(op_ratio0{ir})));
        std_ofp(ir) = std0(ofp_ratio0{ir}(~isnan(ofp_ratio0{ir})));
        std_od(ir) = std0(foci_data{ir}(temp_ind,4)-fake_data{ir}(temp_ind,4));
    end
    hold on
    errorbar(r_size,mean_op,std_op,'go-');
    errorbar(r_size,mean_ofp,std_ofp,'b*-');
    errorbar(r_size,mean_od,std_od,'r--');
%     plot(r_size,mean_op-mean_ofp,'r--')
    title(['Mean absolute enrichment vs integration radius: ',image_folder,', cycle = ',num2str(N_cycle)],'Interpreter','none')
    ylabel(['Protein enrichment (',unit1,')'])
    xlabel(['Radius (pixel)'])
    legend('foci spots','control spots','Difference')






    
function y = std0(x)
    NN = size(x);
    if NN(1) == 1;
        N1 = NN(2);
    else
        N1 = NN(1);
    end
    y = std(x)/sqrt(N1);
        

function y = Hill(beta,x)
 y = beta(4)+beta(3)*x.^beta(1)./(x.^beta(1)+beta(2)^beta(1));
