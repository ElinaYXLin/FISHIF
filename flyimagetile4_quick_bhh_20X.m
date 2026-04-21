function flyimagetile4_quick(lf_name,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Fly embyro image tiling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%clear all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameter setting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
list_sub = 'stacks/';
tile_name = 'tile';
% YY_name = 'Y';
tif_name = 'stack';
single_add = '_new/';
% temp_add = '_temp/';
figure_tail = '.tif';
list_name = 'matchlist.xls';
std_th = 50;

standard_record_confocal = 'Calibration/Results/standard.mat';
% standard_record_60X_confocal = 'Calibration3/Results/standard_60X.mat';
standard_record_60X_confocal = 'Calibration3_01052020/Results/standard_60X.mat';
standard_record_airy = 'Calibration/Results/standard.mat';
standard_record_60X_airy = 'Calibration3_08112019/Results/standard_60X.mat';
standard_record_fast = 'Calibration/Results/standard.mat';
standard_record_60X_fast = 'Calibration3_FA_08112019/Results/standard_60X.mat';
default_obj = '60X';

fast_add = 'F';
airy_add = 'A';
double_add = 'D';
fast_offset = 10000;
delta_fast = [12,0];

overlay_airy = 0.05;
dz_range_airy = 10;
dx_range_airy = 100;
ds_range_airy = 30;
delta_airy = [0,0];

overlay_confocal = 0.05;
dz_range_confocal = 10;
dx_range_confocal = 100;
ds_range_confocal = 30;
delta_confocal = [0,0];

para_name = 'stitch_parameter';
output_tail = '.xls';
mat_tail = '.mat';
image_tail = '*.tif';

old_add = '_new0/';
load_old = false;

if length(varargin) >= 4 && ~isempty(varargin{4})
    p_start = varargin{4}(1);
    p_end = varargin{4}(2);
    if p_start >= p_end
        p_start = 0;
    end
    if p_end > 1
        p_end = 1;
    end
else
    p_start = 0;
    p_end = 1;
end


global standard_data standard_data_60X

if isa(lf_name,'char') && strcmp(lf_name(end-3:end),'.xls')
    list_name0 = lf_name;
    [num_list, folder_list] = xlsread(list_name0);
    folder_list = folder_list(strcmpi('T',folder_list(:,6)),:);
elseif isa(lf_name,'cell')
    folder_list = lf_name;
else
    error('Incorrect input!')
end
[NN1,NN2] = size(folder_list);

for list_I0 = 1:NN1
    folder_name = [folder_list{list_I0,1},list_sub];
    channel_name = eval(folder_list{list_I0,5});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [num_list, file_list] = xlsread([folder_name,list_name]);
    [N1,N2] = size(file_list);
    
    if ~isempty(varargin) && ~isempty(varargin{1})
        list_J_all = varargin{1}{list_I0};
    else
        list_J_all = 1:N1;
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Processing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for list_I = list_J_all%1:N1
        out_name = [file_list{list_I,1}(1:(end-1)),single_add]; %%% output folder name
%         out_name0 = [file_list{list_I,1}(1:(end-1)),temp_add]; %%% temp output folder name
        match_channel = num_list(list_I,4);
        resolution = num_list(list_I,9);
        
        if ~isempty(strfind(file_list{list_I,1},default_obj)) || ~isempty(strfind(file_list{list_I,2},default_obj))
            change_obj = true;
        else
            change_obj = false;
        end
        
        if ~isempty(strfind(file_list{list_I,2},fast_add))
            offset0 = fast_offset;
            delta = delta_fast;
            standard_record = standard_record_fast
            standard_record_60X = standard_record_60X_fast
        elseif ~isempty(strfind(file_list{list_I,2},airy_add))
            offset0 = 0;
            delta = delta_airy;
            standard_record = standard_record_airy
            standard_record_60X = standard_record_60X_airy
        else
            offset0 = 0;
            delta = delta_confocal;
            standard_record = standard_record_confocal
            standard_record_60X = standard_record_60X_confocal
        end
        standard_data = load(standard_record);
        standard_data_60X = load(standard_record_60X);
        
        
        if length(varargin) >= 2 && ~isempty(varargin{2})
            overlayP = varargin{2}{1};
            dz_range = varargin{2}{2};
            dx_range = varargin{2}{3};
            ds_range = varargin{2}{4};
        elseif ~isempty(strfind(file_list{list_I,2},fast_add)) || ~isempty(strfind(file_list{list_I,2},airy_add))
            overlayP = overlay_airy;
            dz_range = dz_range_airy;
            dx_range = dx_range_airy;
            ds_range = ds_range_airy;
        else
            overlayP = overlay_confocal;
            dz_range = dz_range_confocal;
            dx_range = dx_range_confocal;
            ds_range = ds_range_confocal;
        end
        
        if ~isempty(strfind(file_list{list_I,2},double_add))
            Nbin = num_list(list_I,2:3);
        else
            Nbin = ones(1,2);
            Nbin(num_list(list_I,3)) = num_list(list_I,2);
        end
        
        %%%%% Intensity correction mask generation: ---------------------------
        temp = imread([folder_name,file_list{list_I,1},tile_name,'01/',tif_name,'01.tif']);
        corr1_size = size(temp);
        corr1_size(1:2) = corr1_size(1:2)-2*delta;
        imcorr1 = corr_mask(corr1_size,channel_name,resolution,change_obj,1);
        %%%%% -----------------------------------------------------------------
        
        %%%%% Step settings of correlation: -----------------------------------
        Nx = corr1_size(2);
        
        if length(varargin) >= 3 && ~isempty(varargin{3})
            I_tile = varargin{3}(1);
            if length(varargin{3}) > 1
                dxx = varargin{3}(2);
            else
                dxx = [];
            end
            overlayP = find_overlap(imcorr1,[folder_name,file_list{list_I,1},tile_name,num2str(I_tile,'%02u'),'/'],[folder_name,file_list{list_I,1},tile_name,num2str(varargin{3}+1,'%02u'),'/'],delta,match_channel,offset0,dxx)
        end
        
%         dx0 = round(overlayP*max(corr1_size));
        dx0 = round(overlayP*Nx);
        dx = (dx0-dx_range):(dx0+dx_range);
        dx2 = dx-dx(end);
        idx0 = (length(dx)+1)/2;

        Ny = corr1_size(1);
%         dy0 = round(overlayP*max(corr1_size));
        dy0 = round(overlayP*Ny);
        dy = (dy0-dx_range):(dy0+dx_range);
        dy2 = dy-dy(end);

        dz = (-dz_range):dz_range;
        idz0 = (length(dz)+1)/2;
        
        dx_max = dx0*ones(Nbin); dx_max(1) = 0;
        dy_max = dy0*ones(Nbin(1),1); dy_max(1) = 0;
        dz_max = zeros(Nbin);
        dI_max = ones(Nbin);
        dI2_max = ones(Nbin(1),1);
        z_start = zeros(Nbin);
        empty_id = false(Nbin);
        outimage0 = cell(Nbin(1),1);

        %%%%% -----------------------------------------------------------------

%         pool_name = parpool(Nbin(1));

        %%% Horizontal stitching: %%%%%====================================
        I_series = 0;
        
        if load_old && exist([folder_name,file_list{list_I,1}(1:(end-1)),old_add,para_name,mat_tail])
            load([folder_name,file_list{list_I,1}(1:(end-1)),old_add,para_name,mat_tail],'x_start','y_start','z_start')
            x_start0 = x_start;
            y_start0 = y_start;
            z_start0 = z_start;
        end
        
        for I_Y = 1:Nbin(1)
            inimage = cell(1,Nbin(2));
            
            sub_folder1 = [tile_name,num2str(I_series+1,'%02u'),'/'];
            imlist1 = dir([folder_name,file_list{list_I,1},sub_folder1,image_tail]); %%% get the image list from image folder1
            N_Z = length(imlist1);
            Istd1 = zeros(N_Z,1);   %%% standard deviation of the left stack

            for I_layer = 1:N_Z
                im_temp = imread([folder_name,file_list{list_I,1},sub_folder1,imlist1(I_layer).name]);
                %im_temp = im_temp(1+delta(1):end-delta(1),1+delta(2):end-delta(2),:);

                if I_layer == 1
                    inimage{1} = zeros([size(im_temp),N_Z],'uint16');
                end

                %inimage{1}(:,:,:,I_layer) = uint16(double(im_temp-offset0)./imcorr1);
                inimage{1}(:,:,:,I_layer) = uint16(double(im_temp));
            end
                
            %%%% Matching: ----------------------------------------------------
            for I_X = 1:Nbin(2)-1
                I_series = I_series+1;

                sub_folder2 = [tile_name,num2str(I_series+1,'%02u'),'/'];
                for I_layer = 1:N_Z
                    im_temp = imread([folder_name,file_list{list_I,1},sub_folder2,imlist1(I_layer).name]);
                    im_temp = im_temp(1+delta(1):end-delta(1),1+delta(2):end-delta(2),:);

                    if I_layer == 1
                        inimage{I_X+1} = zeros([size(im_temp),N_Z],'uint16');
                    end
                    inimage{I_X+1}(:,:,:,I_layer) = uint16(double(im_temp-offset0)./imcorr1);

                    Istd1(I_layer) = std2(inimage{I_X}(:,Nx-dx(end)+1:Nx,match_channel,I_layer));
%                     Istd2(I_layer) = std2(inimage(:,1:dx(end),match_channel,I_layer,I_X+1));
                end

                [max00,Icontrast1] = max(Istd1);
%                 [~,Icontrast2] = max(Istd2);
                Isel = Icontrast1;

                if length(dx) > 1 || length(dz) > 1
                    im1 = gpuArray(inimage{I_X}(:,Nx-dx(end)+1:Nx,match_channel,:));
                    im2 = gpuArray(inimage{I_X+1}(:,1:dx(end),match_channel,:));
                    Nx2 = size(im1,2);

                    if ~(load_old && exist([folder_name,file_list{list_I,1}(1:(end-1)),old_add,para_name,mat_tail]))
                        corr_tile = gpuArray(-inf(size(dx,2),size(dz,2),'single'));

                        for ix = 1:size(dx,2)
                            for iz = 1:size(dz,2)
                                if Isel+dz(iz) >= 1 && Isel+dz(iz) <= N_Z
                                    corr_tile(ix,iz) = corr2(im1(:,(1-dx2(ix)):Nx2,1,Isel),im2(:,1:(Nx2+dx2(ix)),1,Isel+dz(iz)));
                                end
                            end
                        end
                        if max00 >= std_th
                            [~,ind0] = max(corr_tile(:));
                            [ix0,iz0] = ind2sub(size(corr_tile),ind0);
                        else
                            ix0 = idx0;
                            iz0 = idz0;
                            empty_id(I_Y,I_X+1) = true;
                        end
                        dx_max(I_Y,I_X+1) = dx(ix0);
                        dz_max(I_Y,I_X+1) = dz(iz0);
                    else
                        dx_max(I_Y,I_X+1) = x_start0(I_Y,I_X+1)-1;
                        dz_max(I_Y,I_X+1) = z_start0(I_Y,I_X+1)-z_start0(I_Y,I_X);
                        ix0 = find(dx_max(I_Y,I_X+1) == dx);
                        iz0 = find(dz_max(I_Y,I_X+1) == dz);
    %                     dI_max(I_Y,I_X+1) = gather(sum(sum(inimage{I_X+1}(:,1:dx_max(I_Y,I_X+1),match_channel,Isel+dz_max(I_Y,I_X+1))))/sum(sum(im1(:,(1-dx_max(I_Y,I_X+1)+dx(end)):dx(end),match_channel,Isel))));
                    end
                    
                    dI_max(I_Y,I_X+1) = gather(sum(sum(im2(:,1:(Nx2+dx2(ix0)),1,Isel+dz(iz0))))/sum(sum(im1(:,(1-dx2(ix0)):Nx2,1,Isel))));
                end
            end
            
            I_series = I_series+1;
            %%%% --------------------------------------------------------------
            
            %%%% Stitching: ---------------------------------------------------
            z_max = cumsum(dz_max(I_Y,:));
            I_layer = max(1-z_max):min(N_Z-z_max);
            z_start(I_Y,:) = max(1-z_max)+z_max;
            
            I_layer0 = I_layer;
            outimage0{I_Y} = inimage{1}(:,:,:,I_layer0);
            outimage_20X{I_Y} = inimage{1}(:,:,:,I_layer0);

            for I_X = 2:Nbin(2)
                I_layer0 = I_layer0+dz_max(I_Y,I_X);
                outimage0{I_Y} = cat(2,outimage0{I_Y},inimage{I_X}(:,dx_max(I_Y,I_X)+1:end,:,I_layer0));
            end
            %%%% --------------------------------------------------------------
        end
        clear inimage
        
        
        %%% Vertical stitching: %%%%%======================================
        ds = (-ds_range):ds_range;
        ds2_max = zeros(Nbin(1),1);
        dz2_max = zeros(Nbin(1),1);
        N_xyz = zeros(Nbin(1),3);
%         outimage = uint16(zeros(0));

%         if Nbin(1) > 1
        if (length(dy) > 1 || length(dz) > 1) && Nbin(1) > 1
            for I_Y = 1:Nbin(1)-1
                j1 = max(round(size(outimage0{I_Y},2)*p_start),1);
                j2 = round(size(outimage0{I_Y},2)*p_end);
%                 im1 = gpuArray(outimage0{I_Y}(Ny-dy(end)+1:Ny,j1:j2,match_channel,:));
%                 im2 = gpuArray(outimage0{I_Y+1}(1:dy(end),j1:j2,match_channel,:));
                im1 = gpuArray(outimage0{I_Y}(Ny-dy(end)+1:Ny,:,match_channel,:));
                im2 = gpuArray(outimage0{I_Y+1}(1:dy(end),:,match_channel,:));

                N_Z1 = size(im1,4);
                N_Z2 = size(im2,4);
                Istd1 = gpuArray(zeros(1,N_Z1));   %%% standard deviation of the upper stack
                N_xyz(I_Y,:) = [size(outimage0{I_Y},1),size(outimage0{I_Y},2),size(outimage0{I_Y},4)];

                for I_layer = 1:N_Z1
                    Istd1(I_layer) = std2(im1(:,:,1,I_layer));
%                     Istd1(I_layer) = std2(outimage0{I_Y}(Ny-dy(end)+1:Ny,:,match_channel,I_layer));
                end

                [~,Icontrast1] = max(Istd1);
%                 [~,Icontrast2] = max(Istd2);
                Isel = Icontrast1;
                Ny2 = size(im1,1);

                Ns1 = size(im1,2);
                Ns2 = size(im2,2);
                ds = (-ds_range):ds_range;
                Is1_start = 1-ds.*(ds < 0);
                Is2_start = fliplr(Is1_start);
                Is1_end = min(Ns1,Ns2-ds);
                Is2_end = min(Ns1+ds,Ns2);

                if ~(load_old && exist([folder_name,file_list{list_I,1}(1:(end-1)),old_add,para_name,mat_tail]))
                    corr_tile = gpuArray(-inf(size(dy,2),size(ds,2),size(dz,2),'single'));

    %                 pool_name = parpool(size(ds,2));

                    for is = 1:size(ds,2)
                        for ix = 1:size(dy,2)
                            for iz = 1:size(dz,2)
                                if Isel+dz(iz) >= 1 && Isel+dz(iz) <= N_Z2
                                    corr_tile(ix,is,iz) = corr2(im1((1-dy2(ix)):Ny2,Is1_start(is):Is1_end(is),1,Isel),im2(1:(Ny2+dy2(ix)),Is2_start(is):Is2_end(is),1,Isel+dz(iz)));
                                end
                            end
                        end
                    end

    %                 delete(pool_name)

                    [~,ind0] = max(corr_tile(:));
                    [ix0,is0,iz0] = ind2sub(size(corr_tile),ind0);
                    dy_max(I_Y+1) = dy(ix0);
                    ds2_max(I_Y+1) = ds(is0);
                    dz2_max(I_Y+1) = dz(iz0);
                else
                    dy_max(I_Y+1) = y_start0(I_Y+1,1)-1;
                    ds2_max(I_Y+1) = x_start0(I_Y+1,1)-x_start0(I_Y,1);
                    dz2_max(I_Y+1) = min(z_start0(I_Y+1,:))-min(z_start0(I_Y,:));
                    
                    ix0 = find(dy_max(I_Y+1) == dy);
                    is0 = find(ds2_max(I_Y+1) == ds);
                    iz0 = find(dz2_max(I_Y+1) == dz);
                end
                
                dI2_max(I_Y+1) = gather(sum(sum(im2(1:(Ny2+dy2(ix0)),Is2_start(is0):Is2_end(is0),1,Isel+dz(iz0))))/sum(sum(im1((1-dy2(ix0)):Ny2,Is1_start(is0):Is1_end(is0),1,Isel))));
            end
            
            I_Y = Nbin(1);
            N_xyz(I_Y,:) = [size(outimage0{I_Y},1),size(outimage0{I_Y},2),size(outimage0{I_Y},4)];
            
        else
            for I_Y = 1:Nbin(1)
                N_xyz(I_Y,:) = [size(outimage0{I_Y},1),size(outimage0{I_Y},2),size(outimage0{I_Y},4)];
            end
        end
        
        %%%% --------------------------------------------------------------

        %%%% Stitching: ---------------------------------------------------
        z2_max = cumsum(dz2_max);
        s_max = cumsum(ds2_max);
        I_layer = max(1-z2_max):min(N_xyz(:,3)-z2_max);
        Is_start0 = max(1-s_max);
        Is_end0 = min(N_xyz(:,2)-s_max);
        Is_start = Is_start0+s_max;
        Is_end = Is_end0+s_max;

% % %         I_layer0 = I_layer;
% % %         I_Y = 1;
% % %         outimage = outimage0{I_Y}(:,Is_start(I_Y):Is_end(I_Y),:,I_layer0);
% % % 
% % %         for I_Y = 2:Nbin(1)
% % %             I_layer0 = I_layer0+dz2_max(I_Y);
% % %             outimage = cat(1,outimage,outimage0{I_Y}(dy_max(I_Y)+1:end,Is_start(I_Y):Is_end(I_Y),:,I_layer0));
% % %             outimage0{I_Y} = zeros(0);
% % %         end
% % %         %%%% --------------------------------------------------------------
% % % 
% % %         %%%% Output: ------------------------------------------------------
% % %         for Iz = 1:length(I_layer)
% % %             if exist([folder_name,out_name]) ~= 7
% % %                 mkdir([folder_name,out_name]);
% % %             end
% % %             im_name = [tif_name,num2str(Iz,'%02u'),figure_tail];
% % % %                 imwrite(outimage,[folder_name,out_name,imlist1(I_layer).name]);
% % %             tiffwrite0(outimage(:,:,:,Iz),[folder_name,out_name,im_name]);
% % %         end
% % %             %%%% --------------------------------------------------------------
            
        for Iz = 1:length(I_layer)
            I_layer0 = I_layer(Iz);
            I_Y = 1;
            %outimage = outimage0{I_Y}(:,Is_start(I_Y):Is_end(I_Y),:,I_layer0);
            outimage = outimage_20X{I_Y}(:,:,:,I_layer0);

            for I_Y = 2:Nbin(1)
                I_layer0 = I_layer0+dz2_max(I_Y);
                outimage = cat(1,outimage,outimage0{I_Y}(dy_max(I_Y)+1:end,Is_start(I_Y):Is_end(I_Y),:,I_layer0));
            end
        %%%% --------------------------------------------------------------

        %%%% Output: ------------------------------------------------------
            if exist([folder_name,out_name]) ~= 7
                mkdir([folder_name,out_name]);
            end
            im_name = [tif_name,num2str(Iz,'%02u'),figure_tail];
%                 imwrite(outimage,[folder_name,out_name,imlist1(I_layer).name]);
            
            tiffwrite0(outimage,[folder_name,out_name,im_name]);
        end
        %%%% --------------------------------------------------------------
            
% %             
% %         else
% %             sub_folder1 = [YY_name,'01/'];
% %             movefile([folder_name,out_name0,sub_folder1],[folder_name,out_name]);
% %             
% %             Is_start = 0;
% %             z2_max = 0;
% %         end
    
        %%%% Calculate displacements and output parameters: ---------------
        x_start = dx_max+1; 
        x_start(:,1) = Is_start;
        y_start = repmat(dy_max,1,Nbin(2))+1;
        z_start = z_start+repmat(max(-z2_max)+z2_max,1,Nbin(2));
        save([folder_name,out_name,para_name,mat_tail],'x_start','y_start','z_start','imcorr1','overlayP','delta','offset0','dI_max','dI2_max','empty_id')
        
        if load_old && exist([folder_name,file_list{list_I,1}(1:(end-1)),old_add,para_name,mat_tail]) && (any(x_start(:) ~= x_start0(:)) || any(y_start(:) ~= y_start0(:)) || any(z_start(:) ~= z_start0(:)))
            ['Fail to restitch: ',folder_name,out_name]
        end
        %%%% --------------------------------------------------------------


        file_list(list_I,3) = {out_name}; %%% Output folder name record

        try
            xlswrite([folder_name,list_name],cat(2,file_list,num2cell(num_list)));
        catch
            xlwrite([folder_name,list_name],cat(2,file_list,num2cell(num_list)));
        end
        
%         rmdir([folder_name,out_name0], 's');
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

clear global



function overlayP = find_overlap(imcorr1,folder1,folder2,delta,match_channel,offset0,dxx0)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A function to find the maximum overlap of two tiles %%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

image_tail = '*.tif';
overlayP0 = 0.20;
if isempty(dxx0)
    dxx = 5;
else
    dxx = dxx0;
end
dz_range = 10;
acc0 = 1000; %%% accuracy

Nx = size(imcorr1,2);
dx = dxx:dxx:round(overlayP0*Nx);
dx2 = dx-dx(end);

dz = (-dz_range):dz_range;
        

%% Data loading:%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imlist1 = dir([folder1,image_tail]); %%% get the image list from image folder1
N_Z = length(imlist1);
Istd1 = zeros(N_Z,1);   %%% standard deviation of the left stack
inimage1 = zeros([size(imcorr1,1),size(imcorr1,2),N_Z]);
inimage2 = zeros([size(imcorr1,1),size(imcorr1,2),N_Z]);

for I_layer = 1:N_Z
    im_temp = imread([folder1,imlist1(I_layer).name]);
    inimage1(:,:,I_layer) = im_temp(1+delta(1):end-delta(1),1+delta(2):end-delta(2),match_channel);
    im_temp = imread([folder2,imlist1(I_layer).name]);
    inimage2(:,:,I_layer) = im_temp(1+delta(1):end-delta(1),1+delta(2):end-delta(2),match_channel);
    Istd1(I_layer) = std2(inimage2(:,Nx-dx(end)+1:Nx,I_layer));
end
inimage1 = uint16(double(inimage1-offset0)./repmat(imcorr1(:,:,match_channel),1,1,N_Z));
inimage2 = uint16(double(inimage2-offset0)./repmat(imcorr1(:,:,match_channel),1,1,N_Z));
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Matching: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,Icontrast1] = max(Istd1);
Isel = Icontrast1;

im1 = gpuArray(inimage1(:,Nx-dx(end)+1:Nx,:));
im2 = gpuArray(inimage2(:,1:dx(end),:));
Nx2 = size(im1,2);

corr_tile = gpuArray(-inf(size(dx,2),size(dz,2),'single'));

for ix = 1:size(dx,2)
    for iz = 1:size(dz,2)
        if Isel+dz(iz) >= 1 && Isel+dz(iz) <= N_Z
            corr_tile(ix,iz) = corr2(im1(:,(1-dx2(ix)):Nx2,Isel),im2(:,1:(Nx2+dx2(ix)),Isel+dz(iz)));
        end
    end
end
[~,ind0] = max(corr_tile(:));
[ix0,~] = ind2sub(size(corr_tile),ind0);
dx_max = dx(ix0);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Overlap estimation: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
overlayP = round(dx_max/Nx*acc0)/acc0;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%











