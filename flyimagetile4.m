function flyimagetile4(lf_name)
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
YY_name = 'Y';
tif_name = 'stack';
single_add = '_new/';
temp_add = '_temp/';
figure_tail = '.tif';
list_name = 'matchlist.xls';

standard_record_confocal = 'Calibration/Results/standard.mat';
% standard_record_60X = 'Calibration3/Results/standard_60X.mat';
standard_record_60X_confocal = 'Calibration3_08112019/Results/standard_60X.mat';
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
dx_range_airy = 50;
ds_range_airy = 50;
delta_airy = [0,0];

overlay_confocal = 0;
dz_range_confocal = 0;
dx_range_confocal = 0;
ds_range_confocal = 0;
delta_confocal = [0,0];

para_name = 'stitch_parameter';
output_tail = '.xls';
mat_tail = '.mat';
image_tail = '*.tif';

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Processing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for list_I = 1:N1
        out_name = [file_list{list_I,1}(1:(end-1)),single_add]; %%% output folder name
        out_name0 = [file_list{list_I,1}(1:(end-1)),temp_add]; %%% temp output folder name
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
            standard_record = standard_record_fast;
            standard_record_60X = standard_record_60X_fast;
        elseif ~isempty(strfind(file_list{list_I,2},airy_add))
            offset0 = 0;
            delta = delta_airy;
            standard_record = standard_record_airy;
            standard_record_60X = standard_record_60X_airy;
        else
            offset0 = 0;
            delta = delta_confocal;
            standard_record = standard_record_confocal;
            standard_record_60X = standard_record_60X_confocal;
        end
        standard_data = load(standard_record);
        standard_data_60X = load(standard_record_60X);

        
        if ~isempty(strfind(file_list{list_I,2},fast_add)) || ~isempty(strfind(file_list{list_I,2},airy_add))
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
        dx0 = round(overlayP*max(corr1_size));
%         dx0 = round(overlayP*Nx);
        dx = (dx0-dx_range):(dx0+dx_range);
        dx2 = dx-dx(end);

        Ny = corr1_size(1);
        dy0 = round(overlayP*max(corr1_size));
%         dy0 = round(overlayP*Ny);
        dy = (dy0-dx_range):(dy0+dx_range);
        dy2 = dy-dy(end);

        dz = (-dz_range):dz_range;
        
        dx_max = dx0*ones(Nbin); dx_max(1) = 0;
        dy_max = dy0*ones(Nbin(1),1); dy_max(1) = 0;
        dz_max = zeros(Nbin);
        z_start = zeros(Nbin);

        %%%%% -----------------------------------------------------------------

%         pool_name = parpool(Nbin(1));

        %%% Horizontal stitching: %%%%%====================================
        I_series = 0;
        for I_Y = 1:Nbin(1)
            %%%% Matching: ----------------------------------------------------
            for I_X = 1:Nbin(2)-1
                I_series = I_series+1;
                if length(dx) > 1 || length(dz) > 1
                    
                    sub_folder1 = [tile_name,num2str(I_series,'%02u'),'/'];
                    sub_folder2 = [tile_name,num2str(I_series+1,'%02u'),'/'];
                    imlist1 = dir([folder_name,file_list{list_I,1},sub_folder1,image_tail]); %%% get the image list from image folder1
                    N_Z = length(imlist1);

                    im1 = gpuArray(zeros(0));   %%% left stack
                    im2 = gpuArray(zeros(0));   %%% right stack
                    Istd1 = zeros(0);   %%% standard deviation of the left stack
    %                 Istd2 = zeros(0);   %%% standard deviation of the right stack

                    for I_layer = 1:N_Z
                        temp = double(imread([folder_name,file_list{list_I,1},sub_folder1,imlist1(I_layer).name])-offset0);
                        temp = temp(1+delta(1):end-delta(1),1+delta(2):end-delta(2),:);
                        temp0 = temp(:,Nx-dx(end)+1:Nx,match_channel)./imcorr1(:,Nx-dx(end)+1:Nx,match_channel);
                            im1 = cat(3,im1,temp0);
                            Istd1 = cat(1,Istd1,std2(temp0));
                        temp = double(imread([folder_name,file_list{list_I,1},sub_folder2,imlist1(I_layer).name])-offset0);
                        temp = temp(1+delta(1):end-delta(1),1+delta(2):end-delta(2),:);
                        temp0 = temp(:,1:dx(end),match_channel)./imcorr1(:,1:dx(end),match_channel);
                            im2 = cat(3,im2,temp0);
    %                         Istd2 = cat(1,Istd2,std2(temp0));
                    end

                    [~,Icontrast1] = max(Istd1);
    %                 [~,Icontrast2] = max(Istd2);

% %                     if Icontrast1 >= 1+dz_range && Icontrast1 <= N_Z-dz_range
                        Isel = Icontrast1;
% %                     elseif Icontrast1 >= 1+dz_range
% %                         Isel = N_Z-dz_range;
% %                     else
% %                         Isel = 1+dz_range;
% %                     end

                    corr_tile = gpuArray(-inf(size(dx,2),size(dz,2)));
                    Nx2 = size(im1,2);

                    for ix = 1:size(dx,2)
                        for iz = 1:size(dz,2)
                            if Isel+dz(iz) >= 1 && Isel+dz(iz) <= N_Z
                                corr_tile(ix,iz) = corr2(im1(:,(1-dx2(ix)):Nx2,Isel),im2(:,1:(Nx2+dx2(ix)),Isel+dz(iz)));
                            end
                        end
                    end
                    [~,ind0] = max(corr_tile(:));
                    [ix0,iz0] = ind2sub(size(corr_tile),ind0);
                    dx_max(I_Y,I_X+1) = dx(ix0);
                    dz_max(I_Y,I_X+1) = dz(iz0);
                end
            end
            
            I_series = I_series+1;
            %%%% --------------------------------------------------------------
            
            %%%% Stitching: ---------------------------------------------------
            imlist1 = dir([folder_name,file_list{list_I,1},tile_name,'01/',image_tail]); %%% get the image list from image folder1
            N_Z = length(imlist1);
            z_max = cumsum(dz_max(I_Y,:));
            I_layer = max(1-z_max):min(N_Z-z_max);
            z_start(I_Y,:) = max(1-z_max)+z_max;
            
            for Iz = 1:length(I_layer)
                I_layer0 = I_layer(Iz);
                I_series2 = I_series-Nbin(2)+1;
                sub_folder1 = [tile_name,num2str(I_series2,'%02u'),'/'];

                im_temp = imread([folder_name,file_list{list_I,1},sub_folder1,imlist1(I_layer0).name]);
                im_temp = im_temp(1+delta(1):end-delta(1),1+delta(2):end-delta(2),:);
                outimage = uint16(double(im_temp-offset0)./imcorr1);

                for I_X = 2:Nbin(2)
                    I_layer0 = I_layer0+dz_max(I_Y,I_X);
                    I_series2 = I_series2+1;
                    sub_folder1 = [tile_name,num2str(I_series2,'%02u'),'/'];
                    
                    im_temp = imread([folder_name,file_list{list_I,1},sub_folder1,imlist1(I_layer0).name]);
                    im_temp = im_temp(1+delta(1):end-delta(1),1+delta(2):end-delta(2),:);
                    im_temp = uint16(double(im_temp-offset0)./imcorr1);
                    outimage = cat(2,outimage,im_temp(:,dx_max(I_Y,I_X)+1:end,:));
                end
            %%%% --------------------------------------------------------------
            
            %%%% Output: ------------------------------------------------------
                sub_folder1 = [YY_name,num2str(I_Y,'%02u'),'/'];
                if exist([folder_name,out_name0,sub_folder1]) ~= 7
                    mkdir([folder_name,out_name0,sub_folder1]);
                end
    %                 imwrite(outimage,[folder_name,out_name,imlist1(I_layer).name]);
                tiffwrite0(outimage,[folder_name,out_name0,sub_folder1,imlist1(Iz).name]);
            end
        end
        
        
        %%% Vertical stitching: %%%%%======================================
        ds = (-ds_range):ds_range;
        ds2_max = zeros(Nbin(1),1);
        dz2_max = zeros(Nbin(1),1);
        N_xyz = zeros(Nbin(1),3);

%         if Nbin(1) > 1
            if (length(dy) > 1 || length(dz) > 1) && Nbin(1) > 1
                for I_Y = 1:Nbin(1)-1
                    sub_folder1 = [YY_name,num2str(I_Y,'%02u'),'/'];
                    sub_folder2 = [YY_name,num2str(I_Y+1,'%02u'),'/'];
                    imlist1 = dir([folder_name,out_name0,sub_folder1,image_tail]); %%% get the image list from image folder1
                    N_Z1 = length(imlist1);
                    imlist2 = dir([folder_name,out_name0,sub_folder2,image_tail]); %%% get the image list from image folder1
                    N_Z2 = length(imlist2);

                    im1 = gpuArray(zeros(0));   %%% upper stack
                    im2 = gpuArray(zeros(0));   %%% lower stack
                    Istd1 = zeros(0);   %%% standard deviation of the upper stack
    %                 Istd2 = zeros(0);   %%% standard deviation of the lower stack

                    for I_layer = 1:N_Z1
                        temp = imread([folder_name,out_name0,sub_folder1,imlist1(I_layer).name]);
                        temp0 = temp(Ny-dy(end)+1:Ny,:,match_channel);
                            im1 = cat(3,im1,temp0);
                            Istd1 = cat(1,Istd1,std2(temp0));
                        if I_layer == 1
                            xyz0 = size(temp);
                            N_xyz(I_Y,:) = [xyz0(1:2),N_Z1];
                        end
                    end
                    for I_layer = 1:N_Z2
                        temp = imread([folder_name,out_name0,sub_folder2,imlist2(I_layer).name]);
                        temp0 = temp(1:dy(end),:,match_channel);
                            im2 = cat(3,im2,temp0);
    %                         Istd2 = cat(1,Istd2,std2(temp0));
                        if I_layer == 1
                            xyz0 = size(temp);
                            N_xyz(I_Y+1,:) = [xyz0(1:2),N_Z2];
                        end
                    end

                    [~,Icontrast1] = max(Istd1);
    %                 [~,Icontrast2] = max(Istd2);

% %                     if Icontrast1 >= 1+dz_range && Icontrast1 <= min(N_Z1,N_Z2)-dz_range
                        Isel = Icontrast1;
% %                     elseif Icontrast1 >= 1+dz_range && min(N_Z1,N_Z2)-dz_range > 0
% %                         Isel = min(N_Z1,N_Z2)-dz_range;
% %                     elseif Icontrast1 <= min(N_Z1,N_Z2)-dz_range && 1+dz_range <= min(N_Z1,N_Z2)
% %                         Isel = 1+dz_range;
% %                     else
% %                         Isel = 1;
% %                     end

                    corr_tile = gpuArray(-inf(size(dy,2),size(ds,2),size(dz,2)));
                    Ny2 = size(im1,1);

                    Ns1 = size(im1,2);
                    Ns2 = size(im2,2);
                    ds = (-ds_range):ds_range;
                    Is1_start = 1-ds.*(ds < 0);
                    Is2_start = fliplr(Is1_start);
                    Is1_end = min(Ns1,Ns2-ds);
                    Is2_end = min(Ns1+ds,Ns2);

    %                 pool_name = parpool(size(ds,2));

                    for is = 1:size(ds,2)
                        for ix = 1:size(dy,2)
                            for iz = 1:size(dz,2)
                                if Isel+dz(iz) >= 1 && Isel+dz(iz) <= N_Z2
                                    corr_tile(ix,is,iz) = corr2(im1((1-dy2(ix)):Ny2,Is1_start(is):Is1_end(is),Isel),im2(1:(Ny2+dy2(ix)),Is2_start(is):Is2_end(is),Isel+dz(iz)));
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
                end
            else
                for I_Y = 1:Nbin(1)
                    sub_folder1 = [YY_name,num2str(I_Y,'%02u'),'/'];
                    imlist1 = dir([folder_name,out_name0,sub_folder1,image_tail]); %%% get the image list from image folder1
                    N_Z1 = length(imlist1);

                    temp = imread([folder_name,out_name0,sub_folder1,imlist1(1).name]);
                    xyz0 = size(temp);
                    N_xyz(I_Y,:) = [xyz0(1:2),N_Z1];
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

            for Iz = 1:length(I_layer)
                I_layer0 = I_layer(Iz);
                I_Y = 1;

                sub_folder1 = [YY_name,num2str(I_Y,'%02u'),'/'];
                im_name = [tif_name,num2str(I_layer0,'%02u'),figure_tail];

                im_temp = imread([folder_name,out_name0,sub_folder1,im_name]);
                outimage = im_temp(:,Is_start(I_Y):Is_end(I_Y),:);

                for I_Y = 2:Nbin(1)
                    I_layer0 = I_layer0+dz2_max(I_Y);

                    sub_folder1 = [YY_name,num2str(I_Y,'%02u'),'/'];
                    im_name = [tif_name,num2str(I_layer0,'%02u'),figure_tail];

                    im_temp = imread([folder_name,out_name0,sub_folder1,im_name]);
                    outimage = cat(1,outimage,im_temp(dy_max(I_Y)+1:end,Is_start(I_Y):Is_end(I_Y),:));
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
        save([folder_name,out_name,para_name,mat_tail],'x_start','y_start','z_start','imcorr1','overlayP','delta','offset0')
        %%%% --------------------------------------------------------------


        file_list(list_I,3) = {out_name}; %%% Output folder name record

        try
            xlswrite([folder_name,list_name],cat(2,file_list,num2cell(num_list)));
        catch
            xlwrite([folder_name,list_name],cat(2,file_list,num2cell(num_list)));
        end
        
        rmdir([folder_name,out_name0], 's');
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

clear global