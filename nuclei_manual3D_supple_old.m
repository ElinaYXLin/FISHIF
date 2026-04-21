function nuclei_manual3D_supple(lf_name,varargin)
%clear all
close all
tic

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameter setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
in_folder = 'stacks/';
out_folder = 'masks/';
result_folder = 'Results/';
input_name = 'matchlist.xls';
quick_name = 'quick_mask.mat';
mask_name = 'mask.mat';
mat_tail = '.mat';
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

if length(varargin) >= 2 && ~isempty(varargin{2})
    scale0 = varargin{2};
else
    scale0 = 1;
end
bin0 = round(1/scale0);

if length(varargin) >= 3 && ~isempty(varargin{3})
    area0 = varargin{3};
else
    area0 = 0;
end

if length(varargin) >= 4 && ~isempty(varargin{4})
    z_thresh = varargin{4};
else
    z_thresh = 3;
end

if length(varargin) >= 5 && ~isempty(varargin{5})
    t_fine = varargin{5};
else
    t_fine = true;
end

if length(varargin) >= 6 && ~isempty(varargin{6})
    t_quick = varargin{6};
else
    t_quick = true;
end

if length(varargin) >= 7 && ~isempty(varargin{7})
    t_cut = varargin{7};
else
    t_cut = false;
end

if length(varargin) >= 8 && ~isempty(varargin{8})
    t_em = varargin{8};
else
    t_em = false;
end

for list_I = 1:N1
    [sub_num, sub_list] = xlsread([folder_list{list_I,1},in_folder,input_name]);
    [M1,M2] = size(sub_list);
%     channel_name = eval(folder_list{list_I,5});
    if ~isempty(varargin) && ~isempty(varargin{1})
        list_J_all = varargin{1}{list_I};
    else
        list_J_all = 1:M1;
    end
    
    for list_J = list_J_all%1:M1
        if exist([folder_list{list_I,1},out_folder,sub_list{list_J,3},quick_name]) && t_quick
            load([folder_list{list_I,1},out_folder,sub_list{list_J,3},quick_name]);
        else
            load([folder_list{list_I,1},out_folder,sub_list{list_J,3},mask_name]);
            bw_applied3D_save = logical(mask_stack);
        end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3D mask tiling: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        size0 = size(bw_applied3D_save);

        if area0 > 0
            for ii = 1:size0(3)
% %                 bw_applied3D_save(:,:,ii) = imopen(bw_applied3D_save(:,:,ii),strel('disk',bin0));
                bw_applied3D_save(:,:,ii) = bwareaopen(bw_applied3D_save(:,:,ii),area0);
            end
        end

        if bin0 > 1
            bw_applied3D_save0 = imresize(bw_applied3D_save,1/bin0,'nearest');
            
            if t_cut
                bw_applied3D_save0 = recut(bw_applied3D_save0);
                bw_applied3D_save = imresize(bw_applied3D_save0,bin0,'nearest');
                bw_applied3D_save = bw_applied3D_save(1:size0(1),1:size0(2),:);
                if area0 > 0
                    for ii = 1:size0(3)
                        bw_applied3D_save(:,:,ii) = bwareaopen(bw_applied3D_save(:,:,ii),area0);
                    end
                end
            end
            
            if t_em
                load([folder_list{list_I,1},result_folder,sub_list{list_J,3}(1:end-1),mat_tail],'em_mask')
                if exist('em_mask','var')
                    em_mask0 = imresize(em_mask,1/bin0,'nearest');
                    bw_applied3D_save0 = lim_em(bw_applied3D_save0,em_mask0);
                    clear em_mask em_mask0
                end
            end
            
            mask_stack0 = label3D(bw_applied3D_save0);   %%% mask tiling

            mask_stack1 = imresize(mask_stack0,bin0,'nearest');
            mask_stack1 = mask_stack1(1:size0(1),1:size0(2),:);
            
            if t_fine
                mask_stack = zeros(size0,'uint16');
                clear bw_applied3D_save0 mask_stack0

                Npid = feature('numcores');
                if Npid > 18
                    Npid = 18;
                end
                pid = parpool(Npid);

                parfor ii = 1:size0(3)
                    bw0 = bw_applied3D_save(:,:,ii);
                    bw_conv = bwconvhull(bw0,'objects');

                    D= bwdist(~bw_conv);
                    g2 = imimposemin(-D,bw0);
                    bw_cut0 = imdilate(watershed(g2) == 0,strel('disk',bin0));% & (temp_sel);
                    bw_conv = bw_conv& (~ bw_cut0);

                    bw_label0 = bwlabel(bw_conv);
                    bw_prop = regionprops(bw_label0.*imerode(bw_conv,strel('disk',bin0)),mask_stack1(:,:,ii),'MaxIntensity');   %%% convert bwlabel1 labels to temp_label1 labels
                    bw_temp = [bw_prop.MaxIntensity];
                    bw_temp = uint16([0,bw_temp]);
    %                     ii
                    mask_stack(:,:,ii) = bw_temp(bw_label0+1);
                end
                delete(pid)
                
            else
                mask_stack = mask_stack1;
                clear bw_applied3D_save0 mask_stack0 mask_stack1
                
            end
        else
            if t_cut
                bw_applied3D_save = recut(bw_applied3D_save);
            end
            if t_em
                load([folder_list{list_I,1},result_folder,sub_list{list_J,3}(1:end-1),mat_tail],'em_mask')
                if exist('em_mask','var')
                    bw_applied3D_save = lim_em(bw_applied3D_save,em_mask);
                    clear em_mask
                end
            end
            mask_stack = label3D(bw_applied3D_save,z_thresh);   %%% mask tiling
        end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Segmentation result output: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         if ~exist([folder_list{list_I,1},out_folder,sub_list{list_J,3}],'dir')
%             mkdir([folder_list{list_I,1},out_folder,sub_list{list_J,3}])
%         end
        save([folder_list{list_I,1},out_folder,sub_list{list_J,3},mask_name],'mask_stack','-v7.3')
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        disp(['I = ',num2str(list_I),'J = ',num2str(list_J),', t = ',num2str(toc)])
    
    end
end

toc



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A function to check every layer of the mask and cut joint nuclei
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function bw_out = recut(bw_in)

bw_out = false(size(bw_in));
Npid = feature('numcores');
if Npid > 18
    Npid = 18;
end
pid = parpool(Npid);

parfor ii = 1:size(bw_in,3)
    bw0 = bw_in(:,:,ii);
    bw0  = ellipse_fit(bw0,20);   %%% refine the mask using elliptical shape
    bw0  = ellipse_fit(bw0,20);   %%% refine the mask using elliptical shape
                
    D= -bwdist(~bw0);
    bw_cut = imextendedmin(D,2);
    g2 = imimposemin(D,bw_cut);
    bw_cut = imdilate(watershed(g2) == 0,strel('disk',1));% & (temp_sel);
    bw0 = bw0 & (~ bw_cut);
    
    bw_out(:,:,ii) = bw0
end
delete(pid)



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A function to remove nuclei outside the embryo mask
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function bw_out = lim_em(bw_in,em_mask)

bw_out = false(size(bw_in));
Npid = feature('numcores');
if Npid > 18
    Npid = 18;
end
pid = parpool(Npid);

parfor ii = 1:size(bw_in,3)
    bw0 = bwlabel(bw_in(:,:,ii));
    prop0 = regionprops(bw0,em_mask,'MaxIntensity');
    nu_true = find([prop0.MaxIntensity]);
    bw_out(:,:,ii) = ismember(bw0,nu_true);
end
delete(pid)



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
