function nuclei_manual3D_supple(lf_name,varargin)
%clear all
close all
tic

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A function to label manual segmentation results or refine/doub-check original segmentation results %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% lf_name (input): Data folder for processing, can be in two formats: (a) input as a character array (eg. 'Duallist.xls'): the name of an Excel file containing the list of data folders;
%%                                                                     (b) input as a cell array (eg. {'111/';'222'}): manually list data folders (the delimiter must be ";", not ",");
%% varargin (input): 1. indices of embryos for processing (default: [] for all embryos. Example: {inds for folder 1,inds for folder 2});
%%                   2. Scaling of the segmentation results during processing (default: 1. Example: reduce the size to 1/3: 1/3);
%%                   3. Smallest size of a nucleus (default: 0); 
%%                   4. Smallest z-span of a nucleus (default: 3);
%%                   5. Switch to convex the identified nuclei (default: true);
%%                   6. Switch to load quick_mask.mat (default: true);
%%                   7. Switch to recheck the cutting of overlapped nuclei using watershed (default: false);
%%                   8. Switch to remove the nuclei outside the embryo mask (default: false);
%%                   9. Coordinates of nuclei for removal (To remove a nucleus, include the scaled 2D coordinate of one of its pixel. Only works for labeled mask. E.g., [1,1;2,2].)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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

if length(varargin) >= 9 && ~isempty(varargin{9})
    xy_remove = varargin{9}*bin0;
else
    xy_remove = [];
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
            if ~isempty(xy_remove)
                mask2D = max(mask_stack,[],3);
                in_remove = sub2ind([size(mask_stack,1),size(mask_stack,2)],round(xy_remove(:,2)),round(xy_remove(:,1)));
                ind_remove = mask2D(in_remove);
                ind_nu0 = setdiff(unique(mask_stack),0);
                ind_nu = setdiff(ind_nu0,ind_remove);
                bw_applied3D_save = ismember(mask_stack,ind_nu);
            else
                bw_applied3D_save = logical(mask_stack);
            end
        end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3D mask processing: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        size0 = size(bw_applied3D_save);

        if bin0 > 1
            bw_applied3D_save0 = imresize(bw_applied3D_save,1/bin0,'nearest');
        else
            bw_applied3D_save0 = bw_applied3D_save;
        end
            
        if t_cut
            bw_applied3D_save0 = recut(bw_applied3D_save0);
        end
        
            
        if t_em
            load([folder_list{list_I,1},result_folder,sub_list{list_J,3}(1:end-1),mat_tail],'em_mask')
            if exist('em_mask','var')
                if bin0 > 1
                    em_mask0 = imresize(em_mask,1/bin0,'nearest');
                else
                    em_mask0 = em_mask;
                end
                bw_applied3D_save0 = lim_em(bw_applied3D_save0,em_mask0);
                clear em_mask em_mask0
            end
        end
        
        
        if t_fine
            bw_applied3D_save0 = nu_refine(bw_applied3D_save0);
        end
        
        
        if area0 > 0
            for ii = 1:size0(3)
% %                 bw_applied3D_save(:,:,ii) = imopen(bw_applied3D_save(:,:,ii),strel('disk',bin0));
                bw_applied3D_save0(:,:,ii) = bwareaopen(bw_applied3D_save0(:,:,ii),round(area0));
            end
        end
        
        
        mask_stack0 = label3D(bw_applied3D_save0,z_thresh);   %%% mask tiling
        clear bw_applied3D_save bw_applied3D_save0
        
        
        if bin0 > 1
            mask_stack_temp = imresize(mask_stack0,bin0,'nearest');
            mask_stack = mask_stack_temp(1:size0(1),1:size0(2),:);
        else
            mask_stack = mask_stack0;
        end
        clear mask_stack0 mask_stack_temp
            
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
    
    bw_out(:,:,ii) = bw0;
end
delete(pid)



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A function to check every layer of the mask and cut joint nuclei
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function bw_out = nu_refine(bw_in)

bw_out = false(size(bw_in));
Npid = feature('numcores');
if Npid > 18
    Npid = 18;
end
pid = parpool(Npid);

parfor ii = 1:size(bw_in,3)
    bw0 = bw_in(:,:,ii);
    bw_conv = bwconvhull(bw0,'objects');

    D= bwdist(~bw_conv);
    g2 = imimposemin(-D,bw0);
    bw_cut0 = imdilate(watershed(g2) == 0,strel('disk',1));% & (temp_sel);
    bw_conv = bw_conv& (~ bw_cut0);
    
    bw_out(:,:,ii) = bw_conv;
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
