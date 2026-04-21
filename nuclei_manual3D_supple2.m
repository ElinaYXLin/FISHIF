function nuclei_manual3D_supple2(lf_name,varargin)
%clear all
close all
tic

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameter setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
in_folder = 'stacks/';
out_folder = 'masks/';
input_name = 'matchlist.xls';
quick_name = 'quick_mask.mat';
mask_name = 'mask.mat';
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
        if exist([folder_list{list_I,1},out_folder,sub_list{list_J,3},quick_name])
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
                bw_conv = bwconvhull(bw_applied3D_save(:,:,ii),'objects');
                D= bwdist(~bw_conv);
                g2 = imimposemin(-D,bw_applied3D_save(:,:,ii));
                bw_cut0 = imdilate(watershed(g2) == 0,strel('disk',bin0));% & (temp_sel);
                bw_conv = bw_conv& (~ bw_cut0);
                
                bw_applied3D_save(:,:,ii) = imopen(bw_conv,strel('disk',bin0));
                bw_applied3D_save(:,:,ii) = bwareaopen(bw_applied3D_save(:,:,ii),area0);
            end
        end

        if bin0 > 1
            bw_applied3D_save0 = imresize(bw_applied3D_save,1/bin0,'nearest');
            mask_stack0 = label3D(bw_applied3D_save0);   %%% mask tiling

            mask_stack1 = imresize(mask_stack0,bin0,'nearest');
            mask_stack1 = mask_stack1(1:size0(1),1:size0(2),:);
            mask_stack = zeros(size0,'uint16');

            clear bw_applied3D_save0 mask_stack0

            for ii = 1:size0(3)
                bw_conv = bw_applied3D_save(:,:,ii);

                bw_label0 = bwlabel(bw_conv);
                bw_prop = regionprops(bw_label0.*imerode(bw_conv,strel('disk',bin0)),mask_stack1(:,:,ii),'MaxIntensity');   %%% convert bwlabel1 labels to temp_label1 labels
                bw_temp = [bw_prop.MaxIntensity];
                bw_temp = uint16([0,bw_temp]);
                    ii
                mask_stack(:,:,ii) = bw_temp(bw_label0+1);
            end
        else
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
    end
end

toc



