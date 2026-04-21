function nuclei_overlap_detect(lf_name,varargin)
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
mask_ovname = 'mask_overlap.fig';
mask_imname = 'mask_image.fig';
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
        load([folder_list{list_I,1},out_folder,sub_list{list_J,3},mask_name]);
        load([folder_list{list_I,1},result_folder,sub_list{list_J,3}(1:end-1),mat_tail],'max_image','DAPI_channel');
        if bin0 > 1
            mask_stack = imresize(mask_stack,1/bin0,'nearest');
            max_image = imresize(max_image,1/bin0,'nearest');
        end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% mask overlay with image: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        im111 = zeros(size(mask_stack,1),size(mask_stack,2),3,'uint16');
        figure
            im111(:,:,3) = max_image(:,:,DAPI_channel)*3; 
            im111(:,:,2) = uint16(bwperim(any(logical(mask_stack),3)))*65535; 
            imshow(im111)
        title([folder_list{list_I,1},out_folder,sub_list{list_J,3}])
        clear im111
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Detection result output: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %         save([folder_list{list_I,1},out_folder,sub_list{list_J,3},mask_name],'mask_stack','-v7.3')
        saveas(gcf,[folder_list{list_I,1},out_folder,sub_list{list_J,3},mask_imname])
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% mask overlap plot: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        max_mask = double(max(mask_stack,[],3));
        N_mask = sum(logical(mask_stack),3);
        diff_mask = sum(mask_stack,3)-max_mask.*N_mask ~= 0;
        
        prop0 = regionprops(max_mask,'Centroid');
        nu_xy0 = cell2mat({prop0.Centroid}');
        I_true0 = ~isnan(nu_xy0(:,1)) & ~isnan(nu_xy0(:,2));
        nu_xy = round(nu_xy0(I_true0,:));
        nu_ind = sub2ind(size(max_mask),nu_xy(:,2),nu_xy(:,1));
        nu_center = false(size(max_mask));
        nu_center(nu_ind) = true;
        nu_center = imdilate(nu_center,strel('disk',5));
        
        im000 = zeros([size(mask_stack,1),size(mask_stack,2),3]);
        logi_mask = logical(max_mask);
        im000(:,:,1) = logi_mask-(nu_center & ~diff_mask);
        im000(:,:,2) = logi_mask.*(1-(diff_mask | nu_center));
        im000(:,:,3) = (logi_mask.*(1-diff_mask)) | nu_center;
        
% %         im000 = zeros([size(mask_stack,1),size(mask_stack,2),3]);
% %         im000(:,:,1) = logical(max_mask);
% %         im000(:,:,2) = im000(:,:,1).*(1-diff_mask);
% %         im000(:,:,3) = im000(:,:,2);
        
        figure
            imshow(double(im000))
        title([folder_list{list_I,1},out_folder,sub_list{list_J,3}])
        
        clear im000
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Detection result output: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %         save([folder_list{list_I,1},out_folder,sub_list{list_J,3},mask_name],'mask_stack','-v7.3')
        saveas(gcf,[folder_list{list_I,1},out_folder,sub_list{list_J,3},mask_ovname])
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %         
% %         nu_hollow = sum(abs(diff(logical(mask_stack),1,3)),3)-2 > 0;
% % 
% %         im111 = zeros([size(mask_stack,1),size(mask_stack,2),3]);
% %         im111(:,:,1) = logi_mask-nu_hollow;
% %         im111(:,:,2) = logi_mask;
% %         im111(:,:,3) = logi_mask-nu_hollow;
% %         
% %         figure
% %             imshow(double(im111))
% %         title([folder_list{list_I,1},out_folder,sub_list{list_J,3}])
% %         
% %         clear im111

    end
end

toc



