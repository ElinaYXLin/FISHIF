function mask_mod_temp(lf_name)
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
low_th = 500;
high_th = 50000;
cir_th = 0.8;
th_ratio = 1.0;
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

    for list_J = [1:3,5]%1:M1
        load([folder_list{list_I,1},out_folder,sub_list{list_J,3},mask_name])
        mask_temp = false(size(mask_stack));
%matlabpool
        for image_I = 1:size(mask_stack,3)
                mask_temp(:,:,image_I) = imopen(logical(mask_stack(:,:,image_I)),strel('disk',7));
        end
%matlabpool close

        mask_stack = mask_stack.*mask_temp;
        
        save([folder_list{list_I,1},out_folder,sub_list{list_J,3},mask_name],'mask_stack','-v7.3')
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end
end

toc
