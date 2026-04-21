function image_decross(lf_name,decross_name,varargin)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A program to de-cross multi-channel fly embryo images %%%%%%%%%%%%%%%%%%
%% lf_name (input): Excel file of fly embryo information                 %%
%% dcross_name (input): Mat file of crosstalk information                %%
%% varargin (input): Rename the previous image folder                    %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Parameter setting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
in_folder = 'stacks/';
input_name = 'matchlist.xls';
image_tail = '.tif';
decross_folder = 'Calibration_crosstalk/';

if ~isempty(varargin)
    old_add = varargin{1};
else
    old_add = '_old';
end

if length(varargin) >= 2 && ~isempty(varargin{2})
    inten_ratio = varargin{2}(:);
else
    inten_ratio = [];
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Data loading: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[num_list, folder_list] = xlsread(lf_name);
folder_list = folder_list(strcmpi('T',folder_list(:,6)),:);
[N1,N2] = size(folder_list);

load([decross_folder,decross_name]);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Image processing: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for list_I = 1:N1
    [sub_num, sub_list] = xlsread([folder_list{list_I,1},in_folder,input_name]);
    [M1,M2] = size(sub_list);
    %%% Load channel information:
    channel_name = eval(folder_list{list_I,5});
    ind_color = zeros(size(channel_name));
    for ic = 1:length(ind_color)
        ind_color(ic) = find(strcmp(channel_name{ic},all_color0));
    end
    A0 = CrossA(ind_color,ind_color);
    B0 = CrossB(ind_color');
    
    for list_J = 1:M1
        image_folder = [folder_list{list_I,1},in_folder,sub_list{list_J,3}];
        
        %%% Rename the previous image folder:
        if ~isempty(old_add)
            image_folder_old = [image_folder(1:end-1),old_add,'/'];
            if ~exist(image_folder_old)
                copyfile(image_folder,image_folder_old);
            end
        else
            image_folder_old = image_folder;
        end
        
        %%% Image loading/de-crosstalking/resaving:
        fname = dir([image_folder_old,'*',image_tail]); %%%%% get the image list from image folder
        for I_layer = 1:length(fname)
            %%%%% Image loading:
            im0 = double(imread([image_folder_old,fname(I_layer).name]));
            
            %%%%% De-crosstalk transformation:
            Lxyz = size(im0);
            im0 = reshape(permute(im0,[3,1,2]),Lxyz(3),Lxyz(1)*Lxyz(2));
            if ~isempty(inten_ratio)
                im0 = im0./repmat(inten_ratio,1,Lxyz(1)*Lxyz(2));
            end
            
            if any(B0)
                im1 = A0\(im0-repmat(B0,1,Lxyz(1)*Lxyz(2)));
            else
                im1 = A0\im0;
            end
            
            if ~isempty(inten_ratio)
                im1 = im1.*repmat(inten_ratio,1,Lxyz(1)*Lxyz(2));
            end
            im1 = permute(reshape(im1,Lxyz(3),Lxyz(1),Lxyz(2)),[2,3,1]);
            
            %%%%% Image resaving:
            tiffwrite0(uint16(im1),[image_folder,fname(I_layer).name])
        end
        
    end
end
        
        
        
        
        