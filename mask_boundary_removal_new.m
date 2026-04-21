function mask_boundary_removal_new(filepath,savepath)
%load(filepath);
maskDir  = dir([filepath '*_mask*']); 
sort_nat_name = sort_nat({maskDir.name});
% img = imread([filepath sort_nat_name{1}]);
mkdir(savepath);
try
    OriginalSize=readmatrix([filepath,'bin.xls']);
catch
    [startIndex,~]=regexp(filepath,'\\time');
    filepath_bin=filepath(1:startIndex);
    OriginalSize=readmatrix([filepath_bin,'bin.xls']);
end
for x = 1:length(maskDir)
%     load([filepath sort_nat_name{x}]); 
    masks=imread([filepath sort_nat_name{x}]);
    mask_size = size(masks);
    mask_row = mask_size(1);
    mask_rank = mask_size(2);
    new_masks = masks;
    for i  = 2:mask_row-1   %1:mask_row
        for j = 2:mask_rank-1  %1:mask_rank
            a = masks(i,j);
            if a ~=0
                b = masks(i-1,j);
                c = masks(i+1,j);
                d = masks(i,j+1);
                e = masks(i,j-1);
                f = masks(i-1,j-1);
                g = masks(i-1,j+1);
                h = masks(i+1,j+1);
                k = masks(i+1,j-1);
                if a~=b | a~=c | a~=d | a~=e | a~=f | a~=g | a~=h | a~=k
                    new_masks(i,j) = 0;
                end
            end
        end
    end
    new_masks = imresize(new_masks,OriginalSize,'nearest');
    bw_masks = logical(new_masks);
    save([savepath,'bw_',sort_nat_name{x}(1:end-3),'mat'],'bw_masks')
end
end

                
                
        
        
            
            
    