function cellpose_bwmat(imgPath,savepath)
filename = 'mask.mat';                                % 图像库路径
imgDir  = dir([imgPath '*bw*']); % 遍历所有seg_mat格式文件
sort_nat_name = sort_nat({imgDir.name});
img = cell2mat(struct2cell(load([imgPath sort_nat_name{1}])));
mask_stack = cat(3,img);
for i = 1:length(imgDir) 
    img = cell2mat(struct2cell(load([imgPath sort_nat_name{i}])));
    mask_stack(:,:,i) = img;
end
%mask_stack = double(mask_stack);
mask_stack = label3D_waitbar(mask_stack,3);
mkdir(savepath);
save([savepath,filename],'mask_stack')
end