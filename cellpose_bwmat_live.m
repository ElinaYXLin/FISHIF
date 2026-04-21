function cellpose_bwmat_live(imgPath,savepath)
%imgPath = 'E:\data\Jianing\20210719\stacks\14_6_20210719\tile01\'; 
%savepath = 'E:\data\Jianing\20210719\masks\14_6_20210719_new\'
filename = 'mask_cellpose_stack_bw';                                % 图像库路径
imgDir  = dir([imgPath '*bw_mask*']); % 遍历所有seg_mat格式文件
sort_nat_name = sort_nat({imgDir.name});
file_size = length(imgDir);
total_time = str2num(sort_nat_name{file_size}(13:end-15));
stack_size = file_size/total_time;
%img = cell2mat(struct2cell(load([imgPath sort_nat_name{1}])));
%mask_stack = cat(3,img);
for t = 1:total_time
    for i = stack_size*(t-1)+1:stack_size*t
        img = cell2mat(struct2cell(load([imgPath sort_nat_name{i}])));
        mask_stack = cat(3,img);
        mask_stack(:,:,i) = img;
    end
    %mask_stack = label3D(mask_stack,3);
    save([savepath,filename,num2str(t)],'mask_stack')
    clear mask_stack;
end
end
