function cellpose_bwmat_live_while(imgPath,savepath,thickness)
%imgPath = 'E:\data\Jianing\20210719\stacks\14_6_20210719\tile01\'; 
%savepath = 'E:\data\Jianing\20210719\masks\14_6_20210719_new\'
filename = 'mask.mat';                                % 图像库路径
mkdir(savepath);
imgDir  = dir([imgPath '*bw*']); % 遍历所有seg_mat格式文件
sort_nat_name = sort_nat({imgDir.name});
file_size = length(imgDir);
stack_size = size(dir([imgPath 'bw_time0001*']),1);
if stack_size==0
    stack_size = size(dir([imgPath 'bw_*']),1);
end
total_time = file_size/stack_size;
if stack_size<thickness
    thickness=1;
end
%img = cell2mat(struct2cell(load([imgPath sort_nat_name{1}])));
%mask_stack = cat(3,img);
t = 1;
tic;
h_wait=waitbar(0,'Please wait...','Name','Label3D');
set(h_wait,'unit','normalized','position',[0.4,0.4,0.14,0.08])
while t <= total_time
    exp_img = cell2mat(struct2cell(load([imgPath sort_nat_name{1}])));
    size_exp = size(exp_img);
    mask_stack = false(size_exp(1),size_exp(2),stack_size);
    for i = stack_size*(t-1)+1:stack_size*t
        img = cell2mat(struct2cell(load([imgPath sort_nat_name{i}])));
%         mask_stack = cat(3,img);
        mask_stack(:,:,(i-stack_size*(t-1))) = img;
    end
    mask_stack = label3D(mask_stack,thickness);
    save([savepath,'time',num2str(t,'%04u'),filename],'mask_stack','-v7.3')
    if total_time==1
        save([savepath,filename],'mask_stack','-v7.3')
    end
%     clear mask_stack;
    t = t+1;
    escape_time=toc;
    waitbar((t-1)/total_time,h_wait,['Progress',num2str((t-1)/total_time,'%.2f'),char(13,10)',...
        'Escape time: ',num2str(escape_time/3600,'%.3f'),' hours',char(13,10)',...
        'Remaining time required: ',num2str((total_time-(t-1))*escape_time/3600/(t-1),'%.3f'),' hours'])
end
close(h_wait)
end