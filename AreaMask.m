function AreaMask(threshold_single,threshold_whole,edge,AreaRange,DAPI_channel,maskpath,stackpath)
% threshold=0.3;DAPI_channel=4;
% maskpath='Y:\Imaging\S1\05222020S1_60X\masks\05222020-bcd-488-hb1-TMR-hb7-647_60X_2020_05_22__10_33_53_001_new\';
% stackpath='Y:\Imaging\S1\05222020S1_60X\stacks\05222020-bcd-488-hb1-TMR-hb7-647_60X_2020_05_22__10_33_53_001_new\';
mask_nmae  = dir([maskpath,'*mask.mat']);
time_dim = size(dir([stackpath,'*stack01.tif']),1);
bw_path=[stackpath,'Intermediate\'];
for time_index=1:time_dim
    try
        load([maskpath,'time',num2str(time_index,'%04u'),'mask.mat']);
    end
    stack_name=dir([stackpath,'time',num2str(time_index,'%04u'),'stack*.tif']);
    if size(stack_name,1)==0 && time_dim==1
        stack_name=dir([stackpath,'stack*.tif']);
        load([maskpath,'mask.mat']);
    end
    sgAreaG=[];
    for I_layer = 1:size(mask_stack,3)
        sgArea = regionprops(double(mask_stack(:,:,I_layer)),'Area');
        sgAreaG=[sgAreaG,[sgArea.Area]];
    end
    [sgAreaD,Edges]=histcounts(sgAreaG,100);
    figure
    set(gcf,'Position',[100 100 2000 1200])
    bar(Edges(1:end-1)+(Edges(2)-Edges(1))/2,sgAreaD,'hist')
    ylim([0 1.2*max(sgAreaD(2:end))])
end
end
    