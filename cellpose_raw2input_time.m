function [Runcellpose,time_dim]=cellpose_raw2input_time(rawpath,resize_ratio,channel_num)
mkdir(rawpath,'raw2input');
%rawpath = 'E:\data\'
time_dim=size(dir([rawpath '*stack01.tif']),1);
Runcellpose=cell(1,time_dim);
for timeIndex=1:time_dim
    imgDir  = dir([rawpath,'time',num2str(timeIndex,'%04u'),'stack*.tif']);
    mkdir([rawpath,'raw2input\','time',num2str(timeIndex,'%04u')]);
    savepath = [rawpath,'raw2input\','time',num2str(timeIndex,'%04u'),'\'];
    sort_nat_name = sort_nat({imgDir.name});
    %img = cell2mat(struct2cell(load([imgPath sort_nat_name{1}])));
    for i = 1:length(imgDir) 
        img = imread([rawpath sort_nat_name{i}]);
        img = img(:,:,channel_num);
        img_resize = imresize(img,resize_ratio,'nearest');
        imwrite(img_resize,[savepath,sort_nat_name{i}])
    end
    bin = size(img);
    Runcellpose{1,timeIndex}=['python -m cellpose --dir "',savepath(1:end-1),'"  '];
end
writematrix(bin,[rawpath,'raw2input\','bin.xls']);
end
