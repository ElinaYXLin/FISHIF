function Runcellpose=cellpose_raw2input(rawpath,resize_ratio,channel_num)

close all
tic

mkdir(rawpath,'raw2input');
%rawpath = 'E:\data\'
savepath = [rawpath,'raw2input\'];
%resize_ratio = 1/5
%mkdir(rawpath,'raw2input');
imgDir  = dir([rawpath '*.tif']);
sort_nat_name = sort_nat({imgDir.name});
%img = cell2mat(struct2cell(load([imgPath sort_nat_name{1}])));
for i = 1:length(imgDir) 
    img = imread([rawpath sort_nat_name{i}]);
    img = img(:,:,channel_num);
    img_resize = imresize(img,resize_ratio,'nearest');
    imwrite(img_resize,[savepath,sort_nat_name{i}])
end
bin = size(img);
writematrix(bin,[savepath,'bin.xls']);
Runcellpose=['python -m cellpose --dir "',savepath(1:end-1),'"  '];
end
