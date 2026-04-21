function Runcellpose=cellpose_raw2input_try(rawpath,timeOnOff,resize_ratio,channel_num,Model)
mkdir(rawpath,'raw2input_try');
%rawpath = 'E:\data\'
savepath = [rawpath,'raw2input_try\'];
%resize_ratio = 1/5
%mkdir(rawpath,'raw2input');
time_dim=1;
if timeOnOff
    time_dim=size(dir([rawpath '*stack01.tif']),1);
end
for timeIndex=1:time_dim
    imgDir  = dir([rawpath,'time',num2str(timeIndex,'%04u'),'stack*.tif']);
    if size(imgDir,1)==0
        imgDir  = dir([rawpath 'stack*.tif']);
    end
    img_all=[];
    if strcmp(Model,'All')
        sort_nat_name = sort_nat({imgDir.name});
        for i=1:length(imgDir)
            img = imread([rawpath sort_nat_name{i}]);
            img = img(:,:,channel_num);
            img_all(:,:,i)=im2double(img);
        end
        img_max=im2uint16(max(img_all,[],3));
        img_resize = imresize(img_max,resize_ratio,'nearest');
        imwrite(img_resize,[savepath,'time',num2str(timeIndex,'%04u'),'Try.tif'])
    else
        imgDir=imgDir(ceil(size(imgDir,1)/2),:);
        sort_nat_name = sort_nat({imgDir.name});
        img = imread([rawpath sort_nat_name{1}]);
        img = img(:,:,channel_num);
        img_resize = imresize(img,resize_ratio,'nearest');
        imwrite(img_resize,[savepath,'time',num2str(timeIndex,'%04u'),'Try.tif'])
    end
end
bin = size(img);
writematrix(bin,[savepath,'bin.xls']);
Runcellpose=['python -m cellpose --dir "',savepath(1:end-1),'"  '];
end