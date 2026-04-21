function im_allstack(rawpath,channel_num)
mkdir(rawpath,'allstack');
savepath = [rawpath,'allstack\'];
imgDir  = dir([rawpath 'time0001stack*.tif']);
if size(imgDir,1)==0
    imgDir  = dir([rawpath 'stack*.tif']);
end
img_all=[];
sort_nat_name = sort_nat({imgDir.name});
for i=1:length(imgDir)
    img = imread([rawpath sort_nat_name{i}]);
    img = img(:,:,channel_num);
    img_all(:,:,i)=im2double(img);
end
img_max=im2uint16(max(img_all,[],3));
imwrite(img_max,[savepath,'allstack.tif'])
end