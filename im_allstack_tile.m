function im_allstack_tile(rawpath,channel_num)
TileName  = dir([rawpath 'tile*']);
for ti=1:size(TileName,1)
    tilePath=[rawpath,TileName(ti).name,'\'];
    mkdir(tilePath,'allstack');
    savepath = [tilePath,'allstack\'];
    imgDir  = dir([tilePath 'time0001stack*.tif']);
    if size(imgDir,1)==0
        imgDir  = dir([tilePath 'stack*.tif']);
    end
    img_all=[];
    sort_nat_name = sort_nat({imgDir.name});
    for i=1:length(imgDir)
        img = imread([tilePath sort_nat_name{i}]);
        img = img(:,:,channel_num);
        img_all(:,:,i)=im2double(img);
    end
    img_max=im2uint16(max(img_all,[],3));
    imwrite(img_max,[savepath,'allstack.tif'])
end
end