function im_tile_rotate(rawpath)
TileName  = dir([rawpath 'tile*']);
save0Path = [rawpath,'Rotate\'];mkdir(save0Path);
RotateAngle=[];
for ti=1:size(TileName,1)
    tilePath=[rawpath,TileName(ti).name,'\'];
    rotatePath = [tilePath,'allstack\'];
    Rotate=readstruct([rotatePath,'Rectify.xml']);
    RotateAngle=[RotateAngle Rotate.rotate];
end
AngleMedian=median(RotateAngle,2);
RotateAngleDo=RotateAngle-AngleMedian;

for ti=1:size(TileName,1)
    tilePath=[rawpath,TileName(ti).name,'\'];
    savePath=[save0Path,TileName(ti).name,'\'];mkdir(savePath);
    imgDir  = dir([tilePath 'time0001stack*.tif']);
    if size(imgDir,1)==0
        imgDir  = dir([tilePath 'stack*.tif']);
    end
    img_all=[];
    sort_nat_name = sort_nat({imgDir.name});
    for ii=1:length(imgDir)
        img = imread([tilePath sort_nat_name{ii}]);
        imgR=imrotate(img,RotateAngleDo(ti),"crop");
        imwrite(imgR,[savePath,sort_nat_name{ii}],'Compression','none')
    end
end
end