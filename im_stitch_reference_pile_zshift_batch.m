function im_stitch_reference_pile_zshift_batch(ListPath)
EmbryoList=readcell(ListPath);
EmNum=size(EmbryoList,1)-1;
im1_path0=EmbryoList{2,1};
NucleusChannel1=EmbryoList{2,2};

im1_path0_stitch_reference=[im1_path0(1:end-1),'_new\'];
RefPath = [im1_path0(1:end-1),'_new\'];
FoldRoot=strfind(im1_path0,'\');
FoldNum=['MatchTo ',im1_path0(FoldRoot(end-1)+1:FoldRoot(end)-1)];
FoldZshiftNum=['MatchTo_Zshift_Batch ',im1_path0(FoldRoot(end-1)+1:FoldRoot(end)-1)];
Stitch_reference=load([im1_path0_stitch_reference,'stitch_parameter.mat']);
StitchLoc=size(Stitch_reference.x_start);
X_start=Stitch_reference.x_start';
Y_start=Stitch_reference.y_start';
Z_start=Stitch_reference.z_start';
    
ZShiftTile_start=[];StackNumTile=[];
for Emi=1:EmNum-1
    im2_path0=EmbryoList{2+Emi,1};
    NucleusChannel2=EmbryoList{2+Emi,2};

    rawpath=im2_path0;
    TileName  = dir([rawpath 'tile*']);
    save0Path = [rawpath(1:end-1),'_new\',FoldZshiftNum,'\'];mkdir(save0Path);
    RotateAngle=[];ShiftXY=[];ShiftZ=[];
    %% load shift&rotate
    for ti=1:size(TileName,1)
        tilePath=[rawpath,TileName(ti).name,'\'];
        rotatePath = [tilePath,'allstack\',FoldNum,'\'];
        Rct=readstruct([rotatePath,'Rectify.xml']);
        RotateAngle=[RotateAngle;Rct.rotate];
        ShiftXY=[ShiftXY;Rct.xy];
        ShiftZ=[ShiftZ;Rct.z];
    end
    ZShiftTile_start(:,Emi)=Z_start(:)-ShiftZ;

    imgDirN  = dir([tilePath 'time0001stack*.tif']);
    if size(imgDirN,1)==0
        imgDirN  = dir([tilePath 'stack*.tif']);
    end
    StackNumTile(:,Emi)=length(imgDirN);
end

imgDirN0  = dir([RefPath 'time0001stack*.tif']);
if size(imgDirN0,1)==0
    imgDirN0  = dir([RefPath 'stack*.tif']);
end
StackRefNum=length(imgDirN0);
StackNumTileMin=min(StackNumTile,[],2);

ZShiftTile_startMin=min(ZShiftTile_start,[],'all');
Cat0=0;
if ZShiftTile_startMin<=0
    Cat0=-ZShiftTile_startMin+1;
    ZShiftTile_start=ZShiftTile_start+Cat0;
    StackRefNum=StackRefNum-Cat0;
end
StackNumPile=StackNumTile-max(ZShiftTile_start,[],1);
if min(StackNumPile,[],2)<=StackRefNum
    StackNumPileSame=min(StackNumPile,[],2);
    Cat1=StackRefNum-StackNumPileSame;
else
    StackNumPileSame=StackRefNum;  
end

%reference cat
RefPathCopy=[RefPath,'Reference\'];mkdir(RefPathCopy);
sort_nat_name = sort_nat({imgDirN0.name});
SI=0;
for ii=(Cat0+1):(StackRefNum-Cat1)
    SI=SI+1;
    CopyNum0 = [RefPath sort_nat_name{ii}];
    CopyNum = [RefPathCopy sort_nat_name{SI}];
    copyfile(CopyNum0,CopyNum)
end


for Emi=1:EmNum-1
    im2_path0=EmbryoList{2+Emi,1};
    NucleusChannel2=EmbryoList{2+Emi,2};
    rawpath=im2_path0;
    TileName  = dir([rawpath 'tile*']);
    save0Path = [rawpath(1:end-1),'_new\',FoldZshiftNum,'\'];mkdir(save0Path);
    RotateAngle=[];ShiftXY=[];ShiftZ=[];
    %% load shift&rotate
    for ti=1:size(TileName,1)
        tilePath=[rawpath,TileName(ti).name,'\'];
        rotatePath = [tilePath,'allstack\',FoldNum,'\'];
        Rct=readstruct([rotatePath,'Rectify.xml']);
        RotateAngle=[RotateAngle;Rct.rotate];
        ShiftXY=[ShiftXY;Rct.xy];
        ShiftZ=[ShiftZ;Rct.z];
    end
%check size
img0 = imread([RefPath 'stack01.tif']);
imgtile0 = imread([tilePath 'stack01.tif']);
RefSize=size(img0);
OpenSize=size(imgtile0);
StitchImg=nan([RefSize(1,[1 2])+2*OpenSize(1,[1 2]),OpenSize(3)]);
RotateImgLogUse=ones(OpenSize);
StitchImgLog=~isnan(StitchImg);
StitchImgUse=StitchImg;
%z num
imgDirN  = dir([tilePath 'time0001stack*.tif']);
if size(imgDirN,1)==0
    imgDirN  = dir([tilePath 'stack*.tif']);
end
ImgChannel=cell(length(imgDirN),size(TileName,1));
%
xS=OpenSize(2);yS=OpenSize(1);

% StackNumTileMin=length(imgDirN);

% ZShiftTile_start(ZShiftTile_start<1)=1;
% ZShiftTile_start(ZShiftTile_start>StackNumTileMin)=StackNumTileMin;
% StackNumPile=StackNumTileMin-max(ZShiftTile_start);
StackNumPile=StackNumPileSame;

for ti=1:size(TileName,1)
    tilePath=[rawpath,TileName(ti).name,'\'];
    tilePath0=[im1_path0,TileName(ti).name,'\'];
    % savePath=[save0Path,TileName(ti).name,'\'];mkdir(savePath);
    imgDir  = dir([tilePath 'time0001stack*.tif']);
    if size(imgDir,1)==0
        imgDir  = dir([tilePath 'stack*.tif']);
    end
    img_all=[];
    sort_nat_name = sort_nat({imgDir.name});
    
    TileI=0;
    for ii=ZShiftTile_start(ti,Emi):(ZShiftTile_start(ti,Emi)+StackNumPile)
        TileI=TileI+1;
        img = imread([tilePath sort_nat_name{ii}]);
        TileSize0=size(img);
        TileSize0Half=round(TileSize0/2);
        imgR=imrotate(img,RotateAngle(ti));
        RotateImgLog=imrotate(RotateImgLogUse,RotateAngle(ti))==0;
        imgR=double(imgR);
        imgR(RotateImgLog)=nan;
        TileSizeR=size(imgR);
        % StartAdd=round((TileSizeR-TileSize0)/2);
        % imgRS = imtranslate(imgR,ShiftXY(ti,:),'OutputView','full');
        % ShiftXYCut=min([ShiftXY(ti,:);[0 0]],[],1);
        % StartAdd=StartAdd-[ShiftXYCut ShiftZ(ti)];

        tileCentre=[yS+TileSize0Half(1)-Y_start(ti)...
            xS+TileSize0Half(2)-X_start(ti)];
        UpscaleDirection=ShiftXY(ti,:)>=0;
        switch UpscaleDirection(2)
            case 0
                imgR=cat(2,imgR,nan(TileSizeR(1),2*abs(ShiftXY(ti,2)),TileSizeR(3)));
            case 1
                imgR=cat(2,nan(TileSizeR(1),2*abs(ShiftXY(ti,2)),TileSizeR(3)),imgR);
        end
        TileSizeR=size(imgR);
        switch UpscaleDirection(1)
            case 0
                imgR=cat(1,imgR,nan(2*abs(ShiftXY(ti,1)),TileSizeR(2),TileSizeR(3)));
            case 1
                imgR=cat(1,nan(2*abs(ShiftXY(ti,1)),TileSizeR(2),TileSizeR(3)),imgR);
        end
        TileSizeR=size(imgR);
        TileSizeRHalf=round(TileSizeR/2)-1;
        %imgtile00 = imread([tilePath0 sort_nat_name{ii}]);
        % targetSize = size(imgtile00(:,:,4));
        % imgR0 = imcrop(imgR(:,:,4),centerCropWindow2d(size(imgR(:,:,4)),targetSize));
        % figure;subplot(1,2,1);imshow(imgR0,[]);subplot(1,2,2);imshow(imgtile00(:,:,4),[])
        % figure;imshowpair(imgtile00(:,:,4),imgR0,'falsecolor');
        RowInterval=tileCentre(1)+[(-TileSizeRHalf(1)+1):(TileSizeR(1)-TileSizeRHalf(1))];
        ColInterval=tileCentre(2)+[(-TileSizeRHalf(2)+1):(TileSizeR(2)-TileSizeRHalf(2))];

        % for channel=1:TileSize0(3)
        if ti==1
            StitchImgUse=StitchImg;
        else
            StitchImgUse=ImgChannel{TileI,ti-1};
        end
        StitchImgLog=~isnan(StitchImgUse);

        StitchImgUseN=StitchImgUse;
        StitchImgUseN(RowInterval,ColInterval,:)=imgR(:,:,:);

            StitchImgUseN(StitchImgLog)=StitchImgUse(StitchImgLog);% out of memory?

        StitchImgUse=StitchImgUseN;
        ImgChannel{TileI,ti}=StitchImgUse;
        % end
        
    end
    % Release memory
    if ti>1
        for ii=1:(StackNumPile+1)
            ImgChannel{ii,ti-1}=[];
        end
    end


    if mod(ti+1,StitchLoc(2))==1||StitchLoc(2)==1
        yS=tileCentre(1)+TileSize0Half(1);
        xS=OpenSize(2);
    else
        xS=tileCentre(2)+TileSize0Half(2);
    end
end

try
    ImgShowRef=imread([im1_path0_stitch_reference,'allstack\','allstack.tif']);
catch
    im_allstack(im1_path0_stitch_reference,NucleusChannel1);
    ImgShowRef=imread([im1_path0_stitch_reference,'allstack\','allstack.tif']);
end
ImgShow=[];
for ii=1:(StackNumPile+1)
    ImgOut=ImgChannel{ii,end}(OpenSize(1)+1:end-OpenSize(1),OpenSize(2)+1:end-OpenSize(2),:);
    ImgOut(isnan(ImgOut))=0;
    ImgShow(:,:,ii)=ImgOut(:,:,NucleusChannel2);
    tiffwrite0(uint16(ImgOut),[save0Path,sort_nat_name{ii}])
end
ImgShowRec=max(ImgShow,[],3);
figure;imshowpair(ImgShowRef,ImgShowRec,'falsecolor');
saveas(gcf,[save0Path,'rectify_pair_pile','.fig']);
end