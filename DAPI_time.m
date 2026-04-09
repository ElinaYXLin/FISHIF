%% live dataset for training and testing (three scale) output folder: DAPI'index'_'size', such as DAPI1_75
clear
bhh_path='X:/baohuihan/FISHIF-Time/Confocal/';
[num,txt,raw]=xlsread([bhh_path,'Livedata.xlsx'],'cycle12');%% input Embryo information in Excel
size=75;%% input: cropped image size (nc11:81 301 601; nc12:75 251 401; nc13:71 201 401)
DAPI_th=1.3;%% adjust
for j=1:length(txt(:,1))
    Tstart=num(j,6);
    Tend=num(j,6)+num(j,1)-1;
for time=Tstart:Tend    
    mkdir([bhh_path,char(txt(j,1)),'/DAPI',sprintf('%d',num(j,4)),'_',sprintf('%d',size),'/',sprintf('%02d',time)]); 
Z=num(j,5);
W=num(j,7);
I=cell(Z,1);
se_centroid=cell(Z,1);
for z=1:Z
    I{z}=imread([bhh_path,char(txt(j,1)),'/result',sprintf('%d',num(j,4)),'/result',sprintf('%d',num(j,4)),'_',sprintf('%03d',z),'_',sprintf('%03d',time),'.tiff'],'tiff');%输入tiff
   
    I{z}=I{z}(:,:,1);
end
Iz=zeros(W,W,Z);
for z=1:Z
    Iz(:,:,z)=double(I{z});
end
Iz = Iz/max(Iz(:));
Imax=max(Iz,[],3);
% Imax = Imax/max(Imax(:));

se0=2;
Imax2 = imfilter(Imax,fspecial('gaussian',10,se0),'symmetric','conv');
%Imax2 = Imax2/max(Imax2(:));
Imax2 = (Imax2-min(Imax2(:)))/(max(Imax2(:))-min(Imax2(:)));
figure
imshow(Imax2)

DAPI = graythresh(Imax2);
BW = im2bw(Imax2,DAPI_th*DAPI);%确定BW
BW = imfill(BW,'holes');
se1=5;se2=2;
BW1 = imopen(BW,strel('disk',se1));
BW2 = imclose(BW1,strel('disk',se2));
BW2 = bwareaopen(BW2, 200);
BW3=watershed(BW2);
D= -bwdist(~BW2);
bw_cut = imextendedmin(D,2);
g2 = imimposemin(D,bw_cut);
bw_cut = imdilate(watershed(g2) == 0,strel('disk',1));% & (temp_sel);
BW3 = BW2 & (~ bw_cut);
BW3=BW2;
figure
imshow(BW3)


Inten = zeros(0);
for z = 1:Z

temp = regionprops(BW3,Iz(:,:,z),'MeanIntensity');

Inten(z,:) = [temp.MeanIntensity];

end

[~,zmax] = max(Inten);%确定zmax

se_mask2 = regionprops(BW3,'Centroid');
se_centroid = [se_mask2.Centroid];%Centroid

% se_mask3 = regionprops(BW,'BoundingBox');
% se_BoundingBox = [se_mask3.BoundingBox];%BoundingBox

%write out cut images
wide=(size-1)/2;
for i=1:length(zmax)
    centroid1=round(se_mask2(i).Centroid(1));
    centroid2=round(se_mask2(i).Centroid(2));
    
    if centroid1>wide&&centroid1<(W-wide)&&centroid2>wide&&centroid2<(W-wide)
    Iout=Iz(centroid2-wide:centroid2+wide,centroid1-wide:centroid1+wide,zmax(i));
    Iout = (Iout-min(Iout(:)))/(max(Iout(:))-min(Iout(:)));
    %Iout = Iout/max(Iout(:));
   
    imwrite(Iout,[bhh_path,char(txt(j,1)),'/DAPI',sprintf('%d',num(j,4)),'_',sprintf('%d',size),'/',sprintf('%02d',time),'/DAPI',sprintf('%d',num(j,4)),'_',sprintf('%03d',time),'_',sprintf('%03d',i),'.tiff'],'tiff');
    else
        continue
    end
end
end
end