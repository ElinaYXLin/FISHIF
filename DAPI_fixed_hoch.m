%% fixed data (23651) DNA channel out for relayed training and testing (output folder:DAPI)
clear
figure
bhh_path='X:/baohuihan/FISHIF-Time/Confocal/';
[num,txt,raw]=xlsread([bhh_path,'Fixeddata.xlsx'],'cycle11');%% input Embryo information in Excel
DAPI_th=1.8;%%adjust
T=11;%% inputfor nuclear cycle
if T==11
        size1=70;
        size2=259;
        size3=517;          
elseif T==12
        size1=66;
        size2=221;
        size3=353;            
elseif T==13
        size1=63;
        size2=181;
        size3=361;
end
 
for n=1:length(txt)
    
    load([bhh_path,char(txt(n,1)),'/Results/',char(txt(n,2)),'_new.mat']);
      mkdir([bhh_path,char(txt(n,1)),'/DAPI/DAPI_Hoch_s_w/',char(txt(n,2))]); 
      mkdir([bhh_path,char(txt(n,1)),'/DAPI/DAPI_Hoch_m_w/',char(txt(n,2))]); 
      mkdir([bhh_path,char(txt(n,1)),'/DAPI/DAPI_Hoch_l_w/',char(txt(n,2))]);
      mkdir([bhh_path,char(txt(n,1)),'/DAPI/DAPI_Hoch_s_Aw/',char(txt(n,2))]); 
      mkdir([bhh_path,char(txt(n,1)),'/DAPI/DAPI_Hoch_m_Aw/',char(txt(n,2))]); 
      mkdir([bhh_path,char(txt(n,1)),'/DAPI/DAPI_Hoch_l_Aw/',char(txt(n,2))]);
      mkdir([bhh_path,char(txt(n,1)),'/DAPI/DAPI_Hoch_s_Mw/',char(txt(n,2))]); 
      mkdir([bhh_path,char(txt(n,1)),'/DAPI/DAPI_Hoch_m_Mw/',char(txt(n,2))]); 
      mkdir([bhh_path,char(txt(n,1)),'/DAPI/DAPI_Hoch_l_Mw/',char(txt(n,2))]);
      mkdir([bhh_path,char(txt(n,1)),'/DAPI/DAPI_Hoch_s_Pw/',char(txt(n,2))]); 
      mkdir([bhh_path,char(txt(n,1)),'/DAPI/DAPI_Hoch_m_Pw/',char(txt(n,2))]); 
      mkdir([bhh_path,char(txt(n,1)),'/DAPI/DAPI_Hoch_l_Pw/',char(txt(n,2))]);

    %files = dir(fullfile([bhh_path,char(txt(n,2)),'/stacks/',char(txt(n,2)),'_new/'],'*.tif'));
    files = dir(fullfile([bhh_path,char(txt(n,1)),'/stacks/',char(txt(n,2)),'_new/'],'*.tif'));
    Z = size(files,1);
    clear I Ir Ih Iz Ihz Iz1 Imax Imax2 temp Inten Pixel Std se_centroid
%Z=11;
%W=512;
I=cell(Z,1);
Ir=cell(Z,1);
Ih=cell(Z,1);
se_centroid=cell(Z,1);
 [num2,txt2,raw2]=xlsread([bhh_path,char(txt(n,1)),'/stacks/matchlist.xls']);
for z=1:Z
    I{z}=imread([bhh_path,char(txt(n,1)),'/stacks/',char(txt(n,2)),'_new/stack',sprintf('%02d',z),'.tif'],'tif');%输入tiff
    Ir{z}=I{z}(:,:,1);
    Ih{z}=I{z}(:,:,3);
end
%Iz=zeros(W,W,Z);
for z=1:Z
    Iz(:,:,z)=double(Ir{z});
    Ihz(:,:,z)=double(Ih{z});
end
 em_mask = imresize(em_mask, num2(1,9)/0.132, 'nearest');
 Iz = imresize(Iz, num2(1,9)/0.132, 'nearest');%0.538,0.5
 Ihz = imresize(Ihz, num2(1,9)/0.132, 'nearest');%0.538,0.5

em_mask1=imerode(em_mask,strel('disk',150));
Iz = Iz/max(Iz(:));
Ihz = Ihz/max(Ihz(:));
Iz1=Iz(:,:,:);%focus
Imax=max(Iz1,[],3);
% Imax = Imax/max(Imax(:));
Iz2=Ihz(:,:,:);%focus
Imaxh=max(Iz2,[],3);
%% BW watershed
se0=2;%3
Imax2 = imfilter(Imax,fspecial('gaussian',10,se0),'symmetric','conv');
%Imax2 = Imax2/max(Imax2(:));
Imax2 = (Imax2-min(Imax2(:)))/(max(Imax2(:))-min(Imax2(:)));

DAPI = graythresh(Imax2);
DAPI1=DAPI_th*DAPI;
BW = im2bw(Imax2,DAPI1);%确定BW
BW = imfill(BW,'holes');
se1=4;se2=2;
BW1 = imopen(BW,strel('disk',se1));
BW2 = imclose(BW1,strel('disk',se2));
BW2 = imfill(BW2,'holes');
BW3=BW2;
% D= -bwdist(~BW2);
% bw_cut = imextendedmin(D,2);
% g2 = imimposemin(D,bw_cut);
% bw_cut = imdilate(watershed(g2) == 0,strel('disk',1));% & (temp_sel);
% BW3 = BW2 & (~ bw_cut);
BW3 = bwareaopen(BW3, 300);
imshow(BW3)

Inten = zeros(0);
Std = zeros(0);
Area = zeros(0);
for z = 1:Z

temp = regionprops(BW3,Ihz(:,:,z),'MeanIntensity');%Iz
%temp = regionprops(BW3,Iz(:,:,z),'MeanIntensity');

Inten(z,:) = [temp.MeanIntensity];

Pixel=regionprops(BW3,Ihz(:,:,z),'PixelValues');
%Pixel=regionprops(BW3,Iz(:,:,z),'PixelValues');
for p=1:length(Pixel)
    Pixel(p).pixelvalues=std(Pixel(p).PixelValues);
    Pixel(p).Area=0;
    for j=1:length(Pixel(p).PixelValues)
        if Pixel(p).PixelValues(j)>1*DAPI
            Pixel(p).Area=Pixel(p).Area+1;
        end
    end
    
end
Std(z,:) = [Pixel.pixelvalues];
Area(z,:) = [Pixel.Area];

end

[S,zs]=max(Std);
% for j=1:length(S)
%    if S(j)>1.2*median(S)
%         Inten(zs(j),j)=0;
%    end
% end
[~,zmax] = max(Area);
[MaxI,zmax2] = max(Inten);%确定zmax


se_mask2 = regionprops(BW3,'Centroid');
se_centroid = [se_mask2.Centroid];%Centroid

% se_mask3 = regionprops(BW,'BoundingBox');
% se_BoundingBox = [se_mask3.BoundingBox];%BoundingBox


[W1,W2]=size(Imax);
w1=(size1-1)/2;
w2=(size2-1)/2;
w3=(size3-1)/2;
wide=w3;%31.45 
%% define A M P
EL_info=get_EL(em_mask);
L0=[EL_info(1) EL_info(2)];
L1=[EL_info(3) EL_info(4)];
L2=L0-L1;
La=L1+L2*0.18;
Lm=L1+L2*0.5;
Lp=L1+L2*0.82;
for i=1:length(se_mask2)
   
    centroid1=round(se_mask2(i).Centroid(1));
    centroid2=round(se_mask2(i).Centroid(2));
    %% whole datasets for train
    if centroid1>wide&&centroid1<(W2-wide)&&centroid2>wide&&centroid2<(W1-wide)&&em_mask1(centroid2,centroid1)==1&&MaxI(i)>0.6*median(MaxI)&&abs(zmax(i)-zmax2(i))<2%&&em_mask1(centroid2,centroid1)==1%&&zmax(i)<z
    Iout1=Ihz(centroid2-w1:centroid2+w1,centroid1-w1:centroid1+w1,zmax2(i));
    Iout2=Ihz(centroid2-w2:centroid2+w2,centroid1-w2:centroid1+w2,zmax2(i));
    Iout3=Ihz(centroid2-w3:centroid2+w3,centroid1-w3:centroid1+w3,zmax2(i));
    Iout1 = 1*(Iout1-min(Iout1(:)))/(max(Iout1(:))-min(Iout1(:)));
    Iout2 = 1*(Iout2-min(Iout2(:)))/(max(Iout2(:))-min(Iout2(:)));
    Iout3 = 1*(Iout3-min(Iout3(:)))/(max(Iout3(:))-min(Iout3(:)));

     imwrite(Iout1,[bhh_path,char(txt(n,1)),'/DAPI/DAPI_Hoch_s_w/',char(txt(n,2)),'/DAPI_',sprintf('%03d',i),'.tiff'],'tiff');
      imwrite(Iout2,[bhh_path,char(txt(n,1)),'/DAPI/DAPI_Hoch_m_w/',char(txt(n,2)),'/DAPI_',sprintf('%03d',i),'.tiff'],'tiff');
     imwrite(Iout3,[bhh_path,char(txt(n,1)),'/DAPI/DAPI_Hoch_l_w/',char(txt(n,2)),'/DAPI_',sprintf('%03d',i),'.tiff'],'tiff');

    else
        continue
    end
    %% three region datasets for test
     if centroid1>wide&&centroid1<(W2-wide)&&centroid2>wide&&centroid2<(W1-wide)&&centroid1>(La(1)-Lwide)&&centroid1<(La(1)+Lwide)&&centroid2>(La(2)-Lwide)&&centroid2<(La(2)+Lwide)&&em_mask1(centroid2,centroid1)==1&&MaxI(i)>0.6*median(MaxI)&&abs(zmax(i)-zmax2(i))<2%&&em_mask1(centroid2,centroid1)==1%&&zmax(i)<z
    Iout1=Ihz(centroid2-w1:centroid2+w1,centroid1-w1:centroid1+w1,zmax2(i));
    Iout2=Ihz(centroid2-w2:centroid2+w2,centroid1-w2:centroid1+w2,zmax2(i));
    Iout3=Ihz(centroid2-w3:centroid2+w3,centroid1-w3:centroid1+w3,zmax2(i));
    Iout1 = 1*(Iout1-min(Iout1(:)))/(max(Iout1(:))-min(Iout1(:)));
    Iout2 = 1*(Iout2-min(Iout2(:)))/(max(Iout2(:))-min(Iout2(:)));
    Iout3 = 1*(Iout3-min(Iout3(:)))/(max(Iout3(:))-min(Iout3(:)));
     imwrite(Iout1,[bhh_path,char(txt(n,1)),'/DAPI/DAPI_Hoch_s_Aw/',char(txt(n,2)),'/DAPI_',sprintf('%03d',i),'.tiff'],'tiff');
     imwrite(Iout2,[bhh_path,char(txt(n,1)),'/DAPI/DAPI_Hoch_m_Aw/',char(txt(n,2)),'/DAPI_',sprintf('%03d',i),'.tiff'],'tiff');
     imwrite(Iout3,[bhh_path,char(txt(n,1)),'/DAPI/DAPI_Hoch_l_Aw/',char(txt(n,2)),'/DAPI_',sprintf('%03d',i),'.tiff'],'tiff');
    elseif centroid1>wide&&centroid1<(W2-wide)&&centroid2>wide&&centroid2<(W1-wide)&&centroid1>(Lp(1)-Lwide)&&centroid1<(Lp(1)+Lwide)&&centroid2>(Lp(2)-Lwide)&&centroid2<(Lp(2)+Lwide)&&em_mask1(centroid2,centroid1)==1&&MaxI(i)>0.6*median(MaxI)&&abs(zmax(i)-zmax2(i))<2%&&em_mask1(centroid2,centroid1)==1%&&zmax(i)<z
    Iout1=Ihz(centroid2-w1:centroid2+w1,centroid1-w1:centroid1+w1,zmax2(i));
    Iout2=Ihz(centroid2-w2:centroid2+w2,centroid1-w2:centroid1+w2,zmax2(i));
    Iout3=Ihz(centroid2-w3:centroid2+w3,centroid1-w3:centroid1+w3,zmax2(i));
    Iout1 = 1*(Iout1-min(Iout1(:)))/(max(Iout1(:))-min(Iout1(:)));
    Iout2 = 1*(Iout2-min(Iout2(:)))/(max(Iout2(:))-min(Iout2(:)));
    Iout3 = 1*(Iout3-min(Iout3(:)))/(max(Iout3(:))-min(Iout3(:)));
     imwrite(Iout1,[bhh_path,char(txt(n,1)),'/DAPI/DAPI_Hoch_s_Pw/',char(txt(n,2)),'/DAPI_',sprintf('%03d',i),'.tiff'],'tiff');
     imwrite(Iout2,[bhh_path,char(txt(n,1)),'/DAPI/DAPI_Hoch_m_Pw/',char(txt(n,2)),'/DAPI_',sprintf('%03d',i),'.tiff'],'tiff');
     imwrite(Iout3,[bhh_path,char(txt(n,1)),'/DAPI/DAPI_Hoch_l_Pw/',char(txt(n,2)),'/DAPI_',sprintf('%03d',i),'.tiff'],'tiff');
    elseif centroid1>wide&&centroid1<(W2-wide)&&centroid2>wide&&centroid2<(W1-wide)&&centroid1>(Lm(1)-Lwide)&&centroid1<(Lm(1)+Lwide)&&centroid2>(Lm(2)-Lwide)&&centroid2<(Lm(2)+Lwide)&&em_mask1(centroid2,centroid1)==1&&MaxI(i)>0.6*median(MaxI)&&abs(zmax(i)-zmax2(i))<2%&&em_mask1(centroid2,centroid1)==1%&&zmax(i)<z
    Iout1=Ihz(centroid2-w1:centroid2+w1,centroid1-w1:centroid1+w1,zmax2(i));
    Iout2=Ihz(centroid2-w2:centroid2+w2,centroid1-w2:centroid1+w2,zmax2(i));
    Iout3=Ihz(centroid2-w3:centroid2+w3,centroid1-w3:centroid1+w3,zmax2(i));
    Iout1 = 1*(Iout1-min(Iout1(:)))/(max(Iout1(:))-min(Iout1(:)));
    Iout2 = 1*(Iout2-min(Iout2(:)))/(max(Iout2(:))-min(Iout2(:)));
    Iout3 = 1*(Iout3-min(Iout3(:)))/(max(Iout3(:))-min(Iout3(:)));
     imwrite(Iout1,[bhh_path,char(txt(n,1)),'/DAPI/DAPI_Hoch_s_Mw/',char(txt(n,2)),'/DAPI_',sprintf('%03d',i),'.tiff'],'tiff');
     imwrite(Iout2,[bhh_path,char(txt(n,1)),'/DAPI/DAPI_Hoch_m_Mw/',char(txt(n,2)),'/DAPI_',sprintf('%03d',i),'.tiff'],'tiff');
     imwrite(Iout3,[bhh_path,char(txt(n,1)),'/DAPI/DAPI_Hoch_l_Mw/',char(txt(n,2)),'/DAPI_',sprintf('%03d',i),'.tiff'],'tiff');
    else
        continue
    end
end
end


