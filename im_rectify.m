function [x_shift,y_shift,rotate_angle]=im_rectify(im1_path,im2_path,channel_num)
% uiopen('Y:\Imaging\S1\12232019_S1_60X\stacks\12232019_S1_60X_005_new\stack14.tif',1)
% uiopen('Y:\Imaging\S1\12232019_S1_60X\stacks\12232019_S1_60X_005_new\stack15.tif',1)
% stack1=stack14(:,:,4);
% stack2=stack15(:,:,4);
% stack2_ro=imrotate(stack2,63.4);
% imshow(stack2_ro,[])
% figure;
% imshow(stack2_ro,[])
% figure;imshow(stack1,[])
stack1=imread([im1_path,'allstack\','allstack.tif']);
stack2_ro=imread([im2_path,'allstack\','allstack.tif']);
%% x is row;y is col;
%% cut image of path2 to path1 size
SizeIm1=size(stack1);
SizeIm2=size(stack2_ro);
CatLog=SizeIm1-SizeIm2;
% row
Cat=CatLog(1);
if Cat > 0
    stack2_ro=cat(1,stack2_ro,zeros(Cat,size(stack2_ro,2)));
elseif Cat < 0 
    stack2_ro(end+Cat+1:end,:)=[];
end
%col
Cat=CatLog(2);
if Cat > 0
    stack2_ro=cat(2,stack2_ro,zeros(size(stack2_ro,1),Cat));
elseif Cat < 0 
    stack2_ro(:,end+Cat+1:end)=[];
end
% figure;imshowpair(stack1,stack2_ro,'falsecolor');
%% INI
scan_reg=0.6;
bin_scan=2;
rotate_vector=[-30:2:30];
Resize=0.1;
%% Rough scanning
im_1=imresize(stack1,Resize);
im_2=imresize(stack2_ro,Resize);
[L1,L2]=size(im_1);
[l1,l2]=size(im_2);
cat_size=floor(min([l1,l2])*scan_reg);
cat_size_2=floor(cat_size/2);
im_scan=im_2((1+floor((l1-cat_size)/2)):(1+floor((l1-cat_size)/2)+cat_size),...
    (1+floor((l2-cat_size)/2)):(1+floor((l2-cat_size)/2)+cat_size));
% im_scan=im_2((1+floor(l1*scan_cat)):(end-floor(l1*scan_cat)),...
%     (1+floor(l2*scan_cat)):(end-floor(l2*scan_cat)));
xx_vector=1:bin_scan:(L1-cat_size_2);
yy_vector=1:bin_scan:(L2-cat_size_2);
% time_l=0.1*size(xx_vector,2)*size(yy_vector,2)*size(rotate_vector,2)/3600;
% im_1=im2double(im_1);
ssim_mat=[];
% parpool(100);
parfor angle_index=1:size(rotate_vector,2)
    rotate_angle=rotate_vector(angle_index);
    xx_index=0;
    im_scan_r=imrotate(im_scan,rotate_angle);
    [l1,l2]=size(im_scan_r);
    im_scan_use=im_scan_r((1+floor((l1-cat_size_2)/2)):(1+floor((l1-cat_size_2)/2)+cat_size_2),...
    (1+floor((l2-cat_size_2)/2)):(1+floor((l2-cat_size_2)/2)+cat_size_2));
%     im_scan_use=im2double(im_scan_use);
    ssim_xy=[];
    for xx=xx_vector
        xx_index=xx_index+1;
        yy_index=0;
        for yy=yy_vector
            yy_index=yy_index+1;
            im_match=im_1(xx:(xx+cat_size_2),yy:(yy+cat_size_2));
            ssimval = ssim(im_match,im_scan_use);
            ssim_xy(xx_index,yy_index)=ssimval;
        end 
    end
    ssim_mat(:,:,angle_index)=ssim_xy;
end

[k1,k2,k3]=ind2sub(size(ssim_mat),find(ssim_mat==max(ssim_mat(:))));
xx0=xx_vector(k1(1));yy0=yy_vector(k2(1));rotate_angle0=rotate_vector(k3(1));
%% Micromesh scanning
im_1=stack1;
im_2=stack2_ro;
[L1,L2]=size(im_1);
[l1,l2]=size(im_2);
scan_reg=0.6;cat_size=floor(min([l1,l2])*scan_reg);
cat_size_2=floor(cat_size/2);
im_scan=im_2((1+floor((l1-cat_size)/2)):(1+floor((l1-cat_size)/2)+cat_size),...
    (1+floor((l2-cat_size)/2)):(1+floor((l2-cat_size)/2)+cat_size));

x1=max([floor((xx0-bin_scan)/Resize),1]);
y1=max([floor((yy0-bin_scan)/Resize),1]);
x2=min([floor((xx0+bin_scan)/Resize),L1]);
y2=min([floor((yy0+bin_scan)/Resize),L2]);
rotate_vector=[(rotate_angle0-2):0.1:(rotate_angle0+2)];
xx_vector=x1:1:x2;
yy_vector=y1:1:y2;
% time_l=0.1*size(xx_vector,2)*size(yy_vector,2)*size(rotate_vector,2)/3600;
% im_1=im2double(im_1);
ssim_mat_micr=zeros(size(xx_vector,2),size(yy_vector,2),size(rotate_vector,2));
% parpool(100);
parfor angle_index=1:size(rotate_vector,2)
    rotate_angle=rotate_vector(angle_index);
    xx_index=0;
    im_scan_r=imrotate(im_scan,rotate_angle);
    [l1,l2]=size(im_scan_r);
    im_scan_use=im_scan_r((1+floor((l1-cat_size_2)/2)):(1+floor((l1-cat_size_2)/2)+cat_size_2),...
    (1+floor((l2-cat_size_2)/2)):(1+floor((l2-cat_size_2)/2)+cat_size_2));
%     im_scan_use=im2double(im_scan_use);
    ssim_xy=zeros(size(xx_vector,2),size(yy_vector,2));
    for xx=xx_vector
        xx_index=xx_index+1;
        yy_index=0;
        for yy=yy_vector
            yy_index=yy_index+1;
            try
                im_match=im_1(xx:(xx+cat_size_2),yy:(yy+cat_size_2));
            catch
                continue
            end
            ssimval = ssim(im_match,im_scan_use);
            ssim_xy(xx_index,yy_index)=ssimval;
        end 
    end
    ssim_mat_micr(:,:,angle_index)=ssim_xy;
end

[k1,k2,k3]=ind2sub(size(ssim_mat_micr),find(ssim_mat_micr==max(ssim_mat_micr(:))));
xx=xx_vector(k1(1));yy=yy_vector(k2(1));rotate_angle=rotate_vector(k3(1));
% ssim_mat2=ssim_mat;ssim_mat2(k1,k2,k3)=0;
% [k1_2,k2_2,k3_2]=ind2sub(size(ssim_mat2),find(ssim_mat2==max(ssim_mat2(:))));
% xx2=xx_vector(k1_2);yy2=yy_vector(k2_2);rotate_angle2=rotate_vector(k3_2);
% xx=floor((L1-cat_size_2)/2);yy=floor((L2-cat_size_2)/2);rotate_angle=-63.4;\

im_match=im_1(xx:(xx+cat_size_2),yy:(yy+cat_size_2));
im_scan_r=imrotate(im_scan,rotate_angle);
[l1,l2]=size(im_scan_r);
im_scan_use=im_scan_r((1+floor((l1-cat_size_2)/2)):(1+floor((l1-cat_size_2)/2)+cat_size_2),...
(1+floor((l2-cat_size_2)/2)):(1+floor((l2-cat_size_2)/2)+cat_size_2));
im_scan_use=im2double(im_scan_use);
figure;subplot(1,2,1);imshow(im_match,[]);subplot(1,2,2);imshow(im_scan_use,[]);
saveas(gcf,[im2_path,'allstack\','rectify_scan','.fig']);
close(gcf)
figure;imshowpair(im_match,im_scan_use,'falsecolor');
saveas(gcf,[im2_path,'allstack\','rectify_scan_pair','.fig']);
close(gcf)

im1_center_x=floor((L1-cat_size_2)/2);im1_center_y=floor((L2-cat_size_2)/2);
x_shift=xx-im1_center_x;y_shift=yy-im1_center_y;
im1_rectify=im_1;
im2_rectify=imrotate(im_2,rotate_angle);
im2_rectify=circshift(im2_rectify,[x_shift y_shift]);
[L1_R1,L2_R1]=size(im1_rectify);
[L1_R2,L2_R2]=size(im2_rectify);
if L1_R2<L1_R1
    cat_size1=floor((L1_R1-L1_R2)/2);cat_size2=(L1_R1-L1_R2)-cat_size1;
    im2_rectify=[zeros(cat_size1,L2_R2);im2_rectify;zeros(cat_size2,L2_R2)];
elseif L1_R2>L1_R1
    cat_size1=floor((L1_R2-L1_R1)/2);cat_size2=(L1_R2-L1_R1)-cat_size1;
    im2_rectify(1:cat_size1,:)=[];
    im2_rectify(end-cat_size2+1:end,:)=[];
end
[L1_R2,L2_R2]=size(im2_rectify);
if L2_R2<L2_R1
    cat_size1=floor((L2_R1-L2_R2)/2);cat_size2=(L2_R1-L2_R2)-cat_size1;
    im2_rectify=[zeros(L1_R2,cat_size1),im2_rectify,zeros(L1_R2,cat_size2)];
elseif L2_R2>L2_R1
    cat_size1=floor((L2_R2-L2_R1)/2);cat_size2=(L2_R2-L2_R1)-cat_size1;
    im2_rectify(:,1:cat_size1)=[];
    im2_rectify(:,end-cat_size2+1:end)=[];
end
figure;subplot(1,2,1);imshow(im1_rectify,[]);subplot(1,2,2);imshow(im2_rectify,[]);
saveas(gcf,[im2_path,'allstack\','rectify','.fig']);
close(gcf)
figure;imshowpair(im1_rectify,im2_rectify,'falsecolor');
saveas(gcf,[im2_path,'allstack\','rectify_pair','.fig']);
close(gcf)
%% z rectify
imgCell={};
for ii=1:2%two embryo
%     eval(['im_path=im',num2str(ii),'_path;'])
    if ii==1
        im_path=im1_path;
    else
        im_path=im2_path;
    end
    imgDir  = dir([im_path, 'time0001stack*.tif']);
    if size(imgDir,1)==0
        imgDir  = dir([im_path 'stack*.tif']);
    end
    img_all=[];
    sort_nat_name = sort_nat({imgDir.name});
    parfor i=1:length(imgDir)
        img = imread([im_path sort_nat_name{i}]);
        img = img(:,:,channel_num);
        if ii==1
            im_rectify=img(((1+2*abs(x_shift)):end-2*abs(x_shift)),...
                ((1+2*abs(y_shift)):end-2*abs(y_shift)));
        else
            im2_rectify=imrotate(img,rotate_angle);
            im2_rectify=circshift(im2_rectify,[x_shift y_shift]);
            [L1_R1,L2_R1]=size(im1_rectify);
            [L1_R2,L2_R2]=size(im2_rectify);
            if L1_R2<L1_R1
                cat_size1=floor((L1_R1-L1_R2)/2);cat_size2=(L1_R1-L1_R2)-cat_size1;
                im2_rectify=[zeros(cat_size1,L2_R2);im2_rectify;zeros(cat_size2,L2_R2)];
            elseif L1_R2>L1_R1
                cat_size1=floor((L1_R2-L1_R1)/2);cat_size2=(L1_R2-L1_R1)-cat_size1;
                im2_rectify(1:cat_size1,:)=[];
                im2_rectify(end-cat_size2+1:end,:)=[];
            end
            [L1_R2,L2_R2]=size(im2_rectify);
            if L2_R2<L2_R1
                cat_size1=floor((L2_R1-L2_R2)/2);cat_size2=(L2_R1-L2_R2)-cat_size1;
                im2_rectify=[zeros(L1_R2,cat_size1),im2_rectify,zeros(L1_R2,cat_size2)];
            elseif L2_R2>L2_R1
                cat_size1=floor((L2_R2-L2_R1)/2);cat_size2=(L2_R2-L2_R1)-cat_size1;
                im2_rectify(:,1:cat_size1)=[];
                im2_rectify(:,end-cat_size2+1:end)=[];
            end
            im_rectify=im2_rectify(((1+2*abs(x_shift)):end-2*abs(x_shift)),...
                ((1+2*abs(y_shift)):end-2*abs(y_shift)));
        end
        img_all(:,:,i)=im2double(im_rectify);
    end
    imgCell{ii}=img_all;
end

LayersNum=size(imgCell{1},3);
Zssim=zeros(LayersNum,LayersNum);
parfor z_i1=1:LayersNum
    for z_i2=1:LayersNum
        Zssim(z_i1,z_i2)=ssim(imgCell{1}(:,:,z_i1),imgCell{2}(:,:,z_i2));
    end
end
[~,SsimMaxIndex]=max(Zssim,[],2);
z_shift_vector=-SsimMaxIndex+(1:LayersNum)';
[~,SsimSortIndex]=sort(Zssim,2,'descend');
SsimSortIndexSelect=SsimSortIndex(:,1:3);
ZshiftMatrix=-SsimSortIndexSelect+(1:LayersNum)';
Nshift=histcounts(ZshiftMatrix(:),-LayersNum+1:LayersNum);
ShiftMax=max(Nshift);
ShiftI=find(ShiftMax==Nshift);
Zshift=ShiftI-LayersNum;
[~,ZshiftL]=min(abs(Zshift));
Zshift=Zshift(ZshiftL);
figure
hold on
bar(-LayersNum+1:LayersNum-1,Nshift,'FaceColor',[0.8500 0.3250 0.0980])
plot(Zshift,ShiftMax+0.5,'.','MarkerSize',15,'Color','b')
xlabel('Zshift')
ylabel('Frequency')
title(['Zshift = ',num2str(Zshift)])
saveas(gcf,[im2_path,'allstack\','ZshiftSelect','.fig']);

Index=0;
ZssimShift=nan(LayersNum,2*LayersNum-1);
for shift_i=-LayersNum+1:LayersNum-1
    Index=Index+1;
    ZssimShift(1+abs(shift_i/2)-shift_i/2:LayersNum-abs(shift_i/2)-shift_i/2,Index)=diag(Zssim,shift_i);
end
figure;
BarSsim=bar3(-LayersNum+1:LayersNum-1,log(ZssimShift/(min(ZssimShift,[],'all')))',0.6);
for k = 1:length(BarSsim)
    zdata = BarSsim(k).ZData;
    BarSsim(k).CData = zdata;
    BarSsim(k).FaceColor = 'interp';
end
xlabel('Layers')
ylabel('Z shift')
zlabel('Similarity')
grid off
saveas(gcf,[im2_path,'allstack\','ZshiftBar','.fig']);
close(gcf)


ShiftDiff=nan(2,LayersNum);
ShiftDiff(1,:)=diag(Zssim);
ShiftDiff(2,(1+abs(Zshift/2)-Zshift/2:end-abs(Zshift/2)-Zshift/2))=diag(Zssim,Zshift);
[~,ShowI]=max(-ShiftDiff(1,:)+ShiftDiff(2,:),[],'omitnan');
figure;
subplot(1,2,1)
imshowpair(imgCell{1}(:,:,ShowI),imgCell{2}(:,:,ShowI),'falsecolor');
title('Before Zshift')
subplot(1,2,2)
imshowpair(imgCell{1}(:,:,ShowI),imgCell{2}(:,:,ShowI+Zshift),'falsecolor');
title('After Zshift')
saveas(gcf,[im2_path,'allstack\','ZshiftCompare','.fig']);
close(gcf)

Rectify_S.xy=[x_shift,y_shift];
Rectify_S.z=Zshift;
Rectify_S.rotate=rotate_angle;
Rectify_S.size=[L1_R1,L2_R1];
Rectify_S.rowCat=CatLog(1);
Rectify_S.colCat=CatLog(2);
writestruct(Rectify_S,[im2_path,'allstack\','Rectify.xml']);
end