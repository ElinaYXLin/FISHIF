function [x_shift,y_shift,rotate_angle]=im_rectify_feature(im1_path,im2_path)
% uiopen('Y:\Imaging\S1\12232019_S1_60X\stacks\12232019_S1_60X_005_new\stack14.tif',1)
% uiopen('Y:\Imaging\S1\12232019_S1_60X\stacks\12232019_S1_60X_005_new\stack15.tif',1)
%  stack1=test1(:,:,4);
%  stack2=test2(:,:,4);
% stack2_ro=imrotate(stack2,63.4);
% imshow(stack2_ro,[])
% figure;
% imshow(stack2_ro,[])
% figure;imshow(stack1,[])
 h=openfig('Y:\TestForEyedropper\yxt\b2\stacks\hcr in chip  hb 5-7 b2 org 04_new\allstack\rectify_scan.fig','invisible');
 stack1=h.Children(2).Children.CData;
 stack2_ro=h.Children(1).Children.CData;
 stack2=im2uint16(stack2_ro);
 close(h)
original=stack1;
distorted=stack2_ro;
ptsOriginal  = detectSURFFeatures(original);
ptsDistorted = detectSURFFeatures(distorted);
[featuresOriginal,  validPtsOriginal]  = extractFeatures(original,  ptsOriginal);
[featuresDistorted, validPtsDistorted] = extractFeatures(distorted, ptsDistorted);
indexPairs = matchFeatures(featuresOriginal, featuresDistorted);
matchedOriginal  = validPtsOriginal(indexPairs(:,1));
matchedDistorted = validPtsDistorted(indexPairs(:,2));
figure;
showMatchedFeatures(original,distorted,matchedOriginal,matchedDistorted);
title('Putatively matched points (including outliers)');
[tform, inlierIdx] = estimateGeometricTransform2D(...
    matchedDistorted, matchedOriginal, 'similarity');
inlierDistorted = matchedDistorted(inlierIdx, :);
inlierOriginal  = matchedOriginal(inlierIdx, :);
figure;
showMatchedFeatures(original,distorted,inlierOriginal,inlierDistorted);
title('Matching points (inliers only)');
legend('ptsOriginal','ptsDistorted');
Tinv  = tform.invert.T;
ss = Tinv(2,1);
sc = Tinv(1,1);
scaleRecovered = sqrt(ss*ss + sc*sc)
thetaRecovered = atan2(ss,sc)*180/pi
outputView = imref2d(size(original));
recovered  = imwarp(distorted,tform,'OutputView',outputView);
figure, imshowpair(original,recovered,'montage')









%% INI
scan_reg=0.4;
bin_scan=2;
rotate_vector=[-180:2:180];
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
xx0=xx_vector(k1);yy0=yy_vector(k2);rotate_angle0=rotate_vector(k3);
%% Micromesh scanning
im_1=stack1;
im_2=stack2_ro;
[L1,L2]=size(im_1);
[l1,l2]=size(im_2);
scan_reg=0.4;cat_size=floor(min([l1,l2])*scan_reg);
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
ssim_mat_micr=[];
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
    ssim_mat_micr(:,:,angle_index)=ssim_xy;
end

[k1,k2,k3]=ind2sub(size(ssim_mat_micr),find(ssim_mat_micr==max(ssim_mat_micr(:))));
xx=xx_vector(k1);yy=yy_vector(k2);rotate_angle=rotate_vector(k3);
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
    im2_rectify(end-cat_size2:end,:)=[];
end
[L1_R2,L2_R2]=size(im2_rectify);
if L2_R2<L2_R1
    cat_size1=floor((L2_R1-L2_R2)/2);cat_size2=(L2_R1-L2_R2)-cat_size1;
    im2_rectify=[zeros(L1_R2,cat_size1),im2_rectify,zeros(L1_R2,cat_size2)];
elseif L2_R2>L2_R1
    cat_size1=floor((L2_R2-L2_R1)/2);cat_size2=(L2_R2-L2_R1)-cat_size1;
    im2_rectify(:,1:cat_size1)=[];
    im2_rectify(:,end-cat_size2:end)=[];
end
figure;subplot(1,2,1);imshow(im1_rectify,[]);subplot(1,2,2);imshow(im2_rectify,[]);
saveas(gcf,[im2_path,'allstack\','rectify','.fig']);
close(gcf)
Rectify_S.xy=[x_shift,y_shift];
Rectify_S.rotate=rotate_angle;
writestruct(Rectify_S,[im2_path,'allstack\','Rectify.xml']);
end