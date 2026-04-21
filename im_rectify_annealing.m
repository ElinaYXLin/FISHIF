function [x_shift,y_shift,rotate_angle]=im_rectify_annealing(im1_path,im2_path)
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
%% INI
scan_reg=0.4;
Resize=1;
%% Rough scanning
im_1=imresize(stack1,Resize);
im_2=imresize(stack2_ro,Resize);
[L1,L2]=size(im_1);
[l1,l2]=size(im_2);
cat_size=floor(min([l1,l2])*scan_reg);
cat_size_2=floor(cat_size/2);
im_scan=im_2((1+floor((l1-cat_size)/2)):(1+floor((l1-cat_size)/2)+cat_size),...
    (1+floor((l2-cat_size)/2)):(1+floor((l2-cat_size)/2)+cat_size));

original=stack1;
distorted=stack2_ro;
Size_o=size(original);
Size_d=size(distorted);
ptsOriginal  = detectSURFFeatures(original,'NumScaleLevels',3);
ptsDistorted = detectSURFFeatures(distorted,'NumScaleLevels',3);
Feature_Location_original=ptsOriginal.Location;
Feature_Location_distorted=ptsDistorted.Location;
Cat=@(x,y) x(:,1)>0.1*y(2)&x(:,1)<0.9*y(2)&x(:,2)>0.1*y(1)&x(:,2)<0.9*y(1);
Feature_Location_original=...
    Feature_Location_original(Cat(Feature_Location_original,Size_o),:);
Feature_Location_distorted=...
    Feature_Location_distorted(Cat(Feature_Location_distorted,Size_d),:);
figure;
scatter(Feature_Location_original(:,1),Feature_Location_original(:,2))
hold on
scatter(Feature_Location_distorted(:,1),Feature_Location_distorted(:,2))

Feature_Location_original=Size_o([2 1])-Feature_Location_original;
Feature_Location_distorted=Size_d([2 1])-Feature_Location_distorted;

Feature_Location_original3=[Feature_Location_original,zeros(size(Feature_Location_original,1),1)];
Feature_Location_distorted3=[Feature_Location_distorted,zeros(size(Feature_Location_distorted,1),1)];
[R_m,T_m]=SignalAlignment3D(Feature_Location_original3,Feature_Location_distorted3);
xx=round((L1-cat_size_2)/2-T_m(1));yy=round((L2-cat_size_2)/2+T_m(2));
rotate_angle=-acos(R_m(1,1))/pi*180;
% parpool(100);
% x0=[0 (L1-cat_size_2)/2 (L2-cat_size_2)/2];
x0=[rotate_angle xx yy];
lb=[-10+rotate_angle xx-50 yy-50];
ub=[10+rotate_angle xx+50 yy+50];
options = saoptimset('StallIterLim',100,... % 最高温度
'TolFun',1e-6,... % 最低温度
'InitialTemperature',10);

[Drift,fval,exitflag,output]= simulannealbnd(@(x) ...
        SimilarityJudgment_ssim(x,im_scan,im_1,cat_size_2),x0,lb,ub,options);
disp(['SSIM: ',num2str(1-fval)])
disp([['Rotate: ',num2str(Drift(1))] newline ['Pan: [',num2str(Drift(2:3)),']']])
xx=Drift(2);yy=Drift(3);rotate_angle=Drift(1);

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
figure;imshowpair(im1_rectify,im2_rectify,'falsecolor');
saveas(gcf,[im2_path,'allstack\','rectify_pair','.fig']);
close(gcf)

Rectify_S.xy=[x_shift,y_shift];
Rectify_S.rotate=rotate_angle;
Rectify_S.size=[L1_R1,L2_R1];
writestruct(Rectify_S,[im2_path,'allstack\','Rectify.xml']);
end