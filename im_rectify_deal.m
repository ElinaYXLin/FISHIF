function im_rectify_deal(x_shift,y_shift,rotate_angle,im1_path,im2_path)
imgDir  = dir([im1_path 'stack*.tif']);
mkdir([im1_path,'Restack']);
mkdir([im2_path,'Restack']);
sort_nat_name = sort_nat({imgDir.name});
%img = cell2mat(struct2cell(load([imgPath sort_nat_name{1}])));
for i = 1:size(imgDir,1) 
    img1 = imread([im1_path sort_nat_name{i}]);
    img2 = imread([im2_path sort_nat_name{i}]);
    im1_rectify_allchannel=[];
    im2_rectify_allchannel=[];
    for ii=1:size(img1,3)
        im_1=img1(:,:,ii);
        im_2=img2(:,:,ii);
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
        im1_rectify_allchannel=cat(3,im1_rectify_allchannel,im1_rectify);
        im2_rectify_allchannel=cat(3,im2_rectify_allchannel,im2_rectify);
    end
    tiffwrite0(im1_rectify_allchannel,[im1_path,'Restack\',sort_nat_name{i}])
    tiffwrite0(im2_rectify_allchannel,[im2_path,'Restack\',sort_nat_name{i}])
end
end
