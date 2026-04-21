function EmTPS2(im1_path0,im2_path0,maskPath1,maskPath2,channel_num1,channel_num2,costofnonassignment,dis)
%%% Alignment: warp images(path2) to standard images(path1)
%%% Shihe Zhang 10-12-2024
%%% input sample
    % im1_path0 = 'D:\yxt\data for test\org02\stacks\hcr_in_chip_set1_org02_new\MatchTo_Zshift hcr_in_chip_tll561_gt647_org02\';
    % im2_path0 = 'D:\yxt\data for test\tll561_gt647_org02\stacks\hcr_in_chip_tll561_gt647_org02_new\';
    % maskPath1 = 'D:\yxt\data for test\org02\masks\hcr_in_chip_set1_org02_new\';
    % maskPath2 = 'D:\yxt\data for test\tll561_gt647_org02\masks\hcr_in_chip_tll561_gt647_org02_new\';
    % channel_num1 = 4;
    % channel_num2 = 4;
    
%%% parameter
    % costofnonassignment = 100;
    % If some connections are unreasonable, reduce this value; if connections are insufficient, increase this value.
if ~exist([im1_path0,'\allstack\','allstack.tif'])
    im_allstack(im1_path0,channel_num1)
end
if ~exist([im2_path0,'\allstack\','allstack.tif'])
    im_allstack(im2_path0,channel_num2)
end
stack1=double(imread([im1_path0,'\allstack\','allstack.tif']));
stack2_ro=double(imread([im2_path0,'\allstack\','allstack.tif']));

NewF = [im2_path0,'TPS\'];
NewFS = [NewF,'stacks\'];
if ~exist(NewF)
    mkdir(NewF)
    mkdir(NewFS)
end
% figure;imshowpair(stack2_ro,stack1,'falsecolor');
% [x2,y2] = getpts;%1 3 5 7 -> 2 4 6 8
% Zp = [x2(1:2:end-1),y2(1:2:end-1)]
% Zs = [x2(2:2:end-1),y2(2:2:end-1)]
[Centroid_nucleus1] = GetMaskCentre(channel_num1,{maskPath1,im1_path0});
[Centroid_nucleus2] = GetMaskCentre(channel_num2,{maskPath2,im2_path0});
Centroid_nucleus1 = Centroid_nucleus1(:,[1 2]);
Centroid_nucleus2 = Centroid_nucleus2(:,[1 2]);

costMatrix=[];
for i = size(Centroid_nucleus2, 1):-1:1
    delta = Centroid_nucleus1 - Centroid_nucleus2(i, :);
    costMatrix(i, :) = sqrt(sum(delta .^ 2, 2));
end

[assignments, unassignedTracks, unassignedDetections] = ...
    assignmunkres(costMatrix,costofnonassignment);

Zs = Centroid_nucleus1(assignments(:,2),[2 1]);
Zp = Centroid_nucleus2(assignments(:,1),[2 1]);

disfuc = @(x1,x2) sqrt(sum((x1-x2).^ 2, 2));
D = disfuc(Zs,Zp);
ChooseIndex = D(:,1)>dis(1)&D(:,1)<dis(2);
Zp = Zp(ChooseIndex,:);
Zs = Zs(ChooseIndex,:);
% im = rgb2gray(im);
im = stack2_ro;
shape = size(im);
H = shape(1);
W = shape(2);
outDim = [floor(W), floor(H)];
% the correspondences need at least four points
% Zp = [217, 39; 204, 95;  174, 223; 648, 402] % (x, y) in each row
% Zs = [283, 54; 166, 101; 198, 250; 666, 372]

figure;imshowpair(stack2_ro,stack1,'falsecolor'); hold on;
plot(Zp(:, 1), Zp(:, 2), 'yo')
plot(Zs(:, 1), Zs(:, 2), 'ro')
line([Zp(:, 1)'; Zs(:, 1)'], [Zp(:, 2)'; Zs(:, 2)'],'linewidth',2);
saveas(gcf,[NewF,'TPSpre','.fig']);
% close(gcf)

interp.method = 'nearest';
interp.radius = 5;
interp.power = 1;
% Xw, Yw are the transformed coordinates using warp.
[Xw, Yw, imgw, imgwr, map] = tpswarp(stack2_ro, outDim, Zp, Zs, interp); 

figure;
subplot(2,1,1)
imshowpair(stack2_ro,stack1,'falsecolor');
title('Before')
subplot(2,1,2)
imshowpair(imgw,stack1,'falsecolor');
title('After')
saveas(gcf,[NewF,'TPS','.fig']);
% close(gcf)
save([NewF,'TPS','.mat'], 'Zp', 'Zs', 'assignments', 'unassignedTracks', 'unassignedDetections')

% warp each stack
imgDirN0  = dir([im2_path0 'time0001stack*.tif']);
if size(imgDirN0,1)==0
    imgDirN0  = dir([im2_path0 'stack*.tif']);
end
StackRefNum=length(imgDirN0);
for ii=1:StackRefNum
    img = imread([im2_path0 imgDirN0(ii).name]);
    [~, ~, imgw, ~, ~] = tpswarp(img, outDim, Zp, Zs, interp); 
    tiffwrite0(uint16(imgw),[NewFS,imgDirN0(ii).name])
end

