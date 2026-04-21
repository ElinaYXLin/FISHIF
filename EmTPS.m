function EmTPS(im1_path0,im2_path0,maskPath1,maskPath2,channel_num1,channel_num2)
im1_path0 = 'D:\yxt\data for test\org02\stacks\hcr_in_chip_set1_org02_new\MatchTo_Zshift hcr_in_chip_tll561_gt647_org02\';
im2_path0 = 'D:\yxt\data for test\tll561_gt647_org02\stacks\hcr_in_chip_tll561_gt647_org02_new\';
maskPath1 = 'D:\yxt\data for test\org02\masks\hcr_in_chip_set1_org02_new\';
maskPath2 = 'D:\yxt\data for test\tll561_gt647_org02\masks\hcr_in_chip_tll561_gt647_org02_new\';
channel_num1 = 4;
channel_num2 = 4;
stack1=double(imread([im1_path0,'\allstack\','allstack.tif']));
stack2_ro=double(imread([im2_path0,'\allstack\','allstack.tif']));

NewF = [im2_path0,'TPS\'];
if ~exist(NewF)
    mkdir(NewF)
end
% figure;imshowpair(stack2_ro,stack1,'falsecolor');
% [x2,y2] = getpts;%1 3 5 7 -> 2 4 6 8
% Zp = [x2(1:2:end-1),y2(1:2:end-1)]
% Zs = [x2(2:2:end-1),y2(2:2:end-1)]
[Centroid_nucleus1] = GetMaskCentre(channel_num1,{maskPath1,im1_path0});
[Centroid_nucleus2] = GetMaskCentre(channel_num2,{maskPath2,im2_path0});
[IDX,D] = knnsearch(Centroid_nucleus1(:,[1,2]),Centroid_nucleus2(:,[1,2]),'k',2);
if size(Centroid_nucleus2,1) > size(Centroid_nucleus1,1)
    [B,I] = maxk(D(:,1),size(Centroid_nucleus2,1)-size(Centroid_nucleus1,1));
    IDX(I,:)=[];
    D(I,:)=[];
end

% 找到重复值
[uniqueValues, ~, idx] = unique(IDX(:,1));
% 统计每个值的出现次数
counts = histc(idx, 1:numel(uniqueValues));
% 找到重复值的索引
duplicateIndices = find(counts > 1);
% 初始化矩阵以保存结果
resultMatrix = [];
% 输出重复值及其索引并保存为矩阵
for i = duplicateIndices'
    value = uniqueValues(i);
    indices = find(IDX(:,1) == value);
    % 将值和对应的索引添加到结果矩阵中
    resultMatrix = [resultMatrix; value, indices'];
end
for ii = 1:size(resultMatrix)
    Dm = D(resultMatrix(ii,2:end),:);
    D1 = Dm(1,1)+Dm(2,2);
    D2 = Dm(1,2)+Dm(2,1);
    if D1>D2
        IDX(resultMatrix(ii,3),1) = IDX(resultMatrix(ii,3),2);
        D(resultMatrix(ii,3),1) = D(resultMatrix(ii,3),2);
    else
        IDX(resultMatrix(ii,2),1) = IDX(resultMatrix(ii,2),2);
        D(resultMatrix(ii,2),1) = D(resultMatrix(ii,2),2);
    end
end


Zs = Centroid_nucleus1(IDX(:,1),[2 1]);
Zp = Centroid_nucleus2(:,[2 1]);Zs(I,:)=[];

ChooseIndex = D(:,1)>5&D(:,1)<40;
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
close(gcf)

interp.method = 'nearest';
interp.radius = 10;
interp.power = 2;
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
close(gcf)
save([NewF,'TPS','.mat'], 'Zp', 'Zs')


