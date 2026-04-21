function single_mismatch = mismatch_single(RNA_fit,RNA2_fit,RNA_name,RNA2_name,image_size,image_folder,varargin)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A function to measure z dimension mismatch of two channels: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% single_mismatch: a matrix that record mismatch of all matching foci from two channels. (label of channel 1, label of channel 2, mean x, mean y, x mismatch: x1-x2, y mismatch: y1-y2)
%% RNA_fit: single RNA spot information from channel 1
%% RNA2_fit: single RNA spot information from channel 2
%% RNA_name: RNA channel 1 name
%% RNA2_name: RNA channel 2 name
%% image_size: size of the image
%% image_folder: Image folder name
%% varargin: {dmax}: maximal xy mismatch
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Parameter setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hist_bin = [-5:5];   %%% histogram bin
if isempty(varargin)
    dmax = 4;   %%% maximal xy mismatch
else
    dmax = varargin{1};
end

Inten1 = prod(RNA_fit(:,1:3),2)*2*pi;
Inten2 = prod(RNA2_fit(:,1:3),2)*2*pi;

xyz1 = RNA_fit(:,6:8);   %%% Coordinates of RNA spots in channel 1
xyz2 = RNA2_fit(:,6:8);   %%% Coordinates of RNA spots in channel 2
for Iz = 1:max([xyz1(:,3);xyz2(:,3)])
	dxy{Iz} = pdist2(xyz1(xyz1(:,3) == Iz,1:2),xyz2(xyz2(:,3) == Iz,1:2));   %%% xy distance of RNA spots in two channels
    label1{Iz} = find(xyz1(:,3) == Iz);
    label2{Iz} = find(xyz2(:,3) == Iz);
end
single_mismatch = zeros(0);
xy_mismatch = false([image_size,3]);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Spot matching: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for I1 = 1:size(xyz1,1)
    [dmin,Imin] = min(dxy{xyz1(I1,3)}(label1{xyz1(I1,3)} == I1,:));
    
%     if nnz(label1{xyz1(I1,3)} == I1) ~= 1
%         error([num2str(nnz(label1{xyz1(I1,3)} == I1)),'         ',num2str(I1)])
%     end
%     [~,Imin1] = min(dxy{xyz1(I1,3)}(:,Imin));
%     if (dmin <= dmax) && (label1{xyz1(I1,3)}(Imin1) == I1)
%     size(dmin)
%     size(dxy{xyz1(I1,3)}(:,Imin))
    if (length(dmin) == 1) && (dmin <= dmax) && all(dmin <= dxy{xyz1(I1,3)}(:,Imin))
        I2 = label2{xyz1(I1,3)}(Imin);
        single_mismatch = cat(1,single_mismatch,[I1,I2,(xyz1(I1,1)+xyz2(I2,1))/2,(xyz1(I1,2)+xyz2(I2,2))/2,xyz1(I1,1)-xyz2(I2,1),xyz1(I1,2)-xyz2(I2,2),xyz1(I1,1:3),xyz2(I2,1:3)]);
        xy_mismatch(ceil(xyz1(I1,1)),ceil(xyz1(I1,2)),1) = true;
        xy_mismatch(ceil(xyz2(I2,1)),ceil(xyz2(I2,2)),2) = true;
    end
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Plot z mismatch histogram: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(93)
    plot(Inten1(single_mismatch(:,1)),Inten2(single_mismatch(:,2)),'.')
    xlabel(['Intensity on ',RNA_name])
    ylabel(['Intensity on ',RNA2_name])
    title(['FISH labeling on two channels: ',image_folder,' N1 = ',num2str(size(RNA_fit,1)),' N2 = ',num2str(size(RNA2_fit,1)),' N12 = ',num2str(size(single_mismatch,1))],'Interpreter','none')

figure(94)
    xy_mismatch(:,:,1) = imdilate(xy_mismatch(:,:,1),strel('disk',3));
    xy_mismatch(:,:,2) = imdilate(xy_mismatch(:,:,2),strel('disk',3));
    imshow(double(xy_mismatch))
    title(['Chromatic mismatch on xy dimension (single mRNA): ',image_folder,'(Red: ',RNA_name,', Green: ',RNA2_name,')'],'Interpreter','none')
    
    
    
    