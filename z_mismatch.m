function [foci_mismatch,varargout] = z_mismatch(RNA_fit,RNA2_fit,RNA_name,RNA2_name,image_size,image_folder,varargin)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A function to measure z dimension mismatch of two channels: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% foci_mismatch: a matrix that record mismatch of all matching foci from two channels. (label of channel 1, label of channel 2, mean x, mean y, x mismatch: x1-x2, y mismatch: y1-y2, z mismatch: z1-z2)
%% RNA_fit: RNA foci spot information from channel 1
%% RNA2_fit: RNA foci spot information from channel 2
%% RNA_name: RNA channel 1 name
%% RNA2_name: RNA channel 2 name
%% image_size: size of the image
%% image_folder: Image folder name
%% varargin: {dmax}: maximal xy mismatch
%% varargout: {dxy_min}: xy mismatch
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Parameter setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hist_bin = [-5:5];   %%% histogram bin for z distance
hist_bin2 = 0:0.5:10;   %%% histogram bin for xy distance
if isempty(varargin)
    dmax = 2;   %%% maximal xy mismatch
else
    dmax = varargin{1};
end
xyz1 = RNA_fit(:,6:8);   %%% Coordinates of RNA spots in channel 1
xyz2 = RNA2_fit(:,6:8);   %%% Coordinates of RNA spots in channel 2
dxy = pdist2(xyz1(:,1:2),xyz2(:,1:2));   %%% xy distance of RNA spots in two channels
dxy0 = min(dxy,[],2);
foci_mismatch = zeros(0);
xy_mismatch = false([image_size,3]);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Spot matching: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for I1 = 1:size(xyz1,1)
    [dmin,Imin] = min(dxy(I1,:));
    if dmin <= dmax
        foci_mismatch = cat(1,foci_mismatch,[I1,Imin,(xyz1(I1,1)+xyz2(Imin,1))/2,(xyz1(I1,2)+xyz2(Imin,2))/2,xyz1(I1,1)-xyz2(Imin,1),xyz1(I1,2)-xyz2(Imin,2),xyz1(I1,3)-xyz2(Imin,3),xyz1(I1,1:3),xyz2(Imin,1:3)]);
        xy_mismatch(ceil(xyz1(I1,1)),ceil(xyz1(I1,2)),1) = true;
        xy_mismatch(ceil(xyz2(Imin,1)),ceil(xyz2(Imin,2)),2) = true;
    end
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
varargout = {dxy0};

%% Plot z mismatch histogram: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(91)
subplot(1,2,1)
    hist(foci_mismatch(:,7),hist_bin)
    xlabel(['z(',RNA_name,') - z(',RNA2_name,')'])
    ylabel('#')
    title(['Chromatic mismatch on z dimension: ',image_folder],'Interpreter','none')
subplot(1,2,2)
    hist(dxy0,hist_bin2)
    xlabel(['xy(',RNA_name,') - xy(',RNA2_name,')'])
    ylabel('#')
    title(['Chromatic mismatch on xy dimensions: ',image_folder],'Interpreter','none')

figure(92)
    xy_mismatch(:,:,1) = imdilate(xy_mismatch(:,:,1),strel('disk',3));
    xy_mismatch(:,:,2) = imdilate(xy_mismatch(:,:,2),strel('disk',3));
    imshow(double(xy_mismatch))
    title(['Chromatic mismatch on xy dimension: ',image_folder,'(Red: ',RNA_name,', Green: ',RNA2_name,')'],'Interpreter','none')