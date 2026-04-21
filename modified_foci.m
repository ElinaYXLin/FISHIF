function [foci_bw00,max_image00,SS,foci_layer] = modified_foci(file_name,max_image,RNA_channel,seg_bw,mask_stack,N_cycle,image_folder)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A function to produce input for whole embryo foci analysis using gaussian fitting data:
%%
%% foci_bw00: modified foci_bw;
%% max_image00: modified max_image00;
%% SS: equivalent foci area;
%% foci_layer: z slice of foci;
%% file_name: input Excel file name;
%% max_image: maximal projection image matrix;
%% RNA_channel: RNA channel index;
%% seg_bw: 2D nuclei mask;
%% N_cycle: nuclear cycle;
%% image_folder: image folder name;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r_ratio = 1.2;
max_image00 = zeros(size(max_image));
SS = max_image00(:,:,1);
foci_bw00 = false(size(SS));
foci_layer = zeros(size(seg_bw));
temp0 = false(size(SS));
temp_single = zeros(0);

[foci_list,~,~] = xlsread(file_name);
bwn = bwlabel(seg_bw);
N_nuclei = zeros(1,max(max(bwn)));
%% Foci candidate loading: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii = 1:size(foci_list,1)
    if mask_stack(round(foci_list(ii,6)),round(foci_list(ii,7)),foci_list(ii,8)) && (foci_list(ii,1)*foci_list(ii,2)*foci_list(ii,3)*2*pi > max_image00(round(foci_list(ii,6)),round(foci_list(ii,7)),RNA_channel)*SS(round(foci_list(ii,6)),round(foci_list(ii,7)))) && N_nuclei(bwn(round(foci_list(ii,6)),round(foci_list(ii,7)))) < 4
        foci_bw00(round(foci_list(ii,6)),round(foci_list(ii,7))) = true;
        max_image00(round(foci_list(ii,6)),round(foci_list(ii,7)),RNA_channel) = foci_list(ii,1);
        SS(round(foci_list(ii,6)),round(foci_list(ii,7))) = foci_list(ii,2)*foci_list(ii,3)*2*pi;
        foci_layer(round(foci_list(ii,6)),round(foci_list(ii,7))) = foci_list(ii,8);
        
        add_mask = getnhood(strel('disk',round(r_ratio*sqrt(foci_list(ii,2)*foci_list(ii,3)))));
        centerx = ceil(size(add_mask,1)/2);
        centery = ceil(size(add_mask,2)/2);
        startx = -min((round(foci_list(ii,6))-1),(centerx-1));
        endx = min((size(foci_bw00,1)-round(foci_list(ii,6))),(centerx-1));
        starty = -min((round(foci_list(ii,7))-1),(centery-1));
        endy = min((size(foci_bw00,2)-round(foci_list(ii,7))),(centery-1));
        temp0((round(foci_list(ii,6))+startx):(round(foci_list(ii,6))+endx),(round(foci_list(ii,7))+starty):(round(foci_list(ii,7))+endy)) = temp0((round(foci_list(ii,6))+startx):(round(foci_list(ii,6))+endx),(round(foci_list(ii,7))+starty):(round(foci_list(ii,7))+endy)) | add_mask((centerx+startx):(centerx+endx),(centery+starty):(centery+endy));
        N_nuclei(bwn(round(foci_list(ii,6)),round(foci_list(ii,7)))) = N_nuclei(bwn(round(foci_list(ii,6)),round(foci_list(ii,7))))+1;
    elseif ~mask_stack(round(foci_list(ii,6)),round(foci_list(ii,7)),foci_list(ii,8))
        temp_single = cat(1,temp_single,foci_list(ii,1)*foci_list(ii,2)*foci_list(ii,3)*2*pi);
    end
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Filtering/refinement: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% foci combination: %%% =================================================
fxy_prop = regionprops(bwlabel(temp0).*foci_bw00,max_image00(:,:,RNA_channel).*SS, 'MeanIntensity','Area','WeightedCentroid');
fz_prop = regionprops(bwlabel(temp0).*foci_bw00,foci_layer.*max_image00(:,:,RNA_channel).*SS, 'MeanIntensity');
fs_prop = regionprops(bwlabel(temp0).*foci_bw00,SS, 'MeanIntensity');

foci_inten0 = ([fxy_prop.MeanIntensity]./[fs_prop.MeanIntensity])';
foci_SS = ([fxy_prop.Area].*[fs_prop.MeanIntensity])';
foci_xy = round(cell2mat({fxy_prop.WeightedCentroid}'));
foci_z = round([fz_prop.MeanIntensity]./[fxy_prop.MeanIntensity])';
%%% =======================================================================

%%% foci reloading: %%% ===================================================
max_image00 = zeros(size(max_image));
SS = max_image00(:,:,1);
foci_bw00 = false(size(SS));
foci_layer = zeros(size(seg_bw));

for ii = 1:size(foci_xy,1)
    foci_bw00(foci_xy(ii,2),foci_xy(ii,1)) = true;
    max_image00(foci_xy(ii,2),foci_xy(ii,1),RNA_channel) = foci_inten0(ii);
    SS(foci_xy(ii,2),foci_xy(ii,1)) = foci_SS(ii);
    foci_layer(foci_xy(ii,2),foci_xy(ii,1)) = foci_z(ii);
end
%%% =======================================================================

%%% foci refinement: %%% ==================================================
% %amp_prop = regionprops(foci_bw00,max_image00(:,:,RNA_channel).*SS, 'MaxIntensity');
% if isempty(temp_single)
%     temp_single = 0;
% end
% temp_single = sort(temp_single);
% foci_bw00 = max_image00(:,:,RNA_channel).*SS >  median(temp_single);
% max_image00(:,:,RNA_channel) = max_image00(:,:,RNA_channel).*foci_bw00;
% SS = SS.*foci_bw00;
% foci_layer = foci_layer.*foci_bw00;
%%% =======================================================================

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Image output: %%%======================================================
foci_bw0 = bwmorph(foci_bw00,'thicken',7);
bw_perim_g = bwperim(foci_bw0);
bw_perim_WGA = bwperim(max(mask_stack,[],3));
new_image(:,:,1) = max_image(:,:,RNA_channel);
new_image(:,:,2) = 0;
new_image(:,:,3) = 0;
overlay = imoverlay(new_image, bw_perim_g, [1,1,1]);
overlay = imoverlay(overlay,bw_perim_WGA,[0,0,1]);

figure(4)
imshow(overlay)
title([image_folder,' (cycle = ',num2str(N_cycle),', white: transcription foci recognition, blue: nucleus recognition)'],'Interpreter','none')
%%% =======================================================================
