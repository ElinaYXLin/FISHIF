function [foci_bw3D,max_image00,SS] = modified_foci3D(file_name,max_image,RNA_channel,mask_stack,N_cycle,image_folder,protein_RNA_mismatch,Inten_thresh,varargin)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A function to produce input for whole embryo foci analysis using gaussian fitting data:
%%
%% foci_bw3D: 3D foci mask;
%% max_image00: modified max_image00;
%% SS: equivalent foci area;
%% file_name: input Excel file name;
%% max_image: maximal projection image matrix;
%% RNA_channel: RNA channel index;
%% seg_bw: 2D nuclei mask;
%% N_cycle: nuclear cycle;
%% image_folder: image folder name;
%% protein_RNA_mismatch: z layer mismatch between protein and RNA channels
%% Inten_thresh: foci intensity threshold
%% varargin: 1. Figure number; 2. resolution; 3. whether to use linear representation for max_image00 and SS
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
resolution0 = 0.083;
if length(varargin) < 2 || isempty(varargin{2})
    L_ratio = 1;
else
    L_ratio = varargin{2}/resolution0;
end
r_ratio = 1.2/L_ratio;

if length(varargin) < 3 || isempty(varargin{3})
    use_linear = false;
else
    use_linear = varargin{3};
end

max_image00 = zeros(size(max_image));
SS = max_image00(:,:,1);
foci_bw00 = false(size(SS));
foci_layer = zeros(size(SS));
temp0 = false(size(SS));
temp_single = zeros(0);

[foci_list,~,~] = xlsread(file_name);
foci_list(:,8) = foci_list(:,8)+protein_RNA_mismatch;
%bwn = bwlabel(seg_bw);
N_nuclei = zeros(1,max(mask_stack(:)));
%% Foci candidate loading: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii = 1:size(foci_list,1)
    xx0 = round(foci_list(ii,6));
    yy0 = round(foci_list(ii,7));
    zz0 = foci_list(ii,8); 
    if zz0 <= 1; zz0 = 2; end
    if zz0 >= size(mask_stack,3); zz0 = size(mask_stack,3)-1; end
    if (xx0 > 0) && (xx0 <= size(mask_stack,1)) && (yy0 > 0) && (yy0 <= size(mask_stack,2)) && (zz0 > 0) && (zz0 <= size(mask_stack,3)) && max(max(max(mask_stack(xx0-3:xx0+3,yy0-3:yy0+3,zz0-1:zz0+1)))) && (foci_list(ii,1)*foci_list(ii,2)*foci_list(ii,3)*2*pi >= Inten_thresh) && (foci_list(ii,1)*foci_list(ii,2)*foci_list(ii,3)*2*pi >= max_image00(xx0,yy0,RNA_channel)*SS(xx0,yy0)) && N_nuclei(max(max(max(mask_stack(xx0-3:xx0+3,yy0-3:yy0+3,zz0-1:zz0+1))))) < 4
        foci_bw00(xx0,yy0) = true;
        max_image00(xx0,yy0,RNA_channel) = foci_list(ii,1);
        S0 = foci_list(ii,2)*foci_list(ii,3)*2*pi;
        if foci_list(ii,1) <= 0
            S0 = 1;
        end
        SS(xx0,yy0) = S0;
        foci_layer(xx0,yy0) = zz0;
        
        add_mask = getnhood(strel('disk',round(r_ratio*sqrt(foci_list(ii,2)*foci_list(ii,3)))));
        centerx = ceil(size(add_mask,1)/2);
        centery = ceil(size(add_mask,2)/2);
        startx = -min((xx0-1),(centerx-1));
        endx = min((size(foci_bw00,1)-xx0),(centerx-1));
        starty = -min((yy0-1),(centery-1));
        endy = min((size(foci_bw00,2)-yy0),(centery-1));
        temp0((xx0+startx):(xx0+endx),(yy0+starty):(yy0+endy)) = temp0((xx0+startx):(xx0+endx),(yy0+starty):(yy0+endy)) | add_mask((centerx+startx):(centerx+endx),(centery+starty):(centery+endy));
        N_nuclei(max(max(max(mask_stack(xx0-3:xx0+3,yy0-3:yy0+3,zz0-1:zz0+1))))) = N_nuclei(max(max(max(mask_stack(xx0-3:xx0+3,yy0-3:yy0+3,zz0-1:zz0+1)))))+1;
    elseif (xx0 > 0) && (xx0 <= size(mask_stack,1)) && (yy0 > 0) && (yy0 <= size(mask_stack,2)) && (zz0 > 0) && (zz0 <= size(mask_stack,3)) && ~mask_stack(xx0,yy0,zz0)
        temp_single = cat(1,temp_single,foci_list(ii,1)*foci_list(ii,2)*foci_list(ii,3)*2*pi);
    end
end
% temp0 = temp0.*max(logical(mask_stack),[],3);
temp0 = temp0 & max(logical(mask_stack),[],3);
I0_prop = regionprops(temp0,foci_bw00, 'MaxIntensity');
I_all = find([I0_prop.MaxIntensity]);
temp0 = ismember(bwlabel(temp0), I_all);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Filtering/refinement: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% foci combination: %%% =================================================
fxy_prop = regionprops(bwlabel(temp0).*foci_bw00,max_image00(:,:,RNA_channel).*SS, 'MeanIntensity','Area','WeightedCentroid');
fxy2_prop = regionprops(bwlabel(temp0).*foci_bw00,'Centroid');
fz_prop = regionprops(bwlabel(temp0).*foci_bw00,foci_layer.*max_image00(:,:,RNA_channel).*SS, 'MeanIntensity');
fz2_prop = regionprops(bwlabel(temp0).*foci_bw00,foci_layer, 'MeanIntensity');
fs_prop = regionprops(bwlabel(temp0).*foci_bw00,SS, 'MeanIntensity');

foci_inten0 = ([fxy_prop.MeanIntensity]./[fs_prop.MeanIntensity])';
foci_SS = ([fxy_prop.Area].*[fs_prop.MeanIntensity])';
foci_xy = round(cell2mat({fxy_prop.WeightedCentroid}'));
foci_xy2 = round(cell2mat({fxy2_prop.Centroid}'));
foci_z = round([fz_prop.MeanIntensity]./[fxy_prop.MeanIntensity])';
foci_z2 = round([fz2_prop.MeanIntensity]);
Inan = [fxy_prop.MeanIntensity] == 0;
foci_xy(Inan,:) = foci_xy2(Inan,:);
foci_z(Inan) = foci_z2(Inan);
%%% =======================================================================

%%% foci reloading: %%% ===================================================
foci_bw3D = false(size(mask_stack));
ind_foci = sub2ind(size(mask_stack),foci_xy(:,2),foci_xy(:,1),foci_z(:));
foci_bw3D(ind_foci) = true;

if use_linear
    [~,ind0] = sort(ind_foci);
    max_image00 = foci_inten0(ind0);
    SS = foci_SS(ind0);
else
    max_image00 = zeros(size(mask_stack));
    SS = zeros(size(mask_stack));
    max_image00(ind_foci) = foci_inten0;
    SS(ind_foci) = foci_SS;
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
foci_bw0 = bwmorph(max(foci_bw3D,[],3),'thicken',7);
bw_perim_g = bwperim(foci_bw0);
bw_perim_WGA = bwperim(max(mask_stack,[],3));
new_image(:,:,1) = max_image(:,:,RNA_channel);
new_image(:,:,2) = 0;
new_image(:,:,3) = 0;
overlay = imoverlay(new_image, bw_perim_g, [1,1,1]);
overlay = imoverlay(overlay,bw_perim_WGA,[0,0,1]);

if isempty(varargin) || isempty(varargin{1})
    figure(4)
else
    figure(varargin{1})
end

try
    imshow(overlay)
catch
    imshow_lxt(overlay)
end
title([image_folder,' (cycle = ',num2str(N_cycle),', white: transcription foci recognition, blue: nucleus recognition)'],'Interpreter','none')
%%% =======================================================================
