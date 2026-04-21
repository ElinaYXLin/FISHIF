function [nucleus_profile,DNA_mask,varargout] = DAPI_profile3D(mask_stack,signal_stack,image_folder,N_cycle,varargin)


%% Nucleus DAPI profile calculation: %%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Parameter setting and initiation: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % standard_dorsal = 'Standard protein/Bcd_dorsal.xls';
% % standard_ventral = 'Standard protein/Bcd_ventral.xls';
% % subtp = 'background subtraction';
% % bin2 = 1;
% % Nb2 = 8;
resolution0 = 0.083;
if isempty(varargin) || isempty(varargin{1})
    L_ratio = 1;
else
    L_ratio = varargin{1}/resolution0;
end
if length(varargin) >= 2 && ~isempty(varargin{2})
    r_edge0 = varargin{2};
else
    r_edge0 = 10;
end
r_edge = round(r_edge0/L_ratio);
% % bcd_dxy = xlsread(standard_dorsal);
% % bcd_vxy = xlsread(standard_ventral);
%nucleus_profile = zeros(0);
%cytoplasmic_profile = zeros(0);
Lmin = 0;
Lcbin = 0.02;
Lmax = 1;
Lcwindow = 0.02;
cyto_L = Lmin:Lcbin:Lmax;
cyto_mean = zeros(size(cyto_L));

Lmin = 0;
Lnbin = 0.05;
Lmax = 1;
Lnwindow = 0.05;
nu_L = Lmin:Lnbin:Lmax;
nu_mean = zeros(size(nu_L));
nu_std = zeros(size(nu_L));

start_limit = 0.3;
end_limit = 1-start_limit;



sigma = 1.2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Nucleus region recognition/processing: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
temp = zeros(0);
% % % Inten = zeros(0);
% % % Inten2 = zeros(0);
%mask2 = max(mask_stack,[],3);
% % % seg_bw = max(mask_stack,[],3);

for I_layer = 1:size(mask_stack,3)
%     sg_prop = regionprops(double(imerode(logical(mask_stack(:,:,I_layer)), strel('disk',r_edge))).*double(mask_stack(:,:,I_layer)),double(signal_stack(:,:,I_layer)),'MeanIntensity','PixelValues');
    sg_prop = regionprops(double(imerode(logical(mask_stack(:,:,I_layer)), strel('disk',r_edge))).*double(mask_stack(:,:,I_layer)),double(signal_stack(:,:,I_layer)),'MeanIntensity');
    sg2_prop = regionprops(double(imerode(logical(mask_stack(:,:,I_layer)), strel('disk',r_edge))).*double(mask_stack(:,:,I_layer)),double(signal_stack(:,:,I_layer)).^2,'MeanIntensity');
    if I_layer == 1
%         temp0 = cell(size(mask_stack,3),max(mask_stack(:)));
        temp = zeros(size(mask_stack,3),max(mask_stack(:)));
        temp2 = zeros(size(mask_stack,3),max(mask_stack(:)));
    end
    if ~isempty(sg_prop)
%         temp0(I_layer,1:length(sg2_prop)) = {sg_prop.PixelValues};  %%%%% Binomial partition preparation
        temp(I_layer,1:length(sg2_prop)) = [sg_prop.MeanIntensity];
        temp2(I_layer,1:length(sg2_prop)) = [sg2_prop.MeanIntensity]-[sg_prop.MeanIntensity].^2;
    end
end
temp(isnan(temp)) = 0;
[Inten,max_index] = max(temp,[],1);
Inten2 = temp2(max_index+size(temp,1)*[0:(size(temp,2)-1)]);
nucleus_profile = [Inten',Inten2',max_index'];
varargout = {temp};
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% DNA mask generation: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % DNA_mask = zeros(size(mask_stack));
DNA_mask = zeros(0);
% % for I_nucleus = 1:max(mask_stack(:))
% %     im_temp = signal_stack(mask_stack == I_nucleus);
% %     bw_temp = im2bw(im_temp,graythresh(im_temp));
% %     DNA_mask(mask_stack == I_nucleus) = bw_temp;
% % end
% % %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Figure plot: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(16)
    plot(Inten,Inten2,'o')
    title(['DAPI profile: ',image_folder,', cycle = ',num2str(N_cycle)],'Interpreter','none')
    xlabel('Mean DAPI')
    ylabel('Std DAPI')
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
