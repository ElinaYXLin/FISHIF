clear all

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A program to simulate the local enrichment of Bcd binding %%%%%%%%%%%%%%
%% foci_bw00
%% max_image00
%% SS
%% foci_layer
%% mask_stack
%% signal_stack
%% RNA_stack
%% nucleus_protein_profile
%% image_folder
%% N_cycle
%% resolution
%% resolutionz
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Parameter setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dim_xy = [5000,5000];
N_xy = [20,20];
L_xy = floor(dim_xy(1:2)./(N_xy+1));
lim_L = [0.2,0.8];
resolution = 0.091;
resolutionz = resolution*4;
N_cycle = 12;
image_folder = 'simu/';
nM_con = 6.02e8;
lim_C = [0,100]*1e-9; %%% nM
th_C = 100e-9; %%% nM
h_C = 2;
r_nucleus = 50;

N_site = 2600;
lambda = 2.338;
lim_bind = [0,3];
lim_p = exp(-lambda*lim_bind);
act_bind = 0.5;
act_t = false;

if act_t
    th_foci = 0;
    FISH_name = 'act5C';
else
    th_foci = 1;
    FISH_name = 'hb';
end

psf = fspecial('gaussian', 10, 1.2);
unit_I = 1./resolution./resolution./2./resolutionz./nM_con;
back_I = 0.05*unit_I;
N_layer = 3;

image_folder = 'simu_enrich\Results\';
file_name = ['simu_enrichment_',FISH_name,'_k',num2str(th_C/1e-9),'_multi_psf'];
file_tail = '.mat';
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Mask generation: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nu_mask = false(dim_xy);   %%% nuclei mask
nu_mask(L_xy(1):L_xy(1):N_xy(1)*L_xy(1),L_xy(2):L_xy(2):N_xy(2)*L_xy(2)) = true;
nu_mask = imdilate(nu_mask,strel('disk',r_nucleus));
nucleus_prop = regionprops(nu_mask,'PixelIdxList','Centroid');
image_raw = zeros(size(nu_mask));   %%% background mean protein level
local_raw = zeros(size(nu_mask));   %%% binding site mean protein level
RNA_raw = zeros(size(nu_mask));   %%% FISH level
foci_raw = false(size(nu_mask));   %%% foci mask
nucleus_protein_profile = zeros(length(nucleus_prop),3);

for I_nu = 1:length(nucleus_prop)
    C_nu = lim_C(1)+(lim_C(2)-lim_C(1))*rand(1);
    n_nu =C_nu*resolution*resolution*resolutionz*nM_con;
    ind_nu = nucleus_prop(I_nu).PixelIdxList;
    image_raw(ind_nu) = n_nu;
    
    [~,IX] = sort(rand(size(ind_nu)));
    ind_site = ind_nu(IX(1:N_site));
    str_site = sort(-log(lim_p(1)+(lim_p(2)-lim_p(1))*rand(size(ind_site)))./lambda);
    str_site(end-1:end) = lim_bind(2);
    ind_act = find(str_site >= act_bind,1);
    if act_t
        str_site(ind_act-1:ind_act) = act_bind;
        foci_raw(ind_site(ind_act-1:ind_act)) = true;
        
    else
        foci_raw(ind_site(end-1:end)) = true;
    end
    local_raw(ind_site) = 2*str_site*(C_nu.^h_C./(C_nu.^h_C+th_C.^h_C));
    nucleus_protein_profile(I_nu,:) = [lim_L(1)+nucleus_prop(I_nu).Centroid(1)/dim_xy(2)*(lim_L(2)-lim_L(1)),C_nu,C_nu];
end

%%% Stochastic process:
image_raw = poissrnd(image_raw);
local_raw = poissrnd(local_raw);
RNA_raw = conv2(local_raw.*foci_raw,psf,'same');
foci_raw = foci_raw & (local_raw >= th_foci);
signal_raw = conv2((image_raw+local_raw)*unit_I+back_I*rand(size(image_raw)),psf,'same');

%%% Finalization:
foci_bw00 = foci_raw;
max_image00 = RNA_raw;
SS = double(foci_raw);
foci_layer = double(foci_raw);
mask_stack = repmat(nu_mask,[1,1,N_layer]);
signal_stack = repmat(signal_raw,[1,1,N_layer]);
RNA_stack = repmat(conv2(RNA_raw,psf,'same'),[1,1,N_layer]);

save([image_folder,file_name,file_tail],'foci_bw00','max_image00','SS','foci_layer','mask_stack','signal_stack','RNA_stack','nucleus_protein_profile','image_folder','N_cycle','resolution','resolutionz');
