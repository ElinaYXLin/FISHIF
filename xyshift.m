function dxy = xyshift(XY_raw,protein_RNA_xymismatch,max_image)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A function to calculate the xy compensation for chromatic abberation:
%%
%% dxy: xy compensation
%% XY_raw: original XY coordinates
%% protein_RNA_xymismatch: xy mismatch between protein and RNA channels (M_protein_RNA,C_protein_RNA,x0_protein_RNA,Nbin,Mdim,resolution,resolution_mismatch)
%% max_image: max projection of sample image
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% chromatic aberration compensation: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M_protein_RNA = protein_RNA_xymismatch{1};
C_protein_RNA = protein_RNA_xymismatch{2};
x0_protein_RNA = protein_RNA_xymismatch{3};
Nbin = protein_RNA_xymismatch{4};
Mdim = protein_RNA_xymismatch{5};
resolution = protein_RNA_xymismatch{6};
resolution_mismatch = protein_RNA_xymismatch{7};
dim_div = size(max_image,Mdim)/Nbin;
window_size = [size(max_image,1),size(max_image,2)];
window_size(Mdim) = dim_div;

XY_int = zeros(size(XY_raw));
XY_int(:,Mdim) = floor(XY_raw(:,Mdim)/dim_div)*dim_div;
XY_mod = XY_raw;
XY_mod(:,Mdim) = mod(XY_raw(:,Mdim),dim_div);
XY_center = repmat((1+window_size)/2,size(XY_mod,1),1);
XY_mod_new = (M_protein_RNA*(XY_mod-XY_center)'-repmat(C_protein_RNA,1,size(XY_mod,1))*resolution_mismatch/resolution)'+XY_center;
dxy = round(XY_mod_new+XY_int-XY_raw);
%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

