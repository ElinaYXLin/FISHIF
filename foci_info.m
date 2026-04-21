function [nucleus_RNA_profile0,ind_foci] = foci_info(nucleus_RNA_profile,foci_RNA_profile)

% % % nucleus_profile = [nucleus_distance,[nucleus_prop.MeanIntensity]',Nfoci_nucleus,Ifoci_nucleus,nucleus_background'];
% % % foci_profile = [foci_distance,[N_prop.MaxIntensity]',foci_intensity',foci_intensity2',find(foci_bw3D)];
% foci_RNA_profile = foci_RNA_profile(foci_RNA_profile(:,2) > 0,:);
ind_foci = foci_RNA_profile(:,2);
I_foci = foci_RNA_profile(:,3);

N_nu = size(nucleus_RNA_profile,1);
[n_ind,i_ind] = hist(ind_foci,1:N_nu);
n_add = 2-n_ind';

ind_add = cat(1,find(n_add > 0),find(n_add-1 > 0));
ind_foci = cat(1,ind_foci,ind_add);
I_foci = cat(1,I_foci,zeros(size(ind_add)));

nucleus_RNA_profile0 = nucleus_RNA_profile(ind_foci,:);
nucleus_RNA_profile0(:,4) = I_foci;