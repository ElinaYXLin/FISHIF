function [nucleus_RNA_profile0,ind_foci] = foci_info2(nucleus_RNA_profile,foci_RNA_profile)
%% A function to extract TX for each loci, including the silent ones (New version: exclude nuclei with > 2 foci)
% % % nucleus_profile = [nucleus_distance,[nucleus_prop.MeanIntensity]',Nfoci_nucleus,Ifoci_nucleus,nucleus_background'];
% % % foci_profile = [foci_distance,[N_prop.MaxIntensity]',foci_intensity',foci_intensity2',find(foci_bw3D)];

ind_foci = foci_RNA_profile(:,2);
I_foci = foci_RNA_profile(:,3);

N_nu = size(nucleus_RNA_profile,1);
[n_ind,i_ind] = hist(ind_foci,1:N_nu);
n_add = 2-n_ind';

ind_add = cat(1,find(n_add > 0),find(n_add-1 > 0));
ind_all = find(n_add >= 0);
I2 = ismember(ind_foci,ind_all);
ind_foci = cat(1,ind_foci(I2),ind_add);
I_foci = cat(1,I_foci(I2),zeros(size(ind_add)));

nucleus_RNA_profile0 = nucleus_RNA_profile(ind_foci,:);
nucleus_RNA_profile0(:,4) = I_foci;