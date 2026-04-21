function nucleus_RNA_profile_new = foci_refine(nucleus_RNA_profile,foci_RNA_profile,N_thresh)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A function to revalidate foci based on nascent mRNA number %%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nucleus_RNA_profile_new = nucleus_RNA_profile;
nu0 = zeros(size(nucleus_RNA_profile,1),2);
for ii = 1:size(foci_RNA_profile,1)
    if foci_RNA_profile(ii,3) >= N_thresh && foci_RNA_profile(ii,2) > 0
        nu0(foci_RNA_profile(ii,2),1) = nu0(foci_RNA_profile(ii,2),1)+1;
        nu0(foci_RNA_profile(ii,2),2) = nu0(foci_RNA_profile(ii,2),2)+foci_RNA_profile(ii,3);
    end
end
nucleus_RNA_profile_new(:,[3:4]) = nu0;