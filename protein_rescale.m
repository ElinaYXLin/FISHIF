function [nucleus_protein_profile_ab_out,p_out] = protein_rescale(nucleus_protein_profile,nucleus_protein_profile_ab,quanti_p)

index_true = logical(nucleus_protein_profile(:,4));
p = polyfit(nucleus_protein_profile((nucleus_protein_profile(:,1) <= 0.5) & index_true,2),nucleus_protein_profile((nucleus_protein_profile(:,1) <= 0.5) & index_true,3),1);
p_post = polyfit(nucleus_protein_profile((nucleus_protein_profile(:,1) >= 0.8) & index_true,2),nucleus_protein_profile((nucleus_protein_profile(:,1) >= 0.8) & index_true,3),1);
mean_post = mean(nucleus_protein_profile((nucleus_protein_profile(:,1) >= 0.8) & index_true,2));
p_out = [p(1),mean_post];
% nucleus_protein_profile_ab(:,2) = ((nucleus_protein_profile(:,2)+p_post(2)/p_post(1))*p_post(1)-nucleus_protein_profile(:,3))/(p_post(1)-p(1))/quanti_p(1);
% nucleus_protein_profile_ab(nucleus_protein_profile_ab(:,2)<0,2) = 0;

nucleus_protein_profile_ab_out = nucleus_protein_profile_ab;