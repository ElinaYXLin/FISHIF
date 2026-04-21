function [x_out,y_out] = corr_plot(ind_foci,nucleus_RNA_profile0)

I_nu = ind_foci;
I_nu0 = unique(I_nu);
[SInu,SX] = sort(I_nu);

SInu0 = [SInu+rand(size(SInu));I_nu0];
[~,SX2] = sort(SInu0);
Sn1 = find(SX2 > length(SInu))+1;
Sn2 = Sn1+1;

In1 = SX(SX2(Sn1));
In2 = SX(SX2(Sn2));

x_out = nucleus_RNA_profile0(In1,4);
y_out = nucleus_RNA_profile0(In2,4);