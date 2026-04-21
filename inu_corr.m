function [corr_out,corr_test,corrb_out,corrb_test] = inu_corr(ind_foci,nucleus_RNA_profile0)

I_nu = ind_foci;
I_nu0 = unique(I_nu);
[SInu,SX] = sort(I_nu);

SInu0 = [SInu+rand(size(SInu));I_nu0];
[~,SX2] = sort(SInu0);
Sn1 = find(SX2 > length(SInu))+1;
Sn2 = Sn1+1;

In1 = SX(SX2(Sn1));
In2 = SX(SX2(Sn2));

mix0 = rand(size(In1));
[~,MX] = sort(mix0);
In1 = In1(MX);
In2 = In2(MX);

corr_out = ICC([nucleus_RNA_profile0(In1,4),nucleus_RNA_profile0(In2,4)],'1-1');
corr_test = ICC([nucleus_RNA_profile0(In1,4),circshift(nucleus_RNA_profile0(In2,4),1)],'1-1');

corrb_out = ICC([nucleus_RNA_profile0(In1,4)>0,nucleus_RNA_profile0(In2,4)>0],'1-1');
corrb_test = ICC([nucleus_RNA_profile0(In1,4)>0,circshift(nucleus_RNA_profile0(In2,4),1)>0],'1-1');