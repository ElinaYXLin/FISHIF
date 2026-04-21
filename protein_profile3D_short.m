function out_S = protein_profile3D_short(nucleus_protein_profile,nucleus_protein_profile_ab,mask_stack,resolution,lambda_C,index_true)

N0 = 100;
Nm = ceil(0.02*size(nucleus_protein_profile,1));

rand_I = rand(size(nucleus_protein_profile,1),1);
[~,IX0] = sort(rand_I);
IX0 = IX0(index_true(IX0));
all_S = zeros(N0,1);
IX = IX0(1:N0);

for ii = 1:N0
    I_layer = nucleus_protein_profile(IX(ii),4);
    all_S(ii) = nnz(mask_stack(:,:,I_layer) == IX(ii))*resolution^2;
end
mean_S = mean(all_S);
mean_C = lambda_C(1)+lambda_C(3);

[C_sort,I_sort] = sort(nucleus_protein_profile_ab(:,2),'descend');
C_sort = C_sort(index_true(I_sort));
I_sort = I_sort(index_true(I_sort));

max_layer = nucleus_protein_profile(I_sort(1),4);
max_S = nnz(mask_stack(:,:,max_layer) == I_sort(1))*resolution^2;
max_C = C_sort(1);

all_S = zeros(Nm,1);
for ii = 1:Nm
    I_layer = nucleus_protein_profile(I_sort(ii),4);
    all_S(ii) = nnz(mask_stack(:,:,I_layer) == I_sort(ii))*resolution^2;
end
max_Sm = mean(all_S);
max_Cm = mean(C_sort(1:Nm));

out_S = [mean_S,mean_C,mean_S*mean_C,max_S,max_C,max_S*max_C,max_Sm,max_Cm,max_Sm*max_Cm];