function mean_ratio = protein_profile3D_short2(mask_stack)

mask2D = max(mask_stack > 0,[],3);
em_size = nnz(bwconvhull(mask2D));
nu_size = nnz(mask2D);

mean_ratio = nu_size/(em_size-nu_size);