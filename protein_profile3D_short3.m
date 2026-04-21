function protein_profile3D_short3(nucleus_protein_profile,mask_stack,image_folder,N_cycle)


%% z layer heat map plotting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cmap = [[0,0,0];colormap];
% cmax = max(nucleus_protein_profile(:,4));
% cmin = 0;
cmid = median(nucleus_protein_profile(:,4));
cmax = cmid+5;
cmin = max((cmid-4),1);
c0 = nucleus_protein_profile(:,4);
c0(c0 > cmax) = cmax;
c0(c0 < cmin) = cmin;
cind_profile = ceil((size(cmap,1)-2)*(c0-cmin)/(cmax-cmin))+2;
cind_profile(cind_profile < 2) = 2;
cind_profile(cind_profile > size(cmap,1)) = size(cmap,1);
cind_profile = [1;cind_profile];

nu_cind = cind_profile(mask_stack+1);
nu_cind2D = round(sum(nu_cind-1,3)./(sum(mask_stack > 0,3)+(max(mask_stack,[],3) == 0)))+1;

nu_image = zeros([size(nu_cind2D),3]);
for i_channel = 1:3
    c_temp = cmap(:,i_channel);
    nu_image(:,:,i_channel) = c_temp(nu_cind2D);
end

figure(140)
    clf
    imshow(nu_image)
    
    hc = colorbar;
    ctick = get(hc,'YTick');
    ctick_label = round((ctick-1)/(size(cmap,1)-1)*(cmax-cmin)+cmin);
    set(hc,'YTickLabel',ctick_label)
    title([image_folder,', cycle = ',num2str(N_cycle)],'Interpreter','none')
