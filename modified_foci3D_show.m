function modified_foci3D_show(foci_bw2D,max_image,RNA_channel,mask_stack,N_cycle,image_folder,sub_pos)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A function to plot foci/nuclear mask overlaid on embryo images:
%%
%% foci_bw2D: 2D foci mask;
%% max_image: maximal projection image matrix;
%% RNA_channel: RNA channel index;
%% mask_stack: 3D nuclei mask;
%% N_cycle: nuclear cycle;
%% image_folder: image folder name;
%% sub_plot: subplot coordinate;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Image output: %%%======================================================
foci_bw0 = bwmorph(foci_bw2D,'thicken',7);
bw_perim_g = bwperim(foci_bw0);
bw_perim_WGA = bwperim(max(mask_stack,[],3));
new_image(:,:,1) = max_image(:,:,RNA_channel);
new_image(:,:,2) = 0;
new_image(:,:,3) = 0;
overlay = imoverlay(new_image, bw_perim_g, [1,1,1]);
overlay = imoverlay(overlay,bw_perim_WGA,[0,0,1]);

figure(4)
subaxis(sub_pos(1),sub_pos(2),sub_pos(3), 'Spacing', 0, 'PaddingRight',0.01, 'PaddingLeft',0.02, 'PaddingTop',0.045, 'PaddingBottom',0.04, 'Margin', 0.01);
try
    imshow(overlay)
catch
    imshow_lxt(overlay)
end
title([image_folder,' (cycle = ',num2str(N_cycle),', white: transcription foci recognition, blue: nucleus recognition)'],'Interpreter','none')


figure(40)
try
    imshow(overlay)
catch
    imshow_lxt(overlay)
end
title([image_folder,' (cycle = ',num2str(N_cycle),', white: transcription foci recognition, blue: nucleus recognition)'],'Interpreter','none')

%%% =======================================================================
