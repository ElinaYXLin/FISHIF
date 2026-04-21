function mask_stack_new = mask_dilate(mask_stack,n)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A function to dilate 3D mask stacks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% mask_stack_new (output): Dilated mask stack;
%% mask_stack (input): Original mask stack;
%% n (input): dilation radius;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mask_stack_new = zeros(size(mask_stack));

for ii = 1:size(mask_stack,3)
    mask_temp = bwmorph(logical(mask_stack(:,:,ii)),'thicken',n);
    prop_temp = regionprops(mask_temp,mask_stack(:,:,ii),'MaxIntensity');
    label_temp = [0,[prop_temp.MaxIntensity]];
    mask_stack_new(:,:,ii) = label_temp(bwlabel(mask_temp)+1);
end




