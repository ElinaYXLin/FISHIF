function overlay = DAPI_show(imstack,immask,vlayer,vcontrast,vcontrmin,vmask)

overlay(:,:,2) = double(bwperim(immask(:,:,vlayer)))*vmask;
overlay(:,:,3) = (imstack(:,:,vlayer)-vcontrmin)/(vcontrast-vcontrmin);
overlay(:,:,1) = 0;