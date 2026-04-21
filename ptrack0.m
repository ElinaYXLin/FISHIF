function max_spot = ptrack0(imstack0,mask_out)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A function to extract information of RNA spot candidates %%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% max_spot (output): RNA spot candidate information (I0,sigma_x,sigma_y,sigma_z,A,x0,y0,z0,theta,residual);
%% imstack0 (input): raw image.
%% mask_out (input): RNA spot candidate mask.
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Inten0 = imstack0(mask_out);
ind0 = find(mask_out);
[ii,jj,kk] = ind2sub(size(mask_out),ind0);

max_spot = [Inten0,ones(length(Inten0),3),zeros(length(Inten0),1),ii,jj,kk,zeros(length(Inten0),2)];
