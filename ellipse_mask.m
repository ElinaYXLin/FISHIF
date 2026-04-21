function im_out = ellipse_mask(Lmajor,Lminor,theta)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A function to create an elliptical mask with given size and angle
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% im_out (output): output image;
%% Lmajor (input): half width of the major axis (unit: pixel);
%% Lminor (input): half width of the minor axis (unit: pixel);
%% theta (input): rotational angle of the ellipse (unit: degree);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Initialization: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nmax = max(Lmajor,Lminor);
im_out = false(2*Nmax+1);

XC = Nmax+1;
YC = Nmax+1;
[X00,Y00] = meshgrid(-Nmax:Nmax);
ind00 = sub2ind(size(im_out),Y00(:)+YC,X00(:)+XC);

theta1 = pi*theta/180;
R = [cos(theta1),sin(theta1);-sin(theta1),cos(theta1)];
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Construction of the mask: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
XY11 = [X00(:),Y00(:)]*R;
Itrue = (XY11(:,1)/Lmajor).^2+(XY11(:,2)/Lminor).^2 <= 1;

im_out(ind00) = Itrue;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








