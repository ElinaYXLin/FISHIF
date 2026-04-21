%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A function to create marker for overlapping elliptical areas %%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function bw_out = ellipse_vote(bw_in,fil_im)

r0 = 15;
r1 = 1;
Lmajor0 = 4;
Lminor0 = 2;

Lmajor = 8;
Lminor = 3;
theta0 = 0;

bw_in = imerode(bw_in,strel('disk',3));

bw_out = false(size(bw_in));
im0 = double(fil_im);
im0 = im0/mean(im0(:));
label_in = bwlabel(bw_in);

[gx0,gy0] = gradient(im0);
g0 = sqrt(gx0.^2+gy0.^2);

bw_g0 = im2bw(g0,2.5*graythresh(g0));

bw_bound = bw_g0 & bwperim(bw_in);
[by,bx] = find(bw_bound);
gx = gx0(bw_bound);
gy = gy0(bw_bound);
g = g0(bw_bound);

vx = min(max(bx+round(gx./g*r0),0),size(bw_in,2));
vy = min(max(by+round(gy./g*r0),0),size(bw_in,1));

indv = sub2ind(size(bw_in),vy,vx);
Itrue = label_in(indv) == label_in(bw_bound);
bw_out(indv(Itrue)) = true;

bw_out = imdilate(bw_out,ellipse_mask(Lmajor0,Lminor0,theta0)) & imerode(bw_in,ellipse_mask(Lmajor,Lminor,theta0));
% bw_out = imdilate(bw_out,strel('disk',1)) & imerode(bw_in,ellipse_mask(Lmajor,Lminor,theta0));
% % bw_out = imclose(bw_out,ellipse_mask(Lmajor0,Lminor0,theta0));
bw_out = imfill(bw_out,'hole');