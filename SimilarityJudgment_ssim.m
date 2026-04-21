function ssimval=SimilarityJudgment_ssim(Drift,im_scan,im_1,cat_size_2)
% Drift=[angle x y]
rotate_angle=Drift(1);
xx=round(Drift(2));
yy=round(Drift(3));
im_scan_r=imrotate(im_scan,rotate_angle);
[l1,l2]=size(im_scan_r);
im_scan_use=im_scan_r((1+floor((l1-cat_size_2)/2)):(1+floor((l1-cat_size_2)/2)+cat_size_2),...
(1+floor((l2-cat_size_2)/2)):(1+floor((l2-cat_size_2)/2)+cat_size_2));
im_match=im_1(xx:(xx+cat_size_2),yy:(yy+cat_size_2));
ssimval = 1-ssim(im_match,im_scan_use,'Exponents',[0 0.5 1]);
end