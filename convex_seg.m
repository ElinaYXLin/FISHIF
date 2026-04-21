function convex_image = convex_seg(image0)

convex_image = false(size(image0));
convex_prop = regionprops(logical(image0),'ConvexImage','BoundingBox');
if size(convex_prop(:),1)>0
    for I_temp = 1:size(convex_prop(:),1)
        x1 = uint16(convex_prop(I_temp).BoundingBox(2));
        x2 = uint16(x1+convex_prop(I_temp).BoundingBox(4)-1);
        y1 = uint16(convex_prop(I_temp).BoundingBox(1));
        y2 = uint16(y1+convex_prop(I_temp).BoundingBox(3)-1);
        convex_image(x1:x2,y1:y2) = convex_image(x1:x2,y1:y2) | convex_prop(I_temp).ConvexImage;
    end
end