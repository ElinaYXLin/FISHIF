function [x_start_im,y_start_im,x_center_im,y_center_im] = tile_info_extract(image_folder,Nbin,max_image)

stitch_name = 'stitch_parameter.mat';

%%% Stitch parameter loading: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if exist([image_folder,stitch_name])
    load([image_folder,stitch_name],'x_start','y_start','imcorr1')
    size_tile = size(imcorr1);
    size_im = size(max_image);
else
    x_start = ones(Nbin);
    y_start = ones(Nbin);
    size_im = size(max_image);
    size_tile = size_im(1:2)./Nbin;
end        

x_start_im = ones(size(x_start));
y_start_im = ones(size(x_start));
x_center_im = zeros(size(x_start));
y_center_im = zeros(size(x_start));
for ii = 1:Nbin(1)
    for jj = 1:Nbin(2)
        if jj > 1
            x_start_im(ii,jj) = x_start_im(ii,jj-1)+size_tile(2)-x_start(ii,jj-1)+1;
        end
        if ii > 1
            y_start_im(ii,jj) = y_start_im(ii-1,jj)+size_tile(1)-y_start(ii-1,jj)+1;
        end
        x_center_im(ii,jj) = x_start_im(ii,jj)+(size_tile(2)+1)/2-x_start(ii,jj);
        y_center_im(ii,jj) = y_start_im(ii,jj)+(size_tile(1)+1)/2-y_start(ii,jj);
    end
% %             if x_start_im(ii,end)+size_tile(2)-x_start(ii,end) ~= size_im(2)
% %                 error('dimension 2 mismatch')
% %             end
end
% %         if y_start_im(end,end)+size_tile(1)-y_start(end,end) ~= size_im(1)
% %             error('dimension 1 mismatch')
% %         end
% % % clear x_start y_start imcorr1 size_im size_tile
%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
