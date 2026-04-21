function [x_start_im,y_start_im,z_start_im,x_end_im,y_end_im,z_end_im] = tile_info_extract2(x_start,y_start,z_start,x_end,y_end,z_end)

x_start_im = ones(size(x_start));
y_start_im = ones(size(x_start));
z_start_im = ones(size(z_start));
x_end_im = zeros(size(x_end));
y_end_im = zeros(size(y_end));
z_end_im = zeros(size(z_end));

for ii = 1:size(x_start,1)
    for jj = 1:size(x_start,2)
        if jj > 1
            x_start_im(ii,jj) = x_start_im(ii,jj-1)+x_end(ii,jj-1)-x_start(ii,jj-1)+1;
        end
        x_end_im(ii,jj) = x_start_im(ii,jj)+x_end(ii,jj)-x_start(ii,jj);
        
        if ii > 1
            y_start_im(ii,jj) = y_start_im(ii-1,jj)+y_end(ii-1,jj)-y_start(ii-1,jj)+1;
        end
        y_end_im(ii,jj) = y_start_im(ii,jj)+y_end(ii,jj)-y_start(ii,jj);
        
        z_start_im(ii,jj) = 1;
        z_end_im(ii,jj) = z_end(ii,jj)-z_start(ii,jj)+1;
    end
end
%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
