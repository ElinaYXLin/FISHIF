function y = conv20(x1,x2)
y = zeros(size(x1));
y0 = conv2(x1,x2);
y = y0(1:size(x1,1),1:size(x1,2));
y(:,end) = y(:,end)+sum(y0(1:size(x1,1),size(x1,2)+1:end),2);
y(end,:) = y(end,:)+sum(y0(size(x1,1)+1:end,1:size(x1,2)),1);
y(end,end) = y(end,end)+sum(sum(y0(size(x1,1)+1:end,size(x1,2)+1:end)));