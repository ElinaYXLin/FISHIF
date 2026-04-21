function r2 = cal_r2(x,y,p)

yfit = polyval(p,x);
yresid = y - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(y)-1) * var(y);
r2 = 1 - SSresid/SStotal;
