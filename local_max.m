function ymax0 = local_max(xx0,yy0,Il0,Ir0)

ymax0 = zeros(size(Il0));
for ii = 1:length(Il0)
    Itrue0 = xx0 >= Il0(ii) & xx0 <= Ir0(ii);
    ymax0(ii) = max(yy0(Itrue0));
end