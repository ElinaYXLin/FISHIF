a = 1:100;
b = a.^2;
N = 0;
M = 0;

tic
% distcomp.feature( 'LocalUseMpiexec', false);
pool_name = parpool;
parfor ii = 1:length(a)
    N = N+b(ii);
end
parfor ii = 1:length(a)
    M = M+a(ii);
end
delete(pool_name)
toc