figure
pool_name = parpool;
parfor ii = 1:9
    subplot(3,3,ii)
    plot(rand(1,100),rand(1,100),'.')
end
delete(pool_name)
