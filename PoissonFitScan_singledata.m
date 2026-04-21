function [para,poiss2_dtb_best2] = PoissonFitScan_singledata(ndata2,state02)
%% A program to fit data by poisson distribution

%% RNA2
np1_group=[1:0.1:80];
x_group=1-state02;
P_double=@(k,x,np1)x*poisspdf(k,np1);
ndata2 = ndata2(ndata2 <= 100&ndata2 > 0);
lld_group_poiss2_best2=inf;
para=[0.9 40];
poiss2_dtb=P_double(0:100,para(1),para(2));
poiss2_dtb(1)=0;
poiss2_dtb=poiss2_dtb./sum(poiss2_dtb(:));
poiss2_dtb_best2=poiss2_dtb;
for ii = 1:length(x_group)
    for jj = 1:length(np1_group)
        poiss2_dtb=P_double(0:100,x_group(ii),np1_group(jj));
        poiss2_dtb(1)=0;
        poiss2_dtb=poiss2_dtb./sum(poiss2_dtb(:));
        logPP=log(poiss2_dtb);
        lld_group_poiss2 = -sum(logPP(round(ndata2)+1));
        if lld_group_poiss2 <lld_group_poiss2_best2
            para=[x_group(ii),np1_group(jj)];
            lld_group_poiss2_best2=lld_group_poiss2;
            poiss2_dtb_best2=poiss2_dtb;
        end
    end
end