function [para,poiss2_dtb_best2] = PoissonFitScan_single(ndata2,state02)
%% A program to fit data by poisson distribution

%% RNA2
np1_group=[1:0.1:50];
x_group=0:0.01:1-state02;
% ini_group=[x_group',y_group'];
P_double=@(k,x,y,np1,np2)x*poisspdf(k,np1)+y*poisspdf(k,np2);
% ndata2 = ndata2(ndata2 <= 30);
ndata2 = ndata2(ndata2 <= 50&ndata2 >= 4);
%ndata2 = ndata2(ndata2 ~= 0);
lld_group_poiss2_best2=inf;
% lld_group_poiss2=zeros(length(x_group),length(np1_group),length(np2_group));
% ini_group=repmat(ini_group,length(np1_group)*length(np2_group),1);
% ini_group=[ini_group,repmat(np1_group',length(x_group)*length(np2_group),1),repmat(np2_group',length(x_group)*length(np1_group),1)];
for ii = 1:length(x_group)
    %for iib=1:length(y_group)
    for jj = 1:length(np1_group)
        for kkk =1:length(np2_group)
            poiss2_dtb=P_double(0:50,x_group(ii),y_group(ii),np1_group(jj),np2_group(kkk));
            poiss2_dtb(1:4)=0;
            poiss2_dtb=poiss2_dtb./sum(poiss2_dtb(:));
            logPP=log(poiss2_dtb);
            lld_group_poiss2 = -sum(logPP(round(ndata2)+1));
            if lld_group_poiss2 <lld_group_poiss2_best2
                para=[x_group(ii),y_group(ii),np1_group(jj),np2_group(kkk)];
                lld_group_poiss2_best2=lld_group_poiss2;
                poiss2_dtb_best2=poiss2_dtb;
            end
        end
    end
    %end
end