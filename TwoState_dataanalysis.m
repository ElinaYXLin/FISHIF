function TwoState_dataanalysis
[filename, pathname,~]=uigetfile({'*.xls;*.xlsx'},'Select A Excel');
datapath=fullfile(pathname,filename);
savepath=[pathname,filename(1:end-5),'\'];
mkdir(savepath);
[data_list, ~] = xlsread(datapath);
data_type=size(data_list,2);
Nbin1=2000;
n1=Nbin1;
Ksimu=zeros(data_type,5);
T=[0 1 0];
lb=0*ones(1,5);
ub=100*ones(1,5);
ub(3)=0;lb([1,2])=0.01;ub(5)=1;
ub([1,2])=50;
model='gauss';
for Type_index=1:data_type
    data=data_list(:,Type_index);
    Nbin_data=0:1:200;
    nf = histc(reshape(data,1,length(data)),Nbin_data);nb=Nbin_data;
    nf=nf/sum(nf);
    state0=nf(:,1);
    if strcmp(model,'poission')
        %% poiss
        figure
        bar(nb,nf,'hist');
        ylim([0 0.1])
        xlim([0 100])
        hold on
        % poisson fit
        [para,poiss2_dtb_best2] = PoissonFitScan_singledata(data,state0(1));
        poiss2_dtb_best2(1)=state0(1);
        plot(0:100,poiss2_dtb_best2,'LineWidth',1);
        legend('Data','Poisson')
        saveas(gcf,[savepath,['Dualpioss',num2str(Type_index),'.fig']]);
        saveas(gcf,[savepath,['Dualpioss',num2str(Type_index),'.png']]);
        close(gcf);
    elseif  strcmp(model,'gauss')
        %% gauss
        nf_guass=[0,nf(1,2:end)]';
        fit_1=fit(Nbin_data',nf_guass,'gauss1');
        var_1=(fit_1.c1^2/2);
        mean_1=fit_1.b1;
        para=[fit_1.a1*(2*pi*var_1)^0.5 mean_1];
        lambda_1=var_1./mean_1;
        figure
        bar(Nbin_data',nf_guass,'hist')
        hold on
        plot(fit_1,'r')
        xlim([0 100])
        xlabel('# Num')
        ylabel('# Frequence')
        title('Gauss Fitting')
        saveas(gcf,[savepath,['Gauss',num2str(Type_index),'.fig']]);
        saveas(gcf,[savepath,['Gauss',num2str(Type_index),'.png']]);
        close(gcf);
    end
    %% ini parameter
    k_new1=[0.1*(para(1)),0.1*(1-(para(1))),0,Nbin1*para(2)/n1,0];
    %% FSP: two-state model fit
    [kk_simu,~,~] = simulannealbnd(@(x) ...
        FSP_NSX_alp(x,2,data,T(3),T(1)),k_new1,lb,ub);
    Ksimu(Type_index,:)=kk_simu;
    [~,P_dtb] =  FSP_NSX_alp(kk_simu,2,data,T(3),T(1));
    
    figure
    bar(nb,nf,'hist');
%     ymax=sort(nf);
%     ylim([0 ymax(end)+0.5*ymax(end)])
    ylim([0 0.1])
    xlim([-0.5 80])
    hold on
    plot(0:80,P_dtb(1:81,2),'c','LineWidth',1.5)
    xlabel('# Number')
    ylabel('Frequency')
    axis square
    saveas(gcf,[savepath,'Two_state',num2str(Type_index),'.fig']);
    saveas(gcf,[savepath,'Two_state',num2str(Type_index),'.png']);
    close(gcf);
end
save([savepath,['K','.mat']],'Ksimu')