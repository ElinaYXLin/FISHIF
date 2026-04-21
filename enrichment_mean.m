function enrichment_mean(data73a,data73b,data74,data75,data76a,data76b,data77,foci_data,fake_data,p_ratio0,fp_ratio0,op_ratio0,ofp_ratio0,data_legend0,data_legend,data_legend1,group_name,t_absolute,r_size)

%% Initialization: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%h = fspecial('gaussian', 7, 2);   %%% Local mask for foci protein concentration
% Lmin = 0.15;
% Lmax = 0.85;
N_bin = 15;
rbin = 2;
N_bin0 = floor(N_bin/rbin);
Nmean = 5;
N_fit = 40;

if t_absolute
    unit1 = '#';
    unit2 = 'M';
else
    unit1 = 'A.U.';
    unit2 = 'A.U.';
end

%r_size = 0:(size(data73a,2)/2-1);
h = cell(size(r_size));
for ir = 1:length(r_size)
    h{ir} = double(getnhood(strel('disk',r_size(ir))))./sum(sum(double(getnhood(strel('disk',r_size(ir))))));
end

ppbin = data74(1,:);   %%% binning of the ratio histogram
rrbin = data77(1,:);
mark0 = '-';
mark1 = '--';
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Figure plotting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%color_code = 'rgbcmykrgbcmyk';
color_RGB0 = [1,0,0;0,1,0;0,0,1;0,1,1;1,0,1;1,1,0];
color_RGB = [color_RGB0;[0,0,0];color_RGB0/2;[0.5,0.5,0.5]];
if length(r_size) > size(color_RGB0,1)
    LL2 = ceil(length(r_size)/size(color_RGB0,1));
    color_RGB  = repmat(color_RGB,LL2,1);
end    

figure(73)
hold on
%clf
subplot(1,2,1)
legend_text0 = cell(0);
for ir = 1:length(h)
%    plot(p_ratio0(:,ir),foci_data(:,1),'o','color',color_RGB(ir,:))
%    hold on
    [xtemp,IX] = sort(data73a(:,(ir*2-1)));
    ytemp = data73a(IX,ir*2);
%     N_bin0 = floor(N_bin/rbin);
    N_temp = floor(length(xtemp)/N_bin0);
    x_f = zeros(1,N_bin0);
    y_f = zeros(1,N_bin0);
    x_ferr = zeros(1,N_bin0);
    y_ferr = zeros(1,N_bin0);
    for I_bin = 1:N_bin0
        x_f(I_bin) = mean(xtemp(((I_bin-1)*N_temp+1):I_bin*N_temp));
        x_ferr(I_bin) = std0(xtemp(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y_f(I_bin) = mean(ytemp(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y_ferr(I_bin) = std0(ytemp(((I_bin-1)*N_temp+1):I_bin*N_temp));
    end
    errorbarxy(x_f,y_f,x_ferr,y_ferr,x_ferr,y_ferr,color_RGB(ir,:),color_RGB(ir,:),mark0)
    legend_text0 = cat(2,legend_text0,{['Mean curve over meaned data, r = ',num2str(r_size(ir))]});

    [xtemp,IX] = sort(p_ratio0{ir});
    ytemp = foci_data{ir}(IX,1);
    %N_bin0 = floor(N_bin/rbin);
    N_temp = floor(length(xtemp)/N_bin0);
    x_f = zeros(1,N_bin0);
    y_f = zeros(1,N_bin0);
    x_ferr = zeros(1,N_bin0);
    y_ferr = zeros(1,N_bin0);
    for I_bin = 1:N_bin0
        x_f(I_bin) = mean(xtemp(((I_bin-1)*N_temp+1):I_bin*N_temp));
        x_ferr(I_bin) = std0(xtemp(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y_f(I_bin) = mean(ytemp(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y_ferr(I_bin) = std0(ytemp(((I_bin-1)*N_temp+1):I_bin*N_temp));
    end
    errorbarxy(x_f,y_f,x_ferr,y_ferr,x_ferr,y_ferr,color_RGB(ir,:),color_RGB(ir,:),mark1)
    legend_text0 = cat(2,legend_text0,{['Mean curve over raw data, r = ',num2str(r_size(ir))]});
end
title(['Foci intensity vs foci local protein enrichment: ',group_name])
xlabel('Protein enrichment')
ylabel('Foci intensity (A.U.)')
data_legend0 = cat(2,data_legend0,legend_text0);
legend(data_legend0)
legend('off')


subplot(1,2,2)
hold on
legend_text = cell(0);
for ir = 1:length(h)
%     plot(foci_data(:,2),p_ratio0(:,ir),'*','color',color_RGB(2*ir-1,:))
%     hold on
%     plot(fake_data(:,2),fp_ratio0(:,ir),'.','color',color_RGB(2*ir,:))
    
    [xtemp1,IX1] = sort(data73b(:,(ir*6-5)));
    [xtemp2,IX2] = sort(data73b(:,(ir*6-3)));
    ytemp1 = data73b(IX1,(ir*6-4));
    ytemp2 = data73b(IX2,(ir*6-2));
    temp_ind = (~isnan(ytemp1)) & (~isnan(ytemp2));
    xtemp1 = xtemp1(temp_ind);
    xtemp2 = xtemp2(temp_ind);
    ytemp1 = ytemp1(temp_ind);
    ytemp2 = ytemp2(temp_ind);
    N_temp = floor(length(xtemp1)/N_bin);
    x1_f = zeros(1,N_bin);
    x2_f = zeros(1,N_bin);
    y1_f = zeros(1,N_bin);
    y2_f = zeros(1,N_bin);
    x1_ferr = zeros(1,N_bin);
    x2_ferr = zeros(1,N_bin);
    y1_ferr = zeros(1,N_bin);
    y2_ferr = zeros(1,N_bin);

    for I_bin = 1:N_bin
        x1_f(I_bin) = mean(xtemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
        x1_ferr(I_bin) = std0(xtemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
        x2_f(I_bin) = mean(xtemp2(((I_bin-1)*N_temp+1):I_bin*N_temp));
        x2_ferr(I_bin) = std0(xtemp2(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y1_f(I_bin) = mean(ytemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y1_ferr(I_bin) = std0(ytemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y2_f(I_bin) = mean(ytemp2(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y2_ferr(I_bin) = std0(ytemp2(((I_bin-1)*N_temp+1):I_bin*N_temp));
    end
    errorbarxy(x1_f,y1_f,x1_ferr,y1_ferr,x1_ferr,y1_ferr,color_RGB(2*ir-1,:),color_RGB(2*ir-1,:),mark0)
    errorbarxy(x2_f,y2_f,x2_ferr,y2_ferr,x2_ferr,y2_ferr,color_RGB(2*ir,:),color_RGB(2*ir,:),mark0)
    plot(x1_f,y1_f-y2_f,'color',mean(color_RGB((2*ir-1):(2*ir),:)),'LineStyle',mark0)
    
    %%% Nonlinear fitting 1: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [~,I_middle] = min(abs(y1_f-(mean(y1_f((end-Nmean+1):end))+mean(y1_f(1:Nmean)))/2));
    beta0 = [2,x1_f(I_middle),(mean(y1_f((end-Nmean+1):end))-mean(y1_f(1:Nmean))),mean(y1_f(1:Nmean))];
    try
        [beta1,r,~,~,~] = nlinfit(x1_f,y1_f,@Hill,beta0);
        x_fit = x1_f(1):(x1_f(end)-x1_f(1))/N_fit:x1_f(end);
        plot(x_fit,Hill(beta1,x_fit),'--','color',color_RGB(2*ir-1,:))
    catch err
        plot([x1_f(1),x1_f(end)],[0,0],'--','color',color_RGB(2*ir-1,:))
        beta1 = [nan,nan,nan,nan];
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%% Nonlinear fitting 2: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    try
        [beta2,r,~,~,~] = nlinfit(xtemp1,ytemp1,@Hill,beta0);
        x_fit = xtemp1(1):(xtemp1(end)-xtemp1(1))/N_fit:xtemp1(end);
        plot(x_fit,Hill(beta2,x_fit),'-','color',color_RGB(2*ir-1,:))
    catch err
        plot([xtemp1(1),xtemp1(end)],[0,0],'-','color',color_RGB(2*ir-1,:))
        beta2 = [nan,nan,nan,nan];
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    legend_text = cat(2,legend_text,{['Mean of binned foci data, r = ',num2str(r_size(ir))],['Mean of binned control data, r = ',num2str(r_size(ir))],['Difference of mean of binned, r = ',num2str(r_size(ir))],['Fitting of binned mean data, r = ',num2str(r_size(ir)),', E1 = ',num2str(beta1(4)),', E2 = ',num2str(beta1(3)+beta1(4)),', th = ',num2str(beta1(2)),', h = ',num2str(beta1(1))],['Fitting of mean data, r = ',num2str(r_size(ir)),', E1 = ',num2str(beta2(4)),', E2 = ',num2str(beta2(3)+beta2(4)),', th = ',num2str(beta2(2)),', h = ',num2str(beta2(1))]});


%%% Raw data plotting:
    [xtemp1,IX1] = sort(foci_data{ir}(:,2));
    [xtemp2,IX2] = sort(fake_data{ir}(:,2));
    ytemp1 = p_ratio0{ir}(IX1);
    ytemp2 = fp_ratio0{ir}(IX2);
    temp_ind = (~isnan(ytemp1)) & (~isnan(ytemp2));
    xtemp1 = xtemp1(temp_ind);
    xtemp2 = xtemp2(temp_ind);
    ytemp1 = ytemp1(temp_ind);
    ytemp2 = ytemp2(temp_ind);
    N_temp = floor(length(xtemp1)/N_bin);
    x1_f = zeros(1,N_bin);
    x2_f = zeros(1,N_bin);
    y1_f = zeros(1,N_bin);
    y2_f = zeros(1,N_bin);
    x1_ferr = zeros(1,N_bin);
    x2_ferr = zeros(1,N_bin);
    y1_ferr = zeros(1,N_bin);
    y2_ferr = zeros(1,N_bin);

    for I_bin = 1:N_bin
        x1_f(I_bin) = mean(xtemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
        x1_ferr(I_bin) = std0(xtemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
        x2_f(I_bin) = mean(xtemp2(((I_bin-1)*N_temp+1):I_bin*N_temp));
        x2_ferr(I_bin) = std0(xtemp2(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y1_f(I_bin) = mean(ytemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y1_ferr(I_bin) = std0(ytemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y2_f(I_bin) = mean(ytemp2(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y2_ferr(I_bin) = std0(ytemp2(((I_bin-1)*N_temp+1):I_bin*N_temp));
    end
    errorbarxy(x1_f,y1_f,x1_ferr,y1_ferr,x1_ferr,y1_ferr,color_RGB(2*ir-1,:),color_RGB(2*ir-1,:),mark1)
    errorbarxy(x2_f,y2_f,x2_ferr,y2_ferr,x2_ferr,y2_ferr,color_RGB(2*ir,:),color_RGB(2*ir,:),mark1)
    plot(x1_f,y1_f-y2_f,'color',mean(color_RGB((2*ir-1):(2*ir),:)),'LineStyle',mark1)
    
    %%% Nonlinear fitting 1: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [~,I_middle] = min(abs(y1_f-(mean(y1_f((end-Nmean+1):end))+mean(y1_f(1:Nmean)))/2));
    beta0 = [2,x1_f(I_middle),(mean(y1_f((end-Nmean+1):end))-mean(y1_f(1:Nmean))),mean(y1_f(1:Nmean))];
    try
        [beta1,r,~,~,~] = nlinfit(x1_f,y1_f,@Hill,beta0);
        x_fit = x1_f(1):(x1_f(end)-x1_f(1))/N_fit:x1_f(end);
        plot(x_fit,Hill(beta1,x_fit),'-.','color',color_RGB(2*ir-1,:))
    catch err
        plot([x1_f(1),x1_f(end)],[0,0],'-.','color',color_RGB(2*ir-1,:))
        beta1 = [nan,nan,nan,nan];
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%% Nonlinear fitting 2: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    try
        [beta2,r,~,~,~] = nlinfit(xtemp1,ytemp1,@Hill,beta0);
        x_fit = xtemp1(1):(xtemp1(end)-xtemp1(1))/N_fit:xtemp1(end);
        plot(x_fit,Hill(beta2,x_fit),':','color',color_RGB(2*ir-1,:))
    catch err
        plot([xtemp1(1),xtemp1(end)],[0,0],':','color',color_RGB(2*ir-1,:))
        beta2 = [nan,nan,nan,nan];
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    legend_text = cat(2,legend_text,{['Mean of raw foci data, r = ',num2str(r_size(ir))],['Mean of raw control data, r = ',num2str(r_size(ir))],['Difference of mean of raw, r = ',num2str(r_size(ir))],['Fitting of binned raw data, r = ',num2str(r_size(ir)),', E1 = ',num2str(beta1(4)),', E2 = ',num2str(beta1(3)+beta1(4)),', th = ',num2str(beta1(2)),', h = ',num2str(beta1(1))],['Fitting of raw data, r = ',num2str(r_size(ir)),', E1 = ',num2str(beta2(4)),', E2 = ',num2str(beta2(3)+beta2(4)),', th = ',num2str(beta2(2)),', h = ',num2str(beta2(1))]});

end
title(['Foci local protein enrichment vs nuclear mean protein intensity: ',group_name])
ylabel('Protein enrichment')
xlabel(['P_m_e_a_n (',unit2,')'])
data_legend = cat(2,data_legend,legend_text);
legend(data_legend)
legend('off')


figure(74)
% clf
data_legend1 = cat(2,data_legend1,{'Mean foci data','Mean control data','Mean raw foci data','Mean raw control data'});
for ir = 1:length(h)
    subplot(2,ceil(length(h)/2),ir)
    hold on
    plot(data74(ir*4-3,:),mean(data74((ir*4-2):(length(h)*4):end,:)),mark0,data74(ir*4-1,:),mean(data74((ir*4):(length(h)*4):end,:)),mark0);

    ppbin1 = data74(ir*4-3,:);
    ppbin2 = data74(ir*4-1,:);
    n_foci = hist(p_ratio0{ir},ppbin1);
    n_fake = hist(fp_ratio0{ir},ppbin2);
    plot(ppbin1,(n_foci/sum(n_foci)*100),mark1,ppbin2,(n_fake/sum(n_fake)*100),mark1);
    title(['Foci local protein enrichment histogram (r = ',num2str(r_size(ir)),'): ',group_name])
    ylabel('%')
    xlabel('Protein enrichment')
    legend(data_legend1)
    legend('off')
end



figure(75)
% clf
% legend_text = cell(0);
hold on
for ir = 1:length(h)
%     plot(op_ratio0(:,ir),foci_data(:,3+length(h)+ir),'o','color',color_RGB(ir,:),mark0)
    
    [xtemp,IX] = sort(data75(:,(2*ir-1)));
    ytemp = data75(IX,2*ir);
%     N_bin0 = floor(N_bin/rbin);
    N_temp = floor(length(xtemp)/N_bin0);
    x_f = zeros(1,N_bin0);
    y_f = zeros(1,N_bin0);
    x_ferr = zeros(1,N_bin0);
    y_ferr = zeros(1,N_bin0);
    for I_bin = 1:N_bin0
        x_f(I_bin) = mean(xtemp(((I_bin-1)*N_temp+1):I_bin*N_temp));
        x_ferr(I_bin) = std0(xtemp(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y_f(I_bin) = mean(ytemp(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y_ferr(I_bin) = std0(ytemp(((I_bin-1)*N_temp+1):I_bin*N_temp));
    end
    errorbarxy(x_f,y_f,x_ferr,y_ferr,x_ferr,y_ferr,color_RGB(ir,:),color_RGB(ir,:),mark0)
    
    [xtemp,IX] = sort(op_ratio0{ir});
    ytemp = foci_data{ir}(IX,5);
    N_bin0 = floor(N_bin/rbin);
    N_temp = floor(length(xtemp)/N_bin0);
    x_f = zeros(1,N_bin0);
    y_f = zeros(1,N_bin0);
    x_ferr = zeros(1,N_bin0);
    y_ferr = zeros(1,N_bin0);
    for I_bin = 1:N_bin0
        x_f(I_bin) = mean(xtemp(((I_bin-1)*N_temp+1):I_bin*N_temp));
        x_ferr(I_bin) = std0(xtemp(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y_f(I_bin) = mean(ytemp(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y_ferr(I_bin) = std0(ytemp(((I_bin-1)*N_temp+1):I_bin*N_temp));
    end
    errorbarxy(x_f,y_f,x_ferr,y_ferr,x_ferr,y_ferr,color_RGB(ir,:),color_RGB(ir,:),mark1)
    
    
%     legend_text = cat(2,legend_text,{['foci spots, r = ',num2str(r_size(ir))],['binned data, r = ',num2str(r_size(ir))]});
end
title(['Foci local protein extra fluorescence vs RNA fluorescence: ',group_name])
xlabel('P_e_x_t_r_a')
ylabel('Foci fluorescence (A.U.)')
legend(data_legend0)
legend('off')


figure(76)
% clf
subplot(1,2,1)
hold on
% legend_text = cell(0);
for ir = 1:length(h)
%     plot(op_ratio0(:,ir),foci_data(:,1),'o','color',color_RGB(ir,:))
%     hold on
    
    [xtemp,IX] = sort(data76a(:,(2*ir-1)));
    ytemp = data76a(IX,2*ir);
%     N_bin0 = floor(N_bin/rbin);
    N_temp = floor(length(xtemp)/N_bin0);
    x_f = zeros(1,N_bin0);
    y_f = zeros(1,N_bin0);
    x_ferr = zeros(1,N_bin0);
    y_ferr = zeros(1,N_bin0);
    for I_bin = 1:N_bin0
        x_f(I_bin) = mean(xtemp(((I_bin-1)*N_temp+1):I_bin*N_temp));
        x_ferr(I_bin) = std0(xtemp(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y_f(I_bin) = mean(ytemp(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y_ferr(I_bin) = std0(ytemp(((I_bin-1)*N_temp+1):I_bin*N_temp));
    end
    errorbarxy(x_f,y_f,x_ferr,y_ferr,x_ferr,y_ferr,color_RGB(ir,:),color_RGB(ir,:),mark0)
%     legend_text = cat(2,legend_text,{['local intensity, r = ',num2str(r_size(ir))],['binned data, r = ',num2str(r_size(ir))]});

    [xtemp,IX] = sort(op_ratio0{ir});
    ytemp = foci_data{ir}(IX,1);
    N_bin0 = floor(N_bin/rbin);
    N_temp = floor(length(xtemp)/N_bin0);
    x_f = zeros(1,N_bin0);
    y_f = zeros(1,N_bin0);
    x_ferr = zeros(1,N_bin0);
    y_ferr = zeros(1,N_bin0);
    for I_bin = 1:N_bin0
        x_f(I_bin) = mean(xtemp(((I_bin-1)*N_temp+1):I_bin*N_temp));
        x_ferr(I_bin) = std0(xtemp(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y_f(I_bin) = mean(ytemp(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y_ferr(I_bin) = std0(ytemp(((I_bin-1)*N_temp+1):I_bin*N_temp));
    end
    errorbarxy(x_f,y_f,x_ferr,y_ferr,x_ferr,y_ferr,color_RGB(ir,:),color_RGB(ir,:),mark1)

end
title(['Foci local protein intensity vs foci intensity: ',group_name])
xlabel(['Protein enrichment (',unit1,')'])
ylabel('Foci intensity (A.U.)')
legend(data_legend0)
legend('off')

    
subplot(1,2,2)
% legend_text = cell(0);
hold on
for ir = 1:length(h)
%     plot(foci_data(:,2),op_ratio0(:,ir),'*','color',color_RGB(2*ir-1,:))
%     hold on
%     plot(fake_data(:,2),ofp_ratio0(:,ir),'.','color',color_RGB(2*ir,:))

    [xtemp1,IX1] = sort(data76b(:,(ir*6-5)));
    [xtemp2,IX2] = sort(data76b(:,(ir*6-3)));
    ytemp1 = data76b(IX1,(ir*6-4));
    ytemp2 = data76b(IX2,(ir*6-2));
    temp_ind = (~isnan(ytemp1)) & (~isnan(ytemp2));
    xtemp1 = xtemp1(temp_ind);
    xtemp2 = xtemp2(temp_ind);
    ytemp1 = ytemp1(temp_ind);
    ytemp2 = ytemp2(temp_ind);
    N_temp = floor(length(xtemp1)/N_bin);
    x1_f = zeros(1,N_bin);
    x2_f = zeros(1,N_bin);
    y1_f = zeros(1,N_bin);
    y2_f = zeros(1,N_bin);
    x1_ferr = zeros(1,N_bin);
    x2_ferr = zeros(1,N_bin);
    y1_ferr = zeros(1,N_bin);
    y2_ferr = zeros(1,N_bin);
    
    for I_bin = 1:N_bin
        x1_f(I_bin) = mean(xtemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
        x1_ferr(I_bin) = std0(xtemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
        x2_f(I_bin) = mean(xtemp2(((I_bin-1)*N_temp+1):I_bin*N_temp));
        x2_ferr(I_bin) = std0(xtemp2(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y1_f(I_bin) = mean(ytemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y1_ferr(I_bin) = std0(ytemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y2_f(I_bin) = mean(ytemp2(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y2_ferr(I_bin) = std0(ytemp2(((I_bin-1)*N_temp+1):I_bin*N_temp));
    end
    errorbarxy(x1_f,y1_f,x1_ferr,y1_ferr,x1_ferr,y1_ferr,color_RGB(2*ir-1,:),color_RGB(2*ir-1,:),mark0)
    errorbarxy(x2_f,y2_f,x2_ferr,y2_ferr,x2_ferr,y2_ferr,color_RGB(2*ir,:),color_RGB(2*ir,:),mark0)
    plot(x1_f,y1_f-y2_f,'color',mean(color_RGB((2*ir-1):(2*ir),:)),'LineStyle',mark0)
   
%     legend_text = cat(2,legend_text,{['foci spots, r = ',num2str(r_size(ir))],['control spots, r = ',num2str(r_size(ir))],['binned data, r = ',num2str(r_size(ir))],['binned control data, r = ',num2str(r_size(ir))],['difference of binned data, r = ',num2str(r_size(ir))]});

    %%% Nonlinear fitting 1: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [~,I_middle] = min(abs(y1_f-(mean(y1_f((end-Nmean+1):end))+mean(y1_f(1:Nmean)))/2));
    beta0 = [2,x1_f(I_middle),(mean(y1_f((end-Nmean+1):end))-mean(y1_f(1:Nmean))),mean(y1_f(1:Nmean))];
    try
        [beta1,r,~,~,~] = nlinfit(x1_f,y1_f,@Hill,beta0);
        x_fit = x1_f(1):(x1_f(end)-x1_f(1))/N_fit:x1_f(end);
        plot(x_fit,Hill(beta1,x_fit),'--','color',color_RGB(2*ir-1,:))
    catch err
        plot([x1_f(1),x1_f(end)],[0,0],'--','color',color_RGB(2*ir-1,:))
        beta1 = [nan,nan,nan,nan];
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%% Nonlinear fitting 2: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    try
        [beta2,r,~,~,~] = nlinfit(xtemp1,ytemp1,@Hill,beta0);
        x_fit = xtemp1(1):(xtemp1(end)-xtemp1(1))/N_fit:xtemp1(end);
        plot(x_fit,Hill(beta2,x_fit),'-','color',color_RGB(2*ir-1,:))
    catch err
        plot([xtemp1(1),xtemp1(end)],[0,0],'-','color',color_RGB(2*ir-1,:))
        beta2 = [nan,nan,nan,nan];
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% Raw data plotting:
    [xtemp1,IX1] = sort(foci_data{ir}(:,2));
    [xtemp2,IX2] = sort(fake_data{ir}(:,2));
    ytemp1 = op_ratio0{ir}(IX1);
    ytemp2 = ofp_ratio0{ir}(IX2);
    temp_ind = (~isnan(ytemp1)) & (~isnan(ytemp2));
    xtemp1 = xtemp1(temp_ind);
    xtemp2 = xtemp2(temp_ind);
    ytemp1 = ytemp1(temp_ind);
    ytemp2 = ytemp2(temp_ind);
    N_temp = floor(length(xtemp1)/N_bin);
    x1_f = zeros(1,N_bin);
    x2_f = zeros(1,N_bin);
    y1_f = zeros(1,N_bin);
    y2_f = zeros(1,N_bin);
    x1_ferr = zeros(1,N_bin);
    x2_ferr = zeros(1,N_bin);
    y1_ferr = zeros(1,N_bin);
    y2_ferr = zeros(1,N_bin);

    for I_bin = 1:N_bin
        x1_f(I_bin) = mean(xtemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
        x1_ferr(I_bin) = std0(xtemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
        x2_f(I_bin) = mean(xtemp2(((I_bin-1)*N_temp+1):I_bin*N_temp));
        x2_ferr(I_bin) = std0(xtemp2(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y1_f(I_bin) = mean(ytemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y1_ferr(I_bin) = std0(ytemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y2_f(I_bin) = mean(ytemp2(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y2_ferr(I_bin) = std0(ytemp2(((I_bin-1)*N_temp+1):I_bin*N_temp));
    end
    errorbarxy(x1_f,y1_f,x1_ferr,y1_ferr,x1_ferr,y1_ferr,color_RGB(2*ir-1,:),color_RGB(2*ir-1,:),mark1)
    errorbarxy(x2_f,y2_f,x2_ferr,y2_ferr,x2_ferr,y2_ferr,color_RGB(2*ir,:),color_RGB(2*ir,:),mark1)
    plot(x1_f,y1_f-y2_f,'color',mean(color_RGB((2*ir-1):(2*ir),:)),'LineStyle',mark1)
    
    %%% Nonlinear fitting 1: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [~,I_middle] = min(abs(y1_f-(mean(y1_f((end-Nmean+1):end))+mean(y1_f(1:Nmean)))/2));
    beta0 = [2,x1_f(I_middle),(mean(y1_f((end-Nmean+1):end))-mean(y1_f(1:Nmean))),mean(y1_f(1:Nmean))];
    try
        [beta1,r,~,~,~] = nlinfit(x1_f,y1_f,@Hill,beta0);
        x_fit = x1_f(1):(x1_f(end)-x1_f(1))/N_fit:x1_f(end);
        plot(x_fit,Hill(beta1,x_fit),'-.','color',color_RGB(2*ir-1,:))
    catch err
        plot([x1_f(1),x1_f(end)],[0,0],'-.','color',color_RGB(2*ir-1,:))
        beta1 = [nan,nan,nan,nan];
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%% Nonlinear fitting 2: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    try
        [beta2,r,~,~,~] = nlinfit(xtemp1,ytemp1,@Hill,beta0);
        x_fit = xtemp1(1):(xtemp1(end)-xtemp1(1))/N_fit:xtemp1(end);
        plot(x_fit,Hill(beta2,x_fit),':','color',color_RGB(2*ir-1,:))
    catch err
        plot([xtemp1(1),xtemp1(end)],[0,0],':','color',color_RGB(2*ir-1,:))
        beta2 = [nan,nan,nan,nan];
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
title(['Foci local protein intensity vs nuclear mean protein intensity: ',group_name])
ylabel(['Protein enrichment (',unit1,')'])
xlabel(['P_m_e_a_n (',unit2,')'])
legend(data_legend)
legend('off')


figure(77)
% clf
for ir = 1:length(h)
    subplot(2,ceil(length(h)/2),ir)
    hold on
    plot(data77(ir*4-3,:),mean(data77((ir*4-2):(length(h)*4):end,:)),mark0,data77(ir*4-1,:),mean(data77((ir*4):(length(h)*4):end,:)),mark0);
    Npp1 = data77(ir*4-3,:);
    Npp2 = data77(ir*4-1,:);
    [n_foci,xout_foci] = hist(op_ratio0{ir},Npp1);
    [n_fake,xout_fake] = hist(ofp_ratio0{ir},Npp2);
    plot(xout_foci,(n_foci/sum(n_foci)*100),mark1,xout_fake,(n_fake/sum(n_fake)*100),mark1);
    
    title(['Foci local protein intensity histogram (r = ',num2str(r_size(ir)),'): ',group_name])
    ylabel('%')
    xlabel(['Protein enrichment (',unit1,')'])
    legend(data_legend1)
    legend('off')
end









    
function y = std0(x)
    NN = size(x);
    if NN(1) == 1;
        N1 = NN(2);
    else
        N1 = NN(1);
    end
    y = std(x)/sqrt(N1);
        
function y = Hill(beta,x)
 y = beta(4)+beta(3)*x.^beta(1)./(x.^beta(1)+beta(2)^beta(1));
