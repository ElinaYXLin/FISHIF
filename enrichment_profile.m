function [out73a,out73b,out74,out75,out76a,out76b,out77,p_ratio0,fp_ratio0,op_ratio0,ofp_ratio0,r_size,legend_text0,legend_text,legend_text1] = enrichment_profile(foci_data,fake_data,h,r_size,image_folder,mark0)


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A function to plot the local transcription factor concentration at transcription sites %%
%% foci_data: foci intensity/mean protein intensity profile/local protein intensities/mean protein intensities at foci layers/local RNA intensities
%% fake_data: fake foci intensity/mean protein intensity profile/local protein intensities/local RNA intensities

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Initialization: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%h = fspecial('gaussian', 7, 2);   %%% Local mask for foci protein concentration
% Lmin = 0.15;
% Lmax = 0.85;
N_bin = 15;
rbin = 2;
N_bin0 = floor(N_bin/rbin);

area_h = zeros(size(h));
%r_size = zeros(size(h));
for ir = 1:length(h)
    area_h(ir) = sum(sum(h{ir}));
%    r_size(ir) = floor(size(h{ir},1)/2);
end


ppbin = [-1:0.05:1];   %%% binning of the ratio histogram
rrbin = [-3:0.1:3];

out73a = zeros(N_bin0,ir*2);
out73b = zeros(N_bin ,ir*6);
out74  = zeros(ir*4,length(ppbin));
out75  = zeros(N_bin0,ir*2);
out76a = zeros(N_bin0,ir*2);
out76b = zeros(N_bin ,ir*6);
out77  = zeros(ir*4,length(rrbin));

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Figure plotting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p_ratio0 = cell(size(r_size));
fp_ratio0 = cell(size(r_size));
op_ratio0 = cell(size(r_size));
ofp_ratio0 = cell(size(r_size));
for ir = 1:length(h)
    p_ratio0{ir} = foci_data{ir}(:,4)./foci_data{ir}(:,3)-1;
    fp_ratio0{ir} = fake_data{ir}(:,4)./fake_data{ir}(:,3)-1;
    op_ratio0{ir} = foci_data{ir}(:,4)-foci_data{ir}(:,3);
    ofp_ratio0{ir} = fake_data{ir}(:,4)-fake_data{ir}(:,3);
end



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
    [xtemp,IX] = sort(p_ratio0{ir});
    ytemp = foci_data{ir}(IX,1);
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
    legend_text0 = cat(2,legend_text0,{[image_folder,'binned data, r = ',num2str(r_size(ir))]});
    out73a(:,ir*2-1) = x_f';
    out73a(:,ir*2) = y_f';
end
% title(['Foci intensity vs foci local protein enrichment: ',image_folder,', cycle = ',num2str(N_cycle)])
% xlabel('Protein enrichment')
% ylabel('Foci intensity (A.U.)')
% legend(legend_text)

    
subplot(1,2,2)
hold on
legend_text = cell(0);
for ir = 1:length(h)
%     plot(foci_data(:,2),p_ratio0(:,ir),'*','color',color_RGB(2*ir-1,:))
%     hold on
%     plot(fake_data(:,2),fp_ratio0(:,ir),'.','color',color_RGB(2*ir,:))
    
    [xtemp1,IX1] = sort(foci_data{ir}(:,2));
    [xtemp2,IX2] = sort(fake_data{ir}(:,2));
    ytemp1 = p_ratio0{ir}(IX1);
    ytemp2 = fp_ratio0{ir}(IX2);
%     size(ytemp1)
%     size(ytemp2)
%     size(p_ratio0)
%     size(fp_ratio0)
%     size(IX1)
%     size(IX2)
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
    
    out73b(:,ir*6-5) = x1_f';
    out73b(:,ir*6-4) = y1_f';
    out73b(:,ir*6-3) = x2_f';
    out73b(:,ir*6-2) = y2_f';
    out73b(:,ir*6-1) = x1_f';
    out73b(:,ir*6  ) = y1_f'-y2_f';

    
    legend_text = cat(2,legend_text,{[image_folder,'binned data, r = ',num2str(r_size(ir))],[image_folder,'binned control data, r = ',num2str(r_size(ir))],[image_folder,'difference of binned data, r = ',num2str(r_size(ir))]});
end
% title(['Foci local protein enrichment vs nuclear mean protein intensity: ',image_folder,', cycle = ',num2str(N_cycle)])
% ylabel('Protein enrichment')
% xlabel('P_m_e_a_n (A.U.)')
% legend(legend_text)


figure(74)
% clf
legend_text1 = cell(0);
for ir = 1:length(h)
    subplot(2,ceil(length(h)/2),ir)
    hold on
    n_foci = hist(p_ratio0{ir},ppbin);
    n_fake = hist(fp_ratio0{ir},ppbin);
    plot(ppbin,(n_foci/sum(n_foci)*100),mark0,ppbin,(n_fake/sum(n_fake)*100),mark0);
%     title(['Foci local protein enrichment histogram (r = ',num2str(r_size(ir)),'): ',image_folder,', cycle = ',num2str(N_cycle),', foci mean = ',num2str(mean(p_ratio0((~isnan(p_ratio0(:,ir)))&(p_ratio0(:,ir)>0),ir))),', foci std = ',num2str(std(p_ratio0((~isnan(p_ratio0(:,ir)))&(p_ratio0(:,ir)>0),ir))),', fake foci mean = ',num2str(mean(fp_ratio0((~isnan(fp_ratio0(:,ir)))&(fp_ratio0(:,ir)>0),ir))),', fake foci std = ',num2str(std(fp_ratio0((~isnan(fp_ratio0(:,ir)))&(fp_ratio0(:,ir)>0),ir)))])
%     ylabel('%')
%     xlabel('Protein enrichment')
%     legend('foci spots','control spots')
    out74(ir*4-3,:) = ppbin;
    out74(ir*4-2,:) = (n_foci/sum(n_foci)*100);
    out74(ir*4-1,:) = ppbin;
    out74(ir*4  ,:) = (n_fake/sum(n_fake)*100);    
end
legend_text1 = cat(2,legend_text1,{[image_folder,'foci spots'],[image_folder,'control spots']});


figure(75)
% clf
% legend_text = cell(0);
hold on
for ir = 1:length(h)
%     plot(op_ratio0(:,ir),foci_data(:,3+length(h)+ir),'o','color',color_RGB(ir,:),mark0)
    
    [xtemp,IX] = sort(op_ratio0{ir});
    ytemp = foci_data{ir}(IX,5);
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
    out75(:,ir*2-1) = x_f';
    out75(:,ir*2) = y_f';

%     legend_text = cat(2,legend_text,{['foci spots, r = ',num2str(r_size(ir))],['binned data, r = ',num2str(r_size(ir))]});
end
% title(['Foci local protein extra fluorescence vs RNA fluorescence: ',image_folder,', cycle = ',num2str(N_cycle)])
% xlabel('P_e_x_t_r_a (A.U.)')
% ylabel('Foci fluorescence (A.U.)')
% legend(legend_text)


figure(76)
% clf
subplot(1,2,1)
hold on
% legend_text = cell(0);
for ir = 1:length(h)
%     plot(op_ratio0(:,ir),foci_data(:,1),'o','color',color_RGB(ir,:))
%     hold on
    
    [xtemp,IX] = sort(op_ratio0{ir});
    ytemp = foci_data{ir}(IX,1);
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
    out76a(:,ir*2-1) = x_f';
    out76a(:,ir*2) = y_f';

%     legend_text = cat(2,legend_text,{['local intensity, r = ',num2str(r_size(ir))],['binned data, r = ',num2str(r_size(ir))]});
end
% title(['Foci local protein intensity vs foci intensity: ',image_folder,', cycle = ',num2str(N_cycle)])
% xlabel('Protein intensity increasement (A.U.)')
% ylabel('Foci intensity (A.U.)')
% legend(legend_text)

    
subplot(1,2,2)
% legend_text = cell(0);
hold on
for ir = 1:length(h)
%     plot(foci_data(:,2),op_ratio0(:,ir),'*','color',color_RGB(2*ir-1,:))
%     hold on
%     plot(fake_data(:,2),ofp_ratio0(:,ir),'.','color',color_RGB(2*ir,:))

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
    errorbarxy(x1_f,y1_f,x1_ferr,y1_ferr,x1_ferr,y1_ferr,color_RGB(2*ir-1,:),color_RGB(2*ir-1,:),mark0)
    errorbarxy(x2_f,y2_f,x2_ferr,y2_ferr,x2_ferr,y2_ferr,color_RGB(2*ir,:),color_RGB(2*ir,:),mark0)
    plot(x1_f,y1_f-y2_f,'color',mean(color_RGB((2*ir-1):(2*ir),:)),'LineStyle',mark0)

    out76b(:,ir*6-5) = x1_f';
    out76b(:,ir*6-4) = y1_f';
    out76b(:,ir*6-3) = x2_f';
    out76b(:,ir*6-2) = y2_f';
    out76b(:,ir*6-1) = x1_f';
    out76b(:,ir*6  ) = y1_f'-y2_f';

    
%     legend_text = cat(2,legend_text,{['foci spots, r = ',num2str(r_size(ir))],['control spots, r = ',num2str(r_size(ir))],['binned data, r = ',num2str(r_size(ir))],['binned control data, r = ',num2str(r_size(ir))],['difference of binned data, r = ',num2str(r_size(ir))]});
end
% title(['Foci local protein intensity vs nuclear mean protein intensity: ',image_folder,', cycle = ',num2str(N_cycle)])
% ylabel('Protein intensity increasement (A.U.)')
% xlabel('P_m_e_a_n (A.U.)')
% legend(legend_text)


figure(77)
% clf
for ir = 1:length(h)
    subplot(2,ceil(length(h)/2),ir)
    hold on
    n_foci = hist(op_ratio0{ir},rrbin);
    n_fake = hist(ofp_ratio0{ir},rrbin);
    plot(rrbin,(n_foci/sum(n_foci)*100),mark0,rrbin,(n_fake/sum(n_fake)*100),mark0);
%     title(['Foci local protein intensity histogram (r = ',num2str(r_size(ir)),'): ',image_folder,', cycle = ',num2str(N_cycle),', foci mean = ',num2str(mean(op_ratio0((~isnan(op_ratio0(:,ir)))&(op_ratio0(:,ir)>0),ir))),', foci std = ',num2str(std(op_ratio0((~isnan(op_ratio0(:,ir)))&(op_ratio0(:,ir)>0),ir))),', fake foci mean = ',num2str(mean(ofp_ratio0((~isnan(ofp_ratio0(:,ir)))&(ofp_ratio0(:,ir)>0),ir))),', fake foci std = ',num2str(std(ofp_ratio0((~isnan(ofp_ratio0(:,ir)))&(ofp_ratio0(:,ir)>0),ir)))])
%     ylabel('%')
%     xlabel('Protein intensity increasement')
%     legend('foci spots','control spots')
    out77(ir*4-3,:) = rrbin;
    out77(ir*4-2,:) = (n_foci/sum(n_foci)*100);
    out77(ir*4-1,:) = rrbin;
    out77(ir*4  ,:) = (n_fake/sum(n_fake)*100);    
end









    
function y = std0(x)
    NN = size(x);
    if NN(1) == 1;
        N1 = NN(2);
    else
        N1 = NN(1);
    end
    y = std(x)/sqrt(N1);
        
