function [foci_data,fake_data,h,t_absolute,varargout] = dual_local_local3D(foci_bw00,max_image00,mask_stack,signal_stack,RNA_stack,nucleus_protein_profile,image_folder,N_cycle,varargin)


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A function to plot the local transcription factor concentration at transcription sites %%
%% foci_data: foci intensity/mean protein intensity profile/local protein intensities/mean protein intensities at foci layers/local RNA intensities
%% fake_data: fake foci intensity/mean protein intensity profile/local protein intensities/local RNA intensities

%% foci_bw00: transcription foci center 3D mask
%% max_image00: transcription foci peak intensity (3D)
%% SS: transcription foci intensity area (3D)
%% mask_stack: 3D nuclei mask
%% signal_stack: 3D transcription factor image stack
%% nucleus_protein_profile: nucleus protein mean intensity
%% image_folder: image folder name
%% N_cycle: embryo cycle
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Initialization: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%h = fspecial('gaussian', 7, 2);   %%% Local mask for foci protein concentration
Lmin = 0.15;
Lmax = 0.85;
N_bin = 15;
rbin = 2;
Nmean = 5;
N_fit = 40;
r_edge = 10;

N_fb = 20;
w_fb = 0.05;
N_pb = 20;
w_pb = 0.05;
ppbin = [-1:0.05:1];   %%% binning of the ratio histogram
Npp = length(ppbin);


r_size = 0:1:9;
varargout = {r_size};
h = cell(size(r_size));
mr_size = 3*ones(size(r_size));
mh = cell(size(r_size));
area_h = zeros(size(h));
if isempty(varargin)
    resolution = 1;
    resolutionz = 1/2;
    nM_con = 1;
    unit1 = 'A.U.';
    unit2 = 'A.U.';
    t_absolute = false;
else
    resolution = varargin{1};
    resolutionz = varargin{2};
    nM_con = 6.02e8;
    unit1 = '#';
    unit2 = 'M';
    t_absolute = true;
end

foci_data = cell(size(r_size));   %%% foci data matrix initialization
fake_data = cell(size(r_size));   %%% control data matrix initialization
p_ratio0 = cell(size(r_size));
fp_ratio0 = cell(size(r_size));
op_ratio0 = cell(size(r_size));
ofp_ratio0 = cell(size(r_size));
nucleus_protein_profile((nucleus_protein_profile(:,1) < Lmin)|(nucleus_protein_profile(:,1) > Lmax),:) = nan(sum((nucleus_protein_profile(:,1) < Lmin)|(nucleus_protein_profile(:,1) > Lmax)),size(nucleus_protein_profile,2));
nucleus_protein_profile = [nan(1,size(nucleus_protein_profile,2));nucleus_protein_profile];

foci_raw = foci_bw00;

for ir = 1:length(r_size)
    h{ir} = double(getnhood(strel('disk',r_size(ir))));
    mh{ir} = double(getnhood(strel('disk',r_size(ir)+mr_size(ir))));
    hxy = size(h{ir});
    mhxy = size(mh{ir});
    mh{ir}((mhxy(1)-hxy(1))/2+1:(mhxy(1)+hxy(1))/2,(mhxy(2)-hxy(2))/2+1:(mhxy(2)+hxy(2))/2) = mh{ir}((mhxy(1)-hxy(1))/2+1:(mhxy(1)+hxy(1))/2,(mhxy(2)-hxy(2))/2+1:(mhxy(2)+hxy(2))/2)-h{ir};
%    mh{ir} = mh{ir}*sum(h{ir}(:))/sum(mh{ir}(:));
    h{ir} = h{ir}*resolution*resolution*2*resolutionz*nM_con;
    mh{ir} = mh{ir}*resolution*resolution*2*resolutionz*nM_con;
    
    %h{ir} = h{ir}/sum(sum(h{ir}));
    area_h(ir) = sum(sum(h{ir}));

    %%% Fake spot mask generation: %%% ========================================
    foci_bw00 = foci_raw;
    foci_bw = max(foci_bw00,[],3);   %%% 2D foci mask
    fake_bw00 = false(size(foci_bw00));
    bwi = logical(mask_stack);
    bw0 = imerode(bwi,strel('disk',r_edge+r_size(ir))).*mask_stack;
    bwf = bw0.*foci_bw00;
    foci_bw00 = logical(bwf);
    foci_bw_neg = (~imdilate(foci_bw,strel('disk',2*r_size(ir))));

    f_index = find(bwf);
    for I_f = 1:length(f_index)
        [~,~,foci_layer] = ind2sub(size(mask_stack),f_index(I_f));
        nucleus_prop = (bw0(:,:,foci_layer) == bwf(f_index(I_f))) & foci_bw_neg;
        nucleus_prop_Area = sum(nucleus_prop(:));
        nucleus_prop_PixelIdxList = find(nucleus_prop);
        try 
            if nucleus_prop_Area > 1
                tcontinue = true;
                while tcontinue
                    random_I = ceil(rand(1)*nucleus_prop_Area);
                    I_temp = nucleus_prop_PixelIdxList(random_I);
                    if ~fake_bw00([ind2sub(size(foci_bw),I_temp),foci_layer])
                        fake_bw00([ind2sub(size(foci_bw),I_temp),foci_layer]) = true;
                        tcontinue = false;
                    end
                end
           else
                %random_I = nan;
                foci_bw00(f_index(I_f)) = false;
            end
        catch
            foci_bw00(f_index(I_f)) = false;
        end
        
    end
    foci_data{ir} = zeros(sum(foci_bw00(:)),5);   %%% foci data matrix initialization
    fake_data{ir} = zeros(sum(fake_bw00(:)),5);   %%% control data matrix initialization
    %%% =======================================================================


    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Calculate local protein concentration: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    bw01 = imerode(bwi,strel('disk',r_edge));
    foci_bw00 = logical(foci_bw00);
    fake_bw00 = logical(fake_bw00);
    ind_all = find(foci_bw00 | fake_bw00);
    [x_all,y_all,z_all] = ind2sub(size(mask_stack),ind_all);
    x_h = floor(size(h{ir},1)/2);
    y_h = floor(size(h{ir},2)/2);
    f_signal = zeros(size(signal_stack));   %%% Filtered signal
    r_signal = zeros(size(signal_stack));   %%% Filtered RNA signal
    f_mean = zeros(size(signal_stack));   %%% Filtered signal
    mx_h = floor(size(mh{ir},1)/2);
    my_h = floor(size(mh{ir},2)/2);

    for I_ind = 1:length(x_all)
        f_signal(x_all(I_ind),y_all(I_ind),z_all(I_ind)) = sum(sum(signal_stack(x_all(I_ind)-x_h:x_all(I_ind)+x_h,y_all(I_ind)-y_h:y_all(I_ind)+y_h,z_all(I_ind)).*h{ir}));
        r_signal(x_all(I_ind),y_all(I_ind),z_all(I_ind)) = sum(sum(RNA_stack   (x_all(I_ind)-x_h:x_all(I_ind)+x_h,y_all(I_ind)-y_h:y_all(I_ind)+y_h,z_all(I_ind)).*h{ir}));
        f_ratio = sum(sum(bw01(x_all(I_ind)-x_h:x_all(I_ind)+x_h,y_all(I_ind)-y_h:y_all(I_ind)+y_h,z_all(I_ind)).*h{ir}))/sum(sum(bw01(x_all(I_ind)-mx_h:x_all(I_ind)+mx_h,y_all(I_ind)-my_h:y_all(I_ind)+my_h,z_all(I_ind)).*mh{ir}));
        f_mean(x_all(I_ind),y_all(I_ind),z_all(I_ind)) = sum(sum(signal_stack(x_all(I_ind)-mx_h:x_all(I_ind)+mx_h,y_all(I_ind)-my_h:y_all(I_ind)+my_h,z_all(I_ind)).*mh{ir}))*f_ratio;
    end

    foci_data{ir}(:,4) = f_signal(foci_bw00);
    fake_data{ir}(:,4) = f_signal(fake_bw00);
    foci_data{ir}(:,5) = r_signal(foci_bw00);
    fake_data{ir}(:,5) = r_signal(fake_bw00);
    foci_data{ir}(:,3) = f_mean(foci_bw00);
    fake_data{ir}(:,3) = f_mean(fake_bw00);

    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Calculate mean protein intensity at different plains: %%%%%%%%%%%%%%%%%%
%     mean3_nu = zeros(size(foci_bw00));
%     for I_layer = 1:size(signal_stack,3)
%         nu_prop = regionprops(imerode(bwi, strel('disk',r_size(ir)+r_edge)).*bwn,conv2(signal_stack(:,:,I_layer),h{ir},'same'),'MeanIntensity');
%         temp_nu = [0,[nu_prop.MeanIntensity]];
%         mean3_nu((foci_layer == I_layer)|(fake_layer == I_layer)) = temp_nu(bwn((foci_layer == I_layer)|(fake_layer == I_layer))+1);
%     end
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Collect foci data (foci intensity, foci local protein concentration, mean protein concentration)
    foci_data{ir}(:,1) = max_image00(foci_bw00);
    foci_data{ir}(:,2) = nucleus_protein_profile((mask_stack(foci_bw00)+1),2);

    fake_data{ir}(:,1) = max_image00(fake_bw00);
    fake_data{ir}(:,2) = nucleus_protein_profile((mask_stack(fake_bw00)+1),2);
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    foci_data{ir} = foci_data{ir}(~isnan(foci_data{ir}(:,2)),:);
    fake_data{ir} = fake_data{ir}(~isnan(fake_data{ir}(:,2)),:);
    
    p_ratio0{ir} = foci_data{ir}(:,4)./foci_data{ir}(:,3)-1;
    fp_ratio0{ir} = fake_data{ir}(:,4)./fake_data{ir}(:,3)-1;
    op_ratio0{ir} = foci_data{ir}(:,4)-foci_data{ir}(:,3);
    ofp_ratio0{ir} = fake_data{ir}(:,4)-fake_data{ir}(:,3);
end

%% Figure plotting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
color_RGB0 = [1,0,0;0,1,0;0,0,1;0,1,1;1,0,1;1,1,0];
color_RGB = [color_RGB0;[0,0,0];color_RGB0/2;[0.5,0.5,0.5]];
if length(r_size) > size(color_RGB0,1)
    LL2 = ceil(length(r_size)/size(color_RGB0,1));
    color_RGB  = repmat(color_RGB,LL2,1);
end    
figure(73)
clf
subplot(1,2,1)
legend_text = cell(0);
for ir = 1:length(r_size)
    plot(p_ratio0{ir},foci_data{ir}(:,1),'o','color',color_RGB(ir,:))
    hold on
    [xtemp,IX] = sort(p_ratio0{ir});
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
    errorbarxy(x_f,y_f,x_ferr,y_ferr,x_ferr,y_ferr,color_RGB(ir,:),color_RGB(ir,:))
    legend_text = cat(2,legend_text,{['local intensity, r = ',num2str(r_size(ir))],['binned data, r = ',num2str(r_size(ir))]});
end
title(['Foci intensity vs foci local protein enrichment: ',image_folder,', cycle = ',num2str(N_cycle)],'Interpreter','none')
xlabel('Protein enrichment')
ylabel('Foci intensity (A.U.)')
legend(legend_text)

    
subplot(1,2,2)
legend_text = cell(0);
for ir = 1:length(r_size)
    plot(foci_data{ir}(:,2),p_ratio0{ir},'*','color',color_RGB(2*ir-1,:))
    hold on
    plot(fake_data{ir}(:,2),fp_ratio0{ir},'.','color',color_RGB(2*ir,:))
    
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
    errorbarxy(x1_f,y1_f,x1_ferr,y1_ferr,x1_ferr,y1_ferr,color_RGB(2*ir-1,:),color_RGB(2*ir-1,:),'o')
    errorbarxy(x2_f,y2_f,x2_ferr,y2_ferr,x2_ferr,y2_ferr,color_RGB(2*ir,:),color_RGB(2*ir,:),'*')
    plot(x1_f,y1_f-y2_f,'--','color',mean(color_RGB((2*ir-1):(2*ir),:)))
    
    x_limit = xlim;
    %%% Nonlinear fitting 1: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [~,I_middle] = min(abs(y1_f-(mean(y1_f((end-Nmean+1):end))+mean(y1_f(1:Nmean)))/2));
    beta0 = [2,x1_f(I_middle),(mean(y1_f((end-Nmean+1):end))-mean(y1_f(1:Nmean))),mean(y1_f(1:Nmean))];
    try
        [beta1,r,~,~,~] = nlinfit(x1_f,y1_f,@Hill,beta0);
        x_fit = x1_f(1):(x1_f(end)-x1_f(1))/N_fit:x1_f(end);
        plot(x_fit,Hill(beta1,x_fit),'--','color',color_RGB(2*ir-1,:))
    catch err
        plot(x_limit,[0,0],'--','color',color_RGB(2*ir-1,:))
        beta1 = [nan,nan,nan,nan];
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%% Nonlinear fitting 2: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    try
        [beta2,r,~,~,~] = nlinfit(xtemp1,ytemp1,@Hill,beta0);
        x_fit = xtemp1(1):(xtemp1(end)-xtemp1(1))/N_fit:xtemp1(end);
        plot(x_fit,Hill(beta2,x_fit),'-','color',color_RGB(2*ir-1,:))
    catch err
        plot(x_limit,[0,0],'-','color',color_RGB(2*ir-1,:))
        beta2 = [nan,nan,nan,nan];
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    legend_text = cat(2,legend_text,{['foci spots, r = ',num2str(r_size(ir))],['control spots, r = ',num2str(r_size(ir))],['binned data, r = ',num2str(r_size(ir))],['binned control data, r = ',num2str(r_size(ir))],['difference of binned data, r = ',num2str(r_size(ir))],['Fitting of binned data, r = ',num2str(r_size(ir)),', E1 = ',num2str(beta1(4)),', E2 = ',num2str(beta1(3)+beta1(4)),', th = ',num2str(beta1(2)),', h = ',num2str(beta1(1))],['Fitting of raw data, r = ',num2str(r_size(ir)),', E1 = ',num2str(beta2(4)),', E2 = ',num2str(beta2(3)+beta2(4)),', th = ',num2str(beta2(2)),', h = ',num2str(beta2(1))]});
end
title(['Foci local protein enrichment vs nuclear mean protein intensity: ',image_folder,', cycle = ',num2str(N_cycle)],'Interpreter','none')
ylabel('Protein enrichment')
xlabel(['P_m_e_a_n (',unit2,')'])
legend(legend_text)


figure(74)
clf
for ir = 1:length(r_size)
    subplot(2,ceil(length(r_size)/2),ir)
    n_foci = hist(p_ratio0{ir},ppbin);
    n_fake = hist(fp_ratio0{ir},ppbin);
    plot(ppbin,(n_foci/sum(n_foci)*100),ppbin,(n_fake/sum(n_fake)*100));
    title(['Foci local protein enrichment histogram (r = ',num2str(r_size(ir)),'): ',image_folder,', cycle = ',num2str(N_cycle),', foci mean = ',num2str(mean(p_ratio0{ir}((~isnan(p_ratio0{ir}))&(p_ratio0{ir}>0)))),', foci std = ',num2str(std(p_ratio0{ir}((~isnan(p_ratio0{ir}))&(p_ratio0{ir}>0)))),', fake foci mean = ',num2str(mean(fp_ratio0{ir}((~isnan(fp_ratio0{ir}))&(fp_ratio0{ir}>0)))),', fake foci std = ',num2str(std(fp_ratio0{ir}((~isnan(fp_ratio0{ir}))&(fp_ratio0{ir}>0))))],'Interpreter','none')
    ylabel('%')
    xlabel('Protein enrichment')
    legend('foci spots','control spots')
end


figure(75)
clf
legend_text = cell(0);
for ir = 1:length(r_size)
    plot(op_ratio0{ir},foci_data{ir}(:,5),'o','color',color_RGB(ir,:))
    hold on
    
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
    errorbarxy(x_f,y_f,x_ferr,y_ferr,x_ferr,y_ferr,color_RGB(ir,:),color_RGB(ir,:))
    
    legend_text = cat(2,legend_text,{['foci spots, r = ',num2str(r_size(ir))],['binned data, r = ',num2str(r_size(ir))]});
end
title(['Foci local protein extra fluorescence vs RNA fluorescence: ',image_folder,', cycle = ',num2str(N_cycle)],'Interpreter','none')
xlabel('P_e_x_t_r_a')
ylabel('Foci fluorescence (A.U.)')
legend(legend_text)


figure(76)
clf
subplot(1,2,1)
legend_text = cell(0);
for ir = 1:length(r_size)
    plot(op_ratio0{ir},foci_data{ir}(:,1),'o','color',color_RGB(ir,:))
    hold on
    
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
    errorbarxy(x_f,y_f,x_ferr,y_ferr,x_ferr,y_ferr,color_RGB(ir,:),color_RGB(ir,:))

    legend_text = cat(2,legend_text,{['local intensity, r = ',num2str(r_size(ir))],['binned data, r = ',num2str(r_size(ir))]});
end
title(['Foci local protein intensity vs foci intensity: ',image_folder,', cycle = ',num2str(N_cycle)],'Interpreter','none')
xlabel(['Protein enrichment (',unit1,')'])
ylabel('Foci intensity (A.U.)')
legend(legend_text)

    
subplot(1,2,2)
legend_text = cell(0);
for ir = 1:length(r_size)
    plot(foci_data{ir}(:,2),op_ratio0{ir},'*','color',color_RGB(2*ir-1,:))
    hold on
    plot(fake_data{ir}(:,2),ofp_ratio0{ir},'.','color',color_RGB(2*ir,:))

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
    errorbarxy(x1_f,y1_f,x1_ferr,y1_ferr,x1_ferr,y1_ferr,color_RGB(2*ir-1,:),color_RGB(2*ir-1,:),'o')
    errorbarxy(x2_f,y2_f,x2_ferr,y2_ferr,x2_ferr,y2_ferr,color_RGB(2*ir,:),color_RGB(2*ir,:),'*')
    plot(x1_f,y1_f-y2_f,'--','color',mean(color_RGB((2*ir-1):(2*ir),:)))

    x_limit = xlim;
    %%% Nonlinear fitting 1: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [~,I_middle] = min(abs(y1_f-(mean(y1_f((end-Nmean+1):end))+mean(y1_f(1:Nmean)))/2));
    beta0 = [2,x1_f(I_middle),(mean(y1_f((end-Nmean+1):end))-mean(y1_f(1:Nmean))),mean(y1_f(1:Nmean))];
    try
        [beta1,r,~,~,~] = nlinfit(x1_f,y1_f,@Hill,beta0);
        x_fit = x1_f(1):(x1_f(end)-x1_f(1))/N_fit:x1_f(end);
        plot(x_fit,Hill(beta1,x_fit),'--','color',color_RGB(2*ir-1,:))
    catch err
        plot(x_limit,[0,0],'--','color',color_RGB(2*ir-1,:))
        beta1 = [nan,nan,nan,nan];
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%% Nonlinear fitting 2: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    try
        [beta2,r,~,~,~] = nlinfit(xtemp1,ytemp1,@Hill,beta0);
        x_fit = xtemp1(1):(xtemp1(end)-xtemp1(1))/N_fit:xtemp1(end);
        plot(x_fit,Hill(beta2,x_fit),'-','color',color_RGB(2*ir-1,:))
    catch err
        plot(x_limit,[0,0],'-','color',color_RGB(2*ir-1,:))
        beta2 = [nan,nan,nan,nan];
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    legend_text = cat(2,legend_text,{['foci spots, r = ',num2str(r_size(ir))],['control spots, r = ',num2str(r_size(ir))],['binned data, r = ',num2str(r_size(ir))],['binned control data, r = ',num2str(r_size(ir))],['difference of binned data, r = ',num2str(r_size(ir))],['Fitting of binned data, r = ',num2str(r_size(ir)),', E1 = ',num2str(beta1(4)),', E2 = ',num2str(beta1(3)+beta1(4)),', th = ',num2str(beta1(2)),', h = ',num2str(beta1(1))],['Fitting of raw data, r = ',num2str(r_size(ir)),', E1 = ',num2str(beta2(4)),', E2 = ',num2str(beta2(3)+beta2(4)),', th = ',num2str(beta2(2)),', h = ',num2str(beta2(1))]});
end
title(['Foci local protein intensity vs nuclear mean protein intensity: ',image_folder,', cycle = ',num2str(N_cycle)],'Interpreter','none')
ylabel(['Protein enrichment (',unit1,')'])
xlabel(['P_m_e_a_n (',unit2,')'])
legend(legend_text)


figure(77)
clf
for ir = 1:length(r_size)
    subplot(2,ceil(length(r_size)/2),ir)
    [n_foci,xout_foci] = hist(op_ratio0{ir},Npp);
    [n_fake,xout_fake] = hist(ofp_ratio0{ir},Npp);
    plot(xout_foci,(n_foci/sum(n_foci)*100),xout_fake,(n_fake/sum(n_fake)*100));
    title(['Foci local protein intensity histogram (r = ',num2str(r_size(ir)),'): ',image_folder,', cycle = ',num2str(N_cycle),', foci mean = ',num2str(mean(op_ratio0{ir}((~isnan(op_ratio0{ir}))&(op_ratio0{ir}>0)))),', foci std = ',num2str(std(op_ratio0{ir}((~isnan(op_ratio0{ir}))&(op_ratio0{ir}>0)))),', fake foci mean = ',num2str(mean(ofp_ratio0{ir}((~isnan(ofp_ratio0{ir}))&(ofp_ratio0{ir}>0)))),', fake foci std = ',num2str(std(ofp_ratio0{ir}((~isnan(ofp_ratio0{ir}))&(ofp_ratio0{ir}>0))))],'Interpreter','none')
    ylabel('%')
    xlabel(['Protein enrichment (',unit1,')'])
    legend('foci spots','control spots')
end

figure(80)
clf
    mean_p = zeros(size(r_size));
    mean_fp = zeros(size(r_size));
    std_p = zeros(size(r_size));
    std_fp = zeros(size(r_size));
    for ir = 1:length(r_size)
%         [~,IX1] = sort(foci_data{ir}(:,2));
%         [~,IX2] = sort(fake_data{ir}(:,2));
%         ytemp1 = p_ratio0{ir}(IX1);
%         ytemp2 = fp_ratio0{ir}(IX2);
%         temp_ind = (~isnan(ytemp1)) & (~isnan(ytemp2));
%         mean_p(ir) = mean(ytemp1(temp_ind));
%         mean_fp(ir) = mean(ytemp2(temp_ind));
        mean_p(ir) = mean(p_ratio0{ir}(~isnan(p_ratio0{ir})));
        mean_fp(ir) = mean(fp_ratio0{ir}(~isnan(fp_ratio0{ir})));
        std_p(ir) = std0(p_ratio0{ir}(~isnan(p_ratio0{ir})));
        std_fp(ir) = std0(fp_ratio0{ir}(~isnan(fp_ratio0{ir})));
    end
    hold on
    errorbar(r_size,mean_p,std_p,'go-');
    errorbar(r_size,mean_fp,std_fp,'b*-');
    plot(r_size,mean_p-mean_fp,'r--')
    title(['Mean relative enrichment vs integration radius: ',image_folder,', cycle = ',num2str(N_cycle)],'Interpreter','none')
    ylabel('Protein enrichment')
    xlabel(['Radius (pixel)'])
    legend('foci spots','control spots','Difference')

figure(81)
clf
    mean_op = zeros(size(r_size));
    mean_ofp = zeros(size(r_size));
    std_op = zeros(size(r_size));
    std_ofp = zeros(size(r_size));
    for ir = 1:length(r_size)
%         [~,IX1] = sort(foci_data{ir}(:,2));
%         [~,IX2] = sort(fake_data{ir}(:,2));
%         ytemp1 = op_ratio0{ir}(IX1);
%         ytemp2 = ofp_ratio0{ir}(IX2);
%         temp_ind = (~isnan(ytemp1)) & (~isnan(ytemp2));
%         mean_op(ir) = mean(ytemp1(temp_ind));
%         mean_ofp(ir) = mean(ytemp2(temp_ind));
        mean_op(ir) = mean(op_ratio0{ir}(~isnan(op_ratio0{ir})));
        mean_ofp(ir) = mean(ofp_ratio0{ir}(~isnan(ofp_ratio0{ir})));
        std_op(ir) = std0(op_ratio0{ir}(~isnan(op_ratio0{ir})));
        std_ofp(ir) = std0(ofp_ratio0{ir}(~isnan(ofp_ratio0{ir})));
    end
    hold on
    errorbar(r_size,mean_op,std_p,'go-');
    errorbar(r_size,mean_ofp,std_fp,'b*-');
    plot(r_size,mean_op-mean_ofp,'r--')
    title(['Mean absolute enrichment vs integration radius: ',image_folder,', cycle = ',num2str(N_cycle)],'Interpreter','none')
    ylabel(['Protein enrichment (',unit1,')'])
    xlabel(['Radius (pixel)'])
    legend('foci spots','control spots','Difference')






    
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
