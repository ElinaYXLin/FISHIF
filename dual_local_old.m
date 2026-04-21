function [foci_data,fake_data] = dual_local(foci_bw00,max_image00,SS,foci_layer,mask_stack,signal_stack,RNA_stack,nucleus_protein_profile,image_folder,N_cycle)


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A function to plot the local transcription factor concentration at transcription sites %%
%% foci_data: foci intensity/mean protein intensity profile/local protein intensities/mean protein intensities at foci layers/local RNA intensities
%% fake_data: fake foci intensity/mean protein intensity profile/local protein intensities/local RNA intensities

%% foci_bw00: transcription foci center mask
%% max_image00: transcription foci peak intensity
%% SS: transcription foci intensity area
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

r_size = 0:2;
h = cell(size(r_size));
for ir = r_size
    h{ir+1} = double(getnhood(strel('disk',ir)))./sum(sum(double(getnhood(strel('disk',ir)))));
end
foci_data = zeros(sum(sum(foci_bw00)),3+2*length(h));   %%% foci data matrix

% fake_bw00 = (rand(size(foci_bw00)) <= 0.01) & (~imdilate(foci_bw00,strel('disk',5))) & imerode(logical(max(mask_stack,[],3)),strel('disk',5));   %%% control data mask
% cum_nuclei = cumsum(mask_stack,3);
% bottom_nuclei = max((cum_nuclei == 1).*repmat(permute([1:size(signal_stack,3)],[1,3,2]),[size(foci_bw00),1]),[],3);
% fake_layer = (floor(rand(size(foci_bw00)).*max(cum_nuclei,[],3))+bottom_nuclei).*fake_bw00;   %%% control data layer
% fake_data = zeros(sum(sum(fake_bw00)),4);   %%% control data matrix

%%% Fake spot marsk generation: %%% =======================================
fake_bw00 = false(size(foci_bw00));
fake_layer = zeros(size(foci_bw00));
bwi = logical(max(mask_stack,[],3));
bw0 = bwlabel(imerode(bwi,strel('disk',5)));
nucleus_prop = regionprops((~imdilate(foci_bw00,strel('disk',5))).*bw0,'Area','PixelIdxList');
bwf = bw0.*foci_bw00;
f_index = find(bwf);
for I_f = 1:length(f_index)
    I_temp = nucleus_prop(bwf(f_index(I_f))).PixelIdxList(ceil(rand(1)*nucleus_prop(bwf(f_index(I_f))).Area));
    fake_bw00(I_temp) = true;
    fake_layer(I_temp) = foci_layer(f_index(I_f));
end
fake_data = zeros(sum(sum(fake_bw00)),3+2*length(h));   %%% control data matrix
%%% =======================================================================

bwn = bwlabel(bwi);
nucleus_protein_profile((nucleus_protein_profile(:,1) < Lmin)|(nucleus_protein_profile(:,1) > Lmax),:) = nan(sum((nucleus_protein_profile(:,1) < Lmin)|(nucleus_protein_profile(:,1) > Lmax)),size(nucleus_protein_profile,2));
nucleus_protein_profile = [nan(1,size(nucleus_protein_profile,2));nucleus_protein_profile];

N_fb = 20;
w_fb = 0.05;
N_pb = 20;
w_pb = 0.05;
ppbin = [-1:0.05:1];   %%% binning of the ratio histogram
Npp = length(ppbin);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Calculate local protein concentration: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
foci_bw00 = logical(foci_bw00);
fake_bw00 = logical(fake_bw00);
for ir = 1:length(h)
    f_signal = zeros(size(signal_stack));   %%% Filtered signal
    r_signal = zeros(size(signal_stack));   %%% Filtered RNA signal
    for I_layer = 1:size(signal_stack,3)
        f_signal(:,:,I_layer) = conv2(signal_stack(:,:,I_layer),h{ir},'same').*((foci_layer == I_layer)|(fake_layer == I_layer));
        r_signal(:,:,I_layer) = conv2(RNA_stack(:,:,I_layer),h{ir},'same').*((foci_layer == I_layer)|(fake_layer == I_layer));
    end
    f_signal = f_signal.*mask_stack;
    f_signal0 = max(f_signal,[],3);
    r_signal = r_signal.*mask_stack;
    r_signal0 = max(r_signal,[],3);
    foci_data(:,3+ir) = f_signal0(foci_bw00);
    fake_data(:,3+ir) = f_signal0(fake_bw00);
    foci_data(:,3+length(h)+ir) = r_signal0(foci_bw00);
    fake_data(:,3+length(h)+ir) = r_signal0(fake_bw00);

end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Calculate mean protein intensity at different plains: %%%%%%%%%%%%%%%%%%
mean3_nu = zeros(size(foci_bw00));
for I_layer = 1:size(signal_stack,3)
    nu_prop = regionprops(imerode(bwi, strel('disk',5)).*bwn,signal_stack(:,:,I_layer),'MeanIntensity');
    temp_nu = [0,[nu_prop.MeanIntensity]];
    mean3_nu((foci_layer == I_layer)|(fake_layer == I_layer)) = temp_nu(bwn((foci_layer == I_layer)|(fake_layer == I_layer))+1);
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Collect foci data (foci intensity, foci local protein concentration, mean protein concentration)
foci_data(:,1) = max_image00(foci_bw00).*SS(foci_bw00);
%foci_data(:,2) = f_signal0(foci_bw00);
foci_data(:,2) = nucleus_protein_profile((bwn(foci_bw00)+1),2);
%foci_data(:,3) = r_signal0(foci_bw00);
foci_data(:,3) = mean3_nu(foci_bw00);

fake_data(:,1) = max_image00(fake_bw00).*SS(fake_bw00);
%fake_data(:,2) = f_signal0(fake_bw00);
fake_data(:,2) = nucleus_protein_profile((bwn(fake_bw00)+1),2);
%fake_data(:,3) = r_signal0(fake_bw00);
fake_data(:,3) = mean3_nu(fake_bw00);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Figure plotting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p_ratio0 = foci_data(:,4:(3+length(h)))./repmat(foci_data(:,3),1,length(h))-1;
fp_ratio0 = fake_data(:,4:(3+length(h)))./repmat(fake_data(:,3),1,length(h))-1;
op_ratio0 = foci_data(:,4:(3+length(h)))-repmat(foci_data(:,3),1,length(h));
ofp_ratio0 = fake_data(:,4:(3+length(h)))-repmat(fake_data(:,3),1,length(h));
%color_code = 'rgbcmykrgbcmyk';
color_RGB = [1,0,0;0,1,0;0,0,1;0,1,1;1,0,1;1,1,0;0,0,0;0.5,0.5,0.5];

figure(73)
clf
subplot(1,2,1)
fbin = (max(foci_data(:,1))-min(foci_data(:,1)))/N_fb;
fwin = (max(foci_data(:,1))-min(foci_data(:,1)))*w_fb;
foci_bin = [(min(foci_data(:,1))+fbin/2):fbin:(max(foci_data(:,1))-fbin/2)];
legend_text = cell(0);
for ir = 1:length(h)
    plot(foci_data(:,1),p_ratio0(:,ir),'o','color',color_RGB(ir,:))
    hold on
    ratio_f = zeros(size(foci_bin));
    ratio_ferr = zeros(size(foci_bin));
    for I_bin = 1:length(foci_bin)
        ratio_f(I_bin) = mean(p_ratio0((foci_data(:,1)>=(foci_bin(I_bin)-fwin))&(foci_data(:,1)<=(foci_bin(I_bin)+fwin))&~isnan(p_ratio0(:,ir)),ir));
        ratio_ferr(I_bin) = std0(p_ratio0((foci_data(:,1)>=(foci_bin(I_bin)-fwin))&(foci_data(:,1)<=(foci_bin(I_bin)+fwin))&~isnan(p_ratio0(:,ir)),ir));
    end
    errorbar(foci_bin,ratio_f,ratio_ferr,'color',color_RGB(ir,:))
    legend_text = cat(2,legend_text,{['local intensity, r = ',num2str(r_size(ir))],['binned data, r = ',num2str(r_size(ir))]});
end
title(['Foci local protein enrichment vs foci intensity: ',image_folder,', cycle = ',num2str(N_cycle)])
ylabel('Protein enrichment')
xlabel('Foci intensity (A.U.)')
legend(legend_text)

    
subplot(1,2,2)
pbin = (max(foci_data(:,2))-min(foci_data(:,2)))/N_pb;
pwin = (max(foci_data(:,2))-min(foci_data(:,2)))*w_pb;
prot_bin = [(min(foci_data(:,2))+pbin/2):pbin:(max(foci_data(:,2))-pbin/2)];
legend_text = cell(0);
for ir = 1:length(h)
    plot(foci_data(:,2),p_ratio0(:,ir),'*','color',color_RGB(2*ir-1,:))
    hold on
    plot(fake_data(:,2),fp_ratio0(:,ir),'.','color',color_RGB(2*ir,:))
    ratio_p = zeros(size(prot_bin));
    ratio_perr = zeros(size(prot_bin));
    ratio_fp = zeros(size(prot_bin));
    ratio_fperr = zeros(size(prot_bin));
    for I_bin = 1:length(prot_bin)
        ratio_p(I_bin) = mean(p_ratio0((foci_data(:,2)>=(prot_bin(I_bin)-pwin))&(foci_data(:,2)<=(prot_bin(I_bin)+pwin))&~isnan(p_ratio0(:,ir)),ir));
        ratio_perr(I_bin) = std0(p_ratio0((foci_data(:,2)>=(prot_bin(I_bin)-pwin))&(foci_data(:,2)<=(prot_bin(I_bin)+pwin))&~isnan(p_ratio0(:,ir)),ir));
        ratio_fp(I_bin) = mean(fp_ratio0((fake_data(:,2)>=(prot_bin(I_bin)-pwin))&(fake_data(:,2)<=(prot_bin(I_bin)+pwin))&~isnan(fp_ratio0(:,ir)),ir));
        ratio_fperr(I_bin) = std0(fp_ratio0((fake_data(:,2)>=(prot_bin(I_bin)-pwin))&(fake_data(:,2)<=(prot_bin(I_bin)+pwin))&~isnan(fp_ratio0(:,ir)),ir));
    end
    errorbar(prot_bin,ratio_p,ratio_perr,'color',color_RGB(2*ir-1,:))
    errorbar(prot_bin,ratio_fp,ratio_fperr,'color',color_RGB(2*ir,:))
    plot(prot_bin,ratio_p-ratio_fp,'--','color',mean(color_RGB((2*ir-1):(2*ir),:)))
    legend_text = cat(2,legend_text,{['foci spots, r = ',num2str(r_size(ir))],['control spots, r = ',num2str(r_size(ir))],['binned data, r = ',num2str(r_size(ir))],['binned control data, r = ',num2str(r_size(ir))],['difference of binned data, r = ',num2str(r_size(ir))]});
end
title(['Foci local protein enrichment vs nuclear mean protein intensity: ',image_folder,', cycle = ',num2str(N_cycle)])
ylabel('Protein enrichment')
xlabel('P_m_e_a_n (A.U.)')
legend(legend_text)


figure(74)
clf
for ir = 1:length(h)
    subplot(2,ceil(length(h)/2),ir)
    n_foci = hist(p_ratio0(:,ir),ppbin);
    n_fake = hist(fp_ratio0(:,ir),ppbin);
    plot(ppbin,(n_foci/sum(n_foci)*100),ppbin,(n_fake/sum(n_fake)*100));
    title(['Foci local protein enrichment histogram (r = ',num2str(r_size(ir)),'): ',image_folder,', cycle = ',num2str(N_cycle),', foci mean = ',num2str(mean(p_ratio0((~isnan(p_ratio0(:,ir)))&(p_ratio0(:,ir)>0),ir))),', foci std = ',num2str(std(p_ratio0((~isnan(p_ratio0(:,ir)))&(p_ratio0(:,ir)>0),ir))),', fake foci mean = ',num2str(mean(fp_ratio0((~isnan(fp_ratio0(:,ir)))&(fp_ratio0(:,ir)>0),ir))),', fake foci std = ',num2str(std(fp_ratio0((~isnan(fp_ratio0(:,ir)))&(fp_ratio0(:,ir)>0),ir)))])
    ylabel('%')
    xlabel('Protein enrichment')
    legend('foci spots','control spots')
end


figure(75)
clf
legend_text = cell(0);
for ir = 1:length(h)
    plot(foci_data(:,3+length(h)+ir),op_ratio0(:,ir),'o','color',color_RGB(ir,:))
    hold on
    fbin = (max(foci_data(:,3+length(h)+ir))-min(foci_data(:,3+length(h)+ir)))/N_fb;
    fwin = (max(foci_data(:,3+length(h)+ir))-min(foci_data(:,3+length(h)+ir)))*w_fb;
    foci_bin = [(min(foci_data(:,3+length(h)+ir))+fbin/2):fbin:(max(foci_data(:,3+length(h)+ir))-fbin/2)];
    ratio_f = zeros(size(foci_bin));
    ratio_ferr = zeros(size(foci_bin));
    for I_bin = 1:length(foci_bin)
        ratio_f(I_bin) = mean(op_ratio0((foci_data(:,3+length(h)+ir)>=(foci_bin(I_bin)-fwin))&(foci_data(:,3+length(h)+ir)<=(foci_bin(I_bin)+fwin))&~isnan(op_ratio0(:,ir)),ir));
        ratio_ferr(I_bin) = std0(op_ratio0((foci_data(:,3+length(h)+ir)>=(foci_bin(I_bin)-fwin))&(foci_data(:,3+length(h)+ir)<=(foci_bin(I_bin)+fwin))&~isnan(op_ratio0(:,ir)),ir));
    end
    errorbar(foci_bin,ratio_f,ratio_ferr,'color',color_RGB(ir,:))
    legend_text = cat(2,legend_text,{['foci spots, r = ',num2str(r_size(ir))],['binned data, r = ',num2str(r_size(ir))]});
end
title(['Foci local protein extra fluorescence vs RNA fluorescence: ',image_folder,', cycle = ',num2str(N_cycle)])
ylabel('P_e_x_t_r_a (A.U.)')
xlabel('Foci fluorescence (A.U.)')
legend(legend_text)


figure(76)
clf
subplot(1,2,1)
fbin = (max(foci_data(:,1))-min(foci_data(:,1)))/N_fb;
fwin = (max(foci_data(:,1))-min(foci_data(:,1)))*w_fb;
foci_bin = [(min(foci_data(:,1))+fbin/2):fbin:(max(foci_data(:,1))-fbin/2)];
legend_text = cell(0);
for ir = 1:length(h)
    plot(foci_data(:,1),op_ratio0(:,ir),'o','color',color_RGB(ir,:))
    hold on
    ratio_f = zeros(size(foci_bin));
    ratio_ferr = zeros(size(foci_bin));
    for I_bin = 1:length(foci_bin)
        ratio_f(I_bin) = mean(op_ratio0((foci_data(:,1)>=(foci_bin(I_bin)-fwin))&(foci_data(:,1)<=(foci_bin(I_bin)+fwin))&~isnan(op_ratio0(:,ir)),ir));
        ratio_ferr(I_bin) = std0(op_ratio0((foci_data(:,1)>=(foci_bin(I_bin)-fwin))&(foci_data(:,1)<=(foci_bin(I_bin)+fwin))&~isnan(op_ratio0(:,ir)),ir));
    end
    errorbar(foci_bin,ratio_f,ratio_ferr,'color',color_RGB(ir,:))
    legend_text = cat(2,legend_text,{['local intensity, r = ',num2str(r_size(ir))],['binned data, r = ',num2str(r_size(ir))]});
end
title(['Foci local protein intensity vs foci intensity: ',image_folder,', cycle = ',num2str(N_cycle)])
ylabel('Protein intensity increasement (A.U.)')
xlabel('Foci intensity (A.U.)')
legend(legend_text)

    
subplot(1,2,2)
pbin = (max(foci_data(:,2))-min(foci_data(:,2)))/N_pb;
pwin = (max(foci_data(:,2))-min(foci_data(:,2)))*w_pb;
prot_bin = [(min(foci_data(:,2))+pbin/2):pbin:(max(foci_data(:,2))-pbin/2)];
legend_text = cell(0);
for ir = 1:length(h)
    plot(foci_data(:,2),op_ratio0(:,ir),'*','color',color_RGB(2*ir-1,:))
    hold on
    plot(fake_data(:,2),ofp_ratio0(:,ir),'.','color',color_RGB(2*ir,:))
    ratio_p = zeros(size(prot_bin));
    ratio_perr = zeros(size(prot_bin));
    ratio_fp = zeros(size(prot_bin));
    ratio_fperr = zeros(size(prot_bin));
    for I_bin = 1:length(prot_bin)
        ratio_p(I_bin) = mean(op_ratio0((foci_data(:,2)>=(prot_bin(I_bin)-pwin))&(foci_data(:,2)<=(prot_bin(I_bin)+pwin))&~isnan(op_ratio0(:,ir)),ir));
        ratio_perr(I_bin) = std0(op_ratio0((foci_data(:,2)>=(prot_bin(I_bin)-pwin))&(foci_data(:,2)<=(prot_bin(I_bin)+pwin))&~isnan(op_ratio0(:,ir)),ir));
        ratio_fp(I_bin) = mean(ofp_ratio0((fake_data(:,2)>=(prot_bin(I_bin)-pwin))&(fake_data(:,2)<=(prot_bin(I_bin)+pwin))&~isnan(ofp_ratio0(:,ir)),ir));
        ratio_fperr(I_bin) = std0(ofp_ratio0((fake_data(:,2)>=(prot_bin(I_bin)-pwin))&(fake_data(:,2)<=(prot_bin(I_bin)+pwin))&~isnan(ofp_ratio0(:,ir)),ir));
    end
    errorbar(prot_bin,ratio_p,ratio_perr,'color',color_RGB(2*ir-1,:))
    errorbar(prot_bin,ratio_fp,ratio_fperr,'color',color_RGB(2*ir,:))
    plot(prot_bin,ratio_p-ratio_fp,'--','color',mean(color_RGB((2*ir-1):(2*ir),:)))
    legend_text = cat(2,legend_text,{['foci spots, r = ',num2str(r_size(ir))],['control spots, r = ',num2str(r_size(ir))],['binned data, r = ',num2str(r_size(ir))],['binned control data, r = ',num2str(r_size(ir))],['difference of binned data, r = ',num2str(r_size(ir))]});
end
title(['Foci local protein intensity vs nuclear mean protein intensity: ',image_folder,', cycle = ',num2str(N_cycle)])
ylabel('Protein intensity increasement (A.U.)')
xlabel('P_m_e_a_n (A.U.)')
legend(legend_text)


figure(77)
clf
for ir = 1:length(h)
    subplot(2,ceil(length(h)/2),ir)
    [n_foci,xout_foci] = hist(op_ratio0(:,ir),Npp);
    [n_fake,xout_fake] = hist(ofp_ratio0(:,ir),Npp);
    plot(xout_foci,(n_foci/sum(n_foci)*100),xout_fake,(n_fake/sum(n_fake)*100));
    title(['Foci local protein intensity histogram (r = ',num2str(r_size(ir)),'): ',image_folder,', cycle = ',num2str(N_cycle),', foci mean = ',num2str(mean(op_ratio0((~isnan(op_ratio0(:,ir)))&(op_ratio0(:,ir)>0),ir))),', foci std = ',num2str(std(op_ratio0((~isnan(op_ratio0(:,ir)))&(op_ratio0(:,ir)>0),ir))),', fake foci mean = ',num2str(mean(ofp_ratio0((~isnan(ofp_ratio0(:,ir)))&(ofp_ratio0(:,ir)>0),ir))),', fake foci std = ',num2str(std(ofp_ratio0((~isnan(ofp_ratio0(:,ir)))&(ofp_ratio0(:,ir)>0),ir)))])
    ylabel('%')
    xlabel('Protein intensity increasement')
    legend('foci spots','control spots')
end








    
function y = std0(x)
    NN = size(x);
    if NN(1) == 1;
        N1 = NN(2);
    else
        N1 = NN(1);
    end
    y = std(x)/sqrt(N1);
        
