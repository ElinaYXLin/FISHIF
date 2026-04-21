function [out_data,out_data2,out_data0,out_data20,out_var,out_var2,out_var0,out_var20,out_TX,out_TX2,out_TX0,out_TX20,out_eTX,out_eTX2,out_eTX0,out_eTX20,out_EL,out_EL2,out_EL0,out_EL20,N_foci] = enrichment_compare(nucleus_protein_profile,nucleus_RNA_profile,foci_data,fake_data,h,r_size,Inten_thresh,foci_data2,fake_data2,control_name,image_folder,N_cycle,sub_pos,k_DAPI,N_bin,varargin)

% N_bin = 25;
Lmin = 0;%0.15;
Lmax = 1;%0.85;
Lbin = 0.05;
EL_all = Lmin:Lbin:Lmax;
cover_lim = 0;
r_optimize = 4;
N_grid = 1000;
nM_con = 1;%6.02e8;

if isempty(varargin)
    unit1 = 'A.U.';
    unit2 = 'A.U.';
else
    unit1 = '#';
    unit2 = 'M';
    if length(varargin) >= 2
        r_optimize = varargin{2};
    end
    if length(varargin) >= 3
        z_size = varargin{3};
    end
    if length(varargin) >= 4
        real_name = varargin{4};
    end
    if length(varargin) >= 5
        sigmax = varargin{5}(1);
        sigmay = varargin{5}(2);
        sigmaz = varargin{5}(3);
        sigmax2 = varargin{5}(4);
        sigmay2 = varargin{5}(5);
        sigmaz2 = varargin{5}(6);
    else
        sigmax = 1.35;
        sigmay = 1.35;
        sigmaz = 0.8;
        sigmax2 = 1.8;
        sigmay2 = 1.8;
        sigmaz2 = 1.0;
    end
end

p_ratio0 = cell(size(r_size));
fp_ratio0 = cell(size(r_size));
op_ratio0 = cell(size(r_size));
ofp_ratio0 = cell(size(r_size));

p_ratio02 = cell(size(r_size));
fp_ratio02 = cell(size(r_size));
op_ratio02 = cell(size(r_size));
ofp_ratio02 = cell(size(r_size));

Gm = zeros(size(r_size));
Gm2 = zeros(size(r_size));


for ir = 1:length(r_size)
    ind0 = find(h{ir});
    [xx,yy,zz] = ind2sub(size(h{ir}),ind0);
    G0 = 0;
    G02 = 0;
    for ii = 1:length(xx)
        for jj = 1:length(xx)
            G0 = G0 + exp(-(xx(ii)-xx(jj))^2/4/sigmax^2-(yy(ii)-yy(jj))^2/4/sigmay^2-(zz(ii)-zz(jj))^2/4/sigmaz^2);
            G02 = G02 + exp(-(xx(ii)-xx(jj))^2/4/sigmax2^2-(yy(ii)-yy(jj))^2/4/sigmay2^2-(zz(ii)-zz(jj))^2/4/sigmaz2^2);
        end
    end
    Gm(ir) = G0/length(xx);
    Gm2(ir) = G02/length(xx);
    
    h_area = sum(h{ir}(:));
    foci_data{ir}(:,3) = nM_con*h_area*foci_data{ir}(:,9).*foci_data{ir}(:,12);
    fake_data{ir}(:,3) = nM_con*h_area*fake_data{ir}(:,9).*fake_data{ir}(:,12);
    foci_data2{ir}(:,3) = nM_con*h_area*foci_data2{ir}(:,9).*foci_data2{ir}(:,12);
    fake_data2{ir}(:,3) = nM_con*h_area*fake_data2{ir}(:,9).*fake_data2{ir}(:,12);
    
%     h_area = sum(h{ir}(:));
%     foci_data{ir}(:,3) = nM_con*h_area*foci_data{ir}(:,9).*foci_data{ir}(:,13);
%     fake_data{ir}(:,3) = nM_con*h_area*fake_data{ir}(:,9).*fake_data{ir}(:,13);
%     foci_data2{ir}(:,3) = nM_con*h_area*foci_data2{ir}(:,9).*foci_data2{ir}(:,13);
%     fake_data2{ir}(:,3) = nM_con*h_area*fake_data2{ir}(:,9).*fake_data2{ir}(:,13);
% %     
% %     foci_data{ir}(:,15) = 1./sqrt(2*pi)./sigmaz;
% %     fake_data{ir}(:,15) = 1./sqrt(2*pi)./sigmaz;
% %     foci_data2{ir}(:,15) = 1./sqrt(2*pi)./sigmaz;
% %     fake_data2{ir}(:,15) = 1./sqrt(2*pi)./sigmaz;
%     
% %     foci_data{ir}(:,[3,4]) = foci_data{ir}(:,[3,4])./foci_data{ir}(:,[15,15]);
% %     fake_data{ir}(:,[3,4]) = fake_data{ir}(:,[3,4])./fake_data{ir}(:,[15,15]);
% %     foci_data2{ir}(:,[3,4]) = foci_data2{ir}(:,[3,4])./foci_data2{ir}(:,[15,15]);
% %     fake_data2{ir}(:,[3,4]) = fake_data2{ir}(:,[3,4])./fake_data2{ir}(:,[15,15]);
     

    idx = (foci_data{ir}(:,7)>=Lmin) & (foci_data{ir}(:,7)<=Lmax) & (~isnan(foci_data{ir}(:,4))) & (~isnan(foci_data{ir}(:,3))) & (foci_data{ir}(:,8) > 1) & (foci_data{ir}(:,8) < z_size) & (foci_data{ir}(:,9) >= cover_lim) & (foci_data{ir}(:,1) >= Inten_thresh);
    idxf = (foci_data{ir}(:,7)>=Lmin) & (foci_data{ir}(:,7)<=Lmax) & (~isnan(fake_data{ir}(:,4))) & (~isnan(fake_data{ir}(:,3))) & (fake_data{ir}(:,8) > 1) & (fake_data{ir}(:,8) < z_size) & (foci_data{ir}(:,9) >= cover_lim) & (foci_data{ir}(:,1) >= Inten_thresh);
    idx2 = (foci_data2{ir}(:,7)>=Lmin) & (foci_data2{ir}(:,7)<=Lmax) & (~isnan(foci_data2{ir}(:,4))) & (~isnan(foci_data2{ir}(:,3))) & (foci_data2{ir}(:,8) > 1) & (foci_data2{ir}(:,8) < z_size) & (foci_data2{ir}(:,9) >= cover_lim);
    idx2f = (foci_data2{ir}(:,7)>=Lmin) & (foci_data2{ir}(:,7)<=Lmax) & (~isnan(fake_data2{ir}(:,4))) & (~isnan(fake_data2{ir}(:,3))) & (fake_data2{ir}(:,8) > 1) & (fake_data2{ir}(:,8) < z_size) & (foci_data2{ir}(:,9) >= cover_lim);
    idx = idx & idxf;
    idxf = idx;
    idx2 = idx2 & idx2f;
    idx2f = idx2;
    
    foci_data{ir} = foci_data{ir}(idx,:);
    fake_data{ir} = fake_data{ir}(idxf,:);
    foci_data2{ir} = foci_data2{ir}(idx2,:);
    fake_data2{ir} = fake_data2{ir}(idx2f,:);
    
    %%% align hb and nullo data:
    ind10 = foci_data{ir}(:,6);
    ind20 = foci_data2{ir}(:,6);
    ind1 = false(size(ind10));
    ind2 = false(size(ind20));
    ind22 = zeros(0);
    
    for Id = 1:length(ind10)
        ind_temp = find(ind20.*(~ind2) == ind10(Id),1);
        if ~isempty(ind_temp)
            ind1(Id) = true;
            ind2(ind_temp) = true;
            ind22 = cat(2,ind22,ind_temp);
        end
    end
    foci_data0{ir} = foci_data{ir};
    fake_data0{ir} = fake_data{ir};
    foci_data20{ir} = foci_data2{ir};
    fake_data20{ir} = fake_data2{ir};
    
    foci_data{ir} = foci_data{ir}(ind1,:);
    fake_data{ir} = fake_data{ir}(ind1,:);
    foci_data2{ir} = foci_data2{ir}(ind22,:);
    fake_data2{ir} = fake_data2{ir}(ind22,:);
    
    
    
%     p_ratio0{ir} = foci_data{ir}(:,4)./foci_data{ir}(:,3)-1;
%     fp_ratio0{ir} = fake_data{ir}(:,4)./fake_data{ir}(:,3)-1;
%     op_ratio0{ir} = foci_data{ir}(:,4)-foci_data{ir}(:,3);
    op_ratio00{ir} = foci_data0{ir}(:,4)-foci_data0{ir}(:,3);
    ofp_ratio00{ir} = fake_data0{ir}(:,4)-fake_data0{ir}(:,3);
    
    op_ratio0{ir} = foci_data{ir}(:,4)-foci_data{ir}(:,3);
    ofp_ratio0{ir} = fake_data{ir}(:,4)-fake_data{ir}(:,3);
    
%     p_ratio02{ir} = foci_data2{ir}(:,4)./foci_data2{ir}(:,3)-1;
%     fp_ratio02{ir} = fake_data2{ir}(:,4)./fake_data2{ir}(:,3)-1;
    op_ratio020{ir} = foci_data20{ir}(:,4)-foci_data20{ir}(:,3);
    ofp_ratio020{ir} = fake_data20{ir}(:,4)-fake_data20{ir}(:,3);

    op_ratio02{ir} = foci_data2{ir}(:,4)-foci_data2{ir}(:,3);
    ofp_ratio02{ir} = fake_data2{ir}(:,4)-fake_data2{ir}(:,3);
end

%% Figure plotting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% color_RGB0 = [1,0,0;0,1,0;0,0,1;0,1,1;1,0,1;1,1,0;0.5,0.5,0.5];
% % color_RGB = [color_RGB0;[0,0,0];color_RGB0/2;[0.5,0.5,0.5]];
% color_RGB = [color_RGB0;color_RGB0;color_RGB0;color_RGB0;color_RGB0];
color_RGB = mycolors(7);
color_RGB = [color_RGB;color_RGB;color_RGB];
% if length(r_size) > size(color_RGB0,1)
%     LL2 = ceil(length(r_size)/size(color_RGB0,1))*2;
%     color_RGB  = repmat(color_RGB,LL2,1);
% end    


figure(576)
clf
subplot(1,2,1)
legend_text = cell(0);

%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ir = 1:length(r_size)
%     var_coef = Gm(ir)/4/sqrt(2)/pi/sigmax/sigmay;
    var_coef = Gm(ir)/8/(pi)^1.5/sigmax/sigmay/sigmaz;
    var_coef2 = Gm2(ir)/8/(pi)^1.5/sigmax2/sigmay2/sigmaz2;
    plot(foci_data{ir}(:,2),op_ratio0{ir},'*','color',color_RGB(1,:))
    hold on
    plot(fake_data{ir}(:,2),ofp_ratio0{ir},'.','color',color_RGB(2,:))

    [xtemp1,IX1] = sort(foci_data{ir}(:,2));
    [xtemp2,IX2] = sort(fake_data{ir}(:,2));
    ytemp1 = op_ratio0{ir}(IX1);
    ytemp2 = ofp_ratio0{ir}(IX2);
    ytemp101 = foci_data{ir}(IX1,4);
    ytemp102 = foci_data{ir}(IX1,3)*var_coef;
    ytemp201 = fake_data{ir}(IX2,4);
    ytemp202 = fake_data{ir}(IX2,3)*var_coef;
    ytemp3 = foci_data{ir}(IX1,4)-fake_data{ir}(IX2,4);
    ztemp1 = foci_data{ir}(IX1,7);
%     ytemp3 = ytemp1-ytemp2;
%     ytemp30 = ytemp3;
    temp_ind = (~isnan(ytemp1)) & (~isnan(ytemp2));
    xtemp1 = xtemp1(temp_ind);
    xtemp2 = xtemp2(temp_ind);
    ytemp1 = ytemp1(temp_ind);
    ytemp2 = ytemp2(temp_ind);
    ytemp101 = ytemp101(temp_ind);
    ytemp102 = ytemp102(temp_ind);
    ytemp201 = ytemp201(temp_ind);
    ytemp202 = ytemp202(temp_ind);
    ytemp3 = ytemp3(temp_ind);
    ztemp1 = ztemp1(temp_ind);
    N_temp = floor(length(xtemp1)/N_bin);
    x1_f = zeros(1,N_bin);
    x2_f = zeros(1,N_bin);
    y1_f = zeros(1,N_bin);
    y2_f = zeros(1,N_bin);
    y3_f = zeros(1,N_bin);
    x1_ferr = zeros(1,N_bin);
    x2_ferr = zeros(1,N_bin);
    y1_ferr = zeros(1,N_bin);
    y2_ferr = zeros(1,N_bin);
    y3_ferr = zeros(1,N_bin);
    y1_fvar = zeros(1,N_bin);
    y2_fvar = zeros(1,N_bin);
    y3_fvar = zeros(1,N_bin);
    y1_fverr = zeros(1,N_bin);
    y2_fverr = zeros(1,N_bin);
    y3_fverr = zeros(1,N_bin);
    
    for I_bin = 1:N_bin
        x1_f(I_bin) = mean(xtemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
        x1_ferr(I_bin) = std0(xtemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
        x2_f(I_bin) = mean(xtemp2(((I_bin-1)*N_temp+1):I_bin*N_temp));
        x2_ferr(I_bin) = std0(xtemp2(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y1_f(I_bin) = mean(ytemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y1_ferr(I_bin) = std0(ytemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y1_fvar(I_bin) = var(ytemp101(((I_bin-1)*N_temp+1):I_bin*N_temp))-mean(ytemp102(((I_bin-1)*N_temp+1):I_bin*N_temp));
%         y1_fvar(I_bin) = var(ytemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y1_fverr(I_bin) = y1_fvar(I_bin)/sqrt(N_temp/2);
        y2_f(I_bin) = mean(ytemp2(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y2_ferr(I_bin) = std0(ytemp2(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y2_fvar(I_bin) = var(ytemp201(((I_bin-1)*N_temp+1):I_bin*N_temp))-mean(ytemp202(((I_bin-1)*N_temp+1):I_bin*N_temp));
%         y2_fvar(I_bin) = var(ytemp2(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y2_fverr(I_bin) = y2_fvar(I_bin)/sqrt(N_temp/2);
        y3_f(I_bin) = mean(ytemp3(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y3_ferr(I_bin) = std0(ytemp3(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y3_fvar(I_bin) = var(ytemp101(((I_bin-1)*N_temp+1):I_bin*N_temp))-var(ytemp201(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y3_fverr(I_bin) = y3_fvar(I_bin)/sqrt(N_temp/2);
%         y3_fvar(I_bin) = y1_fvar(I_bin)-y2_fvar(I_bin);
%         y3_fverr(I_bin) = sqrt((y1_fverr(I_bin))^2+(y2_fverr(I_bin))^2);
    end
    errorbarxy(x1_f,y1_f,x1_ferr,y1_ferr,x1_ferr,y1_ferr,color_RGB(1,:),color_RGB(1,:),'o')
    errorbarxy(x2_f,y2_f,x2_ferr,y2_ferr,x2_ferr,y2_ferr,color_RGB(2,:),color_RGB(2,:),'*')
    errorbarxy(x1_f,y3_f,x1_ferr,y3_ferr,x1_ferr,y3_ferr,color_RGB(3,:),color_RGB(3,:),'--')
%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%     if ~isempty(y1_f)
        legend_text = cat(2,legend_text,{[real_name,' - background, r = ',num2str(r_size(ir))],['fake ',real_name,' - background, r = ',num2str(r_size(ir))],['(',real_name,' - background) bin, r = ',num2str(r_size(ir))],['(fake ',real_name,' - background) bin, r = ',num2str(r_size(ir))],['(',real_name,' - fake) bin, r = ',num2str(r_size(ir))]});
%     else
%         legend_text = cat(2,legend_text,{[real_name,'-nullo, r = ',num2str(r_size(ir))],['fake ',real_name,' - background, r = ',num2str(r_size(ir))],['(',real_name,'-nullo) bin, r = ',num2str(r_size(ir))],['(fake ',real_name,' - background) bin, r = ',num2str(r_size(ir))]});
%     end
    if r_size(ir) == r_optimize
        if sub_pos(3) > 0
            out_data2 = [xtemp1,ytemp1,xtemp2,ytemp2,xtemp1,ytemp3];
        end
    end

%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    plot(foci_data2{ir}(:,2),op_ratio02{ir},'*','color',color_RGB(4,:))
    hold on
    plot(fake_data2{ir}(:,2),ofp_ratio02{ir},'.','color',color_RGB(5,:))

    [xtemp1,IX1] = sort(foci_data2{ir}(:,2));
    [xtemp2,IX2] = sort(fake_data2{ir}(:,2));
    ytemp1 = op_ratio02{ir}(IX1);
    ytemp2 = ofp_ratio02{ir}(IX2);
    ytemp111 = foci_data2{ir}(IX1,4);
    ytemp112 = foci_data2{ir}(IX1,3)*var_coef2;
    ytemp211 = fake_data2{ir}(IX2,4);
    ytemp212 = fake_data2{ir}(IX2,3)*var_coef2;
    ytemp3 = foci_data2{ir}(IX1,4)-fake_data2{ir}(IX2,4);
    ytemp4 = foci_data{ir}(IX1,4)-foci_data2{ir}(IX1,4);
%     ytemp3 = ytemp1-ytemp2;
%     ytemp4 = op_ratio0{ir}(IX1)-op_ratio02{ir}(IX1);
%     ytemp4 = ytemp3-ytemp30;
    temp_ind = (~isnan(ytemp1)) & (~isnan(ytemp2)) & (~isnan(ytemp4));
    xtemp1 = xtemp1(temp_ind);
    xtemp2 = xtemp2(temp_ind);
    ytemp1 = ytemp1(temp_ind);
    ytemp2 = ytemp2(temp_ind);
    ytemp3 = ytemp3(temp_ind);
    ytemp4 = ytemp4(temp_ind);
    ytemp111 = ytemp111(temp_ind);
    ytemp112 = ytemp112(temp_ind);
    ytemp211 = ytemp211(temp_ind);
    ytemp212 = ytemp212(temp_ind);
    N_temp = floor(length(xtemp1)/N_bin);
    x21_f = zeros(1,N_bin);
    x22_f = zeros(1,N_bin);
    y21_f = zeros(1,N_bin);
    y22_f = zeros(1,N_bin);
    y23_f = zeros(1,N_bin);
    y24_f = zeros(1,N_bin);
    x21_ferr = zeros(1,N_bin);
    x22_ferr = zeros(1,N_bin);
    y21_ferr = zeros(1,N_bin);
    y22_ferr = zeros(1,N_bin);
    y23_ferr = zeros(1,N_bin);
    y24_ferr = zeros(1,N_bin);
    y21_fvar = zeros(1,N_bin);
    y22_fvar = zeros(1,N_bin);
    y23_fvar = zeros(1,N_bin);
    y24_fvar = zeros(1,N_bin);
    y21_fverr = zeros(1,N_bin);
    y22_fverr = zeros(1,N_bin);
    y23_fverr = zeros(1,N_bin);
    y24_fverr = zeros(1,N_bin);
    
    for I_bin = 1:N_bin
        x21_f(I_bin) = mean(xtemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
        x21_ferr(I_bin) = std0(xtemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
        x22_f(I_bin) = mean(xtemp2(((I_bin-1)*N_temp+1):I_bin*N_temp));
        x22_ferr(I_bin) = std0(xtemp2(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y21_f(I_bin) = mean(ytemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y21_ferr(I_bin) = std0(ytemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y21_fvar(I_bin) = var(ytemp111(((I_bin-1)*N_temp+1):I_bin*N_temp))-mean(ytemp112(((I_bin-1)*N_temp+1):I_bin*N_temp));
%         y21_fvar(I_bin) = var(ytemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y21_fverr(I_bin) = y21_fvar(I_bin)/sqrt(N_temp/2);
        y22_f(I_bin) = mean(ytemp2(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y22_ferr(I_bin) = std0(ytemp2(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y22_fvar(I_bin) = var(ytemp211(((I_bin-1)*N_temp+1):I_bin*N_temp))-mean(ytemp212(((I_bin-1)*N_temp+1):I_bin*N_temp));
%         y22_fvar(I_bin) = var(ytemp2(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y22_fverr(I_bin) = y22_fvar(I_bin)/sqrt(N_temp/2);
        y23_f(I_bin) = mean(ytemp3(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y23_ferr(I_bin) = std0(ytemp3(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y23_fvar(I_bin) = var(ytemp111(((I_bin-1)*N_temp+1):I_bin*N_temp))-var(ytemp211(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y23_fverr(I_bin) = y23_fvar(I_bin)/sqrt(N_temp/2);
%         y23_fvar(I_bin) = y21_fvar(I_bin)-y22_fvar(I_bin);
%         y23_fverr(I_bin) = sqrt((y21_fverr(I_bin))^2+(y22_fverr(I_bin))^2);
        y24_f(I_bin) = mean(ytemp4(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y24_ferr(I_bin) = std0(ytemp4(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y24_fvar(I_bin) = var(ytemp101(((I_bin-1)*N_temp+1):I_bin*N_temp))-var(ytemp111(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y24_fverr(I_bin) = y24_fvar(I_bin)/sqrt(N_temp/2);
%         y24_fvar(I_bin) = y1_fvar(I_bin)-y21_fvar(I_bin);
%         y24_fverr(I_bin) = sqrt((y1_fverr(I_bin))^2+(y21_fverr(I_bin))^2);
    end
    errorbarxy(x21_f,y21_f,x21_ferr,y21_ferr,x21_ferr,y21_ferr,color_RGB(4,:),color_RGB(4,:),'o')
    errorbarxy(x22_f,y22_f,x22_ferr,y22_ferr,x22_ferr,y22_ferr,color_RGB(5,:),color_RGB(5,:),'*')
    errorbarxy(x21_f,y23_f,x21_ferr,y23_ferr,x21_ferr,y23_ferr,color_RGB(6,:),color_RGB(6,:),'--')
    errorbarxy(x21_f,y24_f,x21_ferr,y24_ferr,x21_ferr,y24_ferr,color_RGB(7,:),color_RGB(7,:),'+')
%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
%     if ~isempty(y21_f)
        legend_text = cat(2,legend_text,{[control_name,' - background, r = ',num2str(r_size(ir))],[control_name,' fake - background, r = ',num2str(r_size(ir))],['(',control_name,' - background) bin, r = ',num2str(r_size(ir))],['(',control_name,' fake - background) bin, r = ',num2str(r_size(ir))],['(',control_name,' - fake) bin, r = ',num2str(r_size(ir))],['(',real_name,' - ',control_name,') bin, r = ',num2str(r_size(ir))]});
%     else
%         legend_text = cat(2,legend_text,{[control_name,' foci spots, r = ',num2str(r_size(ir))],[control_name,' control spots, r = ',num2str(r_size(ir))],[control_name,' binned data, r = ',num2str(r_size(ir))],[control_name,' binned control data, r = ',num2str(r_size(ir))],['Real gene vs ',control_name,' difference of binned data, r = ',num2str(r_size(ir))]});
%     end
    
    
    if r_size(ir) == r_optimize
        if sub_pos(3) > 0
            
            out_data = [x1_f',y1_f',x2_f',y2_f',x1_f',y3_f',x21_f',y21_f',x22_f',y22_f',x21_f',y23_f',x21_f',y24_f'];
            out_var  = [x1_f',y1_fvar',x2_f',y2_fvar',x1_f',y3_fvar',x21_f',y21_fvar',x22_f',y22_fvar',x21_f',y23_fvar',x21_f',y24_fvar'];
            out_data2 = cat(2,out_data2,[xtemp1,ytemp1,xtemp2,ytemp2,xtemp1,ytemp3,xtemp1,ytemp4,ztemp1]);
            out_var2 = [xtemp1,ytemp101,ytemp102,xtemp2,ytemp201,ytemp202,xtemp1,ytemp101,ytemp201,xtemp1,ytemp111,ytemp112,xtemp2,ytemp211,ytemp212,xtemp1,ytemp111,ytemp211,xtemp1,ytemp101,ytemp111,ztemp1];

            N_foci = size(foci_data{ir},1);
%             out_data1 = [x21_f',y24_f'];
%             out_data2 = [x1_f',y3_f'];
            
            figure(5076)
            subaxis(sub_pos(1),sub_pos(2),sub_pos(3), 'Spacing', 0, 'PaddingRight',0.01, 'PaddingLeft',0.02, 'PaddingTop',0.045, 'PaddingBottom',0.04, 'Margin', 0.01);
            hold on
            errorbarxy(x1_f,y1_f,x1_ferr,y1_ferr,x1_ferr,y1_ferr,color_RGB(1,:),color_RGB(1,:),'o')
            errorbarxy(x2_f,y2_f,x2_ferr,y2_ferr,x2_ferr,y2_ferr,color_RGB(2,:),color_RGB(2,:),'*')
            errorbarxy(x1_f,y3_f,x1_ferr,y3_ferr,x1_ferr,y3_ferr,color_RGB(3,:),color_RGB(3,:),'--')
            errorbarxy(x21_f,y21_f,x21_ferr,y21_ferr,x21_ferr,y21_ferr,color_RGB(4,:),color_RGB(4,:),'o')
            errorbarxy(x22_f,y22_f,x22_ferr,y22_ferr,x22_ferr,y22_ferr,color_RGB(5,:),color_RGB(5,:),'*')
            errorbarxy(x21_f,y23_f,x21_ferr,y23_ferr,x21_ferr,y23_ferr,color_RGB(6,:),color_RGB(6,:),'--')
            errorbarxy(x21_f,y24_f,x21_ferr,y24_ferr,x21_ferr,y24_ferr,color_RGB(7,:),color_RGB(7,:),'+')
            
            p_fit = polyfit(x21_f,y24_f,1);
            xlim0 = xlim;
            x1_fit = [xlim0(1):(xlim0(2)-xlim0(1))/N_grid:xlim0(2)];
            y1_fit = p_fit(1)*x1_fit+p_fit(2);
            plot(x1_fit,y1_fit,'Color',color_RGB(7,:));
                
%             title([image_folder,char(10),'Foci enrichment vs protein concentration, cycle = ',num2str(N_cycle),char(10),', mean TX = ',num2str(mean(nucleus_RNA_profile(:,4))),', DAPI slope = ',num2str(k_DAPI)],'Interpreter','none')
%             ylabel(['Protein enrichment (',unit1,')'])
%             xlabel(['P_m_e_a_n (',unit2,')'])
% 
%             if ~isempty(y1_f)
                legend_text0 = {['(',real_name,' - background) bin, r = ',num2str(r_size(ir))],['(fake ',real_name,' - background) bin, r = ',num2str(r_size(ir))],['(',real_name,' - fake) bin, r = ',num2str(r_size(ir))]};
%             else
%                 legend_text0 = {['foci spots, r = ',num2str(r_size(ir))],['control spots, r = ',num2str(r_size(ir))]};
%             end
%             if ~isempty(y21_f)
                legend_text0 = cat(2,legend_text0,{['(',control_name,' - background) bin, r = ',num2str(r_size(ir))],['(',control_name,' fake - background) bin, r = ',num2str(r_size(ir))],['(',control_name,' - fake) bin, r = ',num2str(r_size(ir))],['(',real_name,' - ',control_name,') bin, r = ',num2str(r_size(ir))]});
%             else
%                 legend_text0 = cat(2,legend_text0,{[control_name,' foci spots, r = ',num2str(r_size(ir))],[control_name,' control spots, r = ',num2str(r_size(ir))],['Real gene vs ',control_name,' difference of binned data, r = ',num2str(r_size(ir))]});
%             end
            
            legend_text0 = cat(2,legend_text0,{['(',real_name,' - ',control_name,') fit, r = ',num2str(r_size(ir)),': y = ',num2str(p_fit(1)),' * x + ',num2str(p_fit(2))]});
            
%             legend(legend_text0)
%             grid on

            figure(50762)
            subaxis(sub_pos(1),sub_pos(2),sub_pos(3), 'Spacing', 0, 'PaddingRight',0.01, 'PaddingLeft',0.02, 'PaddingTop',0.045, 'PaddingBottom',0.04, 'Margin', 0.01);
            hold on
            errorbarxy(x1_f,y1_fvar,x1_ferr,y1_fverr,x1_ferr,y1_fverr,color_RGB(1,:),color_RGB(1,:),'o')
            errorbarxy(x2_f,y2_fvar,x2_ferr,y2_fverr,x2_ferr,y2_fverr,color_RGB(2,:),color_RGB(2,:),'*')
            errorbarxy(x1_f,y3_fvar,x1_ferr,y3_fverr,x1_ferr,y3_fverr,color_RGB(3,:),color_RGB(3,:),'--')
            errorbarxy(x21_f,y21_fvar,x21_ferr,y21_fverr,x21_ferr,y21_fverr,color_RGB(4,:),color_RGB(4,:),'o')
            errorbarxy(x22_f,y22_fvar,x22_ferr,y22_fverr,x22_ferr,y22_fverr,color_RGB(5,:),color_RGB(5,:),'*')
            errorbarxy(x21_f,y23_fvar,x21_ferr,y23_fverr,x21_ferr,y23_fverr,color_RGB(6,:),color_RGB(6,:),'--')
            errorbarxy(x21_f,y24_fvar,x21_ferr,y24_fverr,x21_ferr,y24_fverr,color_RGB(7,:),color_RGB(7,:),'+')
            
            legend_text0v = {['(',real_name,' - background) bin, r = ',num2str(r_size(ir))],['(fake ',real_name,' - background) bin, r = ',num2str(r_size(ir))],['(',real_name,' - fake) bin, r = ',num2str(r_size(ir))]};
            legend_text0v = cat(2,legend_text0v,{['(',control_name,' - background) bin, r = ',num2str(r_size(ir))],['(',control_name,' fake - background) bin, r = ',num2str(r_size(ir))],['(',control_name,' - fake) bin, r = ',num2str(r_size(ir))],['(',real_name,' - ',control_name,') bin, r = ',num2str(r_size(ir))]});


            figure(5078)
            subaxis(sub_pos(1),sub_pos(2),sub_pos(3), 'Spacing', 0, 'PaddingRight',0.01, 'PaddingLeft',0.02, 'PaddingTop',0.045, 'PaddingBottom',0.04, 'Margin', 0.01);
            hold on
            errorbarxy(x1_f,y1_f./x1_f,x1_ferr,(y1_ferr./y1_f+x1_ferr./x1_f).*y1_f./x1_f,x1_ferr,(y1_ferr./y1_f+x1_ferr./x1_f).*y1_f./x1_f,color_RGB(1,:),color_RGB(1,:),'o')
            errorbarxy(x2_f,y2_f./x2_f,x2_ferr,(y2_ferr./y2_f+x2_ferr./x2_f).*y2_f./x2_f,x2_ferr,(y2_ferr./y2_f+x2_ferr./x2_f).*y2_f./x2_f,color_RGB(2,:),color_RGB(2,:),'*')
            errorbarxy(x1_f,y3_f./x1_f,x1_ferr,(y3_ferr./y3_f+x1_ferr./x1_f).*y3_f./x1_f,x1_ferr,(y3_ferr./y3_f+x1_ferr./x1_f).*y3_f./x1_f,color_RGB(3,:),color_RGB(3,:),'--')
%             errorbarxy(x1_f,y1_f-y2_f,x1_ferr,sqrt(y1_ferr.^2+y2_ferr.^2),x1_ferr,sqrt(y1_ferr.^2+y2_ferr.^2),color_RGB(3,:),color_RGB(3,:),'--')
            errorbarxy(x21_f,y21_f./x21_f,x21_ferr,(y21_ferr./y21_f+x21_ferr./x21_f).*y21_f./x21_f,x21_ferr,(y21_ferr./y21_f+x21_ferr./x21_f).*y21_f./x21_f,color_RGB(4,:),color_RGB(4,:),'o')
            errorbarxy(x22_f,y22_f./x22_f,x22_ferr,(y22_ferr./y22_f+x22_ferr./x22_f).*y22_f./x22_f,x22_ferr,(y22_ferr./y22_f+x22_ferr./x22_f).*y22_f./x22_f,color_RGB(5,:),color_RGB(5,:),'*')
            errorbarxy(x21_f,y23_f./x21_f,x21_ferr,(y23_ferr./y23_f+x21_ferr./x21_f).*y23_f./x21_f,x21_ferr,(y23_ferr./y23_f+x21_ferr./x21_f).*y23_f./x21_f,color_RGB(6,:),color_RGB(6,:),'--')
            errorbarxy(x21_f,y24_f./x21_f,x21_ferr,(y24_ferr./y24_f+x21_ferr./x21_f).*y24_f./x21_f,x21_ferr,(y24_ferr./y24_f+x21_ferr./x21_f).*y24_f./x21_f,color_RGB(7,:),color_RGB(7,:),'+')
%             errorbarxy(x21_f,y21_f-y22_f,x21_ferr,sqrt(y21_ferr.^2+y22_ferr.^2),x21_ferr,sqrt(y21_ferr.^2+y22_ferr.^2),color_RGB(6,:),color_RGB(6,:),'--')
                
%             title([image_folder,char(10),'Foci enrichment/protein concentration vs protein concentration, cycle = ',num2str(N_cycle),char(10),', mean TX = ',num2str(mean(nucleus_RNA_profile(:,4))),', DAPI slope = ',num2str(k_DAPI)],'Interpreter','none')
%             ylabel(['Protein enrichment / P_m_e_a_n (',unit1,'/',unit2,')'])
%             xlabel(['P_m_e_a_n (',unit2,')'])

%             if ~isempty(y1_f)
%                 legend_text0 = cat(2,legend_text0,{['foci spots0, r = ',num2str(r_size(ir))],['control spots0, r = ',num2str(r_size(ir))],['difference of binned data0, r = ',num2str(r_size(ir))]});
%             else
                legend_text01 = {['(',real_name,' - background) bin, r = ',num2str(r_size(ir))],['(fake ',real_name,' - background) bin, r = ',num2str(r_size(ir))],['(',real_name,' - fake) bin, r = ',num2str(r_size(ir))]};
%             end
%             if ~isempty(y21_f)
%                 legend_text0 = cat(2,legend_text0,{[control_name,' foci spots0, r = ',num2str(r_size(ir))],[control_name,' control spots0, r = ',num2str(r_size(ir))],[control_name,' difference of binned data0, r = ',num2str(r_size(ir))]});
%             else
                legend_text01 = cat(2,legend_text01,{['(',control_name,' - background) bin, r = ',num2str(r_size(ir))],['(',control_name,' fake - background) bin, r = ',num2str(r_size(ir))],['(',control_name,' - fake) bin, r = ',num2str(r_size(ir))],['(',real_name,' - ',control_name,') bin, r = ',num2str(r_size(ir))]});
%             end
            

            
            
            
            figure(576)
        end
    end
    

    
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(5761)
    clf
    legend_text_5761 = cell(0);

    plot(op_ratio0{ir},foci_data{ir}(:,1),'*','color',color_RGB(1,:))
    hold on
    plot(foci_data{ir}(:,4)-fake_data{ir}(:,4),foci_data{ir}(:,1),'.','color',color_RGB(3,:))
    plot(op_ratio02{ir},foci_data2{ir}(:,1),'*','color',color_RGB(4,:))
    plot(foci_data2{ir}(:,4)-fake_data2{ir}(:,4),foci_data2{ir}(:,1),'.','color',color_RGB(6,:))
    plot(foci_data{ir}(:,4)-foci_data2{ir}(:,4),foci_data{ir}(:,1),'.','color',color_RGB(7,:))

%     [ytemp1,IX1] = sort(foci_data{ir}(:,1));
    [ztemp1,IX1] = sort(foci_data{ir}(:,2));
    ztemp2 = foci_data{ir}(IX1,7);
    ytemp1 = foci_data{ir}(IX1,1);
    xtemp1 = op_ratio0{ir}(IX1);
    xtemp2 = foci_data{ir}(IX1,4)-fake_data{ir}(IX1,4);
    xtemp3 = foci_data{ir}(IX1,4)-foci_data2{ir}(IX1,4);
    ytemp2 = foci_data2{ir}(IX1,1);
    xtemp21 = op_ratio02{ir}(IX1);
    xtemp22 = foci_data2{ir}(IX1,4)-fake_data2{ir}(IX1,4);
    
    temp_ind = (~isnan(xtemp1)) & (~isnan(xtemp2)) & (~isnan(xtemp3)) & (~isnan(ytemp2)) & (~isnan(xtemp21)) & (~isnan(xtemp22));
    xtemp1 = xtemp1(temp_ind);
    xtemp2 = xtemp2(temp_ind);
    xtemp3 = xtemp3(temp_ind);
    xtemp21 = xtemp21(temp_ind);
    xtemp22 = xtemp22(temp_ind);
    ytemp1 = ytemp1(temp_ind);
    ytemp2 = ytemp2(temp_ind);
    ztemp1 = ztemp1(temp_ind);
    ztemp2 = ztemp2(temp_ind);
    N_temp = floor(length(xtemp1)/N_bin);
    x1_f = zeros(1,N_bin);
    x2_f = zeros(1,N_bin);
    x3_f = zeros(1,N_bin);
    x21_f = zeros(1,N_bin);
    x22_f = zeros(1,N_bin);
    y1_f = zeros(1,N_bin);
    y2_f = zeros(1,N_bin);
    z1_f = zeros(1,N_bin);

    x1_ferr = zeros(1,N_bin);
    x2_ferr = zeros(1,N_bin);
    x3_ferr = zeros(1,N_bin);
    x21_ferr = zeros(1,N_bin);
    x22_ferr = zeros(1,N_bin);
    y1_ferr = zeros(1,N_bin);
    y2_ferr = zeros(1,N_bin);
    z1_ferr = zeros(1,N_bin);
    
    for I_bin = 1:N_bin
        x1_f(I_bin) = mean(xtemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
        x1_ferr(I_bin) = std0(xtemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
        x2_f(I_bin) = mean(xtemp2(((I_bin-1)*N_temp+1):I_bin*N_temp));
        x2_ferr(I_bin) = std0(xtemp2(((I_bin-1)*N_temp+1):I_bin*N_temp));
        x3_f(I_bin) = mean(xtemp3(((I_bin-1)*N_temp+1):I_bin*N_temp));
        x3_ferr(I_bin) = std0(xtemp3(((I_bin-1)*N_temp+1):I_bin*N_temp));
        x21_f(I_bin) = mean(xtemp21(((I_bin-1)*N_temp+1):I_bin*N_temp));
        x21_ferr(I_bin) = std0(xtemp21(((I_bin-1)*N_temp+1):I_bin*N_temp));
        x22_f(I_bin) = mean(xtemp22(((I_bin-1)*N_temp+1):I_bin*N_temp));
        x22_ferr(I_bin) = std0(xtemp22(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y1_f(I_bin) = mean(ytemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y1_ferr(I_bin) = std0(ytemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y2_f(I_bin) = mean(ytemp2(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y2_ferr(I_bin) = std0(ytemp2(((I_bin-1)*N_temp+1):I_bin*N_temp));
        z1_f(I_bin) = mean(ztemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
        z1_ferr(I_bin) = std0(ztemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
    end
    errorbarxy(x1_f,y1_f,x1_ferr,y1_ferr,x1_ferr,y1_ferr,color_RGB(1,:),color_RGB(1,:),'o')
    errorbarxy(x2_f,y1_f,x2_ferr,y1_ferr,x2_ferr,y1_ferr,color_RGB(3,:),color_RGB(3,:),'*')
    errorbarxy(x21_f,y2_f,x21_ferr,y2_ferr,x21_ferr,y2_ferr,color_RGB(4,:),color_RGB(4,:),'o')
    errorbarxy(x22_f,y2_f,x22_ferr,y2_ferr,x22_ferr,y2_ferr,color_RGB(6,:),color_RGB(6,:),'*')
    errorbarxy(x3_f,y1_f,x3_ferr,y1_ferr,x3_ferr,y1_ferr,color_RGB(7,:),color_RGB(7,:),'+')
%    errorbarxy(x1_f,y1_f-y2_f,x1_ferr,sqrt(y1_ferr.^2+y2_ferr.^2),x1_ferr,sqrt(y1_ferr.^2+y2_ferr.^2),color_RGB(3,:),color_RGB(3,:),'--')

    ex1_f = zeros(size(EL_all));
    ex2_f = zeros(size(EL_all));
    ex3_f = zeros(size(EL_all));
    ex21_f = zeros(size(EL_all));
    ex22_f = zeros(size(EL_all));
    ey1_f = zeros(size(EL_all));
    ey2_f = zeros(size(EL_all));
    ez2_f = EL_all;
    ex1_ferr = zeros(size(EL_all));
    ex2_ferr = zeros(size(EL_all));
    ex3_ferr = zeros(size(EL_all));
    ex21_ferr = zeros(size(EL_all));
    ex22_ferr = zeros(size(EL_all));
    ey1_ferr = zeros(size(EL_all));
    ey2_ferr = zeros(size(EL_all));
    ez2_ferr = zeros(size(EL_all));
    
    for I_bin = 1:length(EL_all)
        ex1_f(I_bin) = mean(xtemp1((ztemp2 >= EL_all(I_bin)-Lbin) & (ztemp2 <= EL_all(I_bin)+Lbin)));
        ex1_ferr(I_bin) = std0(xtemp1((ztemp2 >= EL_all(I_bin)-Lbin) & (ztemp2 <= EL_all(I_bin)+Lbin)));
        ex2_f(I_bin) = mean(xtemp2((ztemp2 >= EL_all(I_bin)-Lbin) & (ztemp2 <= EL_all(I_bin)+Lbin)));
        ex2_ferr(I_bin) = std0(xtemp2((ztemp2 >= EL_all(I_bin)-Lbin) & (ztemp2 <= EL_all(I_bin)+Lbin)));
        ex3_f(I_bin) = mean(xtemp3((ztemp2 >= EL_all(I_bin)-Lbin) & (ztemp2 <= EL_all(I_bin)+Lbin)));
        ex3_ferr(I_bin) = std0(xtemp3((ztemp2 >= EL_all(I_bin)-Lbin) & (ztemp2 <= EL_all(I_bin)+Lbin)));
        ex21_f(I_bin) = mean(xtemp21((ztemp2 >= EL_all(I_bin)-Lbin) & (ztemp2 <= EL_all(I_bin)+Lbin)));
        ex21_ferr(I_bin) = std0(xtemp21((ztemp2 >= EL_all(I_bin)-Lbin) & (ztemp2 <= EL_all(I_bin)+Lbin)));
        ex22_f(I_bin) = mean(xtemp22((ztemp2 >= EL_all(I_bin)-Lbin) & (ztemp2 <= EL_all(I_bin)+Lbin)));
        ex22_ferr(I_bin) = std0(xtemp22((ztemp2 >= EL_all(I_bin)-Lbin) & (ztemp2 <= EL_all(I_bin)+Lbin)));
        ey1_f(I_bin) = mean(ytemp1((ztemp2 >= EL_all(I_bin)-Lbin) & (ztemp2 <= EL_all(I_bin)+Lbin)));
        ey1_ferr(I_bin) = std0(ytemp1((ztemp2 >= EL_all(I_bin)-Lbin) & (ztemp2 <= EL_all(I_bin)+Lbin)));
        ey2_f(I_bin) = mean(ytemp2((ztemp2 >= EL_all(I_bin)-Lbin) & (ztemp2 <= EL_all(I_bin)+Lbin)));
        ey2_ferr(I_bin) = std0(ytemp2((ztemp2 >= EL_all(I_bin)-Lbin) & (ztemp2 <= EL_all(I_bin)+Lbin)));
    end
    errorbarxy(ex1_f,ey1_f,ex1_ferr,ey1_ferr,ex1_ferr,ey1_ferr,color_RGB(1,:),color_RGB(1,:),'o')
    errorbarxy(ex2_f,ey1_f,ex2_ferr,ey1_ferr,ex2_ferr,ey1_ferr,color_RGB(3,:),color_RGB(3,:),'*')
    errorbarxy(ex21_f,ey2_f,ex21_ferr,ey2_ferr,ex21_ferr,ey2_ferr,color_RGB(4,:),color_RGB(4,:),'o')
    errorbarxy(ex22_f,ey2_f,ex22_ferr,ey2_ferr,ex22_ferr,ey2_ferr,color_RGB(6,:),color_RGB(6,:),'*')
    errorbarxy(ex3_f,ey1_f,ex3_ferr,ey1_ferr,ex3_ferr,ey1_ferr,color_RGB(7,:),color_RGB(7,:),'+')
%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    legend_text_5761 = cat(2,legend_text_5761,{[real_name,' - background, r = ',num2str(r_size(ir))],[real_name,' - fake, r = ',num2str(r_size(ir))],[control_name,' - background, r = ',num2str(r_size(ir))],[control_name,' - fake, r = ',num2str(r_size(ir))],[real_name,' - ',control_name,', r = ',num2str(r_size(ir))],['(',real_name,' - background) bin, r = ',num2str(r_size(ir))],['(',real_name,' - fake) bin, r = ',num2str(r_size(ir))],['(',control_name,' - background) bin, r = ',num2str(r_size(ir))],['(',control_name,' - fake) bin, r = ',num2str(r_size(ir))],['(',real_name,' - ',control_name,') bin, r = ',num2str(r_size(ir))]});
    legend_text_5761 = cat(2,legend_text_5761,{['(',real_name,' - background) bin (EL), r = ',num2str(r_size(ir))],['(',real_name,' - fake) bin (EL), r = ',num2str(r_size(ir))],['(',control_name,' - background) bin (EL), r = ',num2str(r_size(ir))],['(',control_name,' - fake) bin (EL), r = ',num2str(r_size(ir))],['(',real_name,' - ',control_name,') bin (EL), r = ',num2str(r_size(ir))]});
    

    if r_size(ir) == r_optimize
        if sub_pos(3) > 0
            figure(50761)
            subaxis(sub_pos(1),sub_pos(2),sub_pos(3), 'Spacing', 0, 'PaddingRight',0.01, 'PaddingLeft',0.02, 'PaddingTop',0.045, 'PaddingBottom',0.04, 'Margin', 0.01);
            hold on
            errorbarxy(x1_f,y1_f,x1_ferr,y1_ferr,x1_ferr,y1_ferr,color_RGB(1,:),color_RGB(1,:),'o')
            errorbarxy(x2_f,y1_f,x2_ferr,y1_ferr,x2_ferr,y1_ferr,color_RGB(3,:),color_RGB(3,:),'*')
            errorbarxy(x21_f,y2_f,x21_ferr,y2_ferr,x21_ferr,y2_ferr,color_RGB(4,:),color_RGB(4,:),'o')
            errorbarxy(x22_f,y2_f,x22_ferr,y2_ferr,x22_ferr,y2_ferr,color_RGB(6,:),color_RGB(6,:),'*')
            errorbarxy(x3_f,y1_f,x3_ferr,y1_ferr,x3_ferr,y1_ferr,color_RGB(7,:),color_RGB(7,:),'+')
            
            errorbarxy(ex1_f,ey1_f,ex1_ferr,ey1_ferr,ex1_ferr,ey1_ferr,color_RGB(1,:),color_RGB(1,:),'o')
            errorbarxy(ex2_f,ey1_f,ex2_ferr,ey1_ferr,ex2_ferr,ey1_ferr,color_RGB(3,:),color_RGB(3,:),'*')
            errorbarxy(ex21_f,ey2_f,ex21_ferr,ey2_ferr,ex21_ferr,ey2_ferr,color_RGB(4,:),color_RGB(4,:),'o')
            errorbarxy(ex22_f,ey2_f,ex22_ferr,ey2_ferr,ex22_ferr,ey2_ferr,color_RGB(6,:),color_RGB(6,:),'*')
            errorbarxy(ex3_f,ey1_f,ex3_ferr,ey1_ferr,ex3_ferr,ey1_ferr,color_RGB(7,:),color_RGB(7,:),'+')

            legend_text_57610 = {['(',real_name,' - background) bin, r = ',num2str(r_size(ir))],['(',real_name,' - fake) bin, r = ',num2str(r_size(ir))],['(',control_name,' - background) bin, r = ',num2str(r_size(ir))],['(',control_name,' - fake) bin, r = ',num2str(r_size(ir))],['(',real_name,' - ',control_name,') bin, r = ',num2str(r_size(ir))]};
            legend_text_57610 = cat(2,legend_text_57610,{['(',real_name,' - background) bin (EL), r = ',num2str(r_size(ir))],['(',real_name,' - fake) bin (EL), r = ',num2str(r_size(ir))],['(',control_name,' - background) bin (EL), r = ',num2str(r_size(ir))],['(',control_name,' - fake) bin (EL), r = ',num2str(r_size(ir))],['(',real_name,' - ',control_name,') bin (EL), r = ',num2str(r_size(ir))]});
            
            out_TX = [x1_f',y1_f',z1_f',x2_f',y1_f',z1_f',x21_f',y2_f',z1_f',x22_f',y2_f',z1_f',x3_f',y1_f',z1_f'];
            out_TX2 = [xtemp1,ytemp1,ztemp1,xtemp2,ytemp1,ztemp1,xtemp21,ytemp2,ztemp1,xtemp22,ytemp2,ztemp1,xtemp3,ytemp1,ztemp1];
            out_eTX = [ex1_f',ey1_f',ez2_f',ex2_f',ey1_f',ez2_f',ex21_f',ey2_f',ez2_f',ex22_f',ey2_f',ez2_f',ex3_f',ey1_f',ez2_f'];
            out_eTX2 = [xtemp1,ytemp1,ztemp2,xtemp2,ytemp1,ztemp2,xtemp21,ytemp2,ztemp2,xtemp22,ytemp2,ztemp2,xtemp3,ytemp1,ztemp2];
            
            figure(576)
        end
    end


    
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    plot(foci_data0{ir}(:,2),op_ratio00{ir},'*','color',color_RGB(1,:))
    hold on
    plot(fake_data0{ir}(:,2),ofp_ratio00{ir},'.','color',color_RGB(2,:))

    [xtemp1,IX1] = sort(foci_data0{ir}(:,2));
    [xtemp2,IX2] = sort(fake_data0{ir}(:,2));
    ytemp1 = op_ratio00{ir}(IX1);
    ytemp2 = ofp_ratio00{ir}(IX2);
    ytemp101 = foci_data0{ir}(IX1,4);
    ytemp102 = foci_data0{ir}(IX1,3)*var_coef;
    ytemp201 = fake_data0{ir}(IX2,4);
    ytemp202 = fake_data0{ir}(IX2,3)*var_coef;
    ytemp3 = foci_data0{ir}(IX1,4)-fake_data0{ir}(IX2,4);
    ztemp1 = foci_data0{ir}(IX1,7);
%     ytemp3 = ytemp1-ytemp2;
    temp_ind = (~isnan(ytemp1)) & (~isnan(ytemp2));
    xtemp1 = xtemp1(temp_ind);
    xtemp2 = xtemp2(temp_ind);
    ytemp1 = ytemp1(temp_ind);
    ytemp2 = ytemp2(temp_ind);
    ytemp3 = ytemp3(temp_ind);
    ytemp101 = ytemp101(temp_ind);
    ytemp102 = ytemp102(temp_ind);
    ytemp201 = ytemp201(temp_ind);
    ytemp202 = ytemp202(temp_ind);
    ztemp1 = ztemp1(temp_ind);
    N_temp = floor(length(xtemp1)/N_bin);
    x1_f = zeros(1,N_bin);
    x2_f = zeros(1,N_bin);
    y1_f = zeros(1,N_bin);
    y2_f = zeros(1,N_bin);
    y3_f = zeros(1,N_bin);
    x1_ferr = zeros(1,N_bin);
    x2_ferr = zeros(1,N_bin);
    y1_ferr = zeros(1,N_bin);
    y2_ferr = zeros(1,N_bin);
    y3_ferr = zeros(1,N_bin);
    y1_fvar = zeros(1,N_bin);
    y2_fvar = zeros(1,N_bin);
    y3_fvar = zeros(1,N_bin);
    y1_fverr = zeros(1,N_bin);
    y2_fverr = zeros(1,N_bin);
    y3_fverr = zeros(1,N_bin);
    
    for I_bin = 1:N_bin
        x1_f(I_bin) = mean(xtemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
        x1_ferr(I_bin) = std0(xtemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
        x2_f(I_bin) = mean(xtemp2(((I_bin-1)*N_temp+1):I_bin*N_temp));
        x2_ferr(I_bin) = std0(xtemp2(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y1_f(I_bin) = mean(ytemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y1_ferr(I_bin) = std0(ytemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y1_fvar(I_bin) = var(ytemp101(((I_bin-1)*N_temp+1):I_bin*N_temp))-mean(ytemp102(((I_bin-1)*N_temp+1):I_bin*N_temp));
%         y1_fvar(I_bin) = var(ytemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y1_fverr(I_bin) = y1_fvar(I_bin)/sqrt(N_temp/2);
        y2_f(I_bin) = mean(ytemp2(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y2_ferr(I_bin) = std0(ytemp2(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y2_fvar(I_bin) = var(ytemp201(((I_bin-1)*N_temp+1):I_bin*N_temp))-mean(ytemp202(((I_bin-1)*N_temp+1):I_bin*N_temp));
%         y2_fvar(I_bin) = var(ytemp2(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y2_fverr(I_bin) = y2_fvar(I_bin)/sqrt(N_temp/2);
        y3_f(I_bin) = mean(ytemp3(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y3_ferr(I_bin) = std0(ytemp3(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y3_fvar(I_bin) = var(ytemp101(((I_bin-1)*N_temp+1):I_bin*N_temp))-var(ytemp201(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y3_fverr(I_bin) = y3_fvar(I_bin)/sqrt(N_temp/2);
%         y3_fvar(I_bin) = y1_fvar(I_bin)-y2_fvar(I_bin);
%         y3_fverr(I_bin) = sqrt((y1_fverr(I_bin))^2+(y2_fverr(I_bin))^2);
    end
    errorbarxy(x1_f,y1_f,x1_ferr,y1_ferr,x1_ferr,y1_ferr,color_RGB(1,:),color_RGB(1,:),'o')
    errorbarxy(x2_f,y2_f,x2_ferr,y2_ferr,x2_ferr,y2_ferr,color_RGB(2,:),color_RGB(2,:),'*')
    errorbarxy(x1_f,y3_f,x1_ferr,y3_ferr,x1_ferr,y3_ferr,color_RGB(3,:),color_RGB(3,:),'--')
%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%     if ~isempty(y1_f)
        legend_text = cat(2,legend_text,{[real_name,' - background all, r = ',num2str(r_size(ir))],['fake ',real_name,' - background all, r = ',num2str(r_size(ir))],['(',real_name,' - background) all bin, r = ',num2str(r_size(ir))],['(fake ',real_name,' - background) all bin, r = ',num2str(r_size(ir))],['(',real_name,' - fake) all bin, r = ',num2str(r_size(ir))]});
%     else
%         legend_text = cat(2,legend_text,{[real_name,'-nullo, r = ',num2str(r_size(ir))],['fake ',real_name,' - background, r = ',num2str(r_size(ir))],['(',real_name,'-nullo) bin, r = ',num2str(r_size(ir))],['(fake ',real_name,' - background) bin, r = ',num2str(r_size(ir))]});
%     end
    if r_size(ir) == r_optimize
        if sub_pos(3) > 0
            out_data20 = {[xtemp1,ytemp1,xtemp2,ytemp2,xtemp1,ytemp3,ztemp1]};
            out_var20 = {[xtemp1,ytemp101,ytemp102,xtemp2,ytemp201,ytemp202,xtemp1,ytemp101,ytemp201,ztemp1]};
        end
    end
    

%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    plot(foci_data20{ir}(:,2),op_ratio020{ir},'*','color',color_RGB(4,:))
    hold on
    plot(fake_data20{ir}(:,2),ofp_ratio020{ir},'.','color',color_RGB(5,:))

    [xtemp1,IX1] = sort(foci_data20{ir}(:,2));
    [xtemp2,IX2] = sort(fake_data20{ir}(:,2));
    ytemp1 = op_ratio020{ir}(IX1);
    ytemp2 = ofp_ratio020{ir}(IX2);
    ytemp111 = foci_data20{ir}(IX1,4);
    ytemp112 = foci_data20{ir}(IX1,3)*var_coef2;
    ytemp211 = fake_data20{ir}(IX2,4);
    ytemp212 = fake_data20{ir}(IX2,3)*var_coef2;
    ytemp3 = foci_data20{ir}(IX1,4)-fake_data20{ir}(IX1,4);
    ztemp1 = foci_data20{ir}(IX1,7);
%     ytemp3 = ytemp1-ytemp2;
    temp_ind = (~isnan(ytemp1)) & (~isnan(ytemp2)) & (~isnan(ytemp3));
    xtemp1 = xtemp1(temp_ind);
    xtemp2 = xtemp2(temp_ind);
    ytemp1 = ytemp1(temp_ind);
    ytemp2 = ytemp2(temp_ind);
    ytemp3 = ytemp3(temp_ind);
    ytemp111 = ytemp111(temp_ind);
    ytemp112 = ytemp112(temp_ind);
    ytemp211 = ytemp211(temp_ind);
    ytemp212 = ytemp212(temp_ind);
    ztemp1 = ztemp1(temp_ind);
    N_temp = floor(length(xtemp1)/N_bin);
    x21_f = zeros(1,N_bin);
    x22_f = zeros(1,N_bin);
    y21_f = zeros(1,N_bin);
    y22_f = zeros(1,N_bin);
    y23_f = zeros(1,N_bin);
    x21_ferr = zeros(1,N_bin);
    x22_ferr = zeros(1,N_bin);
    y21_ferr = zeros(1,N_bin);
    y22_ferr = zeros(1,N_bin);
    y23_ferr = zeros(1,N_bin);
    y21_fvar = zeros(1,N_bin);
    y22_fvar = zeros(1,N_bin);
    y23_fvar = zeros(1,N_bin);
    y21_fverr = zeros(1,N_bin);
    y22_fverr = zeros(1,N_bin);
    y23_fverr = zeros(1,N_bin);
    
    for I_bin = 1:N_bin
        x21_f(I_bin) = mean(xtemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
        x21_ferr(I_bin) = std0(xtemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
        x22_f(I_bin) = mean(xtemp2(((I_bin-1)*N_temp+1):I_bin*N_temp));
        x22_ferr(I_bin) = std0(xtemp2(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y21_f(I_bin) = mean(ytemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y21_ferr(I_bin) = std0(ytemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y21_fvar(I_bin) = var(ytemp111(((I_bin-1)*N_temp+1):I_bin*N_temp))-mean(ytemp112(((I_bin-1)*N_temp+1):I_bin*N_temp));
%         y21_fvar(I_bin) = var(ytemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y21_fverr(I_bin) = y21_fvar(I_bin)/sqrt(N_temp/2);
        y22_f(I_bin) = mean(ytemp2(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y22_ferr(I_bin) = std0(ytemp2(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y22_fvar(I_bin) = var(ytemp211(((I_bin-1)*N_temp+1):I_bin*N_temp))-mean(ytemp212(((I_bin-1)*N_temp+1):I_bin*N_temp));
%         y22_fvar(I_bin) = var(ytemp2(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y22_fverr(I_bin) = y22_fvar(I_bin)/sqrt(N_temp/2);
        y23_f(I_bin) = mean(ytemp3(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y23_ferr(I_bin) = std0(ytemp3(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y23_fvar(I_bin) = var(ytemp111(((I_bin-1)*N_temp+1):I_bin*N_temp))-var(ytemp211(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y23_fverr(I_bin) = y23_fvar(I_bin)/sqrt(N_temp/2);
%         y23_fvar(I_bin) = y21_fvar(I_bin)-y22_fvar(I_bin);
%         y23_fverr(I_bin) = sqrt((y21_fverr(I_bin))^2+(y22_fverr(I_bin))^2);
    end
    errorbarxy(x21_f,y21_f,x21_ferr,y21_ferr,x21_ferr,y21_ferr,color_RGB(4,:),color_RGB(4,:),'o')
    errorbarxy(x22_f,y22_f,x22_ferr,y22_ferr,x22_ferr,y22_ferr,color_RGB(5,:),color_RGB(5,:),'*')
    errorbarxy(x21_f,y23_f,x21_ferr,y23_ferr,x21_ferr,y23_ferr,color_RGB(6,:),color_RGB(6,:),'--')
%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
%     if ~isempty(y21_f)
        legend_text = cat(2,legend_text,{[control_name,' - background all, r = ',num2str(r_size(ir))],[control_name,' fake - background all, r = ',num2str(r_size(ir))],['(',control_name,' - background) all bin, r = ',num2str(r_size(ir))],['(',control_name,' fake - background) all bin, r = ',num2str(r_size(ir))],['(',control_name,' - fake) all bin, r = ',num2str(r_size(ir))]});
%     else
%         legend_text = cat(2,legend_text,{[control_name,' foci spots, r = ',num2str(r_size(ir))],[control_name,' control spots, r = ',num2str(r_size(ir))],[control_name,' binned data, r = ',num2str(r_size(ir))],[control_name,' binned control data, r = ',num2str(r_size(ir))],['Real gene vs ',control_name,' difference of binned data, r = ',num2str(r_size(ir))]});
%     end
    
    
    if r_size(ir) == r_optimize
        if sub_pos(3) > 0
            
            out_data0 = [x1_f',y1_f',x2_f',y2_f',x1_f',y3_f',x21_f',y21_f',x22_f',y22_f',x21_f',y23_f'];
            out_var0 = [x1_f',y1_fvar',x2_f',y2_fvar',x1_f',y3_fvar',x21_f',y21_fvar',x22_f',y22_fvar',x21_f',y23_fvar'];
            out_data20 = cat(2,out_data20,{[xtemp1,ytemp1,xtemp2,ytemp2,xtemp1,ytemp3,ztemp1]});
            out_var20 = cat(2,out_var20,{[xtemp1,ytemp111,ytemp112,xtemp2,ytemp211,ytemp212,xtemp1,ytemp111,ytemp211,ztemp1]});

            figure(5076)
            subaxis(sub_pos(1),sub_pos(2),sub_pos(3), 'Spacing', 0, 'PaddingRight',0.01, 'PaddingLeft',0.02, 'PaddingTop',0.045, 'PaddingBottom',0.04, 'Margin', 0.01);
            hold on
            errorbarxy(x1_f,y1_f,x1_ferr,y1_ferr,x1_ferr,y1_ferr,color_RGB(1,:),color_RGB(1,:),'o')
            errorbarxy(x2_f,y2_f,x2_ferr,y2_ferr,x2_ferr,y2_ferr,color_RGB(2,:),color_RGB(2,:),'*')
            errorbarxy(x1_f,y3_f,x1_ferr,y3_ferr,x1_ferr,y3_ferr,color_RGB(3,:),color_RGB(3,:),'--')
            errorbarxy(x21_f,y21_f,x21_ferr,y21_ferr,x21_ferr,y21_ferr,color_RGB(4,:),color_RGB(4,:),'o')
            errorbarxy(x22_f,y22_f,x22_ferr,y22_ferr,x22_ferr,y22_ferr,color_RGB(5,:),color_RGB(5,:),'*')
            errorbarxy(x21_f,y23_f,x21_ferr,y23_ferr,x21_ferr,y23_ferr,color_RGB(6,:),color_RGB(6,:),'--')
                
            title([image_folder,char(10),'Foci enrichment vs protein concentration, cycle = ',num2str(N_cycle),char(10),', mean TX = ',num2str(mean(nucleus_RNA_profile(:,4))),', DAPI slope = ',num2str(k_DAPI),', pair #: ',num2str(size(foci_data{ir},1)),', foci #: ',num2str(size(foci_data0{ir},1)),', ',control_name,' foci #: ',num2str(size(foci_data20{ir},1)),char(10),'foci spots fit, r = ',num2str(r_size(ir)),': y = ',num2str(p_fit(1)),' * x + ',num2str(p_fit(2))],'Interpreter','none')
            ylabel(['Protein enrichment (',unit1,')'])
            xlabel(['P_m_e_a_n (',unit2,')'])

%             if ~isempty(y1_f)
                legend_text0 = cat(2,legend_text0,{['(',real_name,' - background) all bin, r = ',num2str(r_size(ir))],['(fake ',real_name,' - background) all bin, r = ',num2str(r_size(ir))],['(',real_name,' - fake) all bin, r = ',num2str(r_size(ir))]});
%             else
%                 legend_text0 = {['foci spots, r = ',num2str(r_size(ir))],['control spots, r = ',num2str(r_size(ir))]};
%             end
%             if ~isempty(y21_f)
                legend_text0 = cat(2,legend_text0,{['(',control_name,' - background) all bin, r = ',num2str(r_size(ir))],['(',control_name,' fake - background) all bin, r = ',num2str(r_size(ir))],['(',control_name,' - fake) all bin, r = ',num2str(r_size(ir))]});
%             else
%                 legend_text0 = cat(2,legend_text0,{[control_name,' foci spots, r = ',num2str(r_size(ir))],[control_name,' control spots, r = ',num2str(r_size(ir))],['Real gene vs ',control_name,' difference of binned data, r = ',num2str(r_size(ir))]});
%             end
            
            legend(legend_text0)
            legend('hide')
            grid on


            figure(50762)
            subaxis(sub_pos(1),sub_pos(2),sub_pos(3), 'Spacing', 0, 'PaddingRight',0.01, 'PaddingLeft',0.02, 'PaddingTop',0.045, 'PaddingBottom',0.04, 'Margin', 0.01);
            hold on
            errorbarxy(x1_f,y1_fvar,x1_ferr,y1_fverr,x1_ferr,y1_fverr,color_RGB(1,:),color_RGB(1,:),'o')
            errorbarxy(x2_f,y2_fvar,x2_ferr,y2_fverr,x2_ferr,y2_fverr,color_RGB(2,:),color_RGB(2,:),'*')
            errorbarxy(x1_f,y3_fvar,x1_ferr,y3_fverr,x1_ferr,y3_fverr,color_RGB(3,:),color_RGB(3,:),'--')
            errorbarxy(x21_f,y21_fvar,x21_ferr,y21_fverr,x21_ferr,y21_fverr,color_RGB(4,:),color_RGB(4,:),'o')
            errorbarxy(x22_f,y22_fvar,x22_ferr,y22_fverr,x22_ferr,y22_fverr,color_RGB(5,:),color_RGB(5,:),'*')
            errorbarxy(x21_f,y23_fvar,x21_ferr,y23_fverr,x21_ferr,y23_fverr,color_RGB(6,:),color_RGB(6,:),'--')
                
            title([image_folder,char(10),'Foci enrichment variance vs protein concentration, cycle = ',num2str(N_cycle),char(10),', mean TX = ',num2str(mean(nucleus_RNA_profile(:,4))),', DAPI slope = ',num2str(k_DAPI),', pair #: ',num2str(size(foci_data{ir},1)),', foci #: ',num2str(size(foci_data0{ir},1)),', ',control_name,' foci #: ',num2str(size(foci_data20{ir},1)),char(10),'foci spots fit, r = ',num2str(r_size(ir)),': y = ',num2str(p_fit(1)),' * x + ',num2str(p_fit(2))],'Interpreter','none')
            ylabel(['Protein enrichment variance (',unit1,')'])
            xlabel(['P_m_e_a_n (',unit2,')'])

            legend_text0v = cat(2,legend_text0v,{['(',real_name,' - background) all bin, r = ',num2str(r_size(ir))],['(fake ',real_name,' - background) all bin, r = ',num2str(r_size(ir))],['(',real_name,' - fake) all bin, r = ',num2str(r_size(ir))]});
            legend_text0v = cat(2,legend_text0v,{['(',control_name,' - background) all bin, r = ',num2str(r_size(ir))],['(',control_name,' fake - background) all bin, r = ',num2str(r_size(ir))],['(',control_name,' - fake) all bin, r = ',num2str(r_size(ir))]});
           
            legend(legend_text0v)
            legend('hide')
            grid on

            
            
            figure(5078)
            subaxis(sub_pos(1),sub_pos(2),sub_pos(3), 'Spacing', 0, 'PaddingRight',0.01, 'PaddingLeft',0.02, 'PaddingTop',0.045, 'PaddingBottom',0.04, 'Margin', 0.01);
            hold on
            errorbarxy(x1_f,y1_f./x1_f,x1_ferr,(y1_ferr./y1_f+x1_ferr./x1_f).*y1_f./x1_f,x1_ferr,(y1_ferr./y1_f+x1_ferr./x1_f).*y1_f./x1_f,color_RGB(1,:),color_RGB(1,:),'o')
            errorbarxy(x2_f,y2_f./x2_f,x2_ferr,(y2_ferr./y2_f+x2_ferr./x2_f).*y2_f./x2_f,x2_ferr,(y2_ferr./y2_f+x2_ferr./x2_f).*y2_f./x2_f,color_RGB(2,:),color_RGB(2,:),'*')
            errorbarxy(x1_f,y3_f./x1_f,x1_ferr,(y3_ferr./y3_f+x1_ferr./x1_f).*y3_f./x1_f,x1_ferr,(y3_ferr./y3_f+x1_ferr./x1_f).*y3_f./x1_f,color_RGB(3,:),color_RGB(3,:),'--')
%             errorbarxy(x1_f,y1_f-y2_f,x1_ferr,sqrt(y1_ferr.^2+y2_ferr.^2),x1_ferr,sqrt(y1_ferr.^2+y2_ferr.^2),color_RGB(3,:),color_RGB(3,:),'--')
            errorbarxy(x21_f,y21_f./x21_f,x21_ferr,(y21_ferr./y21_f+x21_ferr./x21_f).*y21_f./x21_f,x21_ferr,(y21_ferr./y21_f+x21_ferr./x21_f).*y21_f./x21_f,color_RGB(4,:),color_RGB(4,:),'o')
            errorbarxy(x22_f,y22_f./x22_f,x22_ferr,(y22_ferr./y22_f+x22_ferr./x22_f).*y22_f./x22_f,x22_ferr,(y22_ferr./y22_f+x22_ferr./x22_f).*y22_f./x22_f,color_RGB(5,:),color_RGB(5,:),'*')
            errorbarxy(x21_f,y23_f./x21_f,x21_ferr,(y23_ferr./y23_f+x21_ferr./x21_f).*y23_f./x21_f,x21_ferr,(y23_ferr./y23_f+x21_ferr./x21_f).*y23_f./x21_f,color_RGB(6,:),color_RGB(6,:),'--')
                
            title([image_folder,char(10),'Foci enrichment/protein concentration vs protein concentration, cycle = ',num2str(N_cycle),char(10),', mean TX = ',num2str(mean(nucleus_RNA_profile(:,4))),', DAPI slope = ',num2str(k_DAPI),', pair #: ',num2str(size(foci_data{ir},1)),', foci #: ',num2str(size(foci_data0{ir},1)),', ',control_name,' foci #: ',num2str(size(foci_data20{ir},1))],'Interpreter','none')
            ylabel(['Protein enrichment / P_m_e_a_n (',unit1,'/',unit2,')'])
            xlabel(['P_m_e_a_n (',unit2,')'])

%             if ~isempty(y1_f)
%                 legend_text0 = cat(2,legend_text0,{['foci spots0, r = ',num2str(r_size(ir))],['control spots0, r = ',num2str(r_size(ir))],['difference of binned data0, r = ',num2str(r_size(ir))]});
%             else
                legend_text01 = {['(',real_name,' - background) all bin, r = ',num2str(r_size(ir))],['(fake ',real_name,' - background) all bin, r = ',num2str(r_size(ir))],['(',real_name,' - fake) all bin, r = ',num2str(r_size(ir))]};
%             end
%             if ~isempty(y21_f)
%                 legend_text0 = cat(2,legend_text0,{[control_name,' foci spots0, r = ',num2str(r_size(ir))],[control_name,' control spots0, r = ',num2str(r_size(ir))],[control_name,' difference of binned data0, r = ',num2str(r_size(ir))]});
%             else
                legend_text01 = cat(2,legend_text01,{['(',control_name,' - background) all bin, r = ',num2str(r_size(ir))],['(',control_name,' fake - background) all bin, r = ',num2str(r_size(ir))],['(',control_name,' - fake) all bin, r = ',num2str(r_size(ir))]});
%             end
            
            legend(legend_text01)
            legend('hide')
            grid on
            
            
            figure(576)
        end
    end
    
    
    

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(5761)

    plot(op_ratio00{ir},foci_data0{ir}(:,1),'*','color',color_RGB(1,:))
    hold on
    plot(foci_data0{ir}(:,4)-fake_data0{ir}(:,4),foci_data0{ir}(:,1),'.','color',color_RGB(3,:))
    plot(op_ratio020{ir},foci_data20{ir}(:,1),'*','color',color_RGB(4,:))
    plot(foci_data20{ir}(:,4)-fake_data20{ir}(:,4),foci_data20{ir}(:,1),'.','color',color_RGB(6,:))

%     [ytemp1,IX1] = sort(foci_data0{ir}(:,1));
    [ztemp1,IX1] = sort(foci_data0{ir}(:,2));
    ztemp12 = foci_data0{ir}(IX1,7);
    ytemp1 = foci_data0{ir}(IX1,1);
    xtemp1 = op_ratio00{ir}(IX1);
    xtemp2 = foci_data0{ir}(IX1,4)-fake_data0{ir}(IX1,4);
    [ztemp2,IX2] = sort(foci_data20{ir}(:,2));
    ztemp22 = foci_data20{ir}(IX2,7);
    ytemp2 = foci_data20{ir}(IX2,1);
    xtemp21 = op_ratio020{ir}(IX2);
    xtemp22 = foci_data20{ir}(IX2,4)-fake_data20{ir}(IX2,4);

    temp_ind = (~isnan(xtemp1)) & (~isnan(xtemp2));
    xtemp1 = xtemp1(temp_ind);
    xtemp2 = xtemp2(temp_ind);
    ytemp1 = ytemp1(temp_ind);
    ztemp1 = ztemp1(temp_ind);
    ztemp12 = ztemp12(temp_ind);
    temp_ind = (~isnan(xtemp21)) & (~isnan(xtemp22));
    xtemp21 = xtemp21(temp_ind);
    xtemp22 = xtemp22(temp_ind);
    ytemp2 = ytemp2(temp_ind);
    ztemp2 = ztemp2(temp_ind);
    ztemp22 = ztemp22(temp_ind);
    
    x1_f = zeros(1,N_bin);
    x2_f = zeros(1,N_bin);
    y1_f = zeros(1,N_bin);
    z1_f = zeros(1,N_bin);
    x21_f = zeros(1,N_bin);
    x22_f = zeros(1,N_bin);
    y2_f = zeros(1,N_bin);
    z2_f = zeros(1,N_bin);
    x1_ferr = zeros(1,N_bin);
    x2_ferr = zeros(1,N_bin);
    y1_ferr = zeros(1,N_bin);
    z1_ferr = zeros(1,N_bin);
    x21_ferr = zeros(1,N_bin);
    x22_ferr = zeros(1,N_bin);
    y2_ferr = zeros(1,N_bin);
    z2_ferr = zeros(1,N_bin);
    
    N_temp = floor(length(xtemp1)/N_bin);
    for I_bin = 1:N_bin
        x1_f(I_bin) = mean(xtemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
        x1_ferr(I_bin) = std0(xtemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
        x2_f(I_bin) = mean(xtemp2(((I_bin-1)*N_temp+1):I_bin*N_temp));
        x2_ferr(I_bin) = std0(xtemp2(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y1_f(I_bin) = mean(ytemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y1_ferr(I_bin) = std0(ytemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
        z1_f(I_bin) = mean(ztemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
        z1_ferr(I_bin) = std0(ztemp1(((I_bin-1)*N_temp+1):I_bin*N_temp));
    end
    
    N_temp = floor(length(xtemp21)/N_bin);
    for I_bin = 1:N_bin
        x21_f(I_bin) = mean(xtemp21(((I_bin-1)*N_temp+1):I_bin*N_temp));
        x21_ferr(I_bin) = std0(xtemp21(((I_bin-1)*N_temp+1):I_bin*N_temp));
        x22_f(I_bin) = mean(xtemp22(((I_bin-1)*N_temp+1):I_bin*N_temp));
        x22_ferr(I_bin) = std0(xtemp22(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y2_f(I_bin) = mean(ytemp2(((I_bin-1)*N_temp+1):I_bin*N_temp));
        y2_ferr(I_bin) = std0(ytemp2(((I_bin-1)*N_temp+1):I_bin*N_temp));
        z2_f(I_bin) = mean(ztemp2(((I_bin-1)*N_temp+1):I_bin*N_temp));
        z2_ferr(I_bin) = std0(ztemp2(((I_bin-1)*N_temp+1):I_bin*N_temp));
    end
    errorbarxy(x1_f,y1_f,x1_ferr,y1_ferr,x1_ferr,y1_ferr,color_RGB(1,:),color_RGB(1,:),'o')
    errorbarxy(x2_f,y1_f,x2_ferr,y1_ferr,x2_ferr,y1_ferr,color_RGB(3,:),color_RGB(3,:),'*')
    errorbarxy(x21_f,y2_f,x21_ferr,y2_ferr,x21_ferr,y2_ferr,color_RGB(4,:),color_RGB(4,:),'o')
    errorbarxy(x22_f,y2_f,x22_ferr,y2_ferr,x22_ferr,y2_ferr,color_RGB(6,:),color_RGB(6,:),'*')
%    errorbarxy(x1_f,y1_f-y2_f,x1_ferr,sqrt(y1_ferr.^2+y2_ferr.^2),x1_ferr,sqrt(y1_ferr.^2+y2_ferr.^2),color_RGB(3,:),color_RGB(3,:),'--')

    ex1_f = zeros(size(EL_all));
    ex2_f = zeros(size(EL_all));
    ey1_f = zeros(size(EL_all));
    ez12_f = EL_all;
    ex21_f = zeros(size(EL_all));
    ex22_f = zeros(size(EL_all));
    ey2_f = zeros(size(EL_all));
    ez22_f = EL_all;
    ex1_ferr = zeros(size(EL_all));
    ex2_ferr = zeros(size(EL_all));
    ey1_ferr = zeros(size(EL_all));
    ez12_ferr = zeros(size(EL_all));
    ex21_ferr = zeros(size(EL_all));
    ex22_ferr = zeros(size(EL_all));
    ey2_ferr = zeros(size(EL_all));
    ez22_ferr = zeros(size(EL_all));

    for I_bin = 1:length(EL_all)
        ex1_f(I_bin) = mean(xtemp1((ztemp12 >= EL_all(I_bin)-Lbin) & (ztemp12 <= EL_all(I_bin)+Lbin)));
        ex1_ferr(I_bin) = std0(xtemp1((ztemp12 >= EL_all(I_bin)-Lbin) & (ztemp12 <= EL_all(I_bin)+Lbin)));
        ex2_f(I_bin) = mean(xtemp2((ztemp12 >= EL_all(I_bin)-Lbin) & (ztemp12 <= EL_all(I_bin)+Lbin)));
        ex2_ferr(I_bin) = std0(xtemp2((ztemp12 >= EL_all(I_bin)-Lbin) & (ztemp12 <= EL_all(I_bin)+Lbin)));
        ey1_f(I_bin) = mean(ytemp1((ztemp12 >= EL_all(I_bin)-Lbin) & (ztemp12 <= EL_all(I_bin)+Lbin)));
        ey1_ferr(I_bin) = std0(ytemp1((ztemp12 >= EL_all(I_bin)-Lbin) & (ztemp12 <= EL_all(I_bin)+Lbin)));

        ex21_f(I_bin) = mean(xtemp21((ztemp22 >= EL_all(I_bin)-Lbin) & (ztemp22 <= EL_all(I_bin)+Lbin)));
        ex21_ferr(I_bin) = std0(xtemp21((ztemp22 >= EL_all(I_bin)-Lbin) & (ztemp22 <= EL_all(I_bin)+Lbin)));
        ex22_f(I_bin) = mean(xtemp22((ztemp22 >= EL_all(I_bin)-Lbin) & (ztemp22 <= EL_all(I_bin)+Lbin)));
        ex22_ferr(I_bin) = std0(xtemp22((ztemp22 >= EL_all(I_bin)-Lbin) & (ztemp22 <= EL_all(I_bin)+Lbin)));
        ey2_f(I_bin) = mean(ytemp2((ztemp22 >= EL_all(I_bin)-Lbin) & (ztemp22 <= EL_all(I_bin)+Lbin)));
        ey2_ferr(I_bin) = std0(ytemp2((ztemp22 >= EL_all(I_bin)-Lbin) & (ztemp22 <= EL_all(I_bin)+Lbin)));
    end

    errorbarxy(ex1_f,ey1_f,ex1_ferr,ey1_ferr,ex1_ferr,ey1_ferr,color_RGB(1,:),color_RGB(1,:),'o')
    errorbarxy(ex2_f,ey1_f,ex2_ferr,ey1_ferr,ex2_ferr,ey1_ferr,color_RGB(3,:),color_RGB(3,:),'*')
    errorbarxy(ex21_f,ey2_f,ex21_ferr,ey2_ferr,ex21_ferr,ey2_ferr,color_RGB(4,:),color_RGB(4,:),'o')
    errorbarxy(ex22_f,ey2_f,ex22_ferr,ey2_ferr,ex22_ferr,ey2_ferr,color_RGB(6,:),color_RGB(6,:),'*')
%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        legend_text_5761 = cat(2,legend_text_5761,{[real_name,' - background all, r = ',num2str(r_size(ir))],[real_name,' - fake all, r = ',num2str(r_size(ir))],[control_name,' - background all, r = ',num2str(r_size(ir))],[control_name,' - fake all, r = ',num2str(r_size(ir))],['(',real_name,' - background) all bin, r = ',num2str(r_size(ir))],['(',real_name,' - fake) all bin, r = ',num2str(r_size(ir))],['(',control_name,' - background) all bin, r = ',num2str(r_size(ir))],['(',control_name,' - fake) all bin, r = ',num2str(r_size(ir))]});
        legend_text_5761 = cat(2,legend_text_5761,{['(',real_name,' - background) all bin (EL), r = ',num2str(r_size(ir))],['(',real_name,' - fake) all bin (EL), r = ',num2str(r_size(ir))],['(',control_name,' - background) all bin (EL), r = ',num2str(r_size(ir))],['(',control_name,' - fake) all bin (EL), r = ',num2str(r_size(ir))]});
    
    
    if r_size(ir) == r_optimize
        if sub_pos(3) > 0
            figure(50761)
            subaxis(sub_pos(1),sub_pos(2),sub_pos(3), 'Spacing', 0, 'PaddingRight',0.01, 'PaddingLeft',0.02, 'PaddingTop',0.045, 'PaddingBottom',0.04, 'Margin', 0.01);
            hold on
            errorbarxy(x1_f,y1_f,x1_ferr,y1_ferr,x1_ferr,y1_ferr,color_RGB(1,:),color_RGB(1,:),'o')
            errorbarxy(x2_f,y1_f,x2_ferr,y1_ferr,x2_ferr,y1_ferr,color_RGB(3,:),color_RGB(3,:),'*')
            errorbarxy(x21_f,y2_f,x21_ferr,y2_ferr,x21_ferr,y2_ferr,color_RGB(4,:),color_RGB(4,:),'o')
            errorbarxy(x22_f,y2_f,x22_ferr,y2_ferr,x22_ferr,y2_ferr,color_RGB(6,:),color_RGB(6,:),'*')

            errorbarxy(ex1_f,ey1_f,ex1_ferr,ey1_ferr,ex1_ferr,ey1_ferr,color_RGB(1,:),color_RGB(1,:),'o')
            errorbarxy(ex2_f,ey1_f,ex2_ferr,ey1_ferr,ex2_ferr,ey1_ferr,color_RGB(3,:),color_RGB(3,:),'*')
            errorbarxy(ex21_f,ey2_f,ex21_ferr,ey2_ferr,ex21_ferr,ey2_ferr,color_RGB(4,:),color_RGB(4,:),'o')
            errorbarxy(ex22_f,ey2_f,ex22_ferr,ey2_ferr,ex22_ferr,ey2_ferr,color_RGB(6,:),color_RGB(6,:),'*')

            title([image_folder,char(10),'Foci transcription level vs Foci enrichment, cycle = ',num2str(N_cycle),char(10),', mean TX = ',num2str(mean(nucleus_RNA_profile(:,4))),', DAPI slope = ',num2str(k_DAPI),', pair #: ',num2str(size(foci_data{ir},1)),', foci #: ',num2str(size(foci_data0{ir},1)),', ',control_name,' foci #: ',num2str(size(foci_data20{ir},1))],'Interpreter','none')
            ylabel(['Foci transcription level (',unit1,')'])
            xlabel(['Protein enrichment (',unit1,')'])

            legend_text_57610 = cat(2,legend_text_57610,{['(',real_name,' - background) all bin, r = ',num2str(r_size(ir))],['(',real_name,' - fake) all bin, r = ',num2str(r_size(ir))],['(',control_name,' - background) all bin, r = ',num2str(r_size(ir))],['(',control_name,' - fake) all bin, r = ',num2str(r_size(ir))]});
            legend_text_57610 = cat(2,legend_text_57610,{['(',real_name,' - background) all bin (EL), r = ',num2str(r_size(ir))],['(',real_name,' - fake) all bin (EL), r = ',num2str(r_size(ir))],['(',control_name,' - background) all bin (EL), r = ',num2str(r_size(ir))],['(',control_name,' - fake) all bin (EL), r = ',num2str(r_size(ir))]});
            
            legend(legend_text_57610)
            legend('hide')
            grid on

            out_TX0 = [x1_f',y1_f',z1_f',x2_f',y1_f',z1_f',x21_f',y2_f',z2_f',x22_f',y2_f',z2_f'];
            out_TX20 = {[xtemp1,ytemp1,ztemp1,xtemp2,ytemp1,ztemp1],[xtemp21,ytemp2,ztemp2,xtemp22,ytemp2,ztemp2]};
            out_eTX0 = [ex1_f',ey1_f',ez12_f',ex2_f',ey1_f',ez12_f',ex21_f',ey2_f',ez22_f',ex22_f',ey2_f',ez22_f'];
            out_eTX20 = {[xtemp1,ytemp1,ztemp12,xtemp2,ytemp1,ztemp12],[xtemp21,ytemp2,ztemp22,xtemp22,ytemp2,ztemp22]};

            figure(576)
        end
    end

    
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    
    
    
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
    
end

figure(5761)
title(['Foci transcription level vs Foci protein enrichment: ',image_folder,', cycle = ',num2str(N_cycle)],'Interpreter','none')
ylabel(['Foci transcription level (',unit1,')'])
xlabel(['Protein enrichment (',unit1,')'])
legend(legend_text_5761)
legend('hide')
grid on


figure(576)
title(['Foci protein enrichment vs nuclear mean protein intensity: ',image_folder,', cycle = ',num2str(N_cycle)],'Interpreter','none')
ylabel(['Protein enrichment (',unit1,')'])
xlabel(['P_m_e_a_n (',unit2,')'])
legend(legend_text)
legend('hide')
grid on


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(1,2,2)
legend_text = cell(0);

%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ir = 1:length(r_size)
    bar_mean = zeros(1,13);
    bar_err = zeros(1,13);
    
    ytemp1 = op_ratio0{ir};
    ytemp2 = ofp_ratio0{ir};
%     ytemp3 = foci_data{ir}(:,4)-fake_data{ir}(:,4);
    ytemp3 = ytemp1-ytemp2;
%     ytemp30 = ytemp3;
    temp_ind = (~isnan(ytemp1)) & (~isnan(ytemp2));
    bar_mean(1) = mean(ytemp1(temp_ind));
    bar_mean(2) = mean(ytemp2(temp_ind));
    bar_mean(3) = mean(ytemp3(temp_ind));
    bar_err(1) = std0(ytemp1(temp_ind));
    bar_err(2) = std0(ytemp2(temp_ind));
    bar_err(3) = std0(ytemp3(temp_ind));
    
    ytemp1 = op_ratio02{ir};
    ytemp2 = ofp_ratio02{ir};
%     ytemp3 = foci_data2{ir}(:,4)-fake_data2{ir}(:,4);
    ytemp3 = ytemp1-ytemp2;
%     ytemp4 = ytemp3-ytemp30;
    ytemp4 = op_ratio0{ir}-op_ratio02{ir};
    temp_ind = (~isnan(ytemp1)) & (~isnan(ytemp2)) & (~isnan(ytemp4));
    bar_mean(4) = mean(ytemp1(temp_ind));
    bar_mean(5) = mean(ytemp2(temp_ind));
    bar_mean(6) = mean(ytemp3(temp_ind));
    bar_mean(7) = mean(ytemp4(temp_ind));
    bar_err(4) = std0(ytemp1(temp_ind));
    bar_err(5) = std0(ytemp2(temp_ind));
    bar_err(6) = std0(ytemp3(temp_ind));
    bar_err(7) = std0(ytemp4(temp_ind));

    ytemp1 = op_ratio00{ir};
    ytemp2 = ofp_ratio00{ir};
%     ytemp3 = foci_data0{ir}(:,4)-fake_data0{ir}(:,4);
    ytemp3 = ytemp1-ytemp2;
    temp_ind = (~isnan(ytemp1)) & (~isnan(ytemp2));
    bar_mean(8) = mean(ytemp1(temp_ind));
    bar_mean(9) = mean(ytemp2(temp_ind));
    bar_mean(10) = mean(ytemp3(temp_ind));
    bar_err(8) = std0(ytemp1(temp_ind));
    bar_err(9) = std0(ytemp2(temp_ind));
    bar_err(10) = std0(ytemp3(temp_ind));
    
    ytemp1 = op_ratio020{ir};
    ytemp2 = ofp_ratio020{ir};
%     ytemp3 = foci_data20{ir}(:,4)-fake_data20{ir}(:,4);
    ytemp3 = ytemp1-ytemp2;
    temp_ind = (~isnan(ytemp1)) & (~isnan(ytemp2));
    bar_mean(11) = mean(ytemp1(temp_ind));
    bar_mean(12) = mean(ytemp2(temp_ind));
    bar_mean(13) = mean(ytemp3(temp_ind));
    bar_err(11) = std0(ytemp1(temp_ind));
    bar_err(12) = std0(ytemp2(temp_ind));
    bar_err(13) = std0(ytemp3(temp_ind));

    
    
    bar([1:length(bar_mean)],bar_mean,'FaceColor',color_RGB(ir,:))
    hold on
    errorbar([1:length(bar_mean)],bar_mean,bar_err,'Marker','none','LineStyle','none','color','k')
    legend_text = cat(2,legend_text,{['Mean enrichment (r = ',num2str(r_size(ir)),')']},{['Enrichment error (r = ',num2str(r_size(ir)),')']});
    
    
    
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if r_size(ir) == r_optimize
        if sub_pos(3) > 0
            figure(5077)
%             maximize(5077)
            subaxis(sub_pos(1),sub_pos(2),sub_pos(3), 'Spacing', 0, 'PaddingRight',0.01, 'PaddingLeft',0.02, 'PaddingTop',0.045, 'PaddingBottom',0.04, 'Margin', 0.01);
            hold on
            for I_bar = 1:length(bar_mean)
                bar(I_bar,bar_mean(I_bar),'FaceColor',color_RGB(I_bar,:))
                errorbar(I_bar,bar_mean(I_bar),bar_err(I_bar),'Marker','none','LineStyle','none','color','k')
            end
%                 set(gca,'XTick',[1:13],'XTickLabel',{'Foci-background','Fake-background','Foci-fake',[control_name,char(10),'foci-background'],[control_name,char(10),'fake-background'],[control_name,char(10),'foci-fake'],[real_name,'-',control_name],'Foci-background all','Fake-background all','Foci-fake all',[control_name,char(10),'foci-background all'],[control_name,char(10),'fake-background all'],[control_name,char(10),'foci-fake all']})
                xticklabel_rotate([1:13],90,{'Foci-background','Fake-background','Foci-fake',[control_name,' foci-background'],[control_name,' fake-background'],[control_name,' foci-fake'],[real_name,'-',control_name],'Foci-background all','Fake-background all','Foci-fake all',[control_name,' foci-background all'],[control_name,' fake-background all'],[control_name,' foci-fake all']},'FontSize',7,'FontWeight','bold')
                title([image_folder,char(10),'Mean foci enrichment comparison, cycle = ',num2str(N_cycle),char(10),', mean TX = ',num2str(mean(nucleus_RNA_profile(:,4))),', DAPI slope = ',num2str(k_DAPI),', pair #: ',num2str(size(foci_data{ir},1)),', foci #: ',num2str(size(foci_data0{ir},1)),', ',control_name,' foci #: ',num2str(size(foci_data20{ir},1))],'Interpreter','none')
                ylabel(['P_m_e_a_n (',unit2,')'])
                grid on

            figure(576)
        end
    end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
    
end

title(['Mean foci enrichment comparison: ',image_folder,', cycle = ',num2str(N_cycle)],'Interpreter','none')
ylabel(['P_m_e_a_n (',unit2,')'])
% set(gca,'XTick',[1:13],'XTickLabel',{'Foci-background','Fake-background','Foci-fake',[control_name,char(10),'foci-background'],[control_name,char(10),'fake-background'],[control_name,char(10),'foci-fake'],[real_name,'-',control_name],'Foci-background all','Fake-background all','Foci-fake all',[control_name,char(10),'foci-background all'],[control_name,char(10),'fake-background all'],[control_name,char(10),'foci-fake all']})
xticklabel_rotate([1:13],90,{'Foci-background','Fake-background','Foci-fake',[control_name,' foci-background'],[control_name,' fake-background'],[control_name,' foci-fake'],[real_name,'-',control_name],'Foci-background all','Fake-background all','Foci-fake all',[control_name,' foci-background all'],[control_name,' fake-background all'],[control_name,' foci-fake all']},'FontSize',7,'FontWeight','bold')
legend(legend_text)
legend('hide')
grid on





figure(581)
clf
mean_op = zeros(size(r_size));
mean_ofp = zeros(size(r_size));
mean_odp = zeros(size(r_size));
std_op = zeros(size(r_size));
std_ofp = zeros(size(r_size));
std_odp = zeros(size(r_size));
mean_op2 = zeros(size(r_size));
mean_ofp2 = zeros(size(r_size));
mean_odp2 = zeros(size(r_size));
mean_d = zeros(size(r_size));
std_op2 = zeros(size(r_size));
std_ofp2 = zeros(size(r_size));
std_odp2 = zeros(size(r_size));
std_d = zeros(size(r_size));
for ir = 1:length(r_size)
    mean_op(ir) = mean(op_ratio0{ir});
    mean_ofp(ir) = mean(ofp_ratio0{ir});
%     mean_odp(ir) = mean(foci_data{ir}(:,4)-fake_data{ir}(:,4));
    mean_odp(ir) = mean(op_ratio0{ir}-ofp_ratio0{ir});
    std_op(ir) = std0(op_ratio0{ir});
    std_ofp(ir) = std0(ofp_ratio0{ir});
%     std_odp(ir) = std0(foci_data{ir}(:,4)-fake_data{ir}(:,4));
    std_odp(ir) = std0(op_ratio0{ir}-ofp_ratio0{ir});
    mean_op2(ir) = mean(op_ratio02{ir});
    mean_ofp2(ir) = mean(ofp_ratio02{ir});
%     mean_odp2(ir) = mean(foci_data2{ir}(:,4)-fake_data2{ir}(:,4));
    mean_odp2(ir) = mean(op_ratio02{ir}-ofp_ratio02{ir});
%     mean_d(ir) = mean(foci_data{ir}(:,4)-foci_data2{ir}(:,4));
    mean_d(ir) = mean(op_ratio0{ir}-op_ratio02{ir});
    std_op2(ir) = std0(op_ratio02{ir});
    std_ofp2(ir) = std0(ofp_ratio02{ir});
%     std_odp2(ir) = std0(foci_data2{ir}(:,4)-fake_data2{ir}(:,4));
    std_odp2(ir) = std0(op_ratio0{ir}-op_ratio02{ir});
%     std_d(ir) = std0(foci_data{ir}(:,4)-foci_data2{ir}(:,4));
    std_d(ir) = std0(op_ratio0{ir}-op_ratio02{ir});
end
hold on
errorbar(r_size,mean_op,std_op,'go-');
errorbar(r_size,mean_ofp,std_ofp,'b*-');
errorbar(r_size,mean_odp,std_odp,'r--')
errorbar(r_size,mean_op2,std_op2,'yo-');
errorbar(r_size,mean_ofp2,std_ofp2,'c*-');
errorbar(r_size,mean_odp2,std_odp2,'m--')
errorbar(r_size,mean_d,std_d,'k--')


mean_op0 = zeros(size(r_size));
mean_ofp0 = zeros(size(r_size));
mean_odp0 = zeros(size(r_size));
std_op0 = zeros(size(r_size));
std_ofp0 = zeros(size(r_size));
std_odp0 = zeros(size(r_size));
mean_op20 = zeros(size(r_size));
mean_ofp20 = zeros(size(r_size));
mean_odp20 = zeros(size(r_size));
mean_d0 = zeros(size(r_size));
std_op20 = zeros(size(r_size));
std_ofp20 = zeros(size(r_size));
std_odp20 = zeros(size(r_size));
std_d0 = zeros(size(r_size));
for ir = 1:length(r_size)
    mean_op0(ir) = mean(op_ratio00{ir});
    mean_ofp0(ir) = mean(ofp_ratio00{ir});
%     mean_odp0(ir) = mean(foci_data0{ir}(:,4)-fake_data0{ir}(:,4));
    mean_odp0(ir) = mean(op_ratio00{ir}-ofp_ratio00{ir});
    std_op0(ir) = std0(op_ratio00{ir});
    std_ofp0(ir) = std0(ofp_ratio00{ir});
%     std_odp0(ir) = std0(foci_data0{ir}(:,4)-fake_data0{ir}(:,4));
    std_odp0(ir) = std0(op_ratio00{ir}-ofp_ratio00{ir});
    mean_op20(ir) = mean(op_ratio020{ir});
    mean_ofp20(ir) = mean(ofp_ratio020{ir});
%     mean_odp20(ir) = mean(foci_data20{ir}(:,4)-fake_data20{ir}(:,4));
    mean_odp20(ir) = mean(op_ratio020{ir}-ofp_ratio020{ir});
%     mean_d0(ir) = mean(foci_data0{ir}(:,4))-mean(foci_data20{ir}(:,4));
    mean_d0(ir) = mean(op_ratio00{ir})-mean(op_ratio020{ir});
    std_op20(ir) = std0(op_ratio020{ir});
    std_ofp20(ir) = std0(ofp_ratio020{ir});
%     std_odp20(ir) = std0(foci_data20{ir}(:,4)-fake_data20{ir}(:,4));
    std_odp20(ir) = std0(op_ratio020{ir}-ofp_ratio020{ir});
%     std_d0(ir) = sqrt(std0(foci_data0{ir}(:,4)).^2+std0(foci_data20{ir}(:,4)).^2);
    std_d0(ir) = sqrt(std0(op_ratio020{ir}).^2+std0(ofp_ratio020{ir}).^2);
end
hold on
errorbar(r_size,mean_op0,std_op0,'go-');
errorbar(r_size,mean_ofp0,std_ofp0,'b*-');
errorbar(r_size,mean_odp0,std_odp0,'r--')
errorbar(r_size,mean_op20,std_op20,'yo-');
errorbar(r_size,mean_ofp20,std_ofp20,'c*-');
errorbar(r_size,mean_odp20,std_odp20,'m--')
errorbar(r_size,mean_d0,std_d0,'k--')




title([image_folder,char(10),'Mean absolute enrichment vs integration radius, cycle = ',num2str(N_cycle),char(10),', mean TX = ',num2str(mean(nucleus_RNA_profile(:,4))),', DAPI slope = ',num2str(k_DAPI)],'Interpreter','none')
ylabel(['Protein enrichment (',unit1,')'])
xlabel(['Radius (pixel)'])
legend([real_name,' - background'],['fake ',real_name,' - background'],[real_name,' - fake'],[control_name,' - background'],[control_name,' fake - background'],[control_name,' - fake'],[real_name,' -',control_name],[real_name,' - background all'],'fake - background all',[real_name,' - fake all'],[control_name,' - background all'],[control_name,' - background all'],[control_name,' - fake all'],[real_name,' - ',control_name,' all'])
legend('hide')
grid on

if sub_pos(3) > 0
    figure(5081)
    subaxis(sub_pos(1),sub_pos(2),sub_pos(3), 'Spacing', 0, 'PaddingRight',0.01, 'PaddingLeft',0.02, 'PaddingTop',0.045, 'PaddingBottom',0.04, 'Margin', 0.01);
    hold on
    errorbar(r_size,mean_op,std_op,'go-');
    errorbar(r_size,mean_ofp,std_ofp,'b*-');
    errorbar(r_size,mean_odp,std_odp,'r--')
    errorbar(r_size,mean_op2,std_op2,'yo-');
    errorbar(r_size,mean_ofp2,std_ofp2,'c*-');
    errorbar(r_size,mean_odp2,std_odp2,'m--')
    errorbar(r_size,mean_d,std_d,'k--')

    
    errorbar(r_size,mean_op0,std_op0,'go-');
    errorbar(r_size,mean_ofp0,std_ofp0,'b*-');
    errorbar(r_size,mean_odp0,std_odp0,'r--')
    errorbar(r_size,mean_op20,std_op20,'yo-');
    errorbar(r_size,mean_ofp20,std_ofp20,'c*-');
    errorbar(r_size,mean_odp20,std_odp20,'m--')
    errorbar(r_size,mean_d0,std_d0,'k--')
    
    
    title([image_folder,char(10),'Enrichment vs integration radius: cycle = ',num2str(N_cycle),char(10),', mean TX = ',num2str(mean(nucleus_RNA_profile(:,4))),', DAPI slope = ',num2str(k_DAPI)],'Interpreter','none')
    ylabel(['Protein enrichment (',unit1,')'])
    xlabel(['Radius (pixel)'])
    legend([real_name,' - background'],['fake ',real_name,' - background'],[real_name,' - fake'],[control_name,' - background'],[control_name,' fake - background'],[control_name,' - fake'],[real_name,' - ',control_name],[real_name,' - background all'],'fake - background all',[real_name,' - fake all'],[control_name,' - background all'],[control_name,' - background all'],[control_name,' - fake all'],[real_name,' - ',control_name,' all'])
    legend('hide')
    grid on
end










figure(582)
clf
legend_text = cell(0);
for ir = 1:length(r_size)
    EL_op_mean = zeros(size(EL_all));
    EL_op_std = zeros(size(EL_all));
    EL_ofp_mean = zeros(size(EL_all));
    EL_ofp_std = zeros(size(EL_all));
    EL_op_diff_mean = zeros(size(EL_all));
    EL_op_diff_std = zeros(size(EL_all));
    EL_op2_mean = zeros(size(EL_all));
    EL_op2_std = zeros(size(EL_all));
    EL_ofp2_mean = zeros(size(EL_all));
    EL_ofp2_std = zeros(size(EL_all));
    EL_op2_diff_mean = zeros(size(EL_all));
    EL_op2_diff_std = zeros(size(EL_all));
    EL_diff_mean = zeros(size(EL_all));
    EL_diff_std = zeros(size(EL_all));

    for I_EL = 1:length(EL_all)
        EL_op_mean(I_EL) = mean(op_ratio0{ir}((foci_data{ir}(:,7) >= EL_all(I_EL)-Lbin) & (foci_data{ir}(:,7) <= EL_all(I_EL)+Lbin)));
        EL_op_std(I_EL) = std0(op_ratio0{ir}((foci_data{ir}(:,7) >= EL_all(I_EL)-Lbin) & (foci_data{ir}(:,7) <= EL_all(I_EL)+Lbin)));
        EL_ofp_mean(I_EL) = mean(ofp_ratio0{ir}((fake_data{ir}(:,7) >= EL_all(I_EL)-Lbin) & (fake_data{ir}(:,7) <= EL_all(I_EL)+Lbin)));
        EL_ofp_std(I_EL) = std0(ofp_ratio0{ir}((fake_data{ir}(:,7) >= EL_all(I_EL)-Lbin) & (fake_data{ir}(:,7) <= EL_all(I_EL)+Lbin)));
        EL_op_diff_mean(I_EL) = mean(foci_data{ir}((foci_data{ir}(:,7) >= EL_all(I_EL)-Lbin) & (foci_data{ir}(:,7) <= EL_all(I_EL)+Lbin),4)-fake_data{ir}((foci_data{ir}(:,7) >= EL_all(I_EL)-Lbin) & (foci_data{ir}(:,7) <= EL_all(I_EL)+Lbin),4));
        EL_op_diff_std(I_EL) = std0(foci_data{ir}((foci_data{ir}(:,7) >= EL_all(I_EL)-Lbin) & (foci_data{ir}(:,7) <= EL_all(I_EL)+Lbin),4)-fake_data{ir}((foci_data{ir}(:,7) >= EL_all(I_EL)-Lbin) & (foci_data{ir}(:,7) <= EL_all(I_EL)+Lbin),4));
        EL_op2_mean(I_EL) = mean(op_ratio02{ir}((foci_data{ir}(:,7) >= EL_all(I_EL)-Lbin) & (foci_data{ir}(:,7) <= EL_all(I_EL)+Lbin)));
        EL_op2_std(I_EL) = std0(op_ratio02{ir}((foci_data{ir}(:,7) >= EL_all(I_EL)-Lbin) & (foci_data{ir}(:,7) <= EL_all(I_EL)+Lbin)));
        EL_ofp2_mean(I_EL) = mean(ofp_ratio02{ir}((fake_data{ir}(:,7) >= EL_all(I_EL)-Lbin) & (fake_data{ir}(:,7) <= EL_all(I_EL)+Lbin)));
        EL_ofp2_std(I_EL) = std0(ofp_ratio02{ir}((fake_data{ir}(:,7) >= EL_all(I_EL)-Lbin) & (fake_data{ir}(:,7) <= EL_all(I_EL)+Lbin)));
        EL_op2_diff_mean(I_EL) = mean(foci_data2{ir}((foci_data{ir}(:,7) >= EL_all(I_EL)-Lbin) & (foci_data{ir}(:,7) <= EL_all(I_EL)+Lbin),4)-fake_data2{ir}((foci_data{ir}(:,7) >= EL_all(I_EL)-Lbin) & (foci_data{ir}(:,7) <= EL_all(I_EL)+Lbin),4));
        EL_op2_diff_std(I_EL) = std0(foci_data2{ir}((foci_data{ir}(:,7) >= EL_all(I_EL)-Lbin) & (foci_data{ir}(:,7) <= EL_all(I_EL)+Lbin),4)-fake_data2{ir}((foci_data{ir}(:,7) >= EL_all(I_EL)-Lbin) & (foci_data{ir}(:,7) <= EL_all(I_EL)+Lbin),4));
        EL_diff_mean(I_EL) = mean(foci_data{ir}((foci_data{ir}(:,7) >= EL_all(I_EL)-Lbin) & (foci_data{ir}(:,7) <= EL_all(I_EL)+Lbin),4)-foci_data2{ir}((foci_data{ir}(:,7) >= EL_all(I_EL)-Lbin) & (foci_data{ir}(:,7) <= EL_all(I_EL)+Lbin),4));
        EL_diff_std(I_EL) = std0(foci_data{ir}((foci_data{ir}(:,7) >= EL_all(I_EL)-Lbin) & (foci_data{ir}(:,7) <= EL_all(I_EL)+Lbin),4)-foci_data2{ir}((foci_data{ir}(:,7) >= EL_all(I_EL)-Lbin) & (foci_data{ir}(:,7) <= EL_all(I_EL)+Lbin),4));
    end
    plot(foci_data{ir}(:,7),op_ratio0{ir},'g.')
    hold on
    plot(fake_data{ir}(:,7),ofp_ratio0{ir},'b*')
    plot(foci_data2{ir}(:,7),op_ratio02{ir},'ro')
    plot(fake_data2{ir}(:,7),ofp_ratio02{ir},'y+')
    
    errorbar(EL_all,EL_op_mean,EL_op_std,'g-')
    errorbar(EL_all,EL_ofp_mean,EL_ofp_std,'b-')
    errorbar(EL_all,EL_op_diff_mean,EL_op_diff_std,'r--')
    errorbar(EL_all,EL_op2_mean,EL_op2_std,'y-')
    errorbar(EL_all,EL_ofp2_mean,EL_ofp2_std,'c-')
    errorbar(EL_all,EL_op2_diff_mean,EL_op2_diff_std,'m--')
    errorbar(EL_all,EL_diff_mean,EL_diff_std,'k--')

    legend_text = cat(2,legend_text,{[real_name,' raw enrichment (r = ',num2str(r_size(ir)),')'],[real_name,' fake raw enrichment (r = ',num2str(r_size(ir)),')'],[control_name,' raw enrichment (r = ',num2str(r_size(ir)),')'],[control_name,' fake raw enrichment (r = ',num2str(r_size(ir)),')'],[real_name,' - background (r = ',num2str(r_size(ir)),')'],[real_name,' fake - background (r = ',num2str(r_size(ir)),')'],[real_name,' - fake (r = ',num2str(r_size(ir)),')'],[control_name,' - background (r = ',num2str(r_size(ir)),')'],[control_name,' fake - background (r = ',num2str(r_size(ir)),')'],[control_name,' - fake (r = ',num2str(r_size(ir)),')'],[real_name,' - ',control_name,' (r = ',num2str(r_size(ir)),')']});


    
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    EL_op_mean0 = zeros(size(EL_all));
    EL_op_std0 = zeros(size(EL_all));
    EL_ofp_mean0 = zeros(size(EL_all));
    EL_ofp_std0 = zeros(size(EL_all));
    EL_op_diff_mean0 = zeros(size(EL_all));
    EL_op_diff_std0 = zeros(size(EL_all));
    EL_op2_mean0 = zeros(size(EL_all));
    EL_op2_std0 = zeros(size(EL_all));
    EL_ofp2_mean0 = zeros(size(EL_all));
    EL_ofp2_std0 = zeros(size(EL_all));
    EL_op2_diff_mean0 = zeros(size(EL_all));
    EL_op2_diff_std0 = zeros(size(EL_all));
    EL_diff_mean0 = zeros(size(EL_all));
    EL_diff_std0 = zeros(size(EL_all));

    for I_EL = 1:length(EL_all)
        EL_op_mean0(I_EL) = mean(op_ratio00{ir}((foci_data0{ir}(:,7) >= EL_all(I_EL)-Lbin) & (foci_data0{ir}(:,7) <= EL_all(I_EL)+Lbin)));
        EL_op_std0(I_EL) = std0(op_ratio00{ir}((foci_data0{ir}(:,7) >= EL_all(I_EL)-Lbin) & (foci_data0{ir}(:,7) <= EL_all(I_EL)+Lbin)));
        EL_ofp_mean0(I_EL) = mean(ofp_ratio00{ir}((fake_data0{ir}(:,7) >= EL_all(I_EL)-Lbin) & (fake_data0{ir}(:,7) <= EL_all(I_EL)+Lbin)));
        EL_ofp_std0(I_EL) = std0(ofp_ratio00{ir}((fake_data0{ir}(:,7) >= EL_all(I_EL)-Lbin) & (fake_data0{ir}(:,7) <= EL_all(I_EL)+Lbin)));
        EL_op_diff_mean0(I_EL) = mean(foci_data0{ir}((foci_data0{ir}(:,7) >= EL_all(I_EL)-Lbin) & (foci_data0{ir}(:,7) <= EL_all(I_EL)+Lbin),4)-fake_data0{ir}((foci_data0{ir}(:,7) >= EL_all(I_EL)-Lbin) & (foci_data0{ir}(:,7) <= EL_all(I_EL)+Lbin),4));
        EL_op_diff_std0(I_EL) = std0(foci_data0{ir}((foci_data0{ir}(:,7) >= EL_all(I_EL)-Lbin) & (foci_data0{ir}(:,7) <= EL_all(I_EL)+Lbin),4)-fake_data0{ir}((foci_data0{ir}(:,7) >= EL_all(I_EL)-Lbin) & (foci_data0{ir}(:,7) <= EL_all(I_EL)+Lbin),4));
        EL_op2_mean0(I_EL) = mean(op_ratio020{ir}((foci_data20{ir}(:,7) >= EL_all(I_EL)-Lbin) & (foci_data20{ir}(:,7) <= EL_all(I_EL)+Lbin)));
        EL_op2_std0(I_EL) = std0(op_ratio020{ir}((foci_data20{ir}(:,7) >= EL_all(I_EL)-Lbin) & (foci_data20{ir}(:,7) <= EL_all(I_EL)+Lbin)));
        EL_ofp2_mean0(I_EL) = mean(ofp_ratio020{ir}((fake_data20{ir}(:,7) >= EL_all(I_EL)-Lbin) & (fake_data20{ir}(:,7) <= EL_all(I_EL)+Lbin)));
        EL_ofp2_std0(I_EL) = std0(ofp_ratio020{ir}((fake_data20{ir}(:,7) >= EL_all(I_EL)-Lbin) & (fake_data20{ir}(:,7) <= EL_all(I_EL)+Lbin)));
        EL_op2_diff_mean0(I_EL) = mean(foci_data20{ir}((foci_data20{ir}(:,7) >= EL_all(I_EL)-Lbin) & (foci_data20{ir}(:,7) <= EL_all(I_EL)+Lbin),4)-fake_data20{ir}((foci_data20{ir}(:,7) >= EL_all(I_EL)-Lbin) & (foci_data20{ir}(:,7) <= EL_all(I_EL)+Lbin),4));
        EL_op2_diff_std0(I_EL) = std0(foci_data20{ir}((foci_data20{ir}(:,7) >= EL_all(I_EL)-Lbin) & (foci_data20{ir}(:,7) <= EL_all(I_EL)+Lbin),4)-fake_data20{ir}((foci_data20{ir}(:,7) >= EL_all(I_EL)-Lbin) & (foci_data20{ir}(:,7) <= EL_all(I_EL)+Lbin),4));
        EL_diff_mean0(I_EL) = mean(foci_data0{ir}((foci_data0{ir}(:,7) >= EL_all(I_EL)-Lbin) & (foci_data0{ir}(:,7) <= EL_all(I_EL)+Lbin),4))-mean(foci_data20{ir}((foci_data20{ir}(:,7) >= EL_all(I_EL)-Lbin) & (foci_data20{ir}(:,7) <= EL_all(I_EL)+Lbin),4));
        EL_diff_std0(I_EL) = sqrt(std0(foci_data0{ir}((foci_data0{ir}(:,7) >= EL_all(I_EL)-Lbin) & (foci_data0{ir}(:,7) <= EL_all(I_EL)+Lbin),4)).^2+std0(foci_data20{ir}((foci_data20{ir}(:,7) >= EL_all(I_EL)-Lbin) & (foci_data20{ir}(:,7) <= EL_all(I_EL)+Lbin),4)).^2);
    end
    plot(foci_data0{ir}(:,7),op_ratio00{ir},'g.')
    hold on
    plot(fake_data0{ir}(:,7),ofp_ratio00{ir},'b*')
    plot(foci_data20{ir}(:,7),op_ratio020{ir},'ro')
    plot(fake_data20{ir}(:,7),ofp_ratio020{ir},'y+')
    
    errorbar(EL_all,EL_op_mean0,EL_op_std0,'g-')
    errorbar(EL_all,EL_ofp_mean0,EL_ofp_std0,'b-')
    errorbar(EL_all,EL_op_diff_mean0,EL_op_diff_std0,'r--')
    errorbar(EL_all,EL_op2_mean0,EL_op2_std0,'y-')
    errorbar(EL_all,EL_ofp2_mean0,EL_ofp2_std0,'c-')
    errorbar(EL_all,EL_op2_diff_mean0,EL_op2_diff_std0,'m--')
    errorbar(EL_all,EL_diff_mean0,EL_diff_std0,'k--')

    legend_text = cat(2,legend_text,{[real_name,' all raw enrichment (r = ',num2str(r_size(ir)),')'],[real_name,' fake all raw enrichment (r = ',num2str(r_size(ir)),')'],[control_name,' all raw enrichment (r = ',num2str(r_size(ir)),')'],[control_name,' fake all raw enrichment (r = ',num2str(r_size(ir)),')'],[real_name,' - background all (r = ',num2str(r_size(ir)),')'],[real_name,' fake - background all (r = ',num2str(r_size(ir)),')'],[real_name,' - fake all (r = ',num2str(r_size(ir)),')'],[control_name,' - background all (r = ',num2str(r_size(ir)),')'],[control_name,' fake - background all (r = ',num2str(r_size(ir)),')'],[control_name,' - fake all (r = ',num2str(r_size(ir)),')'],[real_name,' - ',control_name,' all (r = ',num2str(r_size(ir)),')']});
    
    
    
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
    
    if r_size(ir) == r_optimize
        if sub_pos(3) > 0
            out_EL = [EL_all;EL_op_mean;EL_all;EL_ofp_mean;EL_all;EL_op_diff_mean;EL_all;EL_op2_mean;EL_all;EL_ofp2_mean;EL_all;EL_op2_diff_mean;EL_all;EL_diff_mean]';
            out_EL0 = [EL_all;EL_op_mean0;EL_all;EL_ofp_mean0;EL_all;EL_op_diff_mean0;EL_all;EL_op2_mean0;EL_all;EL_ofp2_mean0;EL_all;EL_op2_diff_mean0]';
            out_EL2 = [foci_data{ir}(:,7),op_ratio0{ir},fake_data{ir}(:,7),ofp_ratio0{ir},foci_data{ir}(:,7),op_ratio0{ir}-ofp_ratio0{ir},foci_data2{ir}(:,7),op_ratio02{ir},fake_data2{ir}(:,7),ofp_ratio02{ir},foci_data2{ir}(:,7),op_ratio02{ir}-ofp_ratio02{ir},foci_data{ir}(:,7),op_ratio0{ir}-op_ratio02{ir}];
            out_EL20 = {[foci_data0{ir}(:,7),op_ratio00{ir},fake_data0{ir}(:,7),ofp_ratio00{ir},foci_data0{ir}(:,7),op_ratio00{ir}-ofp_ratio00{ir}],[foci_data20{ir}(:,7),op_ratio020{ir},fake_data20{ir}(:,7),ofp_ratio020{ir},foci_data20{ir}(:,7),op_ratio020{ir}-ofp_ratio020{ir}]};
            
            figure(5082)
            h_axes = subaxis(sub_pos(1),sub_pos(2),sub_pos(3), 'Spacing', 0, 'PaddingRight',0.01, 'PaddingLeft',0.02, 'PaddingTop',0.045, 'PaddingBottom',0.04, 'Margin', 0.01);
            hold on

            plot(foci_data{ir}(:,7),op_ratio0{ir},'g.')
            plot(fake_data{ir}(:,7),ofp_ratio0{ir},'b*')
            plot(foci_data2{ir}(:,7),op_ratio02{ir},'ro')
            plot(fake_data2{ir}(:,7),ofp_ratio02{ir},'y+')

            errorbar(EL_all,EL_op_mean,EL_op_std,'g-')
            errorbar(EL_all,EL_ofp_mean,EL_ofp_std,'b-')
            errorbar(EL_all,EL_op_diff_mean,EL_op_diff_std,'r--')
            errorbar(EL_all,EL_op2_mean,EL_op2_std,'y-')
            errorbar(EL_all,EL_ofp2_mean,EL_ofp2_std,'c-')
            errorbar(EL_all,EL_op2_diff_mean,EL_op2_diff_std,'m--')
            errorbar(EL_all,EL_diff_mean,EL_diff_std,'k--')
           
            plot(foci_data0{ir}(:,7),op_ratio00{ir},'g.')
            plot(fake_data0{ir}(:,7),ofp_ratio00{ir},'b*')
            plot(foci_data20{ir}(:,7),op_ratio020{ir},'ro')
            plot(fake_data20{ir}(:,7),ofp_ratio020{ir},'y+')

            errorbar(EL_all,EL_op_mean0,EL_op_std0,'g-')
            errorbar(EL_all,EL_ofp_mean0,EL_ofp_std0,'b-')
            errorbar(EL_all,EL_op_diff_mean0,EL_op_diff_std0,'r--')
            errorbar(EL_all,EL_op2_mean0,EL_op2_std0,'y-')
            errorbar(EL_all,EL_ofp2_mean0,EL_ofp2_std0,'c-')
            errorbar(EL_all,EL_op2_diff_mean0,EL_op2_diff_std0,'m--')
            errorbar(EL_all,EL_diff_mean0,EL_diff_std0,'k--')
            
            
            
            legend_text0 = {[real_name,' raw enrichment (r = ',num2str(r_size(ir)),')'],[real_name,' fake raw enrichment (r = ',num2str(r_size(ir)),')'],[control_name,' raw enrichment (r = ',num2str(r_size(ir)),')'],[control_name,' fake raw enrichment (r = ',num2str(r_size(ir)),')'],[real_name,' - background (r = ',num2str(r_size(ir)),')'],[real_name,' fake - background (r = ',num2str(r_size(ir)),')'],[real_name,' - fake (r = ',num2str(r_size(ir)),')'],[control_name,' - background (r = ',num2str(r_size(ir)),')'],[control_name,' fake - background (r = ',num2str(r_size(ir)),')'],[control_name,' - fake (r = ',num2str(r_size(ir)),')'],[real_name,' - ',control_name,' (r = ',num2str(r_size(ir)),')']};
            legend_text0 = cat(2,legend_text0,{[real_name,' all raw enrichment (r = ',num2str(r_size(ir)),')'],[real_name,' fake all raw enrichment (r = ',num2str(r_size(ir)),')'],[control_name,' all raw enrichment (r = ',num2str(r_size(ir)),')'],[control_name,' fake all raw enrichment (r = ',num2str(r_size(ir)),')'],[real_name,' - background all (r = ',num2str(r_size(ir)),')'],[real_name,' fake - background all (r = ',num2str(r_size(ir)),')'],[real_name,' - fake all (r = ',num2str(r_size(ir)),')'],[control_name,' - background all (r = ',num2str(r_size(ir)),')'],[control_name,' fake - background all (r = ',num2str(r_size(ir)),')'],[control_name,' - fake all (r = ',num2str(r_size(ir)),')'],[real_name,' - ',control_name,' all (r = ',num2str(r_size(ir)),')']});
            legend(legend_text0)
            legend('hide')
            
            axes
            plot(nucleus_protein_profile(:,1),nucleus_protein_profile(:,2),'Marker','.','LineStyle','none','Color',[0.5,0.5,0.5])
            xlim([Lmin,Lmax])
            set(gca,'YAxisLocation','right','Color','none','Position',get(h_axes,'Position'))
                
            xlim([Lmin,Lmax])
            ylabel(['Protein enrichment (',unit1,')'])
            xlabel(['EL'])
            grid on
            title([image_folder,char(10),'Enrichment vs EL length, cycle = ',num2str(N_cycle),char(10),', mean TX = ',num2str(mean(nucleus_RNA_profile(:,4))),', DAPI slope = ',num2str(k_DAPI),', pair #: ',num2str(size(foci_data{ir},1)),', foci #: ',num2str(size(foci_data0{ir},1)),', ',control_name,' foci #: ',num2str(size(foci_data20{ir},1))],'Interpreter','none')


            figure(582)
        end
    end

    
end

legend_text = cat(2,legend_text,{'Protein profile'});
legend(legend_text)
legend('hide')
xlim([Lmin,Lmax])
ylabel(['Protein enrichment (',unit1,')'])
xlabel(['EL'])

axes
plot(nucleus_protein_profile(:,1),nucleus_protein_profile(:,2),'ko')
set(gca,'YAxisLocation','right','color','none')
xlim([Lmin,Lmax])
    
title(['Enrichment vs EL length: ',image_folder,', cycle = ',num2str(N_cycle)],'Interpreter','none')
grid on






function y = std0(x)
    NN = size(x);
    if NN(1) == 1;
        N1 = NN(2);
    else
        N1 = NN(1);
    end
    y = std(x)/sqrt(N1);
