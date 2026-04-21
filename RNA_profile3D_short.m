function [FISHI_out,FISHN_out,FISHI_out2,FISHN_out2] = RNA_profile3D_short(nucleus_RNA_profile,foci_RNA_profile,cytoplasmic_RNA_profile,image_folder,N_cycle,flip_EL)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A function to re-plot RNA profile %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Parameter setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nucleus_bin = 0:0.05:1;%0.025:0.05:0.975;
average_radius = 0.05;
EL_min = 0.25;
EL_max = 0.75;
intensity_Nbin = 50;
bin_max = min(nucleus_bin+average_radius,1);
bin_min = max(nucleus_bin-average_radius,0);
if flip_EL
    nucleus_RNA_profile(:,1) = 1-nucleus_RNA_profile(:,1);
    cytoplasmic_RNA_profile(:,1) = 1-cytoplasmic_RNA_profile(:,1);
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Re-plot FISH intensity vs EL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nucleus_distance = nucleus_RNA_profile(:,1);
Ifoci_nucleus = nucleus_RNA_profile(:,4);
% figure(5)
% subaxis(sub_pos(1),sub_pos(2),sub_pos(3), 'Spacing', 0, 'PaddingRight',0.01, 'PaddingLeft',0.02, 'PaddingTop',0.045, 'PaddingBottom',0.04, 'Margin', 0.01);
    fi0 = zeros(size(nucleus_bin));
    fi1 = zeros(size(nucleus_bin));
    fi2 = zeros(size(nucleus_bin));
    for I_bin = 1:length(nucleus_bin)
        fi0(I_bin) = mean(nucleus_RNA_profile((nucleus_distance >= bin_min(I_bin))&(nucleus_distance <= bin_max(I_bin)),4));
        fi1(I_bin) = std(nucleus_RNA_profile((nucleus_distance >= bin_min(I_bin))&(nucleus_distance <= bin_max(I_bin)),4));
        fi2(I_bin) = fi1(I_bin)/sqrt(nnz((nucleus_distance >= bin_min(I_bin))&(nucleus_distance <= bin_max(I_bin))));
    end
    FISHI_out = [nucleus_bin',fi0',fi1'];
    good_bin = find(~isnan(fi1));
    
%     plot(nucleus_RNA_profile(:,1),nucleus_RNA_profile(:,2),'bo'); hold on
%     plot(nucleus_RNA_profile(:,1),nucleus_RNA_profile(:,4),'Color',[0.7,0.7,0.7],'Marker','.','LineStyle','none')
%     plot(nucleus_RNA_profile(:,1),nucleus_RNA_profile(:,5),'r^',cytoplasmic_RNA_profile(:,1),cytoplasmic_RNA_profile(:,2),'c*')
%     
%     errorbar(nucleus_bin(good_bin),fi0(good_bin),fi2(good_bin),'k')
%     
%     xpatch = [nucleus_bin(good_bin),nucleus_bin(good_bin(end:-1:1))];
%     ypatch = [fi0(good_bin)-fi1(good_bin),fi0(good_bin(end:-1:1))+fi1(good_bin(end:-1:1))];
%     patch(xpatch,ypatch,'r','FaceAlpha',0.5,'EdgeColor','none');
%     
%     title([image_folder,', cycle = ',num2str(N_cycle)],'Interpreter','none')
%     xlabel('EL')
%     ylabel('# of mRNA')
%     legend('Nucleus mean intensity','Nucleus foci intensity','Nucleus background mean intensity','Cytoplasmic mean intensity','Mean nucleus foci intensity','Std of nucleus foci intensity')
%     legend('hide')
% 
    %%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ELg = nucleus_bin(good_bin & nucleus_bin >= EL_min & nucleus_bin <= EL_max);
    fi0g = fi0(good_bin & nucleus_bin >= EL_min & nucleus_bin <= EL_max);
    fi1g = fi1(good_bin & nucleus_bin >= EL_min & nucleus_bin <= EL_max);
    fi2g = fi2(good_bin & nucleus_bin >= EL_min & nucleus_bin <= EL_max);
    
    [max_FISHI,max_I] = max(fi0g);
    [min_FISHI,min_I] = min(fi0g(max_I:end));
    min_I = max_I+min_I-1;
    max_ELI = ELg(max_I);
    max_stdI = fi1g(max_I);
    max_errI = fi2g(max_I);
    mid_FISHI = (max_FISHI+min_FISHI)/2;
    [~,I0] = min(abs(fi0g(max_I:min_I)-mid_FISHI));
    p0 = polyfit(ELg((I0+max_I-2):(I0+max_I)),fi0g((I0+max_I-2):(I0+max_I)),1);
    mid_ELI = (mid_FISHI-p0(2))/p0(1);
    
    FISHI_out2 = [max_FISHI,max_stdI,max_errI,max_ELI,mid_ELI];
    %%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
% figure(50)
% clf
%     plot(nucleus_RNA_profile(:,1),nucleus_RNA_profile(:,2),'bo'); hold on
%     plot(nucleus_RNA_profile(:,1),nucleus_RNA_profile(:,4),'Color',[0.7,0.7,0.7],'Marker','.','LineStyle','none')
%     plot(nucleus_RNA_profile(:,1),nucleus_RNA_profile(:,5),'r^',cytoplasmic_RNA_profile(:,1),cytoplasmic_RNA_profile(:,2),'c*')
%     errorbar(nucleus_bin(good_bin),fi0(good_bin),fi2(good_bin),'k')
%     patch(xpatch,ypatch,'r','FaceAlpha',0.5,'EdgeColor','none');
%     
%     title([image_folder,', cycle = ',num2str(N_cycle)],'Interpreter','none')
%     xlabel('EL')
%     ylabel('# of mRNA')
%     legend('Nucleus mean intensity','Nucleus foci intensity','Nucleus background mean intensity','Cytoplasmic mean intensity','Mean nucleus foci intensity','Std of nucleus foci intensity')
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Re-plot FISH intensity vs EL new %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(7)
% subaxis(sub_pos(1),sub_pos(2),sub_pos(3), 'Spacing', 0, 'PaddingRight',0.01, 'PaddingLeft',0.02, 'PaddingTop',0.045, 'PaddingBottom',0.04, 'Margin', 0.01);
    fi0 = zeros(size(nucleus_bin));
    fi1 = zeros(size(nucleus_bin));
    fi2 = zeros(size(nucleus_bin));
    for I_bin = 1:length(nucleus_bin)
        fi0(I_bin) = mean(nucleus_RNA_profile((nucleus_distance >= bin_min(I_bin))&(nucleus_distance <= bin_max(I_bin)),3));
        fi1(I_bin) = std(nucleus_RNA_profile((nucleus_distance >= bin_min(I_bin))&(nucleus_distance <= bin_max(I_bin)),3));
        fi2(I_bin) = fi1(I_bin)/sqrt(nnz((nucleus_distance >= bin_min(I_bin))&(nucleus_distance <= bin_max(I_bin))));
    end
    FISHN_out = [nucleus_bin',fi0',fi1'];
    good_bin = find(~isnan(fi1));
    
%     plot(nucleus_RNA_profile(:,1),nucleus_RNA_profile(:,3),'Color',[0.7,0.7,0.7],'Marker','.','LineStyle','none'); hold on
%     
%     errorbar(nucleus_bin(good_bin),fi0(good_bin),fi2(good_bin),'k')
%     
%     xpatch = [nucleus_bin(good_bin),nucleus_bin(good_bin(end:-1:1))];
%     ypatch = [fi0(good_bin)-fi1(good_bin),fi0(good_bin(end:-1:1))+fi1(good_bin(end:-1:1))];
%     patch(xpatch,ypatch,'r','FaceAlpha',0.5,'EdgeColor','none');
%     
%     title([image_folder,', cycle = ',num2str(N_cycle)],'Interpreter','none')
%     xlabel('EL')
%     ylabel('# of foci')
%     legend('Nucleus foci number','Mean nucleus foci number','Std of nucleus foci number')
%     legend('hide')

    %%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ELg = nucleus_bin(good_bin & nucleus_bin >= EL_min & nucleus_bin <= EL_max);
    fi0g = fi0(good_bin & nucleus_bin >= EL_min & nucleus_bin <= EL_max);
    fi1g = fi1(good_bin & nucleus_bin >= EL_min & nucleus_bin <= EL_max);
    fi2g = fi2(good_bin & nucleus_bin >= EL_min & nucleus_bin <= EL_max);
    
    [max_FISHN,max_I] = max(fi0g);
    [min_FISHN,min_I] = min(fi0g(max_I:end));
    min_I = max_I+min_I-1;
    max_ELN = ELg(max_I);
    max_stdN = fi1g(max_I);
    max_errN = fi2g(max_I);
    mid_FISHN = (max_FISHN+min_FISHN)/2;
    [~,I0] = min(abs(fi0g(max_I:min_I)-mid_FISHN));
    p0 = polyfit(ELg((I0+max_I-2):(I0+max_I)),fi0g((I0+max_I-2):(I0+max_I)),1);
    mid_ELN = (mid_FISHN-p0(2))/p0(1);
    
    FISHN_out2 = [max_FISHN,max_stdN,max_errN,max_ELN,mid_ELN];
    %%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
% figure(70)
% clf
%     plot(nucleus_RNA_profile(:,1),nucleus_RNA_profile(:,3),'Color',[0.7,0.7,0.7],'Marker','.','LineStyle','none'); hold on
%     errorbar(nucleus_bin(good_bin),fi0(good_bin),fi2(good_bin),'k')
%     patch(xpatch,ypatch,'r','FaceAlpha',0.5,'EdgeColor','none');
%     
%     title([image_folder,', cycle = ',num2str(N_cycle)],'Interpreter','none')
%     xlabel('EL')
%     ylabel('# of foci')
%     legend('Nucleus foci number','Mean nucleus foci number','Std of nucleus foci number')
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% %% Re-plot FISH spot #/nuclei vs EL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nfoci_nucleus = nucleus_RNA_profile(:,3);
% figure(6)
% subaxis(sub_pos(1),sub_pos(2),sub_pos(3), 'Spacing', 0, 'PaddingRight',0.01, 'PaddingLeft',0.02, 'PaddingTop',0.045, 'PaddingBottom',0.04, 'Margin', 0.01);
%     fn0 = zeros(size(nucleus_bin));
%     fn1 = zeros(size(nucleus_bin));
%     fn2 = zeros(size(nucleus_bin));
%     fn3 = zeros(size(nucleus_bin));
%     for I_bin = 1:length(nucleus_bin)
%         fn0(I_bin) = sum((nucleus_distance(Nfoci_nucleus == 0) >= bin_min(I_bin)).*(nucleus_distance(Nfoci_nucleus == 0) <= bin_max(I_bin)));
%         fn1(I_bin) = sum((nucleus_distance(Nfoci_nucleus == 1) >= bin_min(I_bin)).*(nucleus_distance(Nfoci_nucleus == 1) <= bin_max(I_bin)));
%         fn2(I_bin) = sum((nucleus_distance(Nfoci_nucleus == 2) >= bin_min(I_bin)).*(nucleus_distance(Nfoci_nucleus == 2) <= bin_max(I_bin)));
%         fn3(I_bin) = sum((nucleus_distance(Nfoci_nucleus  > 2) >= bin_min(I_bin)).*(nucleus_distance(Nfoci_nucleus  > 2) <= bin_max(I_bin)));
%     end
%     %fn0 = hist(nucleus_distance(Nfoci_nucleus == 0), nucleus_bin);
%     %fn1 = hist(nucleus_distance(Nfoci_nucleus == 1), nucleus_bin);
%     %fn2 = hist(nucleus_distance(Nfoci_nucleus == 2), nucleus_bin);
%     %fn3 = hist(nucleus_distance(Nfoci_nucleus > 2), nucleus_bin);
%     fn_all = fn0+fn1+fn2+fn3+(fn0+fn1+fn2+fn3 == 0);
%     plot(nucleus_bin,(1-fn0./fn_all)*100,'m',nucleus_bin,fn1./fn_all*100,'b',nucleus_bin,fn2./fn_all*100,'g',nucleus_bin,fn3./fn_all*100,'r')
%     title([image_folder,', cycle = ',num2str(N_cycle)],'Interpreter','none')
%     xlabel('EL')
%     ylabel('%')
%     legend('nuclei with active foci','nuclei with 1 foci','nuclei with 2 foci','nuclei with >3 foci')
%     legend('hide')
% %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% %% Re-plot foci Intensity distribution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% foci_intensity = foci_RNA_profile(:,3);
% figure(8)
% subaxis(sub_pos(1),sub_pos(2),sub_pos(3), 'Spacing', 0, 'PaddingRight',0.01, 'PaddingLeft',0.02, 'PaddingTop',0.045, 'PaddingBottom',0.04, 'Margin', 0.01);
%     [Ifoci_n,xout] = hist(foci_intensity,intensity_Nbin);
%     bar(xout,Ifoci_n/(sum(Ifoci_n)+(sum(Ifoci_n) == 0))*100)
%     title([image_folder,', cycle = ',num2str(N_cycle)],'Interpreter','none')
%     xlabel('Equavalent # of mRNA')
%     ylabel('%')
% %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% %% FISH spot #/nuclei heat map %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(11)
% subaxis(sub_pos(1),sub_pos(2),sub_pos(3), 'Spacing', 0, 'PaddingRight',0.01, 'PaddingLeft',0.02, 'PaddingTop',0.045, 'PaddingBottom',0.04, 'Margin', 0.01);
%     im0 = logical(ismember(mask_stack,find(Nfoci_nucleus == 0)));
%     im1 = logical(ismember(mask_stack,find(Nfoci_nucleus == 1)));
%     imfate(:,:,3) = max(double(im0 | im1),[],3);
%     clear im1
%     im2 = logical(ismember(mask_stack,find(Nfoci_nucleus == 2)));
%     imfate(:,:,2) = max(double(im0 | im2),[],3);
%     clear im2
%     im3 = logical(ismember(mask_stack,find(Nfoci_nucleus > 2)));
%     imfate(:,:,1) = max(double(im0 | im3),[],3);
%     clear im3
%     imshow(imfate);
%     title([image_folder,' (cycle = ',num2str(N_cycle),char(10),', white: 0 foci, blue: 1 foci, green: 2 foci, red: > 2 foci)'],'Interpreter','none');
% 
%     
% figure(110)
%     clf
%     imshow(imfate);
%     title([image_folder,' (cycle = ',num2str(N_cycle),char(10),', white: 0 foci, blue: 1 foci, green: 2 foci, red: > 2 foci)'],'Interpreter','none');
% %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
% 
% 
% %% FISH spot mRNA #/nuclei heat map %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(10)
% subaxis(sub_pos(1),sub_pos(2),sub_pos(3), 'Spacing', 0, 'PaddingRight',0.01, 'PaddingLeft',0.02, 'PaddingTop',0.045, 'PaddingBottom',0.04, 'Margin', 0.01);
%     Ifoci_nucleus = [0;Ifoci_nucleus];
%     imfateI0 = Ifoci_nucleus(mask_stack+1);
%     mask2Dsum = sum(logical(mask_stack),3);
%     imfateI = sum(imfateI0,3)./mask2Dsum;
%     imfateI(~mask2Dsum) = 0;
%     
%     patch([0.5,size(imfateI,2)+0.5,size(imfateI,2)+0.5,0.5],[0.5,0.5,size(imfateI,1)+0.5,size(imfateI,1)+0.5],'k','EdgeColor','none');
%     hold on
%     temp = imagesc(imfateI);
%     alpha(temp,mask2Dsum)
% %     set(gca,'Color','k')
%     axis off
%     axis equal
%     colorbar
%     title([image_folder,' (cycle = ',num2str(N_cycle),')'],'Interpreter','none');
% 
%     
% figure(100)
%     clf
%     patch([0.5,size(imfateI,2)+0.5,size(imfateI,2)+0.5,0.5],[0.5,0.5,size(imfateI,1)+0.5,size(imfateI,1)+0.5],'k','EdgeColor','none');
%     hold on
%     temp = imagesc(imfateI);
%     alpha(temp,mask2Dsum)
% %     set(gca,'Color','k')
%     axis off
%     axis equal
%     colorbar
%     title([image_folder,' (cycle = ',num2str(N_cycle),')'],'Interpreter','none');
% %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% 
