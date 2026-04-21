clear all
close all

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A program to quantify the crosstalk between different channels %%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Parameter setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
list_name = 'RNA2list.xls';
in_folder = 'stacks/';
input_name = 'matchlist.xls';
image_type = '*.tif';
out_folder = 'Results/';
mask_folder = 'masks/';
mask_name = 'mask.mat';
mismatch_name = 'mismatch.xls';
mismatch_name2 = 'mismatch_60X.xls';
xymismatch_name = 'xymismatch.mat';
xymismatch_name2 = 'xymismatch_60X.mat';
hist_folder = 'Histogram/';
hist_folder2 = 'Histogram_RNA2/';
fit_folder = 'Histogram_A/';
fit_folder2 = 'Histogram_A_RNA2/';
N_thresh = 3;
n_dilate = 10;
dmax0 = 7; %%% maximal xy mismatch
fit_add = '_spot_fit';
hist_tail = '_raw.xls';
output_tail = '.xls';
figure_tail = '.fig';
figure_tail2 = '.png';
mat_tail = '.mat';
seg_add = '_seg';
th_add = '_th';
foci_add = '_foci';
fish_add = '_fish';
num_add = '_num';
int_add = '_int';
cmp_add = '_cmp';
sel_add = '_sel';
nu_add = '_nucleus';
cyto_add = '_cytoplasm';
fate_add = '_fate';
RNA_add = '_RNA';
signal2_add = '_RNA2';
D3_add = '_3D';
DAPI_add = '_DAPI';
noise_add = '_noise';
mismatch_add = '_mismatch';
z_add = '_z';
xy_add= '_xy';
fit_add2 = '_fit';
FISH_add = '_FISH';
single_add = '_single';
marker_add = '_marker';
pixel_add = '_pixel';
total_add = '_total';
Nbin0 = 15;

I_old_new = 3;
channel_in = 1;
channel_out = 2;
spot_folder = 'Histogram_RNA2/';
% channel_in = 2;
% channel_out = 3;
% spot_folder = 'Histogram/';
dxyz_RNA = [7,7,0];
dxy_max = 10;
dxy0 = {[0,0],[0,0]};
z_range = 8:20;

sub_ij = [4,6];
corr_add = '_corr';
cross_add = '_cross';
mismatch_add = '_0mis';
intercept_add = '_0int';

sub_dim = [3,4];
xlim0 = [0,6e4];

r_th = [0.8*ones(1,11),2];
r_out = [4,10];
r_Inten = 0.7;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Data loading: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[num_list, folder_list] = xlsread(list_name);
folder_list = folder_list(strcmpi('T',folder_list(:,6)),:);
[N1,N2] = size(folder_list);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ii = 0;
%% Image analysis: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for list_I = 1:N1
    all_foci = cell(0);
    all_foci_xy = cell(0);
    all_single = cell(0);
    [sub_num, sub_list] = xlsread([folder_list{list_I,1},in_folder,input_name]);
    channel_name = eval(folder_list{list_I,5});
    [M1,M2] = size(sub_list);
    
    for list_J = 1:M1
        ii = ii+1;
        image_folder = [folder_list{list_I,1},in_folder,sub_list{list_J,3}];
        result_folder = [folder_list{list_I,1},out_folder];
        load([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mat_tail],'max_image');
        fname1 = dir([image_folder,image_type]);
        Lxyz = [size(max_image,1),size(max_image,2),length(fname1)];
        
        all_color = eval(folder_list{list_I,5});
        in_color = all_color(channel_in);
        out_color = all_color(channel_out);
        
        p_all = cell(length(channel_in),length(channel_out));
        p_all2 = cell(length(channel_in),length(channel_out));
        p_all3 = cell(length(channel_in),length(channel_out));
        p_all4 = cell(length(channel_in),length(channel_out));
        x0_all = zeros(0);
        y0_all = zeros(0);
        x1_all = zeros(0);
        y1_all = zeros(0);
        
        for ic = 1:length(channel_in)
            for jc = 1:length(channel_out)
        %         mask0 = im2bw(max_image(:,:,channel_in),1*graythresh(max_image(:,:,channel_in)));
        %         mask0 = imerode(mask0,strel('disk',2));
                mask0 = im2bw(max_image(:,:,channel_in(ic)),r_th(ii)*graythresh(max_image(:,:,channel_in(ic))));
                mask0 = imopen(mask0,strel('disk',2));
                mask1 = bwlabel(mask0);
                mask2 = dilate3D(mask1,r_out(2))-dilate3D(mask1,r_out(1)); mask2(mask2 < 0) = 0;
                
                props1 = regionprops(mask0,max_image(:,:,channel_in(ic)),'MeanIntensity'); Imean_max = [props1.MeanIntensity];
                props2 = regionprops(mask2,max_image(:,:,channel_in(ic)),'MeanIntensity'); Imean_min = [props2.MeanIntensity];
                Imean_th = (Imean_max-Imean_min)*r_Inten+Imean_min;

                h4 = figure;
                    a(:,:,1) = max_image(:,:,channel_in(ic));
                    a(:,:,2) = bwperim(mask0)*65535;
                    a(:,:,3) = 0;
                    imshow(a)
        
                h0 = figure; maximize(h0); 
                p_all0 = zeros(length(fname1),2);

                for I_layer = 1:length(fname1)
                    im0 = imread([image_folder,fname1(I_layer).name]);
                    
                    I_in = im0(:,:,channel_in(ic)); x0 = double(I_in(mask0));
                    I_out = im0(:,:,channel_out(jc)); y0 = double(I_out(mask0));

                    props3 = regionprops(mask0,I_in,'MeanIntensity','MaxIntensity'); Imean3 = [props3.MeanIntensity]; x1 = Imean3';
                    I3_nosaturate = [props3.MaxIntensity] < 65535;
                    props4 = regionprops(mask0,I_out,'MeanIntensity','MaxIntensity'); Imean4 = [props4.MeanIntensity]; y1 = Imean4';
                    
                    mask3 = ismember(mask1,find(Imean3 >= Imean_th));
                    x0_all = cat(1,x0_all,double(I_in(mask3)));
                    y0_all = cat(1,y0_all,double(I_out(mask3)));
                    x1_all = cat(1,x1_all,x1(I3_nosaturate));
                    y1_all = cat(1,y1_all,y1(I3_nosaturate));
                    
                    p0 = polyfit(x0,y0,1);

        % %             r0 = pca([x0,y0]);
        % %             p0 = [r0(2,1)/r0(1,1),mean(y0)-r0(2,1)/r0(1,1)*mean(x0)];

        % %             p0 = linortfit2(x0,y0);

% % % % % %                     p0 = [std(y0)/std(x0),0];

        % %             cov0 = cov(x0,y0);
        % %             p0 = [cov0(1,2)/std(x0)^2,0];
        
                    p1 = polyfit(x1(I3_nosaturate),y1(I3_nosaturate),1);

                    subplot(sub_dim(1),sub_dim(2),I_layer)
                        plot(x0,y0,'.','Color',0.8*ones(1,3),'DisplayName','Data')
                        hold on
                        plot(x1(I3_nosaturate),y1(I3_nosaturate),'.k','DisplayName','Mean data')
                        plot(xlim0,polyval(p0,xlim0),'r','DisplayName',['Fit: y = ',num2str(p0(1)),' * x + ',num2str(p0(2))])
                        plot(xlim0,polyval(p1,xlim0),'g','DisplayName',['Mean fit: y = ',num2str(p1(1)),' * x + ',num2str(p1(2))])
                        xlabel(in_color{ic})
                        ylabel(out_color{jc})
                        title(['z = ',num2str(I_layer),', <y>/<x> = ',num2str(mean(y0)/mean(x0))])
%                         xlim(xlim0)
                        legend('show')

                    p_all0(I_layer,:) = p0;
                end
                
                p_all{ic,jc} = mean(p_all0,1);
                p_all2{ic,jc} = p_all0;
                set(h0,'Name',[result_folder,sub_list{list_J,3},': k = ',num2str(p_all{ic,jc}(1))]);
                
%                 p0all = [std(y0_all)/std(x0_all),0];
                p0all = polyfit(x0_all,y0_all,1);
                p1all = polyfit(x1_all,y1_all,1);
                
                h5 = figure;
                    plot(x0_all,y0_all,'.','Color',0.8*ones(1,3),'DisplayName','Data')
                    hold on
                    plot(x1_all,y1_all,'.k','DisplayName','Mean data')
                    plot(xlim0,polyval(p0all,xlim0),'r','DisplayName',['Fit: y = ',num2str(p0all(1)),' * x + ',num2str(p0all(2))])
                    plot(xlim0,polyval(p1all,xlim0),'g','DisplayName',['Mean fit: y = ',num2str(p1all(1)),' * x + ',num2str(p1all(2))])
                    xlabel(in_color{ic})
                    ylabel(out_color{jc})
                    title(['All bright spots, <y>/<x> = ',num2str(mean(y0_all)/mean(x0_all))])
%                         xlim(xlim0)
                    legend('show')
                    
                p_all3{ic,jc} = p0all;
                p_all4{ic,jc} = p1all;
%%% Output: %%%============================================================
                saveas(h0,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),'_',in_color{ic},'_',out_color{jc},cross_add,pixel_add,intercept_add,figure_tail]);
                saveas(h0,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),'_',in_color{ic},'_',out_color{jc},cross_add,pixel_add,intercept_add,figure_tail2]);
                saveas(h4,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),'_',in_color{ic},foci_add,intercept_add,figure_tail]);
                saveas(h4,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),'_',in_color{ic},foci_add,intercept_add,figure_tail2]);
                saveas(h5,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),'_',in_color{ic},foci_add,total_add,intercept_add,figure_tail]);
                saveas(h5,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),'_',in_color{ic},foci_add,total_add,intercept_add,figure_tail2]);
% % % % % %         saveas(h2,[result_folder,sub_list{list_J,I_old_new}(1:(length(sub_list{list_J,I_old_new})-1)),xy_add,corr_add,D3_add,figure_tail]);
% % % % % %         saveas(h2,[result_folder,sub_list{list_J,I_old_new}(1:(length(sub_list{list_J,I_old_new})-1)),xy_add,corr_add,D3_add,figure_tail2]);
% % %         saveas(h3,[result_folder,sub_list{list_J,I_old_new}(1:(length(sub_list{list_J,I_old_new})-1)),cross_add,mismatch_add,intercept_add,figure_tail]);
% % %         saveas(h3,[result_folder,sub_list{list_J,I_old_new}(1:(length(sub_list{list_J,I_old_new})-1)),cross_add,mismatch_add,intercept_add,figure_tail2]);
            end
        end
        if exist([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mismatch_add,intercept_add,mat_tail])
            save([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mismatch_add,intercept_add,mat_tail],'p_all','p_all2','p_all3','p_all4','all_color','channel_in','channel_out','-append');
        else
            save([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mismatch_add,intercept_add,mat_tail],'p_all','p_all2','p_all3','p_all4','all_color','channel_in','channel_out');
        end
    end
end
