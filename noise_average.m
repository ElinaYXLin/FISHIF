clear all
close all

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A program to compare the shot noise level of different averaging %%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Parameter setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_folder = '05282011b1/';
im_sub = 'stacks/';
mask_sub = 'masks/';
data_name = {'04272011_9_IF_prewashed_lapse_001_A','04272011_9_IF_prewashed_lapse_002_A'};
name_add = '_new';
im_name = 'stack01.tif';
mask_name = 'mask.mat';
I_channel = 2;
r_edge = 10;

out_folder = 'Results/';
out_name = 'noise_level';
fig_tail = '.fig';
fig_tail2 = '.eps';
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
k_all = zeros(0);
m_all = zeros(0);
for ii = 1:length(data_name)
%% Data loading: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fname = dir([data_folder,im_sub,data_name{ii},'*',name_add]);
    
    %%% load/average IF images:
    im0 = zeros(0);
    ima = zeros(0);
    for jj = 1:length(fname)
        temp = imread([data_folder,im_sub,fname(jj).name,'/',im_name]);
        im0 = cat(3,im0,temp(:,:,I_channel));
        ima = cat(3,ima,mean(im0,3));
    end
    
    %%% load nuclear mask:
    load([data_folder,mask_sub,fname(1).name,'/',mask_name]); 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Data processing: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% calculate mean and variance:
    mask_small = (imerode(logical(mask_stack), strel('disk',r_edge))).*mask_stack;
    mean_IF = nan(size(ima,3),max(mask_small(:)));
    var_IF = nan(size(ima,3),max(mask_small(:)));
    k_IF = zeros(1,size(ima,3));
    m_IF = zeros(1,size(ima,3));
    for I_layer = 1:size(ima,3)
        sg_prop = regionprops(mask_small,ima(:,:,I_layer),'MeanIntensity');
        sg2_prop = regionprops(mask_small,ima(:,:,I_layer).^2,'MeanIntensity');
        mean_IF(I_layer,:) = [sg_prop.MeanIntensity];
        var_IF(I_layer,:) = [sg2_prop.MeanIntensity]-[sg_prop.MeanIntensity].^2;
        I0 = ~isnan(mean_IF(I_layer,:));
        p = polyfit(mean_IF(I_layer,I0),var_IF(I_layer,I0),1);
        k_IF(I_layer) = p(1);
        m_IF(I_layer) = mean(mean_IF(I_layer,I0));
    end
    
    k_all = cat(1,k_all,k_IF/k_IF(end)-1);
    m_all = cat(1,m_all,1-m_IF/m_IF(1));
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Data plot: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(1,length(data_name),ii)
        plot(1:size(ima,3),k_IF/k_IF(end)-1,'r',1:size(ima,3),1-m_IF/m_IF(1),'b')
        xlabel('# average')
        ylabel('\Delta')
        title(data_name{ii})
        ylim([0,1])
        legend('k','mean IF')
end

figure
    errorbar(1:size(ima,3),mean(k_all),std(k_all),'LineStyle','--','Marker','.','Color','r')
    hold on
    errorbar(1:size(ima,3),mean(m_all),std(m_all),'LineStyle','--','Marker','.','Color','b')
    xlim([0,17])
    ylim([0,0.5])
    set(gca,'Units','inches','Position',[1,1,1,1],'FontName','Helvetica','FontSize',7,'XTick',[0:8:16])
    xlabel('# average','FontName','Helvetica','FontSize',8)
    ylabel('\Delta','FontName','Helvetica','FontSize',8)
    legend('k','IF')
    box off
saveas(gcf,[data_folder,out_folder,data_folder(1:end-1),out_name,fig_tail])
export_fig([data_folder,out_folder,data_folder(1:end-1),out_name,fig_tail2],gcf)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    
    
    
    