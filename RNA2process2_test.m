clear

%% Parameter setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
list_name = 'RNA2list.xls';
in_folder = 'stacks/';
input_name = 'matchlist.xls';
image_type = '*.tif';
out_folder = 'Results/';
hist_folder = 'Histogram/';
hist2_folder = 'Histogram_RNA2/';
hist_folder_single = 'Histogram/';
hist2_folder_single = 'Histogram_RNA2/';
fit_tail = '_raw';
output_tail = '.xls';
figure_tail = '.fig';
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
mismatch_add = '_mismatch';
z_add = '_z';
xy_add= '_xy';
fit_add = '_fit';
FISH_add = '_FISH';
single_add = '_single';
Nbin0 = 15;
dmax = 7;
hist_bin0 = [-5:5];   %%% histogram bin for z distance
dye_name = {'DAPI','A488','TMR','A647'};

mismatch_folder = 'mismatch_all/';
xymismatch_name = 'xymismatch_60X_FA';
mismatch_name = 'mismatch_60X_FA';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Data loading: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[num_list, folder_list] = xlsread(list_name);
folder_list = folder_list(strcmpi('T',folder_list(:,6)),:);
[N1,N2] = size(folder_list);

% % mismatch_matrix = cell(length(dye_name)+1);
% % mismatch_matrix(1,2:end) = dye_name;
% % mismatch_matrix(2:end,1) = dye_name';
% % mismatch_matrix(2:end,2:end) = num2cell(zeros(length(dye_name)));
% % 
% % if exist(mismatch_folder) ~= 7
% %     mkdir(mismatch_folder);
% % end
% % 
% % save_I = false;

[~,~,mismatch_matrix] = xlsread([mismatch_folder,mismatch_name,output_tail]);
load([mismatch_folder,xymismatch_name,mat_tail]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Image analysis: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for list_I = 1:N1
    all_foci = cell(0);
    all_single = cell(0);
    all_size = cell(0);
    [sub_num, sub_list] = xlsread([folder_list{list_I,1},in_folder,input_name]);
    channel_name = eval(folder_list{list_I,5});
    [M1,M2] = size(sub_list);
    for list_J = 1:M1
        result_folder = [folder_list{list_I,1},out_folder];
        load([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mat_tail]);
        image_folder = [folder_list{list_I,1},in_folder,sub_list{list_J,3}];
        WGA_channel = sub_num(list_J,6);
        DAPI_channel = sub_num(list_J,7);
% %         RNA_channel = sub_num(list_J,10);
% %         RNA2_channel = sub_num(list_J,8);
        RNA_channel = sub_num(list_J,8);
        RNA2_channel = sub_num(list_J,10);
        resolution = sub_num(list_J,9);
        
        dz = mismatch_matrix{strcmp(channel_name{RNA_channel},mismatch_matrix(1,:)),strcmp(channel_name{RNA2_channel},mismatch_matrix(:,1))};
        M0 = eval([channel_name{RNA_channel},'_',channel_name{RNA2_channel}]);
        C0 = eval([channel_name{RNA_channel},'_',channel_name{RNA2_channel},'_con']);
        size0 = size(max_image);

        RNA_fit = xlsread([folder_list{list_I,1},hist_folder,sub_list{list_J,3}(1:end-1),fit_tail,output_tail]);
        RNA2_fit = xlsread([folder_list{list_I,1},hist2_folder,sub_list{list_J,3}(1:end-1),fit_tail,output_tail]);
        RNA_fit_single = xlsread([folder_list{list_I,1},hist_folder_single,sub_list{list_J,3}(1:end-1),fit_tail,output_tail]);
        RNA2_fit_single = xlsread([folder_list{list_I,1},hist2_folder_single,sub_list{list_J,3}(1:end-1),fit_tail,output_tail]);
        
        x00 = repmat((1+size0(1:2))/2,size(RNA2_fit,1),1);
        RNA2_fit(:,8) = RNA2_fit(:,8)+dz;
        RNA2_fit(:,6:7) = (RNA2_fit(:,6:7)-x00)*M0'+C0'+x00;
        
        x00 = repmat((1+size0(1:2))/2,size(RNA2_fit_single,1),1);
        RNA2_fit_single(:,8) = RNA2_fit_single(:,8)+dz;
        RNA2_fit_single(:,6:7) = (RNA2_fit_single(:,6:7)-x00)*M0'+C0'+x00;
        
        foci_mismatch = z_mismatch(RNA_fit,RNA2_fit,channel_name{RNA_channel},channel_name{RNA2_channel},size(seg_bw),image_folder,dmax);
        all_foci = cat(1,all_foci,{foci_mismatch});
        
        single_mismatch = mismatch_single(RNA_fit_single,RNA2_fit_single,channel_name{RNA_channel},channel_name{RNA2_channel},size(seg_bw),image_folder,dmax);
        all_single = cat(1,all_single,{single_mismatch});
%         [length(RNA_fit_single),length(RNA2_fit_single),length(single_mismatch)]

%         all_size = cat(1,{size(max_image)});
%%% Output: %%%============================================================
% %         if exist(result_folder) ~= 7
% %             mkdir(result_folder);
% %         end
% %         saveas(91,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mismatch_add,z_add,figure_tail]);
% %         saveas(92,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mismatch_add,xy_add,figure_tail]);
% %         saveas(93,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mismatch_add,FISH_add,single_add,figure_tail]);
% %         saveas(94,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mismatch_add,xy_add,single_add,figure_tail]);
% %         save([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mat_tail],'foci_mismatch','single_mismatch','-append');
% %         
% %         if ~isempty(foci_mismatch)
% %             xlswrite([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mismatch_add,z_add,output_tail],foci_mismatch);
% %         end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        
    end


%% Overall foci plot/analysis: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    foci1 = cell2mat(all_foci);
    imdx = zeros(size(seg_bw));
    imdy = zeros(size(seg_bw));
    for I0 = 1:length(foci1)
        imdx(round(foci1(I0,3)),round(foci1(I0,4))) = foci1(I0,5);
        imdy(round(foci1(I0,3)),round(foci1(I0,4))) = foci1(I0,6);
    end

    figure(1)
    subplot(1,2,1)
        Nfil = conv2(double(imdx ~= 0),double(getnhood(strel('disk',40))),'same');
        imfil = conv2(imdx,double(getnhood(strel('disk',40))),'same')./(Nfil+(~Nfil));
        imagesc(imfil')
        colorbar
        axis equal
        axis off
        set(gca,'xtick',[],'ytick',[],'FontName','Arial','FontSize',16,'FontWeight','bold')
%         xlabel('Y','FontName','Arial','FontSize',16,'FontWeight','bold')
%         ylabel('X','FontName','Arial','FontSize',16,'FontWeight','bold')
        title(['X(',channel_name{RNA_channel},') - X(',channel_name{RNA2_channel},'): ',folder_list{list_I,1}],'FontName','Arial','FontSize',16,'FontWeight','bold')

    subplot(1,2,2)
        Nfil = conv2(double(imdy ~= 0),double(getnhood(strel('disk',40))),'same');
        imfil = conv2(imdy,double(getnhood(strel('disk',40))),'same')./(Nfil+(~Nfil));
        imagesc(imfil')
        colorbar
        axis equal
        axis off
        set(gca,'xtick',[],'ytick',[],'FontName','Arial','FontSize',16,'FontWeight','bold')
%         xlabel('Y','FontName','Arial','FontSize',16,'FontWeight','bold')
%         ylabel('X','FontName','Arial','FontSize',16,'FontWeight','bold')
        title(['Y(',channel_name{RNA_channel},') - Y(',channel_name{RNA2_channel},'): ',folder_list{list_I,1}],'Interpreter','none','FontName','Arial','FontSize',16,'FontWeight','bold')

% %     saveas(1,[result_folder,mismatch_add,xy_add,figure_tail]);

    
    figure(2)
    subplot(2,2,1)
%         p0 = polyfit(foci1(foci1(:,3) <= 1354,3),foci1(foci1(:,3) <= 1354,5),1);
        p0 = polyfit(foci1(:,3),foci1(:,5),1);
        lx = [min(foci1(:,3)),max(foci1(:,3))];
        ly = p0(1)*lx+p0(2);
        %%%%%
        [xbin,ybin,xerr,yerr] = equal_bin(foci1(:,3),foci1(:,5),Nbin0);
        pbin = polyfit(xbin,ybin,1);
        lxbin = [min(foci1(:,3)),max(foci1(:,3))];
        lybin = pbin(1)*lx+pbin(2);
        %%%%%
        plot(foci1(:,3),foci1(:,5),'b.',lx,ly,'c-')
        hold on
        errorbarxy(xbin,ybin,xerr,yerr,xerr,yerr,'r','r','.')
        plot(lxbin,lybin,'m-')
        legend('Raw data',['Linear fit: dx = ',num2str(p0(1)),' * x + ',num2str(p0(2))],'binned data',['Linear fit: dx = ',num2str(pbin(1)),' * x + ',num2str(pbin(2))])
        set(gca,'FontName','Arial','FontSize',16,'FontWeight','bold')
        xlabel('X (pixel)','FontName','Arial','FontSize',16,'FontWeight','bold')
        ylabel(['X(',channel_name{RNA_channel},') - X(',channel_name{RNA2_channel},')'],'FontName','Arial','FontSize',16,'FontWeight','bold')
        title(['X(',channel_name{RNA_channel},') - X(',channel_name{RNA2_channel},') vs X : ',folder_list{list_I,1}],'FontName','Arial','FontSize',16,'FontWeight','bold')
    subplot(2,2,2)
        p0 = polyfit(foci1(:,4),foci1(:,5),1);
        lx = [min(foci1(:,4)),max(foci1(:,4))];
        ly = p0(1)*lx+p0(2);
        %%%%%
        [xbin,ybin,xerr,yerr] = equal_bin(foci1(:,4),foci1(:,5),Nbin0);
        pbin = polyfit(xbin,ybin,1);
        lxbin = [min(foci1(:,4)),max(foci1(:,4))];
        lybin = pbin(1)*lx+pbin(2);
        %%%%%
        plot(foci1(:,4),foci1(:,5),'b.',lx,ly,'c-')
        hold on
        errorbarxy(xbin,ybin,xerr,yerr,xerr,yerr,'r','r','.')
        plot(lxbin,lybin,'m-')
        legend('Raw data',['Linear fit: dx = ',num2str(p0(1)),' * y + ',num2str(p0(2))],'binned data',['Linear fit: dx = ',num2str(pbin(1)),' * x + ',num2str(pbin(2))])
        set(gca,'FontName','Arial','FontSize',16,'FontWeight','bold')
        xlabel('Y (pixel)','FontName','Arial','FontSize',16,'FontWeight','bold')
        ylabel(['X(',channel_name{RNA_channel},') - X(',channel_name{RNA2_channel},')'],'FontName','Arial','FontSize',16,'FontWeight','bold')
        title(['X(',channel_name{RNA_channel},') - X(',channel_name{RNA2_channel},') vs Y : ',folder_list{list_I,1}],'FontName','Arial','FontSize',16,'FontWeight','bold')
    subplot(2,2,3)
        p0 = polyfit(foci1(:,3),foci1(:,6),1);
        lx = [min(foci1(:,3)),max(foci1(:,3))];
        ly = p0(1)*lx+p0(2);
        %%%%%
        [xbin,ybin,xerr,yerr] = equal_bin(foci1(:,3),foci1(:,6),Nbin0);
        pbin = polyfit(xbin,ybin,1);
        lxbin = [min(foci1(:,3)),max(foci1(:,3))];
        lybin = pbin(1)*lx+pbin(2);
        %%%%%
        plot(foci1(:,3),foci1(:,6),'b.',lx,ly,'c-')
        hold on
        errorbarxy(xbin,ybin,xerr,yerr,xerr,yerr,'r','r','.')
        plot(lxbin,lybin,'m-')
        legend('Raw data',['Linear fit: dy = ',num2str(p0(1)),' * x + ',num2str(p0(2))],'binned data',['Linear fit: dx = ',num2str(pbin(1)),' * x + ',num2str(pbin(2))])
        set(gca,'FontName','Arial','FontSize',16,'FontWeight','bold')
        xlabel('X (pixel)','FontName','Arial','FontSize',16,'FontWeight','bold')
        ylabel(['Y(',channel_name{RNA_channel},') - Y(',channel_name{RNA2_channel},')'],'FontName','Arial','FontSize',16,'FontWeight','bold')
        title(['Y(',channel_name{RNA_channel},') - Y(',channel_name{RNA2_channel},') vs X : ',folder_list{list_I,1}],'FontName','Arial','FontSize',16,'FontWeight','bold')
    subplot(2,2,4)
        p0 = polyfit(foci1(:,4),foci1(:,6),1);
        lx = [min(foci1(:,4)),max(foci1(:,4))];
        ly = p0(1)*lx+p0(2);
        %%%%%
        [xbin,ybin,xerr,yerr] = equal_bin(foci1(:,4),foci1(:,6),Nbin0);
        pbin = polyfit(xbin,ybin,1);
        lxbin = [min(foci1(:,4)),max(foci1(:,4))];
        lybin = pbin(1)*lx+pbin(2);
        %%%%%
        plot(foci1(:,4),foci1(:,6),'b.',lx,ly,'c-')
        hold on
        errorbarxy(xbin,ybin,xerr,yerr,xerr,yerr,'r','r','.')
        plot(lxbin,lybin,'m-')
        legend('Raw data',['Linear fit: dy = ',num2str(p0(1)),' * y + ',num2str(p0(2))],'binned data',['Linear fit: dx = ',num2str(pbin(1)),' * x + ',num2str(pbin(2))])
        set(gca,'FontName','Arial','FontSize',16,'FontWeight','bold')
        xlabel('Y (pixel)','FontName','Arial','FontSize',16,'FontWeight','bold')
        ylabel(['Y(',channel_name{RNA_channel},') - Y(',channel_name{RNA2_channel},')'],'FontName','Arial','FontSize',16,'FontWeight','bold')
        title(['Y(',channel_name{RNA_channel},') - Y(',channel_name{RNA2_channel},') vs Y : ',folder_list{list_I,1}],'FontName','Arial','FontSize',16,'FontWeight','bold')
% %     saveas(2,[result_folder,mismatch_add,xy_add,fit_add,figure_tail]);
        
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Overall foci plot/analysis: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    foci1 = cell2mat(all_single);
    imdx = zeros(size(seg_bw));
    imdy = zeros(size(seg_bw));
    for I0 = 1:length(foci1)
        imdx(round(foci1(I0,3)),round(foci1(I0,4))) = foci1(I0,5);
        imdy(round(foci1(I0,3)),round(foci1(I0,4))) = foci1(I0,6);
    end

    figure(3)
    subplot(1,2,1)
        Nfil = conv2(double(imdx ~= 0),double(getnhood(strel('disk',40))),'same');
        imfil = conv2(imdx,double(getnhood(strel('disk',40))),'same')./(Nfil+(~Nfil));
        imagesc(imfil')
        colorbar
        axis equal
        axis off
        set(gca,'xtick',[],'ytick',[],'FontName','Arial','FontSize',16,'FontWeight','bold')
%         xlabel('Y','FontName','Arial','FontSize',16,'FontWeight','bold')
%         ylabel('X','FontName','Arial','FontSize',16,'FontWeight','bold')
        title(['X(',channel_name{RNA_channel},') - X(',channel_name{RNA2_channel},'): ',folder_list{list_I,1}],'FontName','Arial','FontSize',16,'FontWeight','bold')

    subplot(1,2,2)
        Nfil = conv2(double(imdy ~= 0),double(getnhood(strel('disk',40))),'same');
        imfil = conv2(imdy,double(getnhood(strel('disk',40))),'same')./(Nfil+(~Nfil));
        imagesc(imfil')
        colorbar
        axis equal
        axis off
        set(gca,'xtick',[],'ytick',[],'FontName','Arial','FontSize',16,'FontWeight','bold')
%         xlabel('Y','FontName','Arial','FontSize',16,'FontWeight','bold')
%         ylabel('X','FontName','Arial','FontSize',16,'FontWeight','bold')
        title(['Y(',channel_name{RNA_channel},') - Y(',channel_name{RNA2_channel},'): ',folder_list{list_I,1}],'Interpreter','none','FontName','Arial','FontSize',16,'FontWeight','bold')

% %     saveas(3,[result_folder,mismatch_add,xy_add,single_add,figure_tail]);

    
    figure(4)
    subplot(2,2,1)
%         p0 = polyfit(foci1(foci1(:,3) <= 1354,3),foci1(foci1(:,3) <= 1354,5),1);
        p0 = polyfit(foci1(:,3),foci1(:,5),1);
        lx = [min(foci1(:,3)),max(foci1(:,3))];
        ly = p0(1)*lx+p0(2);
        %%%%%
        [xbin,ybin,xerr,yerr] = equal_bin(foci1(:,3),foci1(:,5),Nbin0);
        pbin = polyfit(xbin,ybin,1);
        lxbin = [min(foci1(:,3)),max(foci1(:,3))];
        lybin = pbin(1)*lx+pbin(2);
        %%%%%
        plot(foci1(:,3),foci1(:,5),'b.',lx,ly,'c-')
        hold on
        errorbarxy(xbin,ybin,xerr,yerr,xerr,yerr,'r','r','.')
        plot(lxbin,lybin,'m-')
        legend('Raw data',['Linear fit: dx = ',num2str(p0(1)),' * x + ',num2str(p0(2))],'binned data',['Linear fit: dx = ',num2str(pbin(1)),' * x + ',num2str(pbin(2))])
        set(gca,'FontName','Arial','FontSize',16,'FontWeight','bold')
        xlabel('X (pixel)','FontName','Arial','FontSize',16,'FontWeight','bold')
        ylabel(['X(',channel_name{RNA_channel},') - X(',channel_name{RNA2_channel},')'],'FontName','Arial','FontSize',16,'FontWeight','bold')
        title(['X(',channel_name{RNA_channel},') - X(',channel_name{RNA2_channel},') vs X : ',folder_list{list_I,1}],'FontName','Arial','FontSize',16,'FontWeight','bold')
    subplot(2,2,2)
        p0 = polyfit(foci1(:,4),foci1(:,5),1);
        lx = [min(foci1(:,4)),max(foci1(:,4))];
        ly = p0(1)*lx+p0(2);
        %%%%%
        [xbin,ybin,xerr,yerr] = equal_bin(foci1(:,4),foci1(:,5),Nbin0);
        pbin = polyfit(xbin,ybin,1);
        lxbin = [min(foci1(:,4)),max(foci1(:,4))];
        lybin = pbin(1)*lx+pbin(2);
        %%%%%
        plot(foci1(:,4),foci1(:,5),'b.',lx,ly,'c-')
        hold on
        errorbarxy(xbin,ybin,xerr,yerr,xerr,yerr,'r','r','.')
        plot(lxbin,lybin,'m-')
        legend('Raw data',['Linear fit: dx = ',num2str(p0(1)),' * y + ',num2str(p0(2))],'binned data',['Linear fit: dx = ',num2str(pbin(1)),' * x + ',num2str(pbin(2))])
        set(gca,'FontName','Arial','FontSize',16,'FontWeight','bold')
        xlabel('Y (pixel)','FontName','Arial','FontSize',16,'FontWeight','bold')
        ylabel(['X(',channel_name{RNA_channel},') - X(',channel_name{RNA2_channel},')'],'FontName','Arial','FontSize',16,'FontWeight','bold')
        title(['X(',channel_name{RNA_channel},') - X(',channel_name{RNA2_channel},') vs Y : ',folder_list{list_I,1}],'FontName','Arial','FontSize',16,'FontWeight','bold')
    subplot(2,2,3)
        p0 = polyfit(foci1(:,3),foci1(:,6),1);
        lx = [min(foci1(:,3)),max(foci1(:,3))];
        ly = p0(1)*lx+p0(2);
        %%%%%
        [xbin,ybin,xerr,yerr] = equal_bin(foci1(:,3),foci1(:,6),Nbin0);
        pbin = polyfit(xbin,ybin,1);
        lxbin = [min(foci1(:,3)),max(foci1(:,3))];
        lybin = pbin(1)*lx+pbin(2);
        %%%%%
        plot(foci1(:,3),foci1(:,6),'b.',lx,ly,'c-')
        hold on
        errorbarxy(xbin,ybin,xerr,yerr,xerr,yerr,'r','r','.')
        plot(lxbin,lybin,'m-')
        legend('Raw data',['Linear fit: dy = ',num2str(p0(1)),' * x + ',num2str(p0(2))],'binned data',['Linear fit: dx = ',num2str(pbin(1)),' * x + ',num2str(pbin(2))])
        set(gca,'FontName','Arial','FontSize',16,'FontWeight','bold')
        xlabel('X (pixel)','FontName','Arial','FontSize',16,'FontWeight','bold')
        ylabel(['Y(',channel_name{RNA_channel},') - Y(',channel_name{RNA2_channel},')'],'FontName','Arial','FontSize',16,'FontWeight','bold')
        title(['Y(',channel_name{RNA_channel},') - Y(',channel_name{RNA2_channel},') vs X : ',folder_list{list_I,1}],'FontName','Arial','FontSize',16,'FontWeight','bold')
    subplot(2,2,4)
        p0 = polyfit(foci1(:,4),foci1(:,6),1);
        lx = [min(foci1(:,4)),max(foci1(:,4))];
        ly = p0(1)*lx+p0(2);
        %%%%%
        [xbin,ybin,xerr,yerr] = equal_bin(foci1(:,4),foci1(:,6),Nbin0);
        pbin = polyfit(xbin,ybin,1);
        lxbin = [min(foci1(:,4)),max(foci1(:,4))];
        lybin = pbin(1)*lx+pbin(2);
        %%%%%
        plot(foci1(:,4),foci1(:,6),'b.',lx,ly,'c-')
        hold on
        errorbarxy(xbin,ybin,xerr,yerr,xerr,yerr,'r','r','.')
        plot(lxbin,lybin,'m-')
        legend('Raw data',['Linear fit: dy = ',num2str(p0(1)),' * y + ',num2str(p0(2))],'binned data',['Linear fit: dx = ',num2str(pbin(1)),' * x + ',num2str(pbin(2))])
        set(gca,'FontName','Arial','FontSize',16,'FontWeight','bold')
        xlabel('Y (pixel)','FontName','Arial','FontSize',16,'FontWeight','bold')
        ylabel(['Y(',channel_name{RNA_channel},') - Y(',channel_name{RNA2_channel},')'],'FontName','Arial','FontSize',16,'FontWeight','bold')
        title(['Y(',channel_name{RNA_channel},') - Y(',channel_name{RNA2_channel},') vs Y : ',folder_list{list_I,1}],'FontName','Arial','FontSize',16,'FontWeight','bold')
% %     saveas(4,[result_folder,mismatch_add,xy_add,fit_add,single_add,figure_tail]);
        
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    close all
    
% % %% Mismatch estimation/output: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %     xy1 = zeros(0);
% %     xy2 = zeros(0);
% %     z1 = zeros(0);
% %     z2 = zeros(0);
% %     for ii = 1:length(all_foci)
% %         temp = all_foci{ii};
% %         xy1 = cat(1,xy1,temp(:,8:9)-repmat((all_size{ii}(1:2)+[1,1])/2,size(temp,1),1));
% %         xy2 = cat(1,xy2,temp(:,11:12)-repmat((all_size{ii}(1:2)+[1,1])/2,size(temp,1),1));
% %         z1 = cat(1,z1,temp(:,10));
% %         z2 = cat(1,z2,temp(:,13));
% %     end
% %     n0 = hist(z2-z1,hist_bin0);
% %     [~,Imax] = max(n0);
% %     dz = hist_bin0(Imax);    
% %     mismatch_matrix{strcmp(channel_name{RNA2_channel},mismatch_matrix(1,:)),strcmp(channel_name{RNA_channel},mismatch_matrix(:,1))} = dz;
% %     mismatch_matrix{strcmp(channel_name{RNA_channel},mismatch_matrix(1,:)),strcmp(channel_name{RNA2_channel},mismatch_matrix(:,1))} = -dz;
% %     
% %     fit2Dx = polyfitn(xy1,xy2(:,1),1);
% %     fit2Dy = polyfitn(xy1,xy2(:,2),1);
% %     M0 = [fit2Dx.Coefficients;fit2Dy.Coefficients];
% %     
% %     eval([channel_name{RNA2_channel},'_',channel_name{RNA_channel},' = M0(:,1:2);',])
% %     eval([channel_name{RNA2_channel},'_',channel_name{RNA_channel},'_con = M0(:,3);'])
% %     eval([channel_name{RNA2_channel},'_',channel_name{RNA_channel},'_x0 = all_size{1}(1:2)*resolution;'])
% %     
% %     fit2Dx = polyfitn(xy2,xy1(:,1),1);
% %     fit2Dy = polyfitn(xy2,xy1(:,2),1);
% %     M0 = [fit2Dx.Coefficients;fit2Dy.Coefficients];
% %     
% %     eval([channel_name{RNA_channel},'_',channel_name{RNA2_channel},' = M0(:,1:2);',])
% %     eval([channel_name{RNA_channel},'_',channel_name{RNA2_channel},'_con = M0(:,3);'])
% %     eval([channel_name{RNA_channel},'_',channel_name{RNA2_channel},'_x0 = all_size{1}(1:2)*resolution;'])
% %     
% %     if ~save_I
% %         resolution_mismatch = resolution;
% %         save([mismatch_folder,xymismatch_name,mat_tail],'resolution_mismatch');
% %         save_I = true;
% %     end
% % 
% %     save([mismatch_folder,xymismatch_name,mat_tail],[channel_name{RNA2_channel},'_',channel_name{RNA_channel}],[channel_name{RNA2_channel},'_',channel_name{RNA_channel},'_con'],[channel_name{RNA2_channel},'_',channel_name{RNA_channel},'_x0'],'-append');
% %     save([mismatch_folder,xymismatch_name,mat_tail],[channel_name{RNA_channel},'_',channel_name{RNA2_channel}],[channel_name{RNA_channel},'_',channel_name{RNA2_channel},'_con'],[channel_name{RNA_channel},'_',channel_name{RNA2_channel},'_x0'],'-append');
% % 
% % %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

% % xlswrite([mismatch_folder,mismatch_name,output_tail],mismatch_matrix);


