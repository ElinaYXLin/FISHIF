clear all
close all

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
Nbin0 = 15;

hist_binxy = 0:0.5:10;
% sub_pos = [3,3];
% cycle_pos0 = zeros(1,sub_pos(1));
% cycle_range = [11,12,13];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Data loading: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[num_list, folder_list] = xlsread(list_name);
folder_list = folder_list(strcmpi('T',folder_list(:,6)),:);
[N1,N2] = size(folder_list);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Image analysis: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for list_I = 1:N1
    all_foci = cell(0);
    all_foci_xy = cell(0);
    all_single = cell(0);
    [sub_num, sub_list] = xlsread([folder_list{list_I,1},in_folder,input_name]);
    channel_name = eval(folder_list{list_I,5});
    [M1,M2] = size(sub_list);
    
    for list_J = 1:M1
        if isempty(strfind(sub_list{list_J,3},'_60X'))
            [~,~,mismatch_matrix] = xlsread(mismatch_name);
            load(xymismatch_name);
        else
            [~,~,mismatch_matrix] = xlsread(mismatch_name2);
            load(xymismatch_name2);
        end
        
        image_folder = [folder_list{list_I,1},in_folder,sub_list{list_J,3}];
        result_folder = [folder_list{list_I,1},out_folder];
        load([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mat_tail]);
        
        resolution = sub_num(list_J,9);
        resolutionz = sub_num(list_J,11);
        Nbin = sub_num(list_J,2);
        Mdim = sub_num(list_J,3);
        WGA_channel = sub_num(list_J,6);
        DAPI_channel = sub_num(list_J,7);
% %         RNA_channel = sub_num(list_J,10);
% %         RNA2_channel = sub_num(list_J,8);
        RNA_channel = sub_num(list_J,8);
        RNA2_channel = sub_num(list_J,10);
        
        load([folder_list{list_I,1},mask_folder,sub_list{list_J,3},mask_name]);   %%% load 3D mask
        mask_stack = mask_dilate(mask_stack,n_dilate);
        
        z_size = size(mask_stack,3);
        load([folder_list{list_I,1},fit_folder,sub_list{list_J,3}(1:end-1),fit_add,mat_tail]);   %%% load spot intensity fitting result
        Inten_thresh = b*N_thresh;   %%% set foci intensity threshold
        single_Inten = b;
        load([folder_list{list_I,1},fit_folder2,sub_list{list_J,3}(1:end-1),fit_add,mat_tail]);   %%% load spot intensity fitting result
        Inten_thresh2 = b*N_thresh;   %%% set foci intensity threshold
        single_Inten2 = b;
        N_cycle = round(log2(max(mask_stack(:))))+2;   %%% Calculate the nuclear cycle number
        

%%% RNA analysis: %%%======================================================
        all_color = eval(folder_list{list_I,5});
        RNA_color = all_color{RNA_channel};
        RNA2_color = all_color{RNA2_channel};
        RNA_RNA2_mismatch = mismatch_matrix{strcmp(RNA_color,mismatch_matrix(1,:)),strcmp(RNA2_color,mismatch_matrix(:,1))};
        RNA_RNA2_xymismatch = {eval([RNA_color,'_',RNA2_color]),eval([RNA_color,'_',RNA2_color,'_con']),eval([RNA_color,'_',RNA2_color,'_x0']),Nbin,Mdim,resolution,resolution_mismatch};

        [RNA_stack,RNA2_stack,DAPI_stack] = stack3D(imclearborder(seg_bw),RNA_channel,DAPI_channel,RNA2_channel,image_folder,RNA_RNA2_mismatch);   %%% load 3D image stacks
        [nucleus_DAPI_profile,DNA_mask] = DAPI_profile3D(mask_stack,DAPI_stack,image_folder,N_cycle);
        clear DAPI_stack
        
        [foci_bw3D,max_image00,SS] = modified_foci3D([folder_list{list_I,1},hist_folder,sub_list{list_J,3}(1:(end-1)),hist_tail],max_image,RNA_channel,mask_stack,N_cycle,image_folder,0,Inten_thresh);
        foci_bw2D = max(foci_bw3D,[],3);
        max_image00 = max_image00/single_Inten;
        EL_info = get_EL(em_mask);   %%% get EL information (extreme points, EL length) from the embryo mask

        [nucleus_RNA_profile,foci_RNA_profile,cytoplasmic_RNA_profile] = RNA_profile3D(mask_stack,nucleus_DAPI_profile,em_mask,foci_bw3D,max_image00,EL_info,SS,RNA_stack,resolution,image_folder,N_cycle,[]);   %%% Input the foci area matrix from modified_foci.m if gaussian fitting method is used
        RNA_DNA(nucleus_DAPI_profile(:,1),nucleus_RNA_profile(:,[3,4]),image_folder,N_cycle,100);
        TXnoise_plot(nucleus_RNA_profile,N_cycle,image_folder,[101,102]);
        

        [foci_bw3D2,max_image002,SS2] = modified_foci3D([folder_list{list_I,1},hist_folder2,sub_list{list_J,3}(1:(end-1)),hist_tail],max_image,RNA2_channel,mask_stack,N_cycle,image_folder,RNA_RNA2_mismatch,Inten_thresh2,44);
        foci_bw2D2 = max(foci_bw3D2,[],3);
        max_image002 = max_image002/single_Inten2;
        
        [nucleus_signal2_profile,foci_signal2_profile,cytoplasmic_signal2_profile] = RNA_profile3D(mask_stack,nucleus_DAPI_profile,em_mask,foci_bw3D2,max_image002,EL_info,SS2,RNA2_stack,resolution,image_folder,N_cycle,[45,46,47,48,49,411]);   %%% Input the foci area matrix from modified_foci.m if gaussian fitting method is used
        nucleus_signal2_profile(:,1) = nucleus_RNA_profile(:,1);
        
        clear mask_stack RNA_stack RNA2_stack % nucleus_protein_profile_ab
        
%%% RNA match check: %%%===================================================
        RNA_fit = xlsread([folder_list{list_I,1},hist_folder,sub_list{list_J,3}(1:end-1),hist_tail]);
        RNA2_fit = xlsread([folder_list{list_I,1},hist_folder2,sub_list{list_J,3}(1:end-1),hist_tail]);
        RNA_fit_single = xlsread([folder_list{list_I,1},fit_folder,sub_list{list_J,3}(1:end-1),hist_tail]);
        RNA2_fit_single = xlsread([folder_list{list_I,1},fit_folder2,sub_list{list_J,3}(1:end-1),hist_tail]);
        
        [foci_mismatch,dxy0] = z_mismatch(RNA_fit,RNA2_fit,channel_name{RNA_channel},channel_name{RNA2_channel},size(seg_bw),image_folder,dmax0);
        all_foci = cat(1,all_foci,{foci_mismatch});
        all_foci_xy = cat(1,all_foci_xy,{dxy0});
        
        
        single_mismatch = mismatch_single(RNA_fit_single,RNA2_fit_single,channel_name{RNA_channel},channel_name{RNA2_channel},size(seg_bw),image_folder);
        all_single = cat(1,all_single,{single_mismatch});
%         [length(RNA_fit_single),length(RNA2_fit_single),length(single_mismatch)]
        nucleus_DAPI_profile = [nucleus_RNA_profile(:,1),nucleus_DAPI_profile];
        
%%% Fiducial marker identification of nascent mRNAs: %%%===================
        marker_RNA_profile = RNA_marker3D(RNA_fit,RNA2_fit,single_Inten,Inten_thresh,single_Inten2,Inten_thresh2,dmax0,EL_info,image_folder,N_cycle);

%%% Output: %%%============================================================
        if exist(result_folder) ~= 7
            mkdir(result_folder);
        end
        saveas(4,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,seg_add,D3_add,figure_tail]);
        saveas(5,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),fish_add,int_add,D3_add,figure_tail]);
        saveas(6,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),fish_add,num_add,D3_add,figure_tail]);
        saveas(7,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,num_add,D3_add,figure_tail]);
        saveas(8,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,int_add,D3_add,figure_tail]);
        saveas(9,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,int_add,D3_add,cmp_add,figure_tail]);
        saveas(11,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),nu_add,fate_add,D3_add,figure_tail]);

        saveas(44,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,signal2_add,seg_add,D3_add,figure_tail]);
        saveas(45,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),fish_add,signal2_add,int_add,D3_add,figure_tail]);
        saveas(46,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),fish_add,signal2_add,num_add,D3_add,figure_tail]);
        saveas(47,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,signal2_add,num_add,D3_add,figure_tail]);
        saveas(48,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,signal2_add,int_add,D3_add,figure_tail]);
        saveas(49,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,signal2_add,int_add,D3_add,cmp_add,figure_tail]);
        saveas(411,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),nu_add,signal2_add,fate_add,D3_add,figure_tail]);
        
        saveas(16,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),nu_add,DAPI_add,D3_add,figure_tail]);
        
        saveas(100,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),DAPI_add,RNA_add,figure_tail]);
        saveas(101,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),RNA_add,noise_add,figure_tail]);

        saveas(91,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mismatch_add,z_add,figure_tail]);
        saveas(92,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mismatch_add,xy_add,figure_tail]);
        saveas(93,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mismatch_add,FISH_add,single_add,figure_tail]);
        saveas(94,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mismatch_add,xy_add,single_add,figure_tail]);
        saveas(1,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),marker_add,int_add,D3_add,figure_tail]);
        
        save([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mat_tail],'foci_mismatch','single_mismatch','dxy0','marker_RNA_profile','-append');
        save([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mat_tail],'image_folder','image_type','WGA_channel','DAPI_channel','RNA_channel','RNA2_channel','resolution','resolutionz','N_cycle','seg_bw','cyto_bw','max_image','foci_bw2D','foci_bw2D2','nucleus_RNA_profile','foci_RNA_profile','cytoplasmic_RNA_profile','nucleus_signal2_profile','foci_signal2_profile','cytoplasmic_signal2_profile','nucleus_DAPI_profile','-append');
        
        if ~isempty(nucleus_RNA_profile)
           xlswrite([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),RNA_add,nu_add,output_tail],nucleus_RNA_profile);
        end
        if ~isempty(cytoplasmic_RNA_profile)
           xlswrite([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),RNA_add,cyto_add,output_tail],cytoplasmic_RNA_profile);
        end
        if ~isempty(foci_RNA_profile)
           xlswrite([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),RNA_add,foci_add,output_tail],foci_RNA_profile);
        end
        if ~isempty(nucleus_signal2_profile)
           xlswrite([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),signal2_add,nu_add,output_tail],nucleus_signal2_profile);
        end
        if ~isempty(cytoplasmic_signal2_profile)
           xlswrite([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),signal2_add,cyto_add,output_tail],cytoplasmic_signal2_profile);
        end
        if ~isempty(foci_signal2_profile)
           xlswrite([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),signal2_add,foci_add,output_tail],foci_signal2_profile);
        end
        if ~isempty(foci_mismatch)
            xlswrite([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mismatch_add,z_add,output_tail],foci_mismatch);
        end
        
        sub_num(list_J,13) = N_cycle;
        
%         clear ('seg_bw','cyto_bw','max_image','foci_bw','nucleus_RNA_profile','foci_RNA_profile','cytoplasmic_RNA_profile','nucleus_protein_profile','cytoplasmic_protein_profile','N_cycle','foci_data','fake_data','foci_signal2_profile','cytoplasmic_signal2_profile','foci_data2','fake_data2');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        try
            xlswrite([folder_list{list_I,1},in_folder,input_name],cat(2,sub_list,num2cell(sub_num)));
        catch
            xlwrite([folder_list{list_I,1},in_folder,input_name],cat(2,sub_list,num2cell(sub_num)));
        end
        
    end


%% Overall foci plot/analysis: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    foci1 = cell2mat(all_foci);
    foci1xy = cell2mat(all_foci_xy);
    imdx = zeros(size(seg_bw));
    imdy = zeros(size(seg_bw));
    for I0 = 1:length(foci1)
        imdx(round(foci1(I0,3)),round(foci1(I0,4))) = foci1(I0,5);
        imdy(round(foci1(I0,3)),round(foci1(I0,4))) = foci1(I0,6);
    end

    figure(201)
    subplot(1,3,1)
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

    subplot(1,3,2)
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

    subplot(1,3,3)
        np = hist(foci1xy,hist_binxy);
        bar(hist_binxy(1:end-1),np(1:end-1),'hist')
        xlabel('dxy')
        ylabel('#')
        title(['Chromatic mismatch on xy dimensions: ',folder_list{list_I,1}],'Interpreter','none','FontName','Arial','FontSize',16,'FontWeight','bold')

    saveas(201,[result_folder,mismatch_add,xy_add,figure_tail]);

    
    figure(202)
    subplot(2,2,1)
        p0 = polyfit(foci1(foci1(:,3) <= 1354,3),foci1(foci1(:,3) <= 1354,5),1);
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
    saveas(202,[result_folder,mismatch_add,xy_add,fit_add2,figure_tail]);
        
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Overall foci plot/analysis: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    foci1 = cell2mat(all_single);
    imdx = zeros(size(seg_bw));
    imdy = zeros(size(seg_bw));
    for I0 = 1:length(foci1)
        imdx(round(foci1(I0,3)),round(foci1(I0,4))) = foci1(I0,5);
        imdy(round(foci1(I0,3)),round(foci1(I0,4))) = foci1(I0,6);
    end

    figure(203)
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

    saveas(203,[result_folder,mismatch_add,xy_add,single_add,figure_tail]);

    
    figure(204)
    subplot(2,2,1)
        p0 = polyfit(foci1(foci1(:,3) <= 1354,3),foci1(foci1(:,3) <= 1354,5),1);
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
    saveas(204,[result_folder,mismatch_add,xy_add,fit_add2,single_add,figure_tail]);
        
% % % %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




end
