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
Nbin0 = 15;

I_old_new = 3;
channel_in = 3;
channel_out = [2,4];
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
        image_folder = [folder_list{list_I,1},in_folder,sub_list{list_J,I_old_new}];
        result_folder = [folder_list{list_I,1},out_folder];
        load([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mat_tail],'max_image');
        fname1 = dir([image_folder,image_type]);
        Lxyz = [size(max_image,1),size(max_image,2),length(fname1)];
        
        all_color = eval(folder_list{list_I,5});
        in_color = all_color(channel_in);
        out_color = all_color(channel_out);
        
        foci_info = xlsread([folder_list{list_I,1},spot_folder,sub_list{list_J,3}(1:end-1),hist_tail]);
        foci_xyz = round(foci_info(:,6:8));
        im0_i = zeros(0);
        for ii = 1:size(foci_xyz,1)
            [X0,Y0,Z0] = ndgrid(max(foci_xyz(ii,1)-dxyz_RNA(1),1):min(foci_xyz(ii,1)+dxyz_RNA(1),Lxyz(1)),max(foci_xyz(ii,2)-dxyz_RNA(2),1):min(foci_xyz(ii,2)+dxyz_RNA(2),Lxyz(2)),max(foci_xyz(ii,3)-dxyz_RNA(3),1):min(foci_xyz(ii,3)+dxyz_RNA(3),Lxyz(3)));
            im0_i = cat(1,im0_i,sub2ind(Lxyz,X0(:),Y0(:),Z0(:)));
        end
        im0_i = unique(im0_i);
        [im0_x,im0_y,im0_z] = ind2sub(Lxyz,im0_i);
        
        p_all = cell(length(channel_in),length(channel_out));
        if exist(result_folder) ~= 7
            mkdir(result_folder);
        end
% % %         h2 = figure; maximize(h2); set(h2,'Name',[image_folder,': ij displacement'],'PaperPositionMode', 'auto')
        h3 = figure; maximize(h3); set(h3,'Name',[image_folder,': cross-talk'],'PaperPositionMode', 'auto')
        for ic = 1:length(channel_in)
            for jc = 1:length(channel_out)
%%% Crosstalk analysis: %%%================================================
                dmn = zeros(length(fname1),2);
                xx0 = zeros(0);
                yy0 = zeros(0);
                
% % %                 h1 = figure; maximize(h1); set(h1,'Name',[image_folder,': ',in_color{ic},' --> ',out_color{jc}],'PaperPositionMode', 'auto')
                %%%%% xy displacement analysis:
                for list_K = 1:length(fname1)
                    im0 = imread([image_folder,fname1(list_K).name]);
% % %                     r2 = corr2D(im0(:,:,channel_in(ic)),im0(:,:,channel_out(jc)),dxy_max,dxy_max);
% % %                     [m0,n0] = find(r2 == max(r2(:)));
% % %                     dmn(list_K,:) = [m0,n0]-dxy_max-1;
                    
                    is = (im0_z == list_K) & (im0_x-dxy0{ic,jc}(1) >= 1) & (im0_x+dxy0{ic,jc}(1) <= Lxyz(1)) & (im0_y-dxy0{ic,jc}(2) >= 1) & (im0_y+dxy0{ic,jc}(2) <= Lxyz(2));
                    im1_i = sub2ind(size(im0),im0_x(is),im0_y(is),channel_in(ic)*ones(nnz(is),1));
                    im2_i = sub2ind(size(im0),im0_x(is)-dxy0{ic,jc}(1),im0_y(is)-dxy0{ic,jc}(2),channel_out(jc)*ones(nnz(is),1));
                    xx0 = cat(1,xx0,im0(im1_i));
                    yy0 = cat(1,yy0,im0(im2_i));
                    
% % %                     figure(h1)
% % %                     ha1 = subaxis(sub_ij(1),sub_ij(2),list_K, 'Spacing', 0, 'PaddingRight',0.03, 'PaddingLeft',0.05, 'PaddingTop',0.03, 'PaddingBottom',0.03, 'MarginLeft', 0.01, 'MarginRight', 0.001, 'MarginTop', 0.015, 'MarginBottom', 0.025);
% % %                         imagesc([-dxy_max,dxy_max],[-dxy_max,dxy_max],r2)
% % %                         xlabel('j')
% % %                         ylabel('i')
% % %                         title(fname1(list_K).name)
                end
% % %                 saveas(h1,[result_folder,sub_list{list_J,I_old_new}(1:(length(sub_list{list_J,I_old_new})-1)),xy_add,corr_add,'_',in_color{ic},'_',out_color{jc},figure_tail]);
% % %                 saveas(h1,[result_folder,sub_list{list_J,I_old_new}(1:(length(sub_list{list_J,I_old_new})-1)),xy_add,corr_add,'_',in_color{ic},'_',out_color{jc},figure_tail2]);
                
% % %                 figure(h2); maximize(h2);
% % %                 ha2 = subaxis(length(channel_in),length(channel_out),jc,ic, 'Spacing', 0, 'PaddingRight',0.03, 'PaddingLeft',0.05, 'PaddingTop',0.03, 'PaddingBottom',0.03, 'MarginLeft', 0.01, 'MarginRight', 0.001, 'MarginTop', 0.015, 'MarginBottom', 0.025);
% % %                     plot(dmn(:,1),'r')
% % %                     hold on
% % %                     plot(dmn(:,2),'b')
% % %                     legend({'i','j'})
% % %                     xlabel('Slice #')
% % %                     ylabel('dxy (Pixel)')
% % %                     title([in_color{ic},' --> ',out_color{jc}])
                
                %%%%% Intensity correlation analysis:
%                 p0 = polyfit(double(xx0(:)),double(yy0(:)),1);
                p0 = [double(xx0(:))\double(yy0(:)),0];
                [~,ip0] = sort(rand(size(xx0(:))));
                p_all{ic,jc} = p0;
                figure(h3); maximize(h3);
                ha3 = subaxis(length(channel_in),length(channel_out),jc,ic, 'Spacing', 0, 'PaddingRight',0.03, 'PaddingLeft',0.05, 'PaddingTop',0.03, 'PaddingBottom',0.03, 'MarginLeft', 0.01, 'MarginRight', 0.001, 'MarginTop', 0.015, 'MarginBottom', 0.025);
                    plot(xx0(:),yy0(:),'k.')
                    hold on
                    xlim0 = xlim;
                    ylim0 = polyval(p0,xlim0);
                    plot(xlim0,ylim0,'r')
                    legend({'Data',['Fit: k = ',num2str(p0(1))]})
                    xlabel(in_color{ic})
                    ylabel(out_color{jc})
                    title([in_color{ic},' --> ',out_color{jc},': dxy = ',num2str(dxy0{ic,jc}),', r = ',num2str(corr(double(xx0(:)),double(yy0(:)))),', r0 = ',num2str(corr(double(xx0(:)),double(yy0(ip0))))])
            end
        end
%%% Output: %%%============================================================
% % %         saveas(h2,[result_folder,sub_list{list_J,I_old_new}(1:(length(sub_list{list_J,I_old_new})-1)),xy_add,corr_add,D3_add,figure_tail]);
% % %         saveas(h2,[result_folder,sub_list{list_J,I_old_new}(1:(length(sub_list{list_J,I_old_new})-1)),xy_add,corr_add,D3_add,figure_tail2]);
        saveas(h3,[result_folder,sub_list{list_J,I_old_new}(1:(length(sub_list{list_J,I_old_new})-1)),cross_add,mismatch_add,intercept_add,figure_tail]);
        saveas(h3,[result_folder,sub_list{list_J,I_old_new}(1:(length(sub_list{list_J,I_old_new})-1)),cross_add,mismatch_add,intercept_add,figure_tail2]);
        if exist([result_folder,sub_list{list_J,I_old_new}(1:(length(sub_list{list_J,I_old_new})-1)),mismatch_add,intercept_add,mat_tail])
            save([result_folder,sub_list{list_J,I_old_new}(1:(length(sub_list{list_J,I_old_new})-1)),mismatch_add,intercept_add,mat_tail],'p_all','all_color','channel_in','channel_out','dxy0','-append');
        else
            save([result_folder,sub_list{list_J,I_old_new}(1:(length(sub_list{list_J,I_old_new})-1)),mismatch_add,intercept_add,mat_tail],'p_all','all_color','channel_in','channel_out','dxy0');
        end
    end
end
