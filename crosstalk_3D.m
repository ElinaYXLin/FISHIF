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

% channel_in = 3;
% channel_out = [2,4];
channel_in = 2;
channel_out = 3;
dxy_max = 10;
dxy0 = {[0,0],[0,1]};
z_range = 8:20;

sub_ij = [4,6];
corr_add = '_corr';
cross_add = '_cross';
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
        
        all_color = eval(folder_list{list_I,5});
        in_color = all_color(channel_in);
        out_color = all_color(channel_out);

        
        p_all = cell(length(channel_in),length(channel_out));
        if exist(result_folder) ~= 7
            mkdir(result_folder);
        end
        h2 = figure; set(h2,'Name',[image_folder,': ij displacement'])
        h3 = figure; set(h3,'Name',[image_folder,': cross-talk'])
        for ic = 1:length(channel_in)
            for jc = 1:length(channel_out)
%%% Crosstalk analysis: %%%================================================
                fname1 = dir([image_folder,image_type]);
                dmn = zeros(length(fname1),2);
                
                h1 = figure; maximize(h1); set(h1,'Name',[image_folder,': ',in_color{ic},' --> ',out_color{jc}])
                %%%%% xy displacement analysis:
                for list_K = 1:length(fname1)
                    im0 = imread([image_folder,fname1(list_K).name]);
                    r2 = corr2D(im0(:,:,channel_in(ic)),im0(:,:,channel_out(jc)),dxy_max,dxy_max);
                    [m0,n0] = find(r2 == max(r2(:)));
                    dmn(list_K,:) = [m0,n0]-dxy_max-1;
                    
                    if list_K == z_range(1)
                        xx0 = zeros(size(im0,1)-dxy0{ic,jc}(1),size(im0,2)-dxy0{ic,jc}(2),length(z_range));
                        yy0 = xx0;
                    end
                    if any(list_K == z_range)
                        xx0(:,:,list_K) = im0(max(1,dxy0{ic,jc}(1)+1):min(size(im0,1)+dxy0{ic,jc}(1),size(im0,1)),max(1,dxy0{ic,jc}(2)+1):min(size(im0,2)+dxy0{ic,jc}(2),size(im0,2)),channel_in(ic));
                        yy0(:,:,list_K) = im0(max(-dxy0{ic,jc}(1)+1,1):min(size(im0,1),size(im0,1)-dxy0{ic,jc}(1)),max(-dxy0{ic,jc}(2)+1,1):min(size(im0,2),size(im0,2)-dxy0{ic,jc}(2)),channel_out(jc));
                    end
                    
                    figure(h1)
                    ha1 = subaxis(sub_ij(1),sub_ij(2),list_K, 'Spacing', 0, 'PaddingRight',0.03, 'PaddingLeft',0.05, 'PaddingTop',0.03, 'PaddingBottom',0.03, 'MarginLeft', 0.01, 'MarginRight', 0.001, 'MarginTop', 0.015, 'MarginBottom', 0.025);
                        imagesc([-dxy_max,dxy_max],[-dxy_max,dxy_max],r2)
                        xlabel('j')
                        ylabel('i')
                        title(fname1(list_K).name)
                end
                saveas(h1,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),xy_add,corr_add,'_',in_color{ic},'_',out_color{jc},figure_tail]);
                
                figure(h2); maximize(h2);
                ha2 = subaxis(length(channel_in),length(channel_out),jc,ic, 'Spacing', 0, 'PaddingRight',0.03, 'PaddingLeft',0.05, 'PaddingTop',0.03, 'PaddingBottom',0.03, 'MarginLeft', 0.01, 'MarginRight', 0.001, 'MarginTop', 0.015, 'MarginBottom', 0.025);
                    plot(dmn(:,1),'r')
                    hold on
                    plot(dmn(:,2),'b')
                    legend({'i','j'})
                    xlabel('Slice #')
                    ylabel('dxy (Pixel)')
                    title([in_color{ic},' --> ',out_color{jc}])
                
                %%%%% Intensity correlation analysis:
                p0 = polyfit(xx0(:),yy0(:),1);
                p_all{ic,jc} = p0;
                figure(h3); maximize(h3);
                ha3 = subaxis(length(channel_in),length(channel_out),jc,ic, 'Spacing', 0, 'PaddingRight',0.03, 'PaddingLeft',0.05, 'PaddingTop',0.03, 'PaddingBottom',0.03, 'MarginLeft', 0.01, 'MarginRight', 0.001, 'MarginTop', 0.015, 'MarginBottom', 0.025);
                    plot(xx0(:),yy0(:),'k.')
                    hold on
                    xlim0 = xlim;
                    ylim0 = polyval(p0,xlim0);
                    plot(xlim0,ylim0,'r')
                    legend({'Data',['Fit: k = ',num2str(p0(1))]})
                    xlabel('Slice #')
                    ylabel('dxy (Pixel)')
                    title([in_color{ic},' --> ',out_color{jc},': z = ',num2str(z_range(1)),'-',num2str(z_range(end)),', dxy = ',num2str(dxy0{ic,jc})])
            end
        end
%%% Output: %%%============================================================
        saveas(h2,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),xy_add,corr_add,D3_add,figure_tail]);
%         saveas(h3,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),cross_add,figure_tail],'-v7.3');
        save([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mat_tail],'p_all','-append');
    end
end
