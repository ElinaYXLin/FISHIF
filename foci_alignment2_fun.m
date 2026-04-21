function foci_alignment2_fun

clear
% close all

tic
%% Data loading: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imfolder_list = uigetfile_n_dir;

for I_folder = 1:length(imfolder_list)
    imfolder = imfolder_list{I_folder};
    if ~(imfolder(end) == '/' || imfolder(end) == '\')
        imfolder = [imfolder,'/'];
    end
    
    I_cut = strfind(imfolder,'stacks');
    if I_cut
        imname = imfolder(I_cut+7:end);
        imfolder = imfolder(1:I_cut-1);
    else           
        imname = [];
    end
    folder_list = {imfolder};
    foci_alignment_core(folder_list,imname);
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
toc



function foci_alignment_core(folder_list,imname)

xy_mismatch_check = true;
% % % epsilonxy = 0.55;
% % % epsilonxyz = 0.55;
% % % epsilonz = 1;

epsilonxy = 0.4;
epsilonxyz = 0.4;
epsilonz = 2;

%% Parameter setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
list_name = 'Duallist.xls';
in_folder = 'stacks/';
% old_add = '_old';
old_add = '';
input_name = 'matchlist.xls';
mismatch_name = 'mismatch.xls';
mismatch_name2 = 'mismatch_60X.xls';
mismatch_name3 = 'mismatch_60X_FA.xls';
xymismatch_name = 'xymismatch.mat';
xymismatch_name2 = 'xymismatch_60X.mat';
xymismatch_name3 = 'xymismatch_60X_FA.mat';
fast_add = 'F';
airy_add = 'A';
double_add = 'D';
stitch_name = 'stitch_parameter.mat';

flip_list = {'Cad'};
flip0 = false;
image_type = '*.tif';
mask_folder = 'masks/';
mask_name = 'mask.mat';
out_folder = 'Results/';
% out_folder0 = 'Results_3Dz1/';
% out_folder0 = 'Results_decross/';
% out_folder0 = 'Results_original/';
out_folder0 = '';
% % out_folder0 = 'Results0/';
hist_folder = 'Histogram/';
fit_folder = 'Histogram_A/';
hist_folder2 = 'Histogram_RNA2/';
fit_folder2 = 'Histogram_A_RNA2/';
hist_folder_couple = 'Histogram_alignment/';
hist_folder2_couple = 'Histogram_alignment_RNA2/';
fit_add = '_spot_fit';
hist_tail = '_raw.xls';
hist_link_tail = '_link.xls';
N_thresh = 0;
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
protein_add = '_protein';
RNA_add = '_RNA';
signal2_add = '_RNA2';
nu_add = '_nucleus';
cyto_add = '_cytoplasm';
fate_add = '_fate';
reg_add = '_regulation';
fluc_add = '_fluc';
bino_add = '_bino';
local_add = '_local';
hist_add = '_hist';
fake_add = '_fake';
abs_add = '_abs';
rad_add = '_rad';
D3_add = '_3D';
compare_add = '_compare';
DAPI_add = '_DAPI';
noise_add = '_noise';
dist_add = '_dist';
fig_tail = '.fig';
sub_pos = [3,3];
cycle_pos0 = zeros(1,sub_pos(1));
cycle_range = [11,12,13];
dr = 0.05;
r_min = 0:dr:3;
r_max = r_min+dr;
r_range = (r_min+r_max)/2;
z_range = -6:6;
bin3D = {-1:dr:1,-1:dr:1};

rlim0 = [0,1.5];
zlim0 = [-5.5,5.5];
xylim0 = {[-0.4,0.4],[-0.4,0.4]};
% EL_range = [0.2,0.6,0,0.1,0,0.1];
EL_range = [];
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Image analysis: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
list_I = 1;
[sub_num, sub_list] = xlsread([folder_list{list_I,1},in_folder,input_name]);
[M1,M2] = size(sub_list);
% %     if ~isempty(out_folder0)
% %         copyfile([folder_list{list_I,1},out_folder],[folder_list{list_I,1},out_folder0]);
% %     end
if isempty(imname)
    list_J0 = 1:M1;
else
    list_J0 = find(strcmp(sub_list(:,3),imname));
end

for list_J = list_J0
    if isempty(strfind(sub_list{list_J,3},'_60X')) && isempty(strfind(sub_list{list_J,2},'60X')) && isempty(strfind(sub_list{list_J,2},airy_add)) && isempty(strfind(sub_list{list_J,2},fast_add))
        [~,~,mismatch_matrix] = xlsread(mismatch_name);
        load(xymismatch_name)
        mismatch_name
    elseif (~isempty(strfind(sub_list{list_J,3},'_60X')) || ~isempty(strfind(sub_list{list_J,2},'60X'))) && isempty(strfind(sub_list{list_J,2},airy_add)) && isempty(strfind(sub_list{list_J,2},fast_add))
        [~,~,mismatch_matrix] = xlsread(mismatch_name2);
        load(xymismatch_name2)
        mismatch_name2
    elseif (~isempty(strfind(sub_list{list_J,3},'_60X')) || ~isempty(strfind(sub_list{list_J,2},'60X'))) && ~isempty(strfind(sub_list{list_J,2},fast_add))
        [~,~,mismatch_matrix] = xlsread(mismatch_name3);
        load(xymismatch_name3)
        mismatch_name3
    else
        error('No mismatch information')
    end

%         image_folder = [folder_list{list_I,1},in_folder,sub_list{list_J,3}];
    image_folder = [folder_list{list_I,1},in_folder,sub_list{list_J,3}(1:end-1),old_add,'/'];

    if isempty(flip0)
        flip_axis = any(cellfun(@(x) ~isempty(strfind(image_folder,x)),flip_list));   %%% Check whether the protein profile has an opposite direction (to Bcd)
    else
        flip_axis = flip0;
    end

    result_folder = [folder_list{list_I,1},out_folder];
    resolutionz = sub_num(list_J,11);
    signal2_channel = sub_num(list_J,12);

    load([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mat_tail]);
    image_folder = [folder_list{list_I,1},in_folder,sub_list{list_J,3}(1:end-1),old_add,'/'];

    if ~isempty(strfind(sub_list{list_J,2},double_add))
        Nbin = sub_num(list_J,2:3);
        Mdim = 1:2;
    else
        Nbin = ones(1,2);
        Mdim = sub_num(list_J,3);
        Nbin(Mdim) = sub_num(list_J,2);
        Mdim = 1:2;
    end

    %%% Stitch parameter loading: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if exist([image_folder,stitch_name])
        load([image_folder,stitch_name],'x_start','y_start','imcorr1')
        size_tile = size(imcorr1);
        size_im = size(max_image);
    else
        x_start = ones(Nbin);
        y_start = ones(Nbin);
        size_im = size(max_image);
        size_tile = size_im(1:2)./Nbin;
    end        

    x_start_im = ones(size(x_start));
    y_start_im = ones(size(x_start));
    x_center_im = zeros(size(x_start));
    y_center_im = zeros(size(x_start));
    for ii = 1:Nbin(1)
        for jj = 1:Nbin(2)
            if jj > 1
                x_start_im(ii,jj) = x_start_im(ii,jj-1)+size_tile(2)-x_start(ii,jj-1)+1;
            end
            if ii > 1
                y_start_im(ii,jj) = y_start_im(ii-1,jj)+size_tile(1)-y_start(ii-1,jj)+1;
            end
            x_center_im(ii,jj) = x_start_im(ii,jj)+(size_tile(2)+1)/2-x_start(ii,jj);
            y_center_im(ii,jj) = y_start_im(ii,jj)+(size_tile(1)+1)/2-y_start(ii,jj);
        end
% %             if x_start_im(ii,end)+size_tile(2)-x_start(ii,end) ~= size_im(2)
% %                 error('dimension 2 mismatch')
% %             end
    end
% %         if y_start_im(end,end)+size_tile(1)-y_start(end,end) ~= size_im(1)
% %             error('dimension 1 mismatch')
% %         end
    clear x_start y_start imcorr1 size_im size_tile
    %%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % % %     load([folder_list{list_I,1},mask_folder,sub_list{list_J,3},mask_name]);   %%% load 3D mask
% % % %     z_size = size(mask_stack,3);
% %         load([folder_list{list_I,1},fit_folder,sub_list{list_J,3}(1:end-1),fit_add,mat_tail]);   %%% load spot intensity fitting result
    b = 1;
    Inten_thresh = b*N_thresh;   %%% set foci intensity threshold
    single_Inten = b;
% %         load([folder_list{list_I,1},fit_folder2,sub_list{list_J,3}(1:end-1),fit_add,mat_tail]);   %%% load spot intensity fitting result
    Inten_thresh2 = b*N_thresh;   %%% set foci intensity threshold
    single_Inten2 = b;
%         N_cycle = round(log2(max(mask_stack(:))))+2;   %%% Calculate the nuclear cycle number
    N_cycle = sub_num(list_J,13);

% %     all_color = eval(folder_list{list_I,5});
    RNA_color = all_color{RNA_channel};
    signal2_color = all_color{signal2_channel};
    RNA_signal2_mismatch = mismatch_matrix{strcmp(RNA_color,mismatch_matrix(1,:)),strcmp(signal2_color,mismatch_matrix(:,1))};
    if xy_mismatch_check
        RNA_signal2_xymismatch = {eval([RNA_color,'_',signal2_color]),eval([RNA_color,'_',signal2_color,'_con']),eval([RNA_color,'_',signal2_color,'_x0']),Nbin,Mdim,resolution,resolution_mismatch};
    else
        RNA_signal2_xymismatch = [];
    end
    signal2_RNA_mismatch = mismatch_matrix{strcmp(signal2_color,mismatch_matrix(:,1)),strcmp(RNA_color,mismatch_matrix(1,:))};
    if xy_mismatch_check
        signal2_RNA_xymismatch = {eval([signal2_color,'_',RNA_color]),eval([signal2_color,'_',RNA_color,'_con']),eval([signal2_color,'_',RNA_color,'_x0']),Nbin,Mdim,resolution,resolution_mismatch};
    else
        signal2_RNA_xymismatch = [];
    end

    [foci_list,~,~] = xlsread([folder_list{list_I,1},hist_folder,sub_list{list_J,3}(1:(end-1)),hist_tail]);
    [foci_list2,~,~] = xlsread([folder_list{list_I,1},hist_folder2,sub_list{list_J,3}(1:(end-1)),hist_tail]);
    foci_listb = foci_list;
    foci_list2b = foci_list2;
    foci_list2(:,6:8) = position_adjust(foci_list2(:,6:8),max_image,RNA_signal2_mismatch,RNA_signal2_xymismatch,{x_start_im,y_start_im,x_center_im,y_center_im});
    foci_listb(:,6:8) = position_adjust(foci_list(:,6:8),max_image,signal2_RNA_mismatch,signal2_RNA_xymismatch,{x_start_im,y_start_im,x_center_im,y_center_im});

%%% foci coupling:
    xyz1 = foci_list(:,6:8).*repmat([resolution,resolution,resolutionz],size(foci_list,1),1);
    xyz2 = foci_list2(:,6:8).*repmat([resolution,resolution,resolutionz],size(foci_list2,1),1);
        
% % % %     dxyz = pdist2(xyz1,xyz2);
% % % %     [d12min,I2min] = min(dxyz,[],2);
        
    dxyz0 = pdist2(xyz1,xyz2);
    dxy0 = pdist2(xyz1(:,1:2),xyz2(:,1:2));
    dxy = dxy0;
    dxy(dxy > epsilonxy) = 1e8;
    dz = pdist2(foci_list(:,8),foci_list2(:,8));
    dxyz = dxy + 1e8*(dz > epsilonz);

    [I12_couple, I1_single, I2_single] = assignmunkres(dxyz,epsilonxyz);
    I1_couple = I12_couple(:,1);
    I2_couple = I12_couple(:,2);

    foci_list_add = zeros(length(I2_single),size(foci_list2,2));
    foci_list_add(:,6:8) = foci_list2(I2_single,6:8);
    foci_list_new = cat(1,foci_list(I1_couple,:),foci_list(I1_single,:),foci_list_add);
% % %     foci_list_new(:,8) = round(foci_list_new(:,8));
        
    foci_list2_add = zeros(length(I1_single),size(foci_listb,2));
    foci_list2_add(:,6:8) = foci_listb(I1_single,6:8);
    foci_list2_new = cat(1,foci_list2b(I2_couple,:),foci_list2_add,foci_list2b(I2_single,:));
% % %     foci_list2_new(:,8) = round(foci_list2_new(:,8));
        
        
% % % %     [d12min,I2min] = min(dxyz,[],2);
% % % %     I1_couple = find(d12min <= epsilonxyz);
% % % %     I2_couple = I2min(I1_couple);
% % % %     I1_single = find(d12min > epsilonxyz & 2*pi*prod(foci_list(:,1:3),2) >= Inten_thresh);
% % % %         I2_is_single = true(size(foci_list2,1),1);
% % % %         I2_is_single(I2_couple) = false;
% % % %     I2_single = find(I2_is_single & 2*pi*prod(foci_list2(:,1:3),2) >= Inten_thresh2);
% % % % 
% % % %     N2_couple = unique(I2_couple);
% % % %     foci_list_new = zeros(length(N2_couple)+length(I1_single)+length(I2_single),size(foci_list,2));
% % % %     foci_list2_new = zeros(length(N2_couple)+length(I1_single)+length(I2_single),size(foci_list,2));
% % % %     for ii = 1:length(N2_couple)
% % % %         I2_couple0 = N2_couple(ii);
% % % %         I1_couple0 = I1_couple(I2_couple == I2_couple0);
% % % %         if length(I1_couple0) == 1
% % % %             foci_list_new(ii,:) = foci_list(I1_couple0,:);
% % % %         else
% % % %             w1 = prod(foci_list(I1_couple0,1:3),2);
% % % %             w0 = w1/sum(w1);
% % % %             foci_list_new(ii,:) = sum(foci_list(I1_couple0,:).*repmat(w0,1,size(foci_list,2)));
% % % %             foci_list_new(ii,1) = sum(w1)/foci_list_new(ii,2)/foci_list_new(ii,3);
% % % %         end
% % % %         foci_list2_new(ii,:) = foci_list2b(I2_couple0,:);
% % % %     end
% % % %     foci_list_new((length(N2_couple)+1):(length(N2_couple)+length(I1_single)),:) = foci_list(I1_single,:);
% % % %     foci_list_new((length(N2_couple)+length(I1_single)+1):end,:) = 0;
% % % %     foci_list_new((length(N2_couple)+length(I1_single)+1):end,6:8) = foci_list2(I2_single,6:8);
% % % %     foci_list_new(:,8) = round(foci_list_new(:,8));
% % % % 
% % % %     foci_list2_new((length(N2_couple)+1):(length(N2_couple)+length(I1_single)),:) = 0;
% % % %     foci_list2_new((length(N2_couple)+1):(length(N2_couple)+length(I1_single)),6:8) = foci_listb(I1_single,6:8);
% % % %     foci_list2_new((length(N2_couple)+length(I1_single)+1):end,:) = foci_list2b(I2_single,:);
        

    foci_list_link = zeros(length(I1_couple)+length(I1_single)+length(I2_single),size(foci_list,2)+2);
    foci_list_link(1:length(I1_couple),:) = [foci_list(I1_couple,:),I1_couple,I2_couple];
    foci_list_link((length(I1_couple)+1):(length(I1_couple)+length(I1_single)),:) = [foci_list(I1_single,:),I1_single,zeros(size(I1_single))];
    foci_list_link((length(I1_couple)+length(I1_single)+1):end,6:8) = foci_list2(I2_single,6:8);
    foci_list_link((length(I1_couple)+length(I1_single)+1):end,end-1:end) = [zeros(size(I2_single)),I2_single];

    foci_list2_link = zeros(length(I1_couple)+length(I1_single)+length(I2_single),size(foci_list,2)+2);
    foci_list2_link(1:length(I1_couple),:) = [foci_list2b(I2_couple,:),I1_couple,I2_couple];        
    foci_list2_link((length(I1_couple)+1):(length(I1_couple)+length(I1_single)),6:8) = foci_listb(I1_single,6:8);
    foci_list2_link((length(I1_couple)+1):(length(I1_couple)+length(I1_single)),end-1:end) = [I1_single,zeros(size(I1_single))];
    foci_list2_link((length(I1_couple)+length(I1_single)+1):end,:) = [foci_list2b(I2_single,:),zeros(size(I2_single)),I2_single];

%% Output: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Alignment result output:
    figure; maximize(gcf)
    set(gcf,'Name',[folder_list{list_I,1},sub_list{list_J,3}(1:(end-1)),': foci alignment check'])

    subplot(2,3,1)
        [d12min0,I2min0] = min(dxyz0,[],2);
        nd0 = hist(d12min0,r_range);
        nd = nd0./(r_max.^3-r_min.^3)/pi/4*3;
        bar(r_range(1:end-1),nd(1:end-1));
        hold on
        plot(epsilonxyz*ones(1,2),[0,max(nd(1:end-1))],'k--','LineWidth',2)
        xlabel('r_x_y_z (\mum)')
        ylabel('#')
        xlim(rlim0)
        title('Density of raw xyz-distance')

    subplot(2,3,2)
        [d12xymin,I2xymin] = min(dxy0,[],2);
        nd0 = hist(d12xymin,r_range);
        nd = nd0./(r_max.^2-r_min.^2)/pi;
        bar(r_range(1:end-1),nd(1:end-1));
        hold on
        plot(epsilonxy*ones(1,2),[0,max(nd(1:end-1))],'k--','LineWidth',2)
        xlabel('r_x_y (\mum)')
        ylabel('#')
        xlim(rlim0)
        title('Density of xy-distance')

    subplot(2,3,3)
        I1xy_couple = find(d12xymin <= epsilonxy);
        I2xy_couple = I2xymin(I1xy_couple);
        nz = hist(foci_list(I1xy_couple,8)-foci_list2(I2xy_couple,8),z_range);
        bar(z_range,nz);
        hold on
        plot(epsilonz*ones(1,2),[0,max(nz(2:end-1))],'k--','LineWidth',2)
        plot(-epsilonz*ones(1,2),[0,max(nz(2:end-1))],'k--','LineWidth',2)
        xlabel('dz')
        ylabel('#')
        xlim(zlim0)
        title('Distribution of z-distance between xy-paired foci')

    subplot(2,3,4)
        [N0,dxy0] = hist3(xyz1(I1_couple,1:2)-xyz2(I2_couple,1:2),'CdataMode','auto','EdgeColor','none','Ctrs',bin3D);
        imagesc(dxy0{1},dxy0{2},N0');
        axis xy
        colorbar
        xlim(xylim0{1})
        ylim(xylim0{2})
        xlabel('dx (\mum)')
        ylabel('dy (\mum)')
        title('Distribution of xy-displacement between paired foci')

    subplot(2,3,5)
        nxy0 = hist(xyz1(I1_couple,1)-xyz2(I2_couple,1),bin3D{1});
        bar(bin3D{1},nxy0);
        xlim(xylim0{1})
        xlabel('dx (\mum)')
        ylabel('#')
        title('Distribution of x-displacement between paired foci')

    subplot(2,3,6)
        nxy0 = hist(xyz1(I1_couple,2)-xyz2(I2_couple,2),bin3D{2});
        bar(bin3D{2},nxy0);
        xlim(xylim0{2})
        xlabel('dy (\mum)')
        ylabel('#')
        title('Distribution of y-displacement between paired foci')
%% 
            
    result_folder = [folder_list{list_I,1},hist_folder_couple];
    if exist(result_folder) ~= 7
        mkdir(result_folder);
    end
    delete([result_folder,sub_list{list_J,3}(1:(end-1)),hist_tail])
    delete([result_folder,sub_list{list_J,3}(1:(end-1)),hist_link_tail])
    xlswrite([result_folder,sub_list{list_J,3}(1:(end-1)),hist_tail],foci_list_new);
    xlswrite([result_folder,sub_list{list_J,3}(1:(end-1)),hist_link_tail],foci_list_link);
    saveas(gcf,[result_folder,sub_list{list_J,3}(1:(end-1)),dist_add,fig_tail])


    result_folder = [folder_list{list_I,1},hist_folder2_couple];
    if exist(result_folder) ~= 7
        mkdir(result_folder);
    end
    delete([result_folder,sub_list{list_J,3}(1:(end-1)),hist_tail])
    delete([result_folder,sub_list{list_J,3}(1:(end-1)),hist_link_tail])
    xlswrite([result_folder,sub_list{list_J,3}(1:(end-1)),hist_tail],foci_list2_new);
    xlswrite([result_folder,sub_list{list_J,3}(1:(end-1)),hist_link_tail],foci_list2_link);
    saveas(gcf,[result_folder,sub_list{list_J,3}(1:(end-1)),dist_add,fig_tail])
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        

%% Foci pair display with image: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    r0 = 2;
    r1 = 10;
    w0 = 1.8;
    w1 = 4;

    im1 = max_image(:,:,RNA_channel)*15;
    im2 = max_image(:,:,signal2_channel)*15;
    im3 = max_image(:,:,DAPI_channel)*1;
    im_all = cat(3,im1,im2,im3);

    figure; maximize(gcf)
        imshow(im_all)
        hold on
        viscircles(foci_list(:,[7,6]),r0*ones(size(foci_list,1),1),'Color','w','EnhanceVisibility',false,'LineWidth',w0);
        viscircles(foci_list2b(:,[7,6]),r0*ones(size(foci_list2b,1),1),'Color','w','EnhanceVisibility',false,'LineWidth',w0);
        viscircles(foci_list(I1_single,[7,6]),r1*ones(length(I1_single),1),'Color','m','EnhanceVisibility',false,'LineWidth',w0,'LineStyle',':');
        viscircles(foci_list2b(I2_single,[7,6]),r1*ones(length(I2_single),1),'Color','c','EnhanceVisibility',false,'LineWidth',w0,'LineStyle',':');
        rect_show(foci_list(I1_couple,[7,6]),foci_list2b(I2_couple,[7,6]),r1,'Color','y','LineWidth',w0,'LineStyle','-');
%%        
% %     result_folder = [folder_list{list_I,1},hist_folder_couple];
% %         saveas(gcf,[result_folder,sub_list{list_J,3}(1:(end-1)),pair_add,fig_tail])
% %     result_folder = [folder_list{list_I,1},hist_folder2_couple];
% %         saveas(gcf,[result_folder,sub_list{list_J,3}(1:(end-1)),pair_add,fig_tail])
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        
end

toc





function rect_show(xy1,xy2,r0,varargin)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A function to plot rectangles to cover every foci pairs %%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% xy1 (input): the xy coordinates of the first focus in the pair  %%
%% xy1 (input): the xy coordinates of the second focus in the pair %%
%% varargin (input): Line specs                                    %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    c0 = r0*exp(1i*pi/6);

    dxy = xy1-xy2;
    nxy = dxy./sqrt(sum(dxy.^2,2));
    
    comp1 = complex(nxy(:,1),nxy(:,2))*conj(c0);
    d1 = [real(comp1),imag(comp1)];
    comp2 = complex(nxy(:,1),nxy(:,2))*c0;
    d2 = [real(comp2),imag(comp2)];
    
    p1 = xy1+d1;
    p2 = xy1+d2;
    p3 = xy2-d1;
    p4 = xy2-d2;
    
    x_all = [p1(:,1),p2(:,1),p3(:,1),p4(:,1),p1(:,1),nan(size(p1,1),1)]';
    y_all = [p1(:,2),p2(:,2),p3(:,2),p4(:,2),p1(:,2),nan(size(p1,1),1)]';
    
    line(x_all(:),y_all(:),varargin{:})


