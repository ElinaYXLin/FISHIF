clear all
close all

tic
xy_mismatch_check = true;
epsilonxy = 0.55;
epsilonxyz = 0.55;
epsilonz = 1;

%% Parameter setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
list_name = 'Duallist.xls';
in_folder = 'stacks/';
% old_add = '_old';
old_add = '';
input_name = 'matchlist.xls';
mismatch_name = 'mismatch.xls';
mismatch_name2 = 'mismatch_60X.xls';
% % mismatch_name3 = 'mismatch_60X_FA_02142020.xls';
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
% out_folder0 = 'Results0/';
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
r_range = 0:0.02:3;
% EL_range = [0.2,0.6,0,0.1,0,0.1];
EL_range = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Data loading: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[num_list, folder_list] = xlsread(list_name);
% folder_list = folder_list(strcmpi('b',folder_list(:,6)) | strcmpi('ba',folder_list(:,6)),:);
% folder_list = folder_list(strcmpi('b3',folder_list(:,6)),:);
% folder_list = folder_list(strcmpi('h',folder_list(:,6)),:);
% folder_list = folder_list(strcmpi('b3n',folder_list(:,6)),:);
folder_list = folder_list(strcmpi('T',folder_list(:,6)),:);
% folder_list = folder_list(strcmpi('cbgal4',folder_list(:,6)),:);
% folder_list = folder_list(strcmpi('bgal4',folder_list(:,6)),:);
[N1,N2] = size(folder_list);
% run_list = {'[1]','[2]'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Image analysis: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for list_I = 1:N1
    [sub_num, sub_list] = xlsread([folder_list{list_I,1},in_folder,input_name]);
    [M1,M2] = size(sub_list);
    if ~isempty(out_folder0)
        copyfile([folder_list{list_I,1},out_folder],[folder_list{list_I,1},out_folder0]);
    end
    
    for list_J = 1:M1%eval(run_list{list_I})
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
        
        load([folder_list{list_I,1},mask_folder,sub_list{list_J,3},mask_name]);   %%% load 3D mask
        z_size = size(mask_stack,3);
% %         load([folder_list{list_I,1},fit_folder,sub_list{list_J,3}(1:end-1),fit_add,mat_tail]);   %%% load spot intensity fitting result
        b = 1;
        Inten_thresh = b*N_thresh;   %%% set foci intensity threshold
        single_Inten = b;
% %         load([folder_list{list_I,1},fit_folder2,sub_list{list_J,3}(1:end-1),fit_add,mat_tail]);   %%% load spot intensity fitting result
        Inten_thresh2 = b*N_thresh;   %%% set foci intensity threshold
        single_Inten2 = b;
%         N_cycle = round(log2(max(mask_stack(:))))+2;   %%% Calculate the nuclear cycle number
        N_cycle = sub_num(list_J,13);

        all_color = eval(folder_list{list_I,5});
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
        dxyz = pdist2(xyz1,xyz2);
        [d12min,I2min] = min(dxyz,[],2);
        
        I1_couple = find(d12min <= epsilonxyz);
        I2_couple = I2min(I1_couple);
        I1_single = find(d12min > epsilonxyz & 2*pi*prod(foci_list(:,1:3),2) >= Inten_thresh);
            I2_is_single = true(size(foci_list2,1),1);
            I2_is_single(I2_couple) = false;
        I2_single = find(I2_is_single & 2*pi*prod(foci_list2(:,1:3),2) >= Inten_thresh2);
        
        N2_couple = unique(I2_couple);
        foci_list_new = zeros(length(N2_couple)+length(I1_single)+length(I2_single),size(foci_list,2));
        foci_list2_new = zeros(length(N2_couple)+length(I1_single)+length(I2_single),size(foci_list,2));
        for ii = 1:length(N2_couple)
            I2_couple0 = N2_couple(ii);
            I1_couple0 = I1_couple(I2_couple == I2_couple0);
            if length(I1_couple0) == 1
                foci_list_new(ii,:) = foci_list(I1_couple0,:);
            else
                w1 = prod(foci_list(I1_couple0,1:3),2);
                w0 = w1/sum(w1);
                foci_list_new(ii,:) = sum(foci_list(I1_couple0,:).*repmat(w0,1,size(foci_list,2)));
                foci_list_new(ii,1) = sum(w1)/foci_list_new(ii,2)/foci_list_new(ii,3);
            end
            foci_list2_new(ii,:) = foci_list2b(I2_couple0,:);
        end
        foci_list_new((length(N2_couple)+1):(length(N2_couple)+length(I1_single)),:) = foci_list(I1_single,:);
        foci_list_new((length(N2_couple)+length(I1_single)+1):end,:) = 0;
        foci_list_new((length(N2_couple)+length(I1_single)+1):end,6:8) = foci_list2(I2_single,6:8);
        foci_list_new(:,8) = round(foci_list_new(:,8));
       
        foci_list2_new((length(N2_couple)+1):(length(N2_couple)+length(I1_single)),:) = 0;
        foci_list2_new((length(N2_couple)+1):(length(N2_couple)+length(I1_single)),6:8) = foci_listb(I1_single,6:8);
        foci_list2_new((length(N2_couple)+length(I1_single)+1):end,:) = foci_list2b(I2_single,:);
        

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
        
%%% Output: %%%============================================================
        figure(1)
        nd = hist(d12min,r_range);
        bar(r_range(1:end-1),nd(1:end-1)/sum(nd));
        xlabel('r (\mum)')
        ylabel('%')

        result_folder = [folder_list{list_I,1},hist_folder_couple];
        if exist(result_folder) ~= 7
            mkdir(result_folder);
        end
        delete([result_folder,sub_list{list_J,3}(1:(end-1)),hist_tail])
        delete([result_folder,sub_list{list_J,3}(1:(end-1)),hist_link_tail])
        xlswrite([result_folder,sub_list{list_J,3}(1:(end-1)),hist_tail],foci_list_new);
        xlswrite([result_folder,sub_list{list_J,3}(1:(end-1)),hist_link_tail],foci_list_link);
        saveas(1,[result_folder,sub_list{list_J,3}(1:(end-1)),dist_add,fig_tail])
        

        result_folder = [folder_list{list_I,1},hist_folder2_couple];
        if exist(result_folder) ~= 7
            mkdir(result_folder);
        end
        delete([result_folder,sub_list{list_J,3}(1:(end-1)),hist_tail])
        delete([result_folder,sub_list{list_J,3}(1:(end-1)),hist_link_tail])
        xlswrite([result_folder,sub_list{list_J,3}(1:(end-1)),hist_tail],foci_list2_new);
        xlswrite([result_folder,sub_list{list_J,3}(1:(end-1)),hist_link_tail],foci_list2_link);
        saveas(1,[result_folder,sub_list{list_J,3}(1:(end-1)),dist_add,fig_tail])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        
    end
end
toc