clear all
close all

tic
%% Parameter setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
list_name = 'Duallist.xls';
in_folder = 'stacks/';
old_add = '_old';
input_name = 'matchlist.xls';
mismatch_name = 'mismatch.xls';
mismatch_name2 = 'mismatch_60X.xls';
xymismatch_name = 'xymismatch.mat';
xymismatch_name2 = 'xymismatch_60X.mat';
flip_list = {'Cad'};
flip0 = true;
image_type = '*.tif';
mask_folder = 'masks/';
mask_name = 'mask.mat';
out_folder = 'Results/';
% out_folder0 = 'Results_3Dz1/';
% out_folder0 = 'Results_decross/';
out_folder0 = '';
hist_folder = 'Histogram/';
fit_folder = 'Histogram_A/';
hist_folder2 = 'Histogram_RNA2/';
fit_folder2 = 'Histogram_A_RNA2/';
fit_add = '_spot_fit';
hist_tail = '_raw.xls';
N_thresh = 3;
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
cross_add = '_cross';
sub_pos = [3,3];
cycle_pos0 = zeros(1,sub_pos(1));
cycle_range = [11,12,13];
EL_range = [0.2,0.6,0,0.1,0,0.1];

Ith = 3300;
rd = 1:1:5;
Nbin_I = 12;
ccode0 = [0.8,0.8,0.8;0,0,0;1,0,0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Data loading: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[num_list, folder_list] = xlsread(list_name);
folder_list = folder_list(strcmpi('T',folder_list(:,6)),:);
[N1,N2] = size(folder_list);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Image analysis: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for list_I = 1:N1
    [sub_num, sub_list] = xlsread([folder_list{list_I,1},in_folder,input_name]);
    [M1,M2] = size(sub_list);
    if ~isempty(out_folder0)
        copyfile([folder_list{list_I,1},out_folder],[folder_list{list_I,1},out_folder0]);
    end
    
    for list_J = 1:M1%eval(run_list{list_I})
        if isempty(strfind(sub_list{list_J,3},'_60X'))
            [~,~,mismatch_matrix] = xlsread(mismatch_name);
            load(xymismatch_name);
        else
            [~,~,mismatch_matrix] = xlsread(mismatch_name2);
            load(xymismatch_name2);
        end
        
        image_folder = [folder_list{list_I,1},in_folder,sub_list{list_J,3}];
        
        if isempty(flip0)
            flip_axis = any(cellfun(@(x) ~isempty(strfind(image_folder,x)),flip_list));   %%% Check whether the protein profile has an opposite direction (to Bcd)
        else
            flip_axis = flip0;
        end
        
        result_folder = [folder_list{list_I,1},out_folder];
        resolutionz = sub_num(list_J,11);
        signal2_channel = sub_num(list_J,12);
        Nbin = sub_num(list_J,2);
        Mdim = sub_num(list_J,3);
        load([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mat_tail]);

        N_cycle = sub_num(list_J,13);

        all_color = eval(folder_list{list_I,5});
        protein_color = all_color{protein_channel};
        RNA_color = all_color{RNA_channel};
        protein_RNA_mismatch = mismatch_matrix{strcmp(protein_color,mismatch_matrix(1,:)),strcmp(RNA_color,mismatch_matrix(:,1))};
        protein_RNA_xymismatch = {eval([protein_color,'_',RNA_color]),eval([protein_color,'_',RNA_color,'_con']),eval([protein_color,'_',RNA_color,'_x0']),Nbin,Mdim,resolution,resolution_mismatch};
        signal2_color = all_color{signal2_channel};
        protein_signal2_mismatch = mismatch_matrix{strcmp(protein_color,mismatch_matrix(1,:)),strcmp(signal2_color,mismatch_matrix(:,1))};
        protein_signal2_xymismatch = {eval([protein_color,'_',signal2_color]),eval([protein_color,'_',signal2_color,'_con']),eval([protein_color,'_',signal2_color,'_x0']),Nbin,Mdim,resolution,resolution_mismatch};

        [signal_stack0,RNA_stack0,signal2_stack0,~] = stack3D(imclearborder(seg_bw),protein_channel,DAPI_channel,RNA_channel,image_folder,protein_RNA_mismatch,signal2_channel,protein_signal2_mismatch);   %%% load 3D image stacks
        RNA_stack = double(RNA_stack0);
        signal_stack = double(signal_stack0);
        signal2_stack = double(signal2_stack0);
        
        figure('Name',['Channel crosses',image_folder])
        for ir = 1:length(rd)
%%% Extract intensity values in the foci area
            foci_data = xlsread([folder_list{list_I,1},hist_folder,sub_list{list_J,3}(1:end-1),hist_tail]);
            Itrue = foci_data(:,1) >= Ith;
            foci_xyz = round(foci_data(Itrue,6:8));
            foci_ind = sub2ind(size(signal_stack),foci_xyz(:,1),foci_xyz(:,2),foci_xyz(:,3));
            R0 = double(getnhood(strel('disk',rd(ir))));

            for I_layer = 1:size(signal_stack,3)
                RNA_stack(:,:,I_layer) = imfilter(RNA_stack0(:,:,I_layer),R0,'symmetric','conv');
                signal_stack(:,:,I_layer) = imfilter(signal_stack0(:,:,I_layer),R0,'symmetric','conv');
                signal2_stack(:,:,I_layer) = imfilter(signal2_stack0(:,:,I_layer),R0,'symmetric','conv');
            end
            foci_inten(:,1) = RNA_stack(foci_ind);
            foci_inten(:,2) = signal_stack(foci_ind);
            foci_inten(:,3) = signal2_stack(foci_ind);
            x0 = [min(foci_inten(:,1)),max(foci_inten(:,1))];

%%% Plot correlations between different channels:
            subplot(2,length(rd),sub2ind([length(rd),2],ir,1))
                plot(foci_inten(:,1),foci_inten(:,2),'.','Color',ccode0(1,:),'DisplayName','Data')
                %%%%%
                hold on
                [x_N,y_N,xerr_N,yerr_N,~,x_start,x_end] = equal_bin(foci_inten(:,1),foci_inten(:,2),Nbin_I);
                errorbar(x_N,y_N,yerr_N,'LineStyle','none','Marker','.','Color',ccode0(2,:),'DisplayName','Bin')
                %%%%%
                p1 = polyfit(x_N(2:end-1),y_N(2:end-1),1);
                plot(x0,polyval(p1,x0),'Color',ccode0(3,:),'DisplayName',['Fit: k = ',num2str(p1(1))])
                %%%%%
                xlabel(RNA_color)
                ylabel(protein_color)
                title(['r = ',num2str(rd(ir))])
                legend('show')

            subplot(2,length(rd),sub2ind([length(rd),2],ir,2))
                plot(foci_inten(:,1),foci_inten(:,3),'.','Color',ccode0(1,:),'DisplayName','Data')
                %%%%%
                hold on
                [x_N,z_N,xerr_N,zerr_N] = equal_dist(foci_inten(:,1),foci_inten(:,3),x_start,x_end);
                errorbar(x_N,z_N,zerr_N,'LineStyle','none','Marker','.','Color',ccode0(2,:),'DisplayName','Bin')
                %%%%%
                p1 = polyfit(x_N(2:end-1),z_N(2:end-1),1);
                plot(x0,polyval(p1,x0),'Color',ccode0(3,:),'DisplayName',['Fit: k = ',num2str(p1(1))])
                %%%%%
                xlabel(RNA_color)
                ylabel(signal2_color)
                title(['r = ',num2str(rd(ir))])
                legend('show')
        end
        
%%% Output: %%%============================================================
        result_folder = [folder_list{list_I,1},out_folder];
        if exist(result_folder) ~= 7
            mkdir(result_folder);
        end
        saveas(gcf,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),cross_add,figure_tail]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        
    end
end
toc