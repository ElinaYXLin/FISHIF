clear all
close all

tic
%% Parameter setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
folder0 = '\\fly-344e-dt-01.ad.bcm.edu\Heng 64\';
list_name = 'Duallist.xls';
in_folder = 'stacks/';
input_name = 'matchlist.xls';
mismatch_name = 'mismatch.xls';
mismatch_name2 = 'mismatch_60X.xls';
xymismatch_name = 'xymismatch.mat';
xymismatch_name2 = 'xymismatch_60X.mat';
flip_list = {'Cad'};
image_type = '*.tif';
mask_folder = 'masks/';
mask_name = 'mask.mat';
out_folder = 'Results/';
% out_folder0 = 'Results_3Dz1/';
out_folder0 = '';
hist_folder = 'Histogram/';
fit_folder = 'Histogram_A/';
hist_folder2 = 'Histogram_nullo/';
fit_folder2 = 'Histogram_default/';
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
signal2_add = '_nullo';
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
emmask_add = '_emmask';
sub_pos = [4,9];
cycle_pos0 = zeros(1,sub_pos(1));
cycle_range = [11,12,13,14];
EL_range = [0,1];
EL_range2 = [0.8,0.9];
xlim1 = [0,1];
ylim1 = [0,100];
xlim2 = [0,600];
ylim2 = [0,10000];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Data loading: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[num_list, folder_list] = xlsread([folder0,list_name]);
% folder_list = folder_list(strcmpi('b',folder_list(:,6)) | strcmpi('ba',folder_list(:,6)),:);
% folder_list = folder_list(strcmpi('b3',folder_list(:,6)),:);
folder_list = folder_list(strcmpi('h',folder_list(:,6)),:);
% folder_list = folder_list(strcmpi('T',folder_list(:,6)),:);
[N1,N2] = size(folder_list);
% run_list = {'[1]','[2]'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hf1 = figure;
maximize(hf1)
hf2 = figure;
maximize(hf2)
ha1 = zeros(0);
ha2 = zeros(0);
out_txt = cell(0);


%% Image analysis: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for list_I = 1:N1
    [sub_num, sub_list] = xlsread([folder0,folder_list{list_I,1},in_folder,input_name]);
    [M1,M2] = size(sub_list);
    
    for list_J = 1:M1%eval(run_list{list_I})
        image_folder = [folder_list{list_I,1},in_folder,sub_list{list_J,3}];
        result_folder = [folder_list{list_I,1},out_folder];
        load([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mat_tail],'nucleus_protein_profile','nucleus_RNA_profile','cytoplasmic_RNA_profile');
        
        %%% subplot organization: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        N_cycle = sub_num(list_J,13);
        I_cycle = find(N_cycle == cycle_range);
        if N_cycle > max(cycle_range)
            I_cycle = find(max(cycle_range) == cycle_range);
        elseif N_cycle < min(cycle_range)
            I_cycle = find(min(cycle_range) == cycle_range);
        end
        if ~isempty(I_cycle)
            cycle_pos0(I_cycle) = cycle_pos0(I_cycle)+1;
            sub_pos(3) = sub2ind(sub_pos([2,1]),cycle_pos0(I_cycle),I_cycle);
        else
            sub_pos(3) = 0;
        end
        %%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%% Data plot: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        figure(hf1)
        h00 = subaxis(sub_pos(1),sub_pos(2),sub_pos(3), 'Spacing', 0, 'PaddingRight',0.01, 'PaddingLeft',0.02, 'PaddingTop',0.045, 'PaddingBottom',0.04, 'Margin', 0.01);
        ha1 = [ha1,h00];
% % %             plot(nucleus_RNA_profile(:,1),nucleus_RNA_profile(:,4),'.')
% % %             xlabel('EL','FontSize',6)
% % %             ylabel('TX','FontSize',6)
% % %             title(image_folder,'FontSize',8,'Interpreter','none')
        htemp = open([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),fish_add,num_add,D3_add,figure_tail]);
        temp_axes = get(htemp,'Children');
        temp_axes = temp_axes(~strcmp('legend',get(temp_axes,'Tag')));
        copyaxes(temp_axes,h00,true);
        axes(h00);legend('off');
        close(htemp)


        figure(hf2)
        h00 = subaxis(sub_pos(1),sub_pos(2),sub_pos(3), 'Spacing', 0, 'PaddingRight',0.01, 'PaddingLeft',0.02, 'PaddingTop',0.045, 'PaddingBottom',0.04, 'Margin', 0.01);
        ha2 = [ha2,h00];
            Itrue = cytoplasmic_RNA_profile(:,1) >= EL_range(1) & cytoplasmic_RNA_profile(:,1) <= EL_range(2);
            Iposterior = cytoplasmic_RNA_profile(:,1) >= EL_range2(1) & cytoplasmic_RNA_profile(:,1) <= EL_range2(2);
            RNA0 = mean(cytoplasmic_RNA_profile(Iposterior,2));
            protein0 = mean(nucleus_protein_profile(Iposterior,2));
            xx = cytoplasmic_RNA_profile(Itrue,2)-RNA0; xx(xx<0)=0;
            yy = nucleus_protein_profile(Itrue,2)-protein0;yy(yy<0)=0;
            plot(xx,yy,'.')
            xlabel('RNA','FontSize',6)
            ylabel('Protein','FontSize',6)
            title(image_folder,'FontSize',8,'Interpreter','none')
        %%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Output: %%%============================================================
% % %         result_folder = [folder_list{list_I,1},out_folder];
% % %         if exist(result_folder) ~= 7
% % %             mkdir(result_folder);
% % %         end
% % %         saveas(gcf,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),emmask_add,figure_tail]);
% % %         save([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mat_tail],'em_mask','-append');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%         close(1)
        out_txt = cat(1,out_txt,{image_folder,N_cycle,size(nucleus_protein_profile,1)});
    end
end
linkaxes(ha1);
axes(ha1(end));
xlim(xlim1)
ylim(ylim1)

linkaxes(ha2);
axes(ha2(end));
xlim(xlim2)
ylim(ylim2)


toc


