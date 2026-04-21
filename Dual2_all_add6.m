function Dual2_all_add6(total_name)

% clear all
close all

tic
%% Parameter setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
real_name = 'hb';
control_name = 'nullo';
total_folder = 'Regulation_result/';
sub_folder = 'Individual/';
list_name = 'Dual2list.xls';
in_folder = 'stacks/';
input_name = 'matchlist.xls';
mismatch_name = 'mismatch.xls';
mismatch_name2 = 'mismatch_60X.xls';
image_type = '*.tif';
mask_folder = 'masks/';
mask_name = 'mask.mat';
out_folder = 'Results/';
hist_folder = 'Histogram/';
fit_folder = 'Histogram_A/';
fit_folder_protein = 'Histogram_protein_A/';
hist_folder2 = 'Histogram_nullo/';
fit_folder2 = 'Histogram_default/';
fit_add = '_spot_fit';
hist_tail = '_raw.xls';
N_thresh = 10;
output_tail = '.xls';
figure_tail = '.fig';
figure_tail2 = '.jpg';
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
EL_add = '_EL';
D3_add = '_3D';
compare_add = '_compare';
DAPI_add = '_DAPI';
mean_add = '_mean';
slope_add = '_slope';
average_add = '_average';
spot_add = '_spot';
TX_add = '_TX';
var_add = '_var';
abs_add = '_abs';
para_add = '_para';
null_add = '_null';
all_add = '_all';
con_add = '_con';
layer_add= '_layer';
% sub_pos = [3,3];
sub_pos = [4,11];
cycle_pos0 = zeros(1,sub_pos(1));
cycle_range = [11,12];
Nthresh00 = 80;   %%% foci number threshold for enrichment averaging
th_good = 0;%0.5;
cycle_choose = [11,12,13,14];
% cycle_choose = [1:20];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(140);clf
maximize(140)

%% Data loading: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[num_list, folder_list] = xlsread(list_name);
folder_list = folder_list(strcmpi('GB',folder_list(:,6)),:);
[N1,N2] = size(folder_list);

%%%%%%%%%%%%%%%%%%%

out_data = zeros(0);
out_name = cell(0);
N_plus = 0;
flip_EL = false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if total_name(end) == '/' || total_name(end) == '\'
    total_name = total_name(1:end-1);
end
result_folder0 = [total_folder,total_name,'/'];
total_name((total_name == '/') | (total_name == '\')) = '_';

if exist([result_folder0,sub_folder]) ~= 7
    mkdir([result_folder0,sub_folder]);
end



%% Image analysis: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for list_I = 1:N1
    [sub_num, sub_list, sub_all] = xlsread([folder_list{list_I,1},in_folder,input_name]);
    if size(sub_all,2) >= 17
        I_se = ~strcmp('F',sub_all(:,17));
        sub_num = sub_num(I_se,:);
        sub_list = sub_list(I_se,:);
    end
    [M1,M2] = size(sub_list);

    for list_J = 1:M1        
        N_plus = N_plus+1;
        %%% subplot organization: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        N_cycle = sub_num(list_J,13);
        %%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Data loading: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        image_folder = [folder_list{list_I,1},in_folder,sub_list{list_J,3}];
        result_folder = [folder_list{list_I,1},out_folder];
        load([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mat_tail],'nucleus_protein_profile','nucleus_protein_profile_ab','resolution');   %%% load analysis result
        load([folder_list{list_I,1},mask_folder,sub_list{list_J,3},mask_name]);   %%% load 3D mask
        
        protein_profile3D_short3(nucleus_protein_profile,mask_stack,['Sample ',num2str(N_plus)],N_cycle);   %%% protein profile plot (2,22,14,62)

        %%% Image output for individual embryos: %%%%%%%%%%%%%%%%%%%%%%%%%%
        image_name = image_folder(1:end-1);
        image_name((image_name == '/') | (image_name == '\')) = '_';
        saveas(140,[result_folder0,sub_folder,image_name,protein_add,layer_add,D3_add,figure_tail]);
    end
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


toc
