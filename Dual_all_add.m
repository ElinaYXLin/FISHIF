function Dual_all_add(total_name)

% clear all
close all

tic
%% Parameter setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
real_name = 'hb';
control_name = 'nullo';
total_folder = 'Regulation_result/';
sub_folder = 'Individual/';
list_name = 'Duallist.xls';
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
% sub_pos = [3,3];
sub_pos = [3,9];
cycle_pos0 = zeros(1,sub_pos(1));
cycle_range = [11,12,13];
Nthresh00 = 80;   %%% foci number threshold for enrichment averaging
th_good = 0;%0.5;
cycle_choose = [12,13,14];
% cycle_choose = [1:20];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Data loading: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[num_list, folder_list] = xlsread(list_name);
folder_list = folder_list(strcmpi('b',folder_list(:,6)),:);
[N1,N2] = size(folder_list);

%%%%%%%%%%%%%%%%%%%
figure(30);clf

%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%
maximize(30)

%%%%%%%%%%%%%%%%%%%

out_para = zeros(0);
out_name = cell(0);
out_title = {'File name','Cycle','Imax','lambda_I','Imin','NImax','Cmax','lambda_C','Cmin','NCmax','NCmax_par','H_I','C0_I','Amax_I','Amin_I','H_N','C0_N','Amax_N','Amin_N','H_A','C0_A','Pmax_A','Amin_A'};
N_plus = 0;

RNAI_out = zeros(0);
RNAN_out = zeros(0);
null_out = zeros(0);
embryo_info = zeros(0);
Nccode = 6;
color_code = mycolors(Nccode);
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
    [sub_num, sub_list] = xlsread([folder_list{list_I,1},in_folder,input_name]);
    [M1,M2] = size(sub_list);

    for list_J = 1:M1        
        %%% Data loading: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        image_folder = [folder_list{list_I,1},in_folder,sub_list{list_J,3}];
        result_folder = [folder_list{list_I,1},out_folder];
        load([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mat_tail]);   %%% load analysis result
        N_cycle = sub_num(list_J,13);

        %%% Nuclear image plot: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        figure(30)
        im_out = max_image(:,:,DAPI_channel);
        imshow(double(im_out)/max(double(im_out(:))));
        title([image_folder,': cycle ',num2str(N_cycle)],'interpreter','none');
        %%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%% Image output for individual embryos: %%%%%%%%%%%%%%%%%%%%%%%%%%
        image_name = image_folder(1:end-1);
        image_name((image_name == '/') | (image_name == '\')) = '_';
        
        saveas(30,[result_folder0,sub_folder,image_name,DAPI_add,figure_tail]);

    end
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


