clear all
close all

tic
%% Parameter setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
list_name = 'RNA2list_liuliu.xls';
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
sub_pos = [3,3];
cycle_pos0 = zeros(1,sub_pos(1));
cycle_range = [11,12,13];
xls_header = {'N','x','y','z'};

xls_header2 = {'File name','Mean','Std'};
out_name0 = 'foci_N';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Data loading: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[num_list, folder_list] = xlsread(list_name);
folder_list = folder_list(strcmpi('T',folder_list(:,6)),:);
[N1,N2] = size(folder_list);
% run_list = {'[1]','[2]'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Image analysis: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for list_I = 1:N1
    [sub_num, sub_list] = xlsread([folder_list{list_I,1},in_folder,input_name]);
    [M1,M2] = size(sub_list);
    IntenN_all = cell(M1,3);
    
    for list_J = 1:M1%eval(run_list{list_I})
        image_folder = [folder_list{list_I,1},in_folder,sub_list{list_J,3}];
        result_folder = [folder_list{list_I,1},out_folder];
        resolutionz = sub_num(list_J,11);
        Nbin = sub_num(list_J,2);
        Mdim = sub_num(list_J,3);

        load([folder_list{list_I,1},fit_folder,sub_list{list_J,3}(1:end-1),fit_add,mat_tail]);   %%% load spot intensity fitting result
        Inten_thresh = b*N_thresh;   %%% set foci intensity threshold
        single_Inten = b;
        
        raw_data = xlsread([folder_list{list_I,1},hist_folder,sub_list{list_J,3}(1:(end-1)),hist_tail]);
        Inten0 = raw_data(:,1).*raw_data(:,2).*raw_data(:,3)*2*pi;
        Itrue = Inten0 >= Inten_thresh;
        IntenN = Inten0(Itrue)/b;
        IntenXYZ = raw_data(Itrue,6:8);
        
        IntenN_all(list_J,:) = {image_folder,mean(IntenN),std(IntenN)};
        
        figure(1)
        hist(IntenN,0:5:100);
        xlabel('Nascent mRNAs per locus')
        ylabel('#')
        title([folder_list{list_I,1},sub_list{list_J,3}(1:end-1),': mean = ',num2str(mean(IntenN))])
%%% Output: %%%============================================================
        result_folder = [folder_list{list_I,1},out_folder];
        if exist(result_folder) ~= 7
            mkdir(result_folder);
        end
        saveas(1,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,figure_tail]);
        
        if ~isempty(IntenN)
           xlswrite([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,output_tail],cat(1,xls_header,num2cell([IntenN,IntenXYZ])));
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        
    end
    xlswrite([folder_list{list_I,1},out_name0,output_tail],cat(1,xls_header2,IntenN_all));
    
end
toc