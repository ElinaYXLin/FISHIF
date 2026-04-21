function Dual_all_DAPI(total_name)

% clear all
close all

tic
%% Parameter setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
real_name = 'hb';
control_name = 'act5C';
total_folder = 'Regulation_result/';
sub_folder = 'DAPI/';
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
N_thresh = 3;
output_tail = '.xls';
figure_tail = '.fig';
figure_tail2 = '.eps';
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
p_add = '_scale';
null_add = '_null';
all_add = '_all';
con_add = '_con';
% sub_pos = [3,3];
sub_pos = [3,1];
cycle_pos0 = zeros(1,sub_pos(1));
cycle_range = [11,12,13,14];
Nthresh00 = 80;   %%% foci number threshold for enrichment averaging
th_good = 0;%0.5;
cycle_choose = [11,12];
% cycle_choose = [1:20];
Nresize = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Data loading: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[num_list, folder_list] = xlsread(list_name);
folder_list = folder_list(strcmpi('b',folder_list(:,6)) | strcmpi('ba',folder_list(:,6)),:);
[N1,N2] = size(folder_list);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if total_name(end) == '/' || total_name(end) == '\'
    total_name = total_name(1:end-1);
end
result_folder0 = [total_folder,total_name,'/'];
total_name((total_name == '/') | (total_name == '\')) = '_';

if exist([result_folder0,sub_folder]) ~= 7
    mkdir([result_folder0,sub_folder]);
end


N_plus = 0;

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
        %%% Data loading: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        image_folder = [folder_list{list_I,1},in_folder,sub_list{list_J,3}];
        result_folder = [folder_list{list_I,1},out_folder];
        load([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mat_tail],'max_image','DAPI_channel');   %%% load analysis result
        N_plus = N_plus+1;
        %%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%% subplot organization: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        N_cycle = sub_num(list_J,13);
        I_plus = mod(N_plus,sub_pos(1)*sub_pos(2));
        if I_plus <= 0
            I_plus = sub_pos(1)*sub_pos(2);
        elseif I_plus == 1
            figure
            set(gcf,'Units','inches')
            set(gcf,'Position',[4,0.5,8.5*Nresize,11*Nresize])
        end
        %%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%% Nuclear image plot: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        subaxis(sub_pos(1),sub_pos(2),I_plus, 'Spacing', 0, 'PaddingRight',0.04, 'PaddingLeft',0.04, 'PaddingTop',0.03, 'PaddingBottom',0.03, 'MarginLeft', 0.04, 'MarginRight', 0.005, 'MarginTop', 0.015, 'MarginBottom', 0.025);
        im_out = max_image(:,:,DAPI_channel);
        imshow(double(im_out)/max(double(im_out(:))));
        title([image_folder,': Number ',num2str(N_plus),', Cycle ',num2str(N_cycle)],'Interpreter','none','FontName','Arial','FontSize',10,'FontWeight','bold');
        %%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %% Output: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if I_plus == sub_pos(1)*sub_pos(2)
            saveas(gcf,[result_folder0,sub_folder,total_name,DAPI_add,'_',num2str(N_plus/sub_pos(1)/sub_pos(2)+1,'%2d'),figure_tail]);
            set(gcf,'PaperPositionMode','auto')
%             export_fig([result_folder0,sub_folder,total_name,DAPI_add,'_',num2str(N_plus/sub_pos(1)/sub_pos(2)+1,'%2d'),figure_tail2],'-nocrop','-transparent',gcf)
%             saveas(gcf,[result_folder0,sub_folder,total_name,DAPI_add,'_',num2str(N_plus/sub_pos(1)/sub_pos(2)+1,'%2d'),figure_tail2],'epsc');
            print
            close(gcf)
        end
    end
end

if I_plus ~= sub_pos(1)*sub_pos(2)
    saveas(gcf,[result_folder0,sub_folder,total_name,DAPI_add,'_',num2str(ceil(N_plus/sub_pos(1)/sub_pos(2))+1,'%2d'),figure_tail]);
    set(gcf,'PaperPositionMode','auto')
%     export_fig([result_folder0,sub_folder,total_name,DAPI_add,'_',num2str(ceil(N_plus/sub_pos(1)/sub_pos(2))+1,'%2d'),figure_tail2],'-nocrop','-transparent',gcf)
%     saveas(gcf,[result_folder0,sub_folder,total_name,DAPI_add,'_',num2str(ceil(N_plus/sub_pos(1)/sub_pos(2))+1,'%2d'),figure_tail2],'epsc');
    print
end

toc

end


