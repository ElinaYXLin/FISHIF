clear all
tic
%% Parameter setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
list_name = 'Duallist.xls';
in_folder = 'stacks/';
input_name = 'matchlist.xls';
image_type = '*.tif';
out_folder = 'Results/';
hist_folder = 'Histogram/';
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
protein_add = '_protein';
RNA_add = '_RNA';
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
    for list_J = 1:M1
        image_folder = [folder_list{list_I,1},in_folder,sub_list{list_J,3}];
        %WGA_channel = sub_num(list_J,6);
        %DAPI_channel = sub_num(list_J,7);
        %RNA_channel = sub_num(list_J,10);
        %protein_channel = sub_num(list_J,8);
        %resolution = sub_num(list_J,9);
        %[seg_bw,cyto_bw,max_image] = nuclei_seg4(image_folder,image_type,WGA_channel,DAPI_channel,resolution);

        %foci_bw = RNA_seg(seg_bw,max_image,RNA_channel,WGA_channel,resolution,image_folder);
        result_folder = [folder_list{list_I,1},out_folder];
        load([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mat_tail]);

        r_size = 0:5;
        h = cell(size(r_size));
        nM_con = 6.02e8;
        for ir = r_size
            h{ir+1} = double(getnhood(strel('disk',ir))).*resolution.*resolution.*2*4*resolution*nM_con;
        end
        t_absolute = true;
        
        save([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mat_tail],'h','t_absolute','-append');
        
    end
end
toc