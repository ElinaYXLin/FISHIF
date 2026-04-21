clear all
close all

tic

%% Parameter setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
list_name = 'Duallist.xls';
in_folder = 'stacks/';
input_name = 'matchlist.xls';

image_type = '*.tif';
mask_folder = 'masks/';
mask_name = 'mask';
out_folder = 'Results/';
fit_folder_protein = 'Histogram_protein_A/';
fit_folder_protein2 = 'Histogram_protein_M/';
fit_folder_protein3 = 'Histogram_protein_P/';
mask2_add = '_mask_area';
fit_add = '_spot_fit';
hist_tail = '_raw.xls';
N_thresh = 3;
N_thresh2 = 0;
EL_check_name = 'RNA_stack';
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

bin0 = 20;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Data loading: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[num_list, folder_list] = xlsread(list_name);
% folder_list = folder_list(strcmpi('b',folder_list(:,6)) | strcmpi('ba',folder_list(:,6)),:);
% folder_list = folder_list(strcmpi('b3',folder_list(:,6)),:);
% folder_list = folder_list(strcmpi('h',folder_list(:,6)),:);
% folder_list = folder_list(strcmpi('b3n',folder_list(:,6)),:);
folder_list = folder_list(strcmpi('T',folder_list(:,6)),:);
% % folder_list = folder_list(strcmpi('cgal4new',folder_list(:,6)),:);
% folder_list = folder_list(strcmpi('h4k5ac',folder_list(:,6)),:);
% folder_list = folder_list(strcmpi('cbgal4',folder_list(:,6)),:);
% folder_list = folder_list(strcmpi('bgal4',folder_list(:,6)),:);
[N1,N2] = size(folder_list);
% % run_list = {[1,2,4,5]};
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Image analysis: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for list_I = 1:N1
    [sub_num, sub_list] = xlsread([folder_list{list_I,1},in_folder,input_name]);
    [M1,M2] = size(sub_list);
    
    for list_J = 1:M1%run_list{list_I}
        load([folder_list{list_I,1},mask_folder,sub_list{list_J,3},mask_name,mat_tail]);   %%% load 3D mask
        mask1 = any(mask_stack,3);
        
        load([folder_list{list_I,1},fit_folder_protein,sub_list{list_J,3}(1:end-1),mask2_add,mat_tail]);   %%% load 3D mask from the spot fitting result
        mask2 = imerode(any(~embryo_mask,3),strel('disk',2));
        
        load([folder_list{list_I,1},fit_folder_protein2,sub_list{list_J,3}(1:end-1),mask2_add,mat_tail]);   %%% load 3D mask from the spot fitting result
        mask3 = imerode(any(~embryo_mask,3),strel('disk',2));
        
        load([folder_list{list_I,1},fit_folder_protein3,sub_list{list_J,3}(1:end-1),mask2_add,mat_tail]);   %%% load 3D mask from the spot fitting result
        mask4 = imerode(any(~embryo_mask,3),strel('disk',2));
        
        if bin0 > 1
            mask1 = imresize(mask1,1/bin0);
            mask2 = imresize(mask2,1/bin0);
            mask3 = imresize(mask3,1/bin0);
            mask4 = imresize(mask4,1/bin0);
        end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        


%% Find the position of the spot fitting mask in the entire mask: %%%%%%%%%
        a = imfilter(double(mask1),double(mask2),'same','corr');
        [~,I0] = max(a(:));
        [x0,y0] = ind2sub(size(mask1),I0);
        
        b = imfilter(double(mask1),double(mask3),'same','corr');
        [~,I0] = max(b(:));
        [x1,y1] = ind2sub(size(mask1),I0);
        
        c = imfilter(double(mask1),double(mask4),'same','corr');
        [~,I0] = max(c(:));
        [x2,y2] = ind2sub(size(mask1),I0);
        
        disp([folder_list{list_I,1},mask_folder,sub_list{list_J,3}])
        disp(size(mask1)*bin0)
        disp([x0-round(size(mask2,1)/2),x0+round(size(mask2,1)/2),y0-round(size(mask2,2)/2),y0+round(size(mask2,2)/2)]*bin0)
        disp([x1-round(size(mask3,1)/2),x1+round(size(mask3,1)/2),y1-round(size(mask3,2)/2),y1+round(size(mask3,2)/2)]*bin0)
        disp([x2-round(size(mask4,1)/2),x2+round(size(mask4,1)/2),y2-round(size(mask4,2)/2),y2+round(size(mask4,2)/2)]*bin0)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    end
end
toc

