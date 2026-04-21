clear

%% Parameter setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
list_name = 'Duallist.xls';
in_folder = 'stacks/';
input_name = 'matchlist.xls';
image_type = '*.tif';
out_folder = 'Results/';
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
shell_add = '_shell';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Data loading: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[num_list, folder_list] = xlsread(list_name);
[N1,N2] = size(folder_list);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Image analysis: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for list_I = 1:N1
    [sub_num, sub_list] = xlsread([folder_list{list_I,1},in_folder,input_name]);
    [M1,M2] = size(sub_list);
    for list_J = 1:14%M1
        result_folder = [folder_list{list_I,1},out_folder];
        load([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mat_tail]);
        
        embryo_region = false(size(seg_bw));
        if any(any(seg_bw))
            temp_label = 2*seg_bw;
            temp_label(1,1) = 1;
            embryo_prop = regionprops(temp_label,'ConvexImage','BoundingBox');
            x1 = uint16(embryo_prop(2).BoundingBox(2));
            x2 = uint16(x1+embryo_prop(2).BoundingBox(4)-1);
            y1 = uint16(embryo_prop(2).BoundingBox(1));
            y2 = uint16(y1+embryo_prop(2).BoundingBox(3)-1);
            embryo_region(x1:x2,y1:y2) = embryo_prop(2).ConvexImage;
        end
        
        embryo_shell = embryo_region & (~imerode(embryo_region, strel('disk',floor(100))));
        new_prop = regionprops(seg_bw,embryo_shell,'MaxIntensity');
        keepIdx = find([new_prop.MaxIntensity] > 0);
        seg_new = logical(ismember(bwlabel(seg_bw),keepIdx));
        cyto_new = embryo_shell & (~seg_new);
        cyto_new = imerode(cyto_new, strel('disk',floor(5)));
        [nucleus_protein_profile,cytoplasmic_protein_profile] = protein_profile(seg_new,cyto_new,max_image,protein_channel,image_folder,N_cycle);
            
%%% Output: %%%============================================================
        result_folder = [folder_list{list_I,1},out_folder];
        if exist(result_folder) ~= 7
            mkdir(result_folder);
        end
        saveas(2,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,shell_add,figure_tail]);
        saveas(14,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,shell_add,fate_add,figure_tail]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        
    end
end