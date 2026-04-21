clear all
tic
%% Parameter setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
image_folder = 'simu_enrich\';
load_type = '*.mat';
image_type = '*.tif';
out_folder = 'Results\';
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
rad_add = '_rad';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Data loading: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
result_folder = [image_folder,out_folder];
simu_info = dir([result_folder,load_type]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Image analysis: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for list_J = 1:length(simu_info)
    load_name = simu_info(list_J).name;
    image_folder = [result_folder,load_name];
    load([result_folder,load_name]);
    [foci_data,fake_data,h,t_absolute] = dual_local_local(foci_bw00,max_image00,SS,foci_layer,mask_stack,signal_stack,RNA_stack,nucleus_protein_profile,image_folder,N_cycle,resolution,resolutionz);

    saveas(73,[result_folder,load_name(1:(length(load_name)-4)),protein_add,local_add,foci_add,figure_tail]);
    saveas(74,[result_folder,load_name(1:(length(load_name)-4)),protein_add,local_add,hist_add,figure_tail]);
    saveas(75,[result_folder,load_name(1:(length(load_name)-4)),protein_add,local_add,RNA_add,figure_tail]);
    saveas(76,[result_folder,load_name(1:(length(load_name)-4)),protein_add,local_add,abs_add,foci_add,figure_tail]);
    saveas(77,[result_folder,load_name(1:(length(load_name)-4)),protein_add,local_add,abs_add,hist_add,figure_tail]);
    saveas(80,[result_folder,load_name(1:(length(load_name)-4)),protein_add,local_add,rad_add,figure_tail]);
    saveas(81,[result_folder,load_name(1:(length(load_name)-4)),protein_add,local_add,abs_add,rad_add,figure_tail]);

    save([result_folder,load_name],'foci_data','fake_data','h','t_absolute','-append');

    if ~isempty(nucleus_protein_profile)
       xlswrite([result_folder,load_name(1:(length(load_name)-4)),protein_add,nu_add,output_tail],nucleus_protein_profile);
    end
    if ~isempty(foci_data)
        for I_data = 1:length(foci_data)
            if ~isempty(foci_data{I_data})
                xlswrite([result_folder,load_name(1:(length(load_name)-4)),protein_add,local_add,foci_add,output_tail],foci_data{I_data},I_data);
            end
        end
    end
    if ~isempty(fake_data)
        for I_data = 1:length(fake_data)
            if ~isempty(fake_data{I_data})
                xlswrite([result_folder,load_name(1:(length(load_name)-4)),protein_add,local_add,fake_add,foci_add,output_tail],fake_data{I_data},I_data);
            end
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        

end
toc