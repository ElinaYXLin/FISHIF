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
rad_add = '_rad';
scatter_add = '_scatter';
N_theta = 20;
N_rho = 20;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Data loading: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[num_list, folder_list] = xlsread(list_name);
folder_list = folder_list(strcmpi('T',folder_list(:,6)),:);
[N1,N2] = size(folder_list);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Image analysis: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for list_I = 3:N1
    [sub_num, sub_list] = xlsread([folder_list{list_I,1},in_folder,input_name]);
    [M1,M2] = size(sub_list);
    for list_J = 1:M1
        image_folder = [folder_list{list_I,1},in_folder,sub_list{list_J,3}];
        result_folder = [folder_list{list_I,1},out_folder];
        load([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mat_tail]);
        [mask_stack,signal_stack,RNA_stack,DAPI_stack] = mask3D(imclearborder(seg_bw),protein_channel,DAPI_channel,RNA_channel,image_folder);
        [foci_bw00,~,~,~] = modified_foci([folder_list{list_I,1},hist_folder,sub_list{list_J,3}(1:(end-1)),hist_tail],max_image,RNA_channel,imclearborder(seg_bw),mask_stack,N_cycle,image_folder);
        close(4)
        bwn = bwlabel(seg_bw);
        nu_prop = regionprops(logical(bwn),'Centroid');
        nu_xy = cell2mat({nu_prop.Centroid}');
        foci_I = find(foci_bw00);
        [foci_x0,foci_y0] = ind2sub(size(foci_bw00),foci_I);
        foci_x = foci_x0-nu_xy(bwn(foci_I),2);
        foci_y = foci_y0-nu_xy(bwn(foci_I),1);
        [foci_theta,foci_rho] = cart2pol(foci_x,foci_y);
        
        figure(1)
            subplot(1,3,1)
                plot(foci_x,foci_y,'.')
                xlabel('x (pixel)')
                ylabel('y (pixel)')
                title(['Foci scatter plot: ',image_folder],'Interpreter','none')
            subplot(1,3,2)
                [n_theta,x_theta] = hist(foci_theta,N_theta);
                plot(x_theta/pi*180,n_theta./sum(n_theta)*100,'-*')
                xlabel('theta (degree)')
                ylabel('%')
                title(['Foci angular distribution: ',image_folder],'Interpreter','none')            
            subplot(1,3,3)
                [n_rho,x_rho] = hist(foci_rho,N_rho);
                plot(x_rho,n_rho./x_rho./sum(n_rho)*100,'-o')
                xlabel('rho (pixel)')
                ylabel('%/pixel')
                title(['Foci radial distribution: ',image_folder],'Interpreter','none')  
            
%%% Output: %%%============================================================
        result_folder = [folder_list{list_I,1},out_folder];
        if exist(result_folder) ~= 7
            mkdir(result_folder);
        end
        saveas(1,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,scatter_add,figure_tail]);

        
        clear ('seg_bw','cyto_bw','max_image','foci_bw00','foci_bw','nucleus_RNA_profile','foci_RNA_profile','cytoplasmic_RNA_profile','nucleus_protein_profile','cytoplasmic_protein_profile','N_cycle','foci_data','fake_data');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        
    end
end
toc