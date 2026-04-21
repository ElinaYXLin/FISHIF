clear all

%% Parameter setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
list_name = 'Duallist_temp.xls';
in_folder = 'stacks/';
input_name = 'matchlist.xls';
standard_record = 'Calibration/Results/standard.mat';
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
immu_add = '_immu';
fluc_add = '_fluc';
bino_add = '_bino';
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
    channel_name = eval(folder_list{list_I,5});
    for list_J = 1:M1
        image_folder = [folder_list{list_I,1},in_folder,sub_list{list_J,3}];
        Mdim = sub_num(list_J,3);
        Nbin1 = sub_num(list_J,2);
        WGA_channel = sub_num(list_J,6);
        DAPI_channel = sub_num(list_J,7);
        RNA_channel = sub_num(list_J,10);
        protein_channel = sub_num(list_J,8);
        resolution = sub_num(list_J,9);
        resolutionz = sub_num(list_J,11);

        %[seg_bw,cyto_bw,max_image,N_cycle,WGA_th0,DAPI_th0] = nuclei_seg4(image_folder,image_type,WGA_channel,DAPI_channel,resolution);
        [seg_bw,cyto_bw,max_image,N_cycle,WGA_th0,DAPI_th0] = nuclei_segDAPI(image_folder,image_type,DAPI_channel,resolution);
        [mask_stack,signal_stack,~] = mask3D(seg_bw,protein_channel,DAPI_channel,RNA_channel,image_folder);
%         max_size0 = size(max_image);
%         max_size0(Mdim) = size(max_image,Mdim)/Nbin1;
%         imcorr1 = corr_mask(max_size0,channel_name,resolution);
%         Ntile = ones(size(max_size0));
%         Ntile(Mdim) = Nbin1;
%         imcorr = repmat(imcorr1,Ntile);
%         max_image = uint16(double(max_image)./imcorr);
        
        [foci_bw,psize0,g_threshold0] = RNA_seg(seg_bw,max_image,RNA_channel,WGA_channel,resolution,image_folder,N_cycle);
        [nucleus_RNA_profile,foci_RNA_profile,cytoplasmic_RNA_profile] = RNA_profile(seg_bw,cyto_bw,foci_bw,max_image,RNA_channel,resolution,image_folder,N_cycle);
%         [nucleus_protein_profile,cytoplasmic_protein_profile] = protein_profile(seg_bw,cyto_bw,max_image,signal_channel,image_folder,N_cycle)
        [nucleus_protein_profile,cytoplasmic_protein_profile] = protein_profile3(seg_bw,cyto_bw,max_image,protein_channel,mask_stack,signal_stack,image_folder,N_cycle,resolution);
        clear mask_stack signal_stack
        dual_profile(nucleus_protein_profile,nucleus_RNA_profile,foci_RNA_profile,image_folder,N_cycle);

        figure(16)
            temp_image = zeros(size(max_image,1),size(max_image,2),3);
            temp_image(:,:,2) = double(max_image(:,:,protein_channel))./max(max(double(max_image(:,:,protein_channel))));
            imshow(temp_image);
            title(['Immunofluorescence signal: ',image_folder,', cycle = ',num2str(N_cycle)],'Interpreter','none');
        figure(17)
            imshow(double(max_image(:,:,RNA_channel))./max(max(double(max_image(:,:,RNA_channel)))));
            title(['FISH signal: ',image_folder,', cycle = ',num2str(N_cycle)],'Interpreter','none');
        
        
%%% Output: %%%============================================================
        result_folder = [folder_list{list_I,1},out_folder];
        if exist(result_folder) ~= 7
            mkdir(result_folder);
        end
        saveas(51,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),nu_add,seg_add,figure_tail]);
        saveas(2,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,figure_tail]);
        saveas(3,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,th_add,figure_tail]);
        saveas(4,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,seg_add,figure_tail]);
        saveas(5,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),fish_add,int_add,figure_tail]);
        saveas(6,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),fish_add,num_add,figure_tail]);
        saveas(7,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,num_add,figure_tail]);
        saveas(8,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,int_add,figure_tail]);
        saveas(9,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,int_add,cmp_add,figure_tail]);
        saveas(10,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,th_add,sel_add,figure_tail]);
        saveas(11,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),nu_add,fate_add,figure_tail]);
        saveas(12,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),nu_add,reg_add,figure_tail]);
        saveas(13,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,reg_add,figure_tail]);
        saveas(14,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,fate_add,figure_tail]);
        saveas(15,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),nu_add,foci_add,reg_add,figure_tail]);
        saveas(16,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),immu_add,figure_tail]);
        saveas(17,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),fish_add,figure_tail]);
        saveas(62,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,fluc_add,figure_tail]);
        saveas(63,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,bino_add,figure_tail]);
        save([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mat_tail],'image_folder','image_type','WGA_channel','DAPI_channel','RNA_channel','protein_channel','resolution','seg_bw','cyto_bw','max_image','foci_bw','nucleus_RNA_profile','foci_RNA_profile','cytoplasmic_RNA_profile','nucleus_protein_profile','cytoplasmic_protein_profile','N_cycle','psize0','g_threshold0','WGA_th0','DAPI_th0','resolutionz');
        
        if ~isempty(nucleus_RNA_profile)
            xlswrite([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),RNA_add,nu_add,output_tail],nucleus_RNA_profile);
        end
        if ~isempty(cytoplasmic_RNA_profile)
            xlswrite([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),RNA_add,cyto_add,output_tail],cytoplasmic_RNA_profile);
        end
        if ~isempty(foci_RNA_profile)
            xlswrite([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),RNA_add,foci_add,output_tail],foci_RNA_profile);
        end
        if ~isempty(nucleus_protein_profile)
            xlswrite([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,nu_add,output_tail],nucleus_protein_profile);
        end
        if ~isempty(cytoplasmic_protein_profile)
            xlswrite([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,cyto_add,output_tail],cytoplasmic_protein_profile);
        end
        sub_num(list_J,13) = N_cycle;
        
        clear ('seg_bw','cyto_bw','max_image','foci_bw','nucleus_RNA_profile','foci_RNA_profile','cytoplasmic_RNA_profile','nucleus_protein_profile','cytoplasmic_protein_profile','N_cycle','WGA_th0','DAPI_th0');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        
    end
    xlswrite([folder_list{list_I,1},in_folder,input_name],cat(2,sub_list,num2cell(sub_num)));
end

% clear global standard_data