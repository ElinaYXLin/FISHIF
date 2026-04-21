function Dualprocess0_specified(lf_name,sp_name)
%% Parameter setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(lf_name)
    list_name=lf_name;
else
    list_name = 'Duallist.xls';
end
in_folder = 'stacks/';
input_name = 'matchlist.xls';
time_name='_time';
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
emmask_add = '_emmask';
double_add = 'D';
% size0 = 5000;
bin0 = 3;
Smin = 100;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Data loading: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isa(lf_name,'char') && strcmp(lf_name(end-3:end),'.xls')
    [num_list, folder_list] = xlsread(list_name);
    folder_list = folder_list(strcmpi('T',folder_list(:,6)),:);
%     pathname=replace(lf_name,'Duallist.xls','');
%     folder_new=cellfun(@(x) [pathname,x],folder_list(:,1),'UniformOutput',false);
    folder_new=folder_list;
elseif isa(lf_name,'cell')
    folder_list = lf_name;
end
% folder_list(:,1)=folder_new;
[N1,N2] = size(folder_list);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Image analysis: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for list_I = 1:N1
    folder_name=folder_list{list_I,1};
    try
        raw = readcell([folder_name,in_folder,input_name]);
    catch
        folder_name=folder_new{list_I,1};
        raw = readcell([folder_name,in_folder,input_name]);
    end
    sub_num=cell2mat(raw(:,cellfun(@isnumeric,raw(1,:))));
    sub_list=raw;
    sub_list(:,~cellfun(@ischar,raw(1,:)))={''};
    [M1,M2] = size(sub_list);
%     channel_name = eval(folder_list{list_I,5});
    
    list_J_all=1:M1;
    list_J_all=list_J_all(cellfun(@(x) strcmp(x(1:end-1),sp_name(1:end-1)),sub_list(:,3)));
    
    for list_J = list_J_all
        image_folder = [folder_name,in_folder,sub_list{list_J,3}];
        if isempty(sub_list{list_J,3})
            image_folder = [folder_name,in_folder,sub_list{list_J,1},'tile01\'];
        end
        
        if ~isempty(strfind(sub_list{list_J,2},double_add))
            Nbin = sub_num(list_J,2:3);
            Mdim = 1:2;
        else
            Mdim = sub_num(list_J,3);
            Nbin = sub_num(list_J,2);
        end
        
        WGA_channel = sub_num(list_J,6);
        DAPI_channel = sub_num(list_J,7);
        RNA_channel = sub_num(list_J,10);
        protein_channel = sub_num(list_J,8);
        signal2_channel = sub_num(list_J,12);
        resolution = sub_num(list_J,9);
        resolutionz = sub_num(list_J,11);
        try
            all_color = eval(folder_list{list_I,5});
        catch
            all_color=[];
        end
        
        time_num = size(dir([image_folder,['*stack01','*.tif']]),1);
        for time_index=1:time_num
            figure_tail=[time_name,num2str(time_index,'%04u'),'.fig'];
            mat_tail=[time_name,num2str(time_index,'%04u'),'.mat'];
            temp_name = dir([image_folder,'time',num2str(time_index,'%04u'),image_type]);
            if size(temp_name,1)==0
                temp_name = dir([image_folder,image_type]);
            end
            temp0 = imread([image_folder,temp_name(1).name]);
    % %         if size0 > 0
    % %             scale0 = min(size0/max(size(temp0)),1);
    % %         else
    % %             scale0 = 1;
    % %         end
            if bin0 > 1
                scale0 = 1/bin0;
            else
                scale0 = 1;
            end
            clear temp_name temp0

            %[seg_bw,cyto_bw,max_image,N_cycle,WGA_th0,DAPI_th0] = nuclei_seg4(image_folder,image_type,WGA_channel,DAPI_channel,resolution);
            try
                [seg_bw,cyto_bw,max_image,N_cycle,WGA_th0,DAPI_th0] = nuclei_segDAPI(image_folder,['time',num2str(time_index,'%04u'),image_type],DAPI_channel,resolution,scale0,Smin);
            catch
                [seg_bw,cyto_bw,max_image,N_cycle,WGA_th0,DAPI_th0] = nuclei_segDAPI(image_folder,image_type,DAPI_channel,resolution,scale0,Smin);
            end
    %         [mask_stack,signal_stack,~] = mask3D(seg_bw,protein_channel,DAPI_channel,RNA_channel,image_folder);

            em_mask = get_emmask(image_folder,[],DAPI_channel);
            figure(1)
                out_image = double(repmat(bwperim(em_mask),[1,1,3]));
                out_image(:,:,3) = out_image(:,:,3)+double(max_image(:,:,DAPI_channel))/max(max(double(max_image(:,:,DAPI_channel))));
                if scale0 < 1
                    out_image1 = imresize(out_image,scale0);
                else
                    out_image1 = out_image;
                end
                imshow(out_image1)
                title(image_folder,'Interpreter','none')

    %         max_size0 = size(max_image);
    %         max_size0(Mdim) = size(max_image,Mdim)/Nbin1;
    %         imcorr1 = corr_mask(max_size0,channel_name,resolution);
    %         Ntile = ones(size(max_size0));
    %         Ntile(Mdim) = Nbin1;
    %         imcorr = repmat(imcorr1,Ntile);
    %         max_image = uint16(double(max_image)./imcorr);

            [foci_bw,psize0,g_threshold0] = RNA_seg(seg_bw,max_image,RNA_channel,WGA_channel,resolution,image_folder,N_cycle,scale0);
            try
                [nucleus_RNA_profile,foci_RNA_profile,cytoplasmic_RNA_profile] = RNA_profile(seg_bw,cyto_bw,foci_bw,max_image,RNA_channel,resolution,image_folder,N_cycle,[],[],scale0);
            catch
                nucleus_RNA_profile=zeros(1,5);foci_RNA_profile=zeros(1,4);cytoplasmic_RNA_profile=zeros(1,2);
            end
            try
                [nucleus_protein_profile,cytoplasmic_protein_profile] = protein_profile(seg_bw,cyto_bw,max_image,protein_channel,image_folder,N_cycle);
            catch
                nucleus_protein_profile=zeros(1,2);cytoplasmic_protein_profile=zeros(1,2);
            end
                %         [nucleus_protein_profile,cytoplasmic_protein_profile] = protein_profile3(seg_bw,cyto_bw,max_image,protein_channel,mask_stack,signal_stack,image_folder,N_cycle,resolution);
            clear mask_stack signal_stack
            dual_profile(nucleus_protein_profile,nucleus_RNA_profile,foci_RNA_profile,image_folder,N_cycle);

            figure(16)
                temp_image = zeros(size(max_image,1),size(max_image,2),3);
                temp_image(:,:,2) = double(max_image(:,:,protein_channel))./max(max(double(max_image(:,:,protein_channel))));
                if scale0 < 1
                    temp_image1 = imresize(temp_image,scale0);
                else
                    temp_image1 = temp_image;
                end
                imshow(temp_image1);
                title(['Immunofluorescence signal: ',image_folder,', cycle = ',num2str(N_cycle)],'Interpreter','none');
            figure(17)
                if scale0 < 1
                    temp_image1 = imresize(double(max_image(:,:,RNA_channel))./max(max(double(max_image(:,:,RNA_channel)))),scale0);
                else
                    temp_image1 = double(max_image(:,:,RNA_channel))./max(max(double(max_image(:,:,RNA_channel))));
                end
                imshow(temp_image1);
                title(['FISH signal: ',image_folder,', cycle = ',num2str(N_cycle)],'Interpreter','none');


    %%% Output: %%%============================================================
            result_folder = [folder_name,out_folder];
            if exist(result_folder) ~= 7
                mkdir(result_folder);
            end
            saveas(1,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),emmask_add,figure_tail]);
            saveas(51,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),nu_add,seg_add,figure_tail]);
            saveas(3,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,th_add,figure_tail]);
            saveas(4,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,seg_add,figure_tail]);
            saveas(10,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,th_add,sel_add,figure_tail]);
            saveas(12,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),nu_add,reg_add,figure_tail]);
            saveas(13,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,reg_add,figure_tail]);
            saveas(15,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),nu_add,foci_add,reg_add,figure_tail]);
            saveas(16,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),immu_add,figure_tail]);
            saveas(17,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),fish_add,figure_tail]);
            try
                saveas(2,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,figure_tail]);
                saveas(5,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),fish_add,int_add,figure_tail]);
                saveas(6,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),fish_add,num_add,figure_tail]);
                saveas(7,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,num_add,figure_tail]);
                saveas(8,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,int_add,figure_tail]);
                saveas(9,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,int_add,cmp_add,figure_tail]);
                saveas(11,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),nu_add,fate_add,figure_tail]);
                saveas(14,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,fate_add,figure_tail]);
            end
    % %         saveas(62,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,fluc_add,figure_tail]);
    % %         saveas(63,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,bino_add,figure_tail]);
            save([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mat_tail],'image_folder','image_type','WGA_channel','DAPI_channel','RNA_channel','protein_channel','signal2_channel','resolution','seg_bw','cyto_bw','em_mask','max_image','foci_bw','nucleus_RNA_profile','foci_RNA_profile','cytoplasmic_RNA_profile','nucleus_protein_profile','cytoplasmic_protein_profile','N_cycle','psize0','g_threshold0','WGA_th0','DAPI_th0','resolutionz','all_color','Nbin','Mdim');
            if time_num==1
                save([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),'.mat'],'image_folder','image_type','WGA_channel','DAPI_channel','RNA_channel','protein_channel','signal2_channel','resolution','seg_bw','cyto_bw','em_mask','max_image','foci_bw','nucleus_RNA_profile','foci_RNA_profile','cytoplasmic_RNA_profile','nucleus_protein_profile','cytoplasmic_protein_profile','N_cycle','psize0','g_threshold0','WGA_th0','DAPI_th0','resolutionz','all_color','Nbin','Mdim');
            end
    % %         if ~isempty(nucleus_RNA_profile)
    % %             xlswrite([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),RNA_add,nu_add,output_tail],nucleus_RNA_profile);
    % %         end
    % %         if ~isempty(cytoplasmic_RNA_profile)
    % %             xlswrite([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),RNA_add,cyto_add,output_tail],cytoplasmic_RNA_profile);
    % %         end
    % %         if ~isempty(foci_RNA_profile)
    % %             xlswrite([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),RNA_add,foci_add,output_tail],foci_RNA_profile);
    % %         end
    % %         if ~isempty(nucleus_protein_profile)
    % %             xlswrite([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,nu_add,output_tail],nucleus_protein_profile);
    % %         end
    % %         if ~isempty(cytoplasmic_protein_profile)
    % %             xlswrite([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,cyto_add,output_tail],cytoplasmic_protein_profile);
    % %         end
            sub_num(list_J,13) = N_cycle;

            clear ('em_mask','seg_bw','cyto_bw','max_image','foci_bw','nucleus_RNA_profile','foci_RNA_profile','cytoplasmic_RNA_profile','nucleus_protein_profile','cytoplasmic_protein_profile','N_cycle','WGA_th0','DAPI_th0');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
            close all
        end
    end
    sub_list(:,4:4+size(sub_num,2)-1) = num2cell(sub_num);
    xlswrite([folder_name,in_folder,input_name],sub_list);
end
end

% clear global standard_data