clear all
close all

%% Parameter setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
list_name = 'RNA2list.xls';
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
nu_add = '_nucleus';
cyto_add = '_cytoplasm';
fate_add = '_fate';
emmask_add = '_emmask';
double_add = 'D';
% size0 = 5000;
bin0 = 1;
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
        
        if ~isempty(strfind(sub_list{list_J,2},double_add))
            Nbin = sub_num(list_J,2:3);
            Mdim = 1:2;
        else
            Mdim = sub_num(list_J,3);
            Nbin = sub_num(list_J,2);
        end
        
        WGA_channel = sub_num(list_J,6);
        DAPI_channel = sub_num(list_J,7);
        signal_channel = sub_num(list_J,8);
        RNA_channel = signal_channel;
        signal2_channel = sub_num(list_J,10);
        resolution = sub_num(list_J,9);
        resolutionz = sub_num(list_J,11);
        all_color = eval(folder_list{list_I,5});


        temp_name = dir([image_folder,image_type]);
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
        [seg_bw,cyto_bw,max_image,N_cycle,WGA_th0,DAPI_th0] = nuclei_segDAPI(image_folder,image_type,DAPI_channel,resolution,scale0);
%         [mask_stack,signal_stack,~] = mask3D(seg_bw,protein_channel,DAPI_channel,RNA_channel,image_folder);

        em_mask = get_emmask(image_folder,[],DAPI_channel);
        
% %         figure(1)
% %             out_image = double(repmat(bwperim(em_mask),[1,1,3]));
% %             out_image(:,:,3) = out_image(:,:,3)+double(max_image(:,:,DAPI_channel))/max(max(double(max_image(:,:,DAPI_channel))));
% %             if scale0 < 1
% %                 out_image1 = imresize(out_image,scale0);
% %             else
% %                 out_image1 = out_image;
% %             end
% %             imshow(out_image1)
% %             title(image_folder,'Interpreter','none')
% % 
        if max(max(bwlabel(seg_bw))) >= 20
            [foci_bw,psize0,g_threshold0] = RNA_seg(seg_bw,max_image,signal_channel,WGA_channel,resolution,image_folder,N_cycle,scale0);
            [nucleus_profile,foci_profile,cytoplasmic_profile] = RNA_profile(seg_bw,cyto_bw,foci_bw,max_image,signal_channel,resolution,image_folder,N_cycle,[],[],scale0);
        else
            foci_bw = false(size(seg_bw));
            psize0 = 0;
            g_threshold0 = 0;
            nucleus_profile = 0;
            foci_profile = 0;
            cytoplasmic_profile = 0;
            figure(3)
            figure(4)
            figure(5)
            figure(6)
            figure(7)
            figure(8)
            figure(9)
            figure(10)
            figure(11)
        end
        
% %         figure(17)
% %             if scale0 < 1
% %                 temp_image1 = imresize(double(max_image(:,:,RNA_channel))./max(max(double(max_image(:,:,RNA_channel)))),scale0);
% %             else
% %                 temp_image1 = double(max_image(:,:,RNA_channel))./max(max(double(max_image(:,:,RNA_channel))));
% %             end
% %             imshow(temp_image1);
% %             title(['FISH signal: ',image_folder,', cycle = ',num2str(N_cycle)],'Interpreter','none');
            
%%% Output: %%%============================================================
        result_folder = [folder_list{list_I,1},out_folder];
        if exist(result_folder) ~= 7
            mkdir(result_folder);
        end
% %         saveas(1,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),emmask_add,figure_tail]);
% %         saveas(51,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),nu_add,seg_add,figure_tail]);
% %         saveas(3,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,th_add,figure_tail]);
% %         saveas(4,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,seg_add,figure_tail]);
% %         saveas(5,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),fish_add,int_add,figure_tail]);
% %         saveas(6,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),fish_add,num_add,figure_tail]);
% %         saveas(7,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,num_add,figure_tail]);
% %         saveas(8,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,int_add,figure_tail]);
% %         saveas(9,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,int_add,cmp_add,figure_tail]);
% %         saveas(10,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,th_add,sel_add,figure_tail]);
% %         saveas(11,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),nu_add,fate_add,figure_tail]);
% %         saveas(17,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),fish_add,figure_tail]);
        save([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mat_tail],'em_mask','image_folder','image_type','WGA_channel','DAPI_channel','signal_channel','RNA_channel','signal2_channel','resolution','seg_bw','cyto_bw','max_image','foci_bw','nucleus_profile','foci_profile','cytoplasmic_profile','N_cycle','psize0','g_threshold0','WGA_th0','DAPI_th0','resolutionz','all_color','Nbin','Mdim');
        
% %         if ~isempty(nucleus_profile)
% %             xlswrite([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),nu_add,output_tail],nucleus_profile);
% %         end
% %         if ~isempty(cytoplasmic_profile)
% %             xlswrite([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),cyto_add,output_tail],cytoplasmic_profile);
% %         end
% %         if ~isempty(foci_profile)
% %             xlswrite([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,output_tail],foci_profile);
% %         end
        sub_num(list_J,13) = N_cycle;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        
    end
    xlswrite([folder_list{list_I,1},in_folder,input_name],cat(2,sub_list,num2cell(sub_num)));
end

clear all
close all
